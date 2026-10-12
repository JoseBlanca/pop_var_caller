//! **The repeat-tract fitting bench**: draw one motif length's classes from a known model, fit them
//! with each method asked for, and write how close each method came to the truth, how well it
//! predicts loci it was not fitted on, and what it cost (`doc/devel/ng/research/str_fit_research_plan.md`
//! §2, §3.1).
//!
//! One run is one **cell** of the plan's grid — a motif length, a number of samples, a read depth,
//! a truth and a seed. The drawing is the count-level simulator
//! (`ssr_fit::bench_fixtures::cohort_draw`); the methods, the curves and the held-out score are in
//! `ssr_fit::bench_fixtures::fit_bench`, whose documentation says how to add a method.
//!
//! ```text
//! cargo run --release --features bench-fixtures --example ng_str_fit_bench -- \
//!   --period 1 --samples 60 --depth 3 --seed 1 --truth pcr --methods A,B \
//!   --out rows.tsv --truth-out truth.tsv
//! ```
//!
//! # What it writes
//!
//! `--out`: one row per method, class and number, long format, tab-separated:
//!
//! | column | meaning |
//! |---|---|
//! | `cell`, `period`, `samples`, `depth`, `truth_variant`, `seed` | the cell |
//! | `method` | `A_raw`, `A_curves`, `B_raw`, `B_curves`, …; the bare name (`A`) on cost rows |
//! | `repeats` | the class's reference repeat count, or `all` |
//! | `parameter` | below |
//! | `truth`, `estimate`, `standard_error` | `NA` where there is none |
//!
//! Parameters: `level`, `shorter_share`, `fall_off`, `concentration`, `between_locus_share`
//! (`1 / (1 + concentration)`), `share_<offset>` for every allele class; the held-out score
//! `heldout_ln_per_read` and `heldout_ln_per_locus` (truth = the truth's own score on the same loci,
//! standard error = the Monte Carlo error) with `heldout_loci`, `heldout_reads` and
//! `heldout_loci_on_fallback` (loci of classes scored at calling's fallback rung, a flat shape at the
//! stated concentration, for want of their own shares or concentration); and the costs
//! `fit_seconds` (the fit on every locus), `heldout_fit_seconds`, `heldout_score_seconds`,
//! `draw_seconds` and `peak_rss_mib` (Linux only; elsewhere `NA`).
//!
//! **Peak memory is the process's**, so a run with several methods reports the largest of them. To
//! measure each method, run one method a process, as `scripts/research/str_fit_grid.sh` does; it
//! also records the process's peak through `/usr/bin/time`, which works on macOS too.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;
use rayon::prelude::*;

use pop_var_caller::parameter_estimation::joint::census::RECORDED_OFFSET_RANGE;
use pop_var_caller::parameter_estimation::joint::ssr_fit::bench_fixtures::cohort_draw::{
    ClassTruth, CohortModel, HomozygoteExcessTruth, LociPerClass, ParameterCurve, StutterRegime,
    draw_period, write_truth_tsv,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::bench_fixtures::fit_bench::{
    ClassEstimate, ClassNumbers, CountingFitErrorMethod, FoldPlan, HeldOutSettings, PeriodFit,
    Sampler, TheClimb, TheCountingFit, numbers_for_scoring, score_locus, with_the_curves,
    without_fold,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::{
    ALLELE_SPAN, SsrFitConfig, StratumEvidence, TractReads,
};

#[derive(Parser)]
#[command(about = "Score repeat-tract fitting methods on a cohort drawn from a known model")]
struct Cli {
    /// The motif length.
    #[arg(long, default_value_t = 1)]
    period: u8,
    /// Samples in the cohort.
    #[arg(long)]
    samples: Option<usize>,
    /// Mean spanning reads per sample per locus.
    #[arg(long)]
    depth: Option<f64>,
    #[arg(long, default_value_t = 1)]
    seed: u64,
    /// `pcr` or `pcr_free` (the default truth, `CohortModel::tomato_like`), or a TOML model file.
    #[arg(long, default_value = "pcr")]
    truth: String,
    /// Loci per class: `tomato2` (as tomato2 has them) or a number.
    #[arg(long)]
    loci: Option<String>,
    /// At most this many loci per class.
    #[arg(long)]
    loci_at_most: Option<usize>,
    /// Only these reference repeat counts, e.g. `8-16` or `8,10,12`.
    #[arg(long)]
    repeats: Option<String>,
    /// Every sample's homozygote excess, in place of the truth's.
    #[arg(long)]
    excess: Option<f64>,
    /// One stutter multiplier per library, comma-separated.
    #[arg(long, value_delimiter = ',')]
    library_multipliers: Vec<f64>,
    /// Tell the methods the libraries' true multipliers.
    #[arg(long)]
    multipliers_told: bool,
    /// The methods, comma-separated: `A` (the climb), `B` (the counting fit, spread from pairs),
    /// `B_within` (the counting fit, spread from each sample's own two copies).
    #[arg(long, value_delimiter = ',', default_value = "B")]
    methods: Vec<String>,
    /// Folds for the held-out score; 0 switches it off.
    #[arg(long, default_value_t = 5)]
    folds: usize,
    /// At most this many held-out loci scored per class.
    #[arg(long, default_value_t = 100)]
    heldout_loci: usize,
    /// The Monte Carlo error a locus's held-out score is drawn down to.
    #[arg(long, default_value_t = 0.05)]
    target_se: f64,
    /// The most Monte Carlo draws a locus gets.
    #[arg(long, default_value_t = 65_536)]
    max_draws: usize,
    /// The held-out sampler beyond two samples: `automatic`, `dirichlet` or `sequential`.
    #[arg(long, default_value = "automatic")]
    sampler: String,
    /// How many classes the climb fits at once (all the machine's threads by default).
    #[arg(long)]
    climb_strata_at_once: Option<usize>,
    /// Where the rows go.
    #[arg(long)]
    out: PathBuf,
    /// Where the truth goes, one row a class.
    #[arg(long)]
    truth_out: Option<PathBuf>,
    /// The cell's name in the rows; built from the cell's settings by default.
    #[arg(long)]
    cell: Option<String>,
    /// Where every held-out locus's score goes, one row a method variant and locus, with the read
    /// offset most of the locus's reads report and the share the fit gave that offset.
    #[arg(long)]
    locus_out: Option<PathBuf>,
    /// How far either side of the reference alleles are drawn, in place of the truth's (14).
    #[arg(long)]
    allele_reach: Option<i32>,
    /// One between-locus share F for every class, in place of the truth's curve.
    #[arg(long)]
    between_locus_share: Option<f64>,
}

/// The method a name on the command line stands for. **A new method adds a line here.**
fn method_named(name: &str, cli: &Cli) -> Box<dyn PeriodFit> {
    match name {
        "A" => {
            let threads = std::thread::available_parallelism().map_or(1, NonZeroUsize::get);
            Box::new(TheClimb {
                config: SsrFitConfig {
                    strata_at_once: NonZeroUsize::new(cli.climb_strata_at_once.unwrap_or(threads))
                        .unwrap_or(NonZeroUsize::MIN),
                    ..SsrFitConfig::default()
                },
            })
        }
        "B" => Box::new(TheCountingFit::default()),
        "B_within" => Box::new(TheCountingFit {
            pairs: false,
            ..TheCountingFit::default()
        }),
        "B_noerrors" => Box::new(TheCountingFit {
            errors: CountingFitErrorMethod::None,
            ..TheCountingFit::default()
        }),
        "B_noprior" => Box::new(TheCountingFit {
            share_prior_tracts: 0.0,
            ..TheCountingFit::default()
        }),
        "B_bootstrap" => Box::new(TheCountingFit {
            errors: CountingFitErrorMethod::Bootstrap(50),
            ..TheCountingFit::default()
        }),
        other => panic!(
            "no method is called {other}; the methods are A, B, B_within, B_noerrors, \
             B_noprior and B_bootstrap"
        ),
    }
}

fn repeats_wanted(text: &str) -> Vec<u64> {
    if let Some((low, high)) = text.split_once('-') {
        let low: u64 = low.parse().expect("a repeat count");
        let high: u64 = high.parse().expect("a repeat count");
        (low..=high).collect()
    } else {
        text.split(',')
            .map(|value| value.parse().expect("a repeat count"))
            .collect()
    }
}

fn model_of(cli: &Cli) -> CohortModel {
    let mut model = match cli.truth.as_str() {
        "pcr" => CohortModel::tomato_like(cli.period, StutterRegime::Pcr),
        "pcr_free" => CohortModel::tomato_like(cli.period, StutterRegime::PcrFree),
        path => toml::from_str(&std::fs::read_to_string(path).expect("the model file reads"))
            .expect("the model file parses"),
    };
    if let Some(samples) = cli.samples {
        model.samples = samples;
    }
    if let Some(depth) = cli.depth {
        model.depth = depth;
    }
    model.seed = cli.seed;
    model.loci = match (cli.loci.as_deref(), model.loci.clone()) {
        (Some("tomato2"), _) => LociPerClass::Tomato2 {
            at_most: cli.loci_at_most,
        },
        (Some(number), _) => LociPerClass::Each(
            number
                .parse::<usize>()
                .expect("loci a class: tomato2 or a number")
                .min(cli.loci_at_most.unwrap_or(usize::MAX)),
        ),
        (None, LociPerClass::Tomato2 { at_most }) => LociPerClass::Tomato2 {
            at_most: cli.loci_at_most.or(at_most),
        },
        (None, LociPerClass::Each(loci)) => {
            LociPerClass::Each(loci.min(cli.loci_at_most.unwrap_or(usize::MAX)))
        }
        (None, LociPerClass::Listed(listed)) => LociPerClass::Listed(
            listed
                .into_iter()
                .map(|(repeats, loci)| (repeats, loci.min(cli.loci_at_most.unwrap_or(usize::MAX))))
                .collect(),
        ),
    };
    if let Some(repeats) = &cli.repeats {
        model.repeat_counts = Some(repeats_wanted(repeats));
    }
    if let Some(reach) = cli.allele_reach {
        model.true_allele_reach = reach;
    }
    if let Some(share) = cli.between_locus_share {
        model.between_locus_share = ParameterCurve::Constant(share);
    }
    if let Some(excess) = cli.excess {
        model.homozygote_excess = HomozygoteExcessTruth::Each(excess);
    }
    if !cli.library_multipliers.is_empty() {
        model
            .library_multipliers
            .clone_from(&cli.library_multipliers);
    }
    if cli.multipliers_told {
        model.multipliers_told = true;
    }
    model
}

/// The rows a run writes, with the cell's columns in front of each.
struct Rows {
    cell_columns: String,
    lines: Vec<String>,
}

fn number(value: Option<f64>) -> String {
    match value {
        Some(value) if value.is_finite() => format!("{value:.8}"),
        _ => "NA".to_owned(),
    }
}

impl Rows {
    fn push(
        &mut self,
        method: &str,
        repeats: &str,
        parameter: &str,
        truth: Option<f64>,
        estimate: Option<f64>,
        error: Option<f64>,
    ) {
        self.lines.push(format!(
            "{}\t{method}\t{repeats}\t{parameter}\t{}\t{}\t{}",
            self.cell_columns,
            number(truth),
            number(estimate),
            number(error)
        ));
    }
}

/// Every number of every class against its truth.
fn recovery_rows(rows: &mut Rows, method: &str, estimates: &[ClassEstimate], truth: &[ClassTruth]) {
    for (estimate, truth) in estimates.iter().zip(truth) {
        let repeats = truth.stratum.reference_repeats.to_string();
        let slippage = estimate.slippage;
        let errors = &estimate.errors;
        rows.push(
            method,
            &repeats,
            "level",
            Some(truth.slippage.level),
            slippage.map(|s| s.level),
            errors.level,
        );
        rows.push(
            method,
            &repeats,
            "shorter_share",
            Some(truth.slippage.shorter_share),
            slippage.map(|s| s.shorter_share),
            errors.shorter_share,
        );
        rows.push(
            method,
            &repeats,
            "fall_off",
            Some(truth.slippage.fall_off),
            slippage.map(|s| s.fall_off),
            errors.fall_off,
        );
        rows.push(
            method,
            &repeats,
            "concentration",
            Some(truth.concentration),
            estimate.concentration,
            errors.concentration,
        );
        rows.push(
            method,
            &repeats,
            "between_locus_share",
            Some(truth.between_locus_share()),
            estimate.concentration.map(|c| 1.0 / (1.0 + c)),
            errors
                .concentration
                .zip(estimate.concentration)
                .map(|(error, c)| error / ((1.0 + c) * (1.0 + c))),
        );
        let shares = truth.shares_in_window(ALLELE_SPAN);
        for (index, share) in shares.iter().enumerate() {
            let offset = index as i32 - ALLELE_SPAN;
            rows.push(
                method,
                &repeats,
                &format!("share_{offset}"),
                Some(*share),
                estimate
                    .allele_length_variation
                    .as_ref()
                    .and_then(|alv| alv.get(index).copied()),
                errors.allele_length_variation.get(index).copied().flatten(),
            );
        }
    }
}

/// One held-out locus's score under a method: its class, the method's log-likelihood and that
/// score's Monte Carlo variance, the truth's log-likelihood, and its reads.
type HeldOutLocus = (usize, f64, f64, f64, u64);

/// Held-out sums for one method variant, per class.
#[derive(Clone, Default)]
struct HeldOut {
    ln: f64,
    variance: f64,
    truth_ln: f64,
    reads: u64,
    loci: usize,
    /// Of `loci`, those whose class was scored at calling's fallback rung (no concentration or no
    /// shares of its own).
    fallback_loci: usize,
}

fn heldout_rows(rows: &mut Rows, method: &str, per_class: &[HeldOut], truth: &[ClassTruth]) {
    let mut all = HeldOut::default();
    let write = |rows: &mut Rows, repeats: &str, sums: &HeldOut| {
        if sums.loci == 0 {
            return;
        }
        let reads = sums.reads.max(1) as f64;
        let loci = sums.loci as f64;
        let error = sums.variance.sqrt();
        rows.push(
            method,
            repeats,
            "heldout_ln_per_read",
            Some(sums.truth_ln / reads),
            Some(sums.ln / reads),
            Some(error / reads),
        );
        rows.push(
            method,
            repeats,
            "heldout_ln_per_locus",
            Some(sums.truth_ln / loci),
            Some(sums.ln / loci),
            Some(error / loci),
        );
        rows.push(method, repeats, "heldout_loci", None, Some(loci), None);
        rows.push(
            method,
            repeats,
            "heldout_loci_on_fallback",
            None,
            Some(sums.fallback_loci as f64),
            None,
        );
        rows.push(
            method,
            repeats,
            "heldout_reads",
            None,
            Some(sums.reads as f64),
            None,
        );
    };
    for (sums, truth) in per_class.iter().zip(truth) {
        write(rows, &truth.stratum.reference_repeats.to_string(), sums);
        all.ln += sums.ln;
        all.variance += sums.variance;
        all.truth_ln += sums.truth_ln;
        all.reads += sums.reads;
        all.loci += sums.loci;
        all.fallback_loci += sums.fallback_loci;
    }
    write(rows, "all", &all);
}

fn peak_rss_mib() -> Option<f64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    status.lines().find_map(|line| {
        line.strip_prefix("VmHWM:")?
            .split_whitespace()
            .next()?
            .parse::<f64>()
            .ok()
            .map(|kib| kib / 1024.0)
    })
}

/// A locus's Monte Carlo seed, the same for every method so their scores share their draws.
fn locus_seed(settings: &HeldOutSettings, repeats: u64, locus: usize) -> u64 {
    settings.seed ^ (repeats << 40) ^ (locus as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
}

fn main() {
    let cli = Cli::parse();
    let model = model_of(&cli);
    let truth_variant = match cli.truth.as_str() {
        "pcr" | "pcr_free" => cli.truth.clone(),
        path => PathBuf::from(path).file_stem().map_or_else(
            || "file".to_owned(),
            |stem| stem.to_string_lossy().into_owned(),
        ),
    };
    let cell = cli.cell.clone().unwrap_or_else(|| {
        format!(
            "p{}_n{}_d{}_{}",
            model.period, model.samples, model.depth, truth_variant
        )
    });

    let drawing = Instant::now();
    let drawn = draw_period(&model, RECORDED_OFFSET_RANGE);
    let draw_seconds = drawing.elapsed().as_secs_f64();
    let rows_with_reads: usize = drawn
        .strata
        .iter()
        .map(StratumEvidence::rows_with_reads)
        .sum();
    eprintln!(
        "{cell} seed {}: drew {} classes, {} loci, {rows_with_reads} sample-locus rows with reads \
         in {draw_seconds:.1} s",
        model.seed,
        drawn.strata.len(),
        drawn.truth.iter().map(|t| t.loci).sum::<usize>(),
    );
    if let Some(path) = &cli.truth_out {
        let file = File::create(path).expect("the truth file opens");
        write_truth_tsv(&drawn.truth, ALLELE_SPAN, BufWriter::new(file))
            .expect("the truth file writes");
    }

    let mut rows = Rows {
        cell_columns: format!(
            "{cell}\t{}\t{}\t{}\t{truth_variant}\t{}",
            model.period, model.samples, model.depth, model.seed
        ),
        lines: Vec::new(),
    };
    rows.push(
        "draw",
        "all",
        "draw_seconds",
        None,
        Some(draw_seconds),
        None,
    );
    rows.push(
        "draw",
        "all",
        "rows_with_reads",
        None,
        Some(rows_with_reads as f64),
        None,
    );

    let settings = HeldOutSettings {
        folds: cli.folds,
        loci_scored_per_class: cli.heldout_loci,
        target_standard_error: cli.target_se,
        max_draws: cli.max_draws,
        sampler: match cli.sampler.as_str() {
            "automatic" => Sampler::Automatic,
            "dirichlet" => Sampler::DirichletOnly,
            "sequential" => Sampler::SequentialOnly,
            other => panic!("no sampler is called {other}"),
        },
        ..HeldOutSettings::default()
    };
    let curve_config = SsrFitConfig::default();
    let plans: Vec<FoldPlan> = drawn
        .strata
        .iter()
        .map(|evidence| FoldPlan::of(evidence, &settings, evidence.stratum.reference_repeats))
        .collect();
    // Every scored locus, with its fold.
    let scored: Vec<(usize, usize, usize)> = plans
        .iter()
        .enumerate()
        .flat_map(|(class, plan)| {
            plan.scored
                .iter()
                .enumerate()
                .filter(|(_, scored)| **scored)
                .map(move |(locus, _)| (class, locus, plan.fold_of_locus[locus]))
        })
        .collect();
    let score_of = |class: usize, locus: usize, numbers: &ClassNumbers| {
        let evidence = &drawn.strata[class];
        score_locus(
            &evidence.tracts[locus],
            numbers,
            &evidence.libraries,
            &drawn.homozygote_excess,
            evidence.read_span,
            &settings,
            locus_seed(&settings, evidence.stratum.reference_repeats, locus),
        )
    };
    // The truth's own score on every scored locus, once.
    let truth_scores: Vec<f64> = if settings.folds > 0 {
        let started = Instant::now();
        let scores: Vec<(f64, usize, f64)> = scored
            .par_iter()
            .map(|(class, locus, _)| {
                let score = score_of(
                    *class,
                    *locus,
                    &ClassNumbers::of_truth(&drawn.truth[*class], ALLELE_SPAN),
                );
                (score.ln_likelihood, score.draws, score.variance.sqrt())
            })
            .collect();
        let mut draws: Vec<usize> = scores.iter().map(|(_, draws, _)| *draws).collect();
        draws.sort_unstable();
        let short_of_target = scores
            .iter()
            .filter(|(_, _, error)| *error >= settings.target_standard_error)
            .count();
        let quantile = |q: f64| {
            draws
                .get(((draws.len() as f64 - 1.0) * q) as usize)
                .copied()
        };
        eprintln!(
            "{cell}: scored the truth on {} held-out loci in {:.1} s; Monte Carlo draws a locus: \
             median {:?}, 90th percentile {:?}, most {:?}; {short_of_target} loci short of the \
             error target",
            scored.len(),
            started.elapsed().as_secs_f64(),
            quantile(0.5),
            quantile(0.9),
            draws.last(),
        );
        rows.push(
            "truth",
            "all",
            "heldout_median_draws",
            None,
            quantile(0.5).map(|d| d as f64),
            None,
        );
        rows.push(
            "truth",
            "all",
            "heldout_loci_short_of_target",
            None,
            Some(short_of_target as f64),
            None,
        );
        let scores: Vec<f64> = scores.into_iter().map(|(ln, _, _)| ln).collect();
        rows.push(
            "truth",
            "all",
            "heldout_score_seconds",
            None,
            Some(started.elapsed().as_secs_f64()),
            None,
        );
        scores
    } else {
        Vec::new()
    };

    let mut locus_lines: Vec<String> = Vec::new();
    // An empty name is no method: a run that only draws, to measure the draw's own memory.
    for name in cli.methods.iter().filter(|name| !name.is_empty()) {
        let method = method_named(name, &cli);
        let base = method.name();
        let started = Instant::now();
        let raw = method.fit_period(&drawn.strata, &drawn.homozygote_excess, ALLELE_SPAN);
        let fit_seconds = started.elapsed().as_secs_f64();
        eprintln!("{cell}: {base} fitted every locus in {fit_seconds:.1} s");
        rows.push(&base, "all", "fit_seconds", None, Some(fit_seconds), None);
        let variants = |raw: Vec<ClassEstimate>, classes: &[StratumEvidence]| {
            if method.draws_its_own_curves() {
                vec![(base.clone(), raw)]
            } else {
                let smoothed = with_the_curves(&raw, classes, &curve_config);
                vec![
                    (format!("{base}_raw"), raw),
                    (format!("{base}_curves"), smoothed),
                ]
            }
        };
        let full = variants(raw, &drawn.strata);
        for (variant, estimates) in &full {
            recovery_rows(&mut rows, variant, estimates, &drawn.truth);
        }

        if settings.folds > 0 {
            let mut sums = vec![vec![HeldOut::default(); drawn.strata.len()]; full.len()];
            let mut fit_seconds = 0.0;
            let mut score_seconds = 0.0;
            for fold in 0..settings.folds {
                let training: Vec<StratumEvidence> = drawn
                    .strata
                    .iter()
                    .zip(&plans)
                    .map(|(evidence, plan)| without_fold(evidence, plan, fold))
                    .collect();
                let started = Instant::now();
                let raw = method.fit_period(&training, &drawn.homozygote_excess, ALLELE_SPAN);
                let fitted = variants(raw, &training);
                fit_seconds += started.elapsed().as_secs_f64();
                let started = Instant::now();
                let in_fold: Vec<(usize, usize, f64)> = scored
                    .iter()
                    .zip(&truth_scores)
                    .filter(|((_, _, locus_fold), _)| *locus_fold == fold)
                    .map(|((class, locus, _), truth)| (*class, *locus, *truth))
                    .collect();
                for (variant_index, (variant, estimates)) in fitted.iter().enumerate() {
                    // A class without its own shares or concentration is scored at calling's
                    // fallback rung and counted, not skipped.
                    let scored_numbers = numbers_for_scoring(estimates);
                    let numbers: Vec<Option<ClassNumbers>> = scored_numbers
                        .iter()
                        .map(|scored| scored.as_ref().map(|s| s.numbers.clone()))
                        .collect();
                    let on_the_fallback: Vec<bool> = scored_numbers
                        .iter()
                        .map(|scored| scored.as_ref().is_some_and(|s| s.on_the_fallback))
                        .collect();
                    let scores: Vec<Option<HeldOutLocus>> = in_fold
                        .par_iter()
                        .map(|(class, locus, truth)| {
                            let numbers = numbers[*class].as_ref()?;
                            let score = score_of(*class, *locus, numbers);

                            Some((
                                *class,
                                score.ln_likelihood,
                                score.variance,
                                *truth,
                                score.reads,
                            ))
                        })
                        .collect();
                    if cli.locus_out.is_some() {
                        for ((class, locus, _), score) in in_fold.iter().zip(&scores) {
                            let Some((_, ln, _, truth, reads)) = score else {
                                continue;
                            };
                            let (offset, share) = most_reported_offset(
                                &drawn.strata[*class].tracts[*locus],
                                drawn.strata[*class].read_span,
                                numbers[*class].as_ref(),
                            );
                            locus_lines.push(format!(
                                "{variant}\t{}\t{locus}\t{ln}\t{truth}\t{reads}\t{offset}\t{share:e}",
                                drawn.strata[*class].stratum.reference_repeats
                            ));
                        }
                    }
                    for (class, ln, variance, truth, reads) in scores.into_iter().flatten() {
                        let into = &mut sums[variant_index][class];
                        into.ln += ln;
                        into.variance += variance;
                        into.truth_ln += truth;
                        into.reads += reads;
                        into.loci += 1;
                        into.fallback_loci += usize::from(on_the_fallback[class]);
                    }
                }
                score_seconds += started.elapsed().as_secs_f64();
            }
            eprintln!(
                "{cell}: {base} held out {} folds: fits {fit_seconds:.1} s, scores \
                 {score_seconds:.1} s",
                settings.folds
            );
            rows.push(
                &base,
                "all",
                "heldout_fit_seconds",
                None,
                Some(fit_seconds),
                None,
            );
            rows.push(
                &base,
                "all",
                "heldout_score_seconds",
                None,
                Some(score_seconds),
                None,
            );
            for ((variant, _), per_class) in full.iter().zip(&sums) {
                heldout_rows(&mut rows, variant, per_class, &drawn.truth);
            }
        }
    }
    if let Some(path) = &cli.locus_out {
        let mut out = BufWriter::new(File::create(path).expect("the locus file opens"));
        writeln!(
            out,
            "variant\trepeats\tlocus\theldout_ln\ttruth_ln\treads\tmost_reported_offset\tfitted_share"
        )
        .expect("the locus file writes");
        for line in &locus_lines {
            writeln!(out, "{line}").expect("the locus file writes");
        }
    }
    rows.push("process", "all", "peak_rss_mib", None, peak_rss_mib(), None);

    let partial = cli.out.with_extension("partial");
    {
        let mut out = BufWriter::new(File::create(&partial).expect("the rows file opens"));
        writeln!(
            out,
            "cell\tperiod\tsamples\tdepth\ttruth_variant\tseed\tmethod\trepeats\tparameter\ttruth\t\
             estimate\tstandard_error"
        )
        .expect("the rows file writes");
        for line in &rows.lines {
            writeln!(out, "{line}").expect("the rows file writes");
        }
        out.flush().expect("the rows file writes");
    }
    std::fs::rename(&partial, &cli.out).expect("the rows file moves into place");
}

/// **The read offset most of a locus's reads report**, summed over its samples, and the share the
/// fit gave that offset's allele class; at an end bucket, the smaller of the two classes it can
/// hold (`±8` and `±9`), since either could be the allele.
fn most_reported_offset(
    tract: &TractReads,
    read_span: i32,
    numbers: Option<&ClassNumbers>,
) -> (i32, f64) {
    let mut totals = vec![0_u64; (2 * read_span + 1) as usize];
    for sample in &tract.samples {
        for (_, counts) in &sample.by_library {
            for (into, count) in totals.iter_mut().zip(counts) {
                *into += u64::from(*count);
            }
        }
    }
    let bucket = totals
        .iter()
        .enumerate()
        .max_by_key(|(_, count)| **count)
        .map_or(0, |(index, _)| index);
    let offset = bucket as i32 - read_span;
    let share = numbers.map_or(f64::NAN, |numbers| {
        let span = (numbers.allele_length_variation.len() as i32 - 1) / 2;
        let at = |offset: i32| numbers.allele_length_variation[(offset + span) as usize];
        if offset.abs() == read_span && span > read_span {
            at(offset).min(at(offset + offset.signum()))
        } else {
            at(offset)
        }
    });
    (offset, share)
}
