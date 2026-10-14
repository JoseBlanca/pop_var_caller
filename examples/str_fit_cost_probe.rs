//! **Where the repeat-tract fits' time goes, a class at a time, on real or drawn evidence** — for
//! the cost study `doc/devel/ng/research/str_fit/cost_2026-10-14.md`, which asks why the fitting
//! bench costs more a sample-locus row than the counting fit's real runs.
//!
//! Loads one cohort's classes, either a real one from its censuses (`real`, the arguments of
//! `ng_str_counting_fit`) or a cell of the bench's tomato-like truth (`drawn`), and prints a row a
//! class for the counting fit (B), timed as the real driver runs it: classes one after another,
//! each spreading its loci over the thread pool. Every time is **process CPU time** (every thread's,
//! from the process CPU clock), with wall time beside it; CPU time is the measure the study uses, because
//! other jobs share the machine.
//!
//! Each class's row gives what the evidence holds (rows with reads, reads a row, read-length
//! buckets with reads a row) and what one unit of each kind of work costs on it:
//!
//! - `estep_cpu`: one round of counting (the E-step over every row) at the fit's start;
//! - `pairs_cpu`: one reading of the spread off pairs of samples (`fit_spread_by_pairs`: the
//!   pairs' table built once, then the 61-point search) at the start;
//! - `fit_cpu`: the whole fit (`fit_class_by_counting_with_pairs`, what the driver runs), with
//!   `rounds_in_all` (rounds of counting) and `readings` (readings of the spread);
//! - `errors_cpu`: the sandwich errors the bench adds (`--errors`).
//!
//! `--methods D,D_tract_prior` also fits each motif length by penalised counting, as the bench
//! does, and prints one row a motif length with its own split of the time.
//!
//! ```text
//! cargo run --release --features bench-fixtures --example str_fit_cost_probe -- \
//!     drawn --period 1 --samples 2000 --depth 3 --seed 1 --loci-at-most 100
//! cargo run --release --features bench-fixtures --example str_fit_cost_probe -- \
//!     real --reference REF --catalog CATALOG --psp PSPS --output unused --homozygote-excess 0.98
//! ```

use std::time::Instant;

use clap::{Parser, Subcommand};

use pop_var_caller::cli::estimate_parameters::{
    EstimateParametersArgs, OpenedCensuses, contig_id_of, open_the_censuses,
};
use pop_var_caller::parameter_estimation::joint::census::RECORDED_OFFSET_RANGE;
use pop_var_caller::parameter_estimation::joint::ssr_fit::bench_fixtures::cohort_draw::{
    CohortModel, LociPerClass, StutterRegime, draw_period,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::composite_fit::{
    expected_counts, fit_class_by_counting_with_pairs, fit_spread_by_pairs,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::counting_fit_errors::sandwich_errors;
use pop_var_caller::parameter_estimation::joint::ssr_fit::penalised_fit::{
    PenalisedFitConfig, fit_period_by_penalised_counting,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::{
    DEFAULT_REFUSAL_FLOOR, Slippage, SsrFitConfig, StratumEvidence,
};
use pop_var_caller::run::{
    CohortTractStrata, every_read_group_pooled, the_tract_strata_of_a_cohort,
};

#[derive(Parser)]
#[command(about = "Where the repeat-tract fits' time goes, a class at a time")]
struct Cli {
    #[command(subcommand)]
    source: Source,
    /// How far either side of the reference length alleles may sit (the shipped ±9 by default).
    #[arg(long, global = true)]
    allele_span: Option<i32>,
    /// Also time the sandwich errors the bench adds to B.
    #[arg(long, global = true)]
    errors: bool,
    /// Skip the whole fit (only the unit costs).
    #[arg(long, global = true)]
    units_only: bool,
    /// Skip the unit costs (one round, one reading of F): only the whole fit and its errors.
    #[arg(long, global = true)]
    no_units: bool,
    /// Comma-separated: B (default), D, D_tract_prior.
    #[arg(long, global = true, default_value = "B")]
    methods: String,
    /// Only the classes of this motif length.
    #[arg(long, global = true)]
    only_period: Option<u8>,
    /// The counting fit's tolerance: it stops when no number moves by more than this in a round.
    #[arg(long, global = true, default_value_t = 1e-7)]
    tolerance: f64,
    /// Write every fitted number and error at full precision here, a line each, so that two
    /// builds' fits can be compared number by number.
    #[arg(long, global = true)]
    numbers: Option<std::path::PathBuf>,
}

/// Lines `method period repeats name value`, every value in Rust's shortest exact form.
#[derive(Default)]
struct Numbers(Vec<String>);

impl Numbers {
    fn push(&mut self, method: &str, evidence: &StratumEvidence, name: &str, value: Option<f64>) {
        self.0.push(format!(
            "{method}\t{}\t{}\t{name}\t{}",
            evidence.stratum.period,
            evidence.stratum.reference_repeats,
            value.map_or_else(|| "NA".to_owned(), |v| format!("{v:?}"))
        ));
    }
}

#[derive(Subcommand)]
enum Source {
    /// A real cohort, read from its censuses as `ng_str_counting_fit` reads it.
    Real {
        #[command(flatten)]
        args: EstimateParametersArgs,
        /// One homozygote excess for every sample.
        #[arg(long, default_value_t = 0.98)]
        homozygote_excess: f64,
    },
    /// One cell of the bench's tomato-like PCR truth.
    Drawn {
        #[arg(long)]
        period: u8,
        #[arg(long)]
        samples: usize,
        #[arg(long)]
        depth: f64,
        #[arg(long, default_value_t = 1)]
        seed: u64,
        #[arg(long, default_value_t = 100)]
        loci_at_most: usize,
    },
}

/// This process's CPU time so far, every thread's, in seconds.
fn cpu_seconds() -> f64 {
    let time = rustix::time::clock_gettime(rustix::time::ClockId::ProcessCPUTime);
    time.tv_sec as f64 + time.tv_nsec as f64 * 1e-9
}

/// CPU and wall seconds a closure took, and what it returned.
fn timed<T>(work: impl FnOnce() -> T) -> (T, f64, f64) {
    let (cpu, wall) = (cpu_seconds(), Instant::now());
    let out = work();
    (out, cpu_seconds() - cpu, wall.elapsed().as_secs_f64())
}

/// Reads a row holds and read-length buckets with at least one read, summed over the class; then
/// how many distinct rows of reads the class holds (two samples whose reads fall in the same
/// buckets in the same numbers are one), over the whole class and summed tract by tract.
fn read_shape(evidence: &StratumEvidence) -> (u64, u64, usize, usize) {
    let mut reads = 0_u64;
    let mut buckets = 0_u64;
    let mut in_class = std::collections::HashSet::new();
    let mut by_tract = 0;
    for tract in &evidence.tracts {
        let mut in_tract = std::collections::HashSet::new();
        for sample in &tract.samples {
            let width = sample.by_library.first().map_or(0, |(_, r)| r.len());
            let mut pooled = vec![0_u32; width];
            for (_, library) in &sample.by_library {
                for (into, count) in pooled.iter_mut().zip(library) {
                    *into += count;
                }
            }
            reads += pooled.iter().map(|count| u64::from(*count)).sum::<u64>();
            buckets += pooled.iter().filter(|count| **count > 0).count() as u64;
            in_tract.insert(pooled.clone());
            in_class.insert(pooled);
        }
        by_tract += in_tract.len();
    }
    (reads, buckets, in_class.len(), by_tract)
}

fn main() {
    let cli = Cli::parse();
    let (load, load_cpu, load_wall) = timed(|| load(&cli.source));
    let (strata, excess, label) = load;
    eprintln!(
        "{label}: {} classes loaded in {load_wall:.1} s wall, {load_cpu:.1} s CPU; {} rayon threads",
        strata.len(),
        rayon::current_num_threads()
    );
    let mut config = SsrFitConfig::default();
    if let Some(span) = cli.allele_span {
        config.allele_span = span;
    }
    let start = Slippage {
        level: 0.05,
        shorter_share: 0.5,
        fall_off: 0.3,
    };
    let methods: Vec<&str> = cli.methods.split(',').collect();
    let classes = (2 * config.allele_span + 1) as usize;
    let flat = vec![1.0 / classes as f64; classes];
    let mut numbers = Numbers::default();
    let wanted = |evidence: &StratumEvidence| {
        cli.only_period
            .is_none_or(|period| evidence.stratum.period == period)
    };

    if methods.contains(&"B") {
        println!(
            "source\tspan\tperiod\trepeats\ttracts\trows\treads_per_row\tbuckets_per_row\testep_cpu\testep_wall\tpairs_cpu\tpairs_wall\tfit_cpu\tfit_wall\trounds_in_all\treadings\terrors_cpu\tlevel\tconcentration\tdistinct_rows_in_class\tdistinct_rows_by_tract"
        );
        let mut totals = (0_usize, 0.0, 0.0, 0_usize, 0_usize);
        let mut error_totals = (0.0, 0.0);
        for evidence in strata.iter().filter(|e| wanted(e)) {
            if evidence.tracts_with_reads() < DEFAULT_REFUSAL_FLOOR {
                continue;
            }
            let rows = evidence.rows_with_reads();
            let (reads, buckets, in_class, by_tract) = read_shape(evidence);
            let (estep_cpu, estep_wall, pairs_cpu, pairs_wall) = if cli.no_units {
                (f64::NAN, f64::NAN, f64::NAN, f64::NAN)
            } else {
                let (_, estep_cpu, estep_wall) =
                    timed(|| expected_counts(evidence, start, &flat, 1.0, &excess, &config));
                let (_, pairs_cpu, pairs_wall) =
                    timed(|| fit_spread_by_pairs(evidence, start, &flat, &excess, &config));
                (estep_cpu, estep_wall, pairs_cpu, pairs_wall)
            };
            let (fit, fit_cpu, fit_wall) = if cli.units_only {
                (None, f64::NAN, f64::NAN)
            } else {
                let (fit, cpu, wall) = timed(|| {
                    fit_class_by_counting_with_pairs(
                        evidence,
                        start,
                        &excess,
                        &config,
                        2000,
                        cli.tolerance,
                        10,
                    )
                });
                (Some(fit), cpu, wall)
            };
            let errors_cpu = match (&fit, cli.errors) {
                (Some(fit), true) => {
                    let (errors, cpu, wall) =
                        timed(|| sandwich_errors(evidence, fit, &excess, &config));
                    error_totals.0 += cpu;
                    error_totals.1 += wall;
                    let errors = errors.sandwich;
                    numbers.push("B", evidence, "error_level", errors.level.value());
                    numbers.push("B", evidence, "error_shorter", errors.shorter_share.value());
                    numbers.push("B", evidence, "error_fall_off", errors.fall_off.value());
                    numbers.push(
                        "B",
                        evidence,
                        "error_concentration",
                        errors.concentration.value(),
                    );
                    for (index, error) in errors.allele_length_variation.iter().enumerate() {
                        numbers.push(
                            "B",
                            evidence,
                            &format!("error_share_{index}"),
                            error.value(),
                        );
                    }
                    cpu
                }
                _ => f64::NAN,
            };
            if let Some(fit) = &fit {
                numbers.push("B", evidence, "level", Some(fit.slippage.level));
                numbers.push("B", evidence, "shorter", Some(fit.slippage.shorter_share));
                numbers.push("B", evidence, "fall_off", Some(fit.slippage.fall_off));
                numbers.push("B", evidence, "concentration", fit.concentration);
                numbers.push(
                    "B",
                    evidence,
                    "rounds_in_all",
                    Some(fit.rounds_in_all as f64),
                );
                for (index, share) in fit.allele_length_variation.iter().enumerate() {
                    numbers.push("B", evidence, &format!("share_{index}"), Some(*share));
                }
                totals.0 += rows;
                totals.1 += fit_cpu;
                totals.2 += fit_wall;
                totals.3 += fit.rounds_in_all * rows;
                totals.4 += fit.spread_readings * rows;
            }
            println!(
                "{label}\t{}\t{}\t{}\t{}\t{rows}\t{:.3}\t{:.3}\t{estep_cpu:.3}\t{estep_wall:.3}\t{pairs_cpu:.3}\t{pairs_wall:.3}\t{fit_cpu:.3}\t{fit_wall:.3}\t{}\t{}\t{errors_cpu:.3}\t{}\t{}\t{in_class}\t{by_tract}",
                config.allele_span,
                evidence.stratum.period,
                evidence.stratum.reference_repeats,
                evidence.tracts_with_reads(),
                reads as f64 / rows.max(1) as f64,
                buckets as f64 / rows.max(1) as f64,
                fit.as_ref().map_or(0, |f| f.rounds_in_all),
                fit.as_ref().map_or(0, |f| f.spread_readings),
                fit.as_ref()
                    .map_or_else(|| "NA".to_owned(), |f| format!("{:.6}", f.slippage.level)),
                fit.as_ref()
                    .and_then(|f| f.concentration)
                    .map_or_else(|| "NA".to_owned(), |c| format!("{c:.6}")),
            );
        }
        eprintln!(
            "B over {} rows: {:.1} s CPU, {:.1} s wall; {:.0} row-rounds, {:.0} row-readings; \
             errors {:.1} s CPU, {:.1} s wall",
            totals.0,
            totals.1,
            totals.2,
            totals.3 as f64,
            totals.4 as f64,
            error_totals.0,
            error_totals.1
        );
    }

    for method in methods.iter().filter(|m| m.starts_with('D')) {
        let fit_config = PenalisedFitConfig {
            start,
            max_rounds: 2000,
            tolerance: cli.tolerance,
            pair_alternations: 10,
            standard_errors: true,
            allele_prior_tracts: if *method == "D_tract_prior" { 1.0 } else { 0.0 },
            ..PenalisedFitConfig::default()
        };
        let mut periods: Vec<u8> = strata.iter().map(|e| e.stratum.period).collect();
        periods.sort_unstable();
        periods.dedup();
        for period in periods {
            if cli.only_period.is_some_and(|only| only != period) {
                continue;
            }
            let classes: Vec<StratumEvidence> = strata
                .iter()
                .filter(|e| e.stratum.period == period)
                .cloned()
                .collect();
            let rows: usize = classes.iter().map(StratumEvidence::rows_with_reads).sum();
            let (fit, cpu, wall) =
                timed(|| fit_period_by_penalised_counting(&classes, &excess, &config, &fit_config));
            for (evidence, class) in classes.iter().zip(&fit.classes) {
                let Some(class) = class else { continue };
                numbers.push(method, evidence, "level", Some(class.slippage.level));
                numbers.push(
                    method,
                    evidence,
                    "shorter",
                    Some(class.slippage.shorter_share),
                );
                numbers.push(method, evidence, "fall_off", Some(class.slippage.fall_off));
                numbers.push(method, evidence, "concentration", class.concentration);
                for (index, share) in class.allele_length_variation.iter().enumerate() {
                    numbers.push(method, evidence, &format!("share_{index}"), Some(*share));
                }
                if let Some(errors) = &class.errors {
                    numbers.push(method, evidence, "error_level", errors.level);
                    numbers.push(
                        method,
                        evidence,
                        "error_concentration",
                        errors.concentration,
                    );
                    numbers.push(
                        method,
                        evidence,
                        "error_between",
                        errors.between_locus_share,
                    );
                    for (index, error) in errors.allele_length_variation.iter().enumerate() {
                        numbers.push(
                            method,
                            evidence,
                            &format!("error_share_{index}"),
                            Some(*error),
                        );
                    }
                }
            }
            println!(
                "{label}\t{method}\tspan {}\tperiod {period}\trows {rows}\tcpu {cpu:.1}\twall {wall:.1}\trounds {}\tclass_rounds {}\texpectation_s {:.1}\tstrength_s {:.1}\tpairs_s {:.1}\tpair_alternations {}\terrors_s {:.1}",
                config.allele_span,
                fit.rounds,
                fit.class_rounds_counted,
                fit.expectation_seconds,
                fit.strength_search_seconds,
                fit.pair_seconds,
                fit.pair_alternations,
                fit.error_seconds,
            );
        }
    }
    if let Some(path) = &cli.numbers {
        std::fs::write(path, numbers.0.join("\n") + "\n").expect("the numbers file writes");
    }
}

/// The classes, each sample's homozygote excess, and a label for the rows.
fn load(source: &Source) -> (Vec<StratumEvidence>, Vec<f64>, String) {
    match source {
        Source::Real {
            args,
            homozygote_excess,
        } => {
            let OpenedCensuses {
                mut evidence, plan, ..
            } = open_the_censuses(args).unwrap_or_else(|error| panic!("{error:?}"));
            let samples = evidence.sample_names().count();
            let contig_of = contig_id_of(&plan);
            let pooled = every_read_group_pooled(&evidence);
            let CohortTractStrata { strata, .. } =
                the_tract_strata_of_a_cohort(&mut evidence, &plan.loci, &contig_of, &pooled)
                    .unwrap_or_else(|error| panic!("{error:?}"));
            (
                strata,
                vec![*homozygote_excess; samples],
                format!("real_n{samples}"),
            )
        }
        Source::Drawn {
            period,
            samples,
            depth,
            seed,
            loci_at_most,
        } => {
            let mut model = CohortModel::tomato_like(*period, StutterRegime::Pcr);
            model.samples = *samples;
            model.depth = *depth;
            model.seed = *seed;
            model.loci = LociPerClass::Tomato2 {
                at_most: Some(*loci_at_most),
            };
            let drawn = draw_period(&model, RECORDED_OFFSET_RANGE);
            (
                drawn.strata,
                drawn.homozygote_excess,
                format!("drawn_p{period}_n{samples}_d{depth}_s{seed}_l{loci_at_most}"),
            )
        }
    }
}
