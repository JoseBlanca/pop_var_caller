//! **Can each library's stutter multiplier be recovered?** A prototype of
//! `doc/devel/ng/spec/str_library_stutter.md`, run on a cohort of the simulated stuttering
//! benchmark, whose libraries stutter at a known multiplier each.
//!
//! It reads the cohort's psp files as `estimate-parameters` does, takes every repeat class's other
//! numbers — the direction and fall-off of a slip, the allele-length frequencies and their spread —
//! from the pooled fit the benchmark already made, and each sample's homozygote excess from the same
//! file. Then it fits the class rates and the libraries' multipliers together, with those numbers
//! held (`ssr_fit::library_stutter`), and compares both with the simulator's truth:
//!
//! - each library's multiplier against its true one, both scaled to a read-weighted geometric mean of
//!   one, with how many land within two standard errors;
//! - each class rate against the simulator's rate for that class times the true multipliers' mean,
//!   beside the pooled fit's rate against the same.
//!
//! ```text
//! cargo run --release --example ng_str_library_multiplier_prototype -- \
//!   --reference RUN/reference.fa --catalog RUN/reference.fa.repeats.parquet --psp RUN/psps \
//!   --output RUN/unused --fitted RUN/ng_fit.parameters.toml \
//!   --truth-libraries RUN/libraries.tsv --truth-rates RUN/stutter_rates.tsv
//! ```

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;

use pop_var_caller::calling::parameters_file::ParametersFile;
use pop_var_caller::cli::estimate_parameters::{
    EstimateParametersArgs, OpenedCensuses, contig_id_of, open_the_censuses,
};
use pop_var_caller::float;
use pop_var_caller::parameter_estimation::joint::ssr_fit::library_stutter::{
    HeldClass, fit_rates_and_multipliers,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::{
    Slippage, SsrFitConfig, StratumEvidence,
};
use pop_var_caller::run::{
    CohortTractStrata, every_read_group_pooled, the_tract_strata_of_a_cohort,
};

#[derive(Parser)]
#[command(
    about = "Fit each library's stutter multiplier on a simulated cohort and compare with truth"
)]
struct Cli {
    #[command(flatten)]
    args: EstimateParametersArgs,
    /// The pooled fit's parameters file, whose other numbers are held.
    #[arg(long)]
    fitted: PathBuf,
    /// The simulator's `libraries.tsv`: sample, level, multiplier, shorter share.
    #[arg(long)]
    truth_libraries: PathBuf,
    /// The simulator's `stutter_rates.tsv`: level, period, repeats, rate.
    #[arg(long)]
    truth_rates: PathBuf,
    /// Pseudo-reads pulling a multiplier toward one (spec §3.3).
    #[arg(long, default_value_t = 20.0)]
    pull: f64,
    /// Rounds of the alternation at most.
    #[arg(long, default_value_t = 50)]
    rounds: usize,
}

fn read_tsv(path: &PathBuf) -> Vec<Vec<String>> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()))
        .lines()
        .skip(1)
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect()
}

fn median(mut values: Vec<f64>) -> f64 {
    values.sort_by(f64::total_cmp);
    values[values.len() / 2]
}

/// The geometric mean of `values` weighted by `weights`.
fn weighted_geometric_mean(values: &[f64], weights: &[f64]) -> f64 {
    let total: f64 = weights.iter().sum();
    float::exp(
        values
            .iter()
            .zip(weights)
            .map(|(v, w)| w * float::ln(*v))
            .sum::<f64>()
            / total,
    )
}

fn main() {
    let cli = Cli::parse();
    let started = Instant::now();
    let OpenedCensuses {
        mut evidence, plan, ..
    } = open_the_censuses(&cli.args).unwrap_or_else(|error| panic!("{error:?}"));
    let names: Vec<String> = evidence.sample_names().map(str::to_owned).collect();
    let contig_of = contig_id_of(&plan);
    let pooled = every_read_group_pooled(&evidence);
    let CohortTractStrata { strata, .. } =
        the_tract_strata_of_a_cohort(&mut evidence, &plan.loci, &contig_of, &pooled)
            .unwrap_or_else(|error| panic!("{error}"));
    eprintln!(
        "gathered {} classes from {} samples in {:.1} s",
        strata.len(),
        names.len(),
        started.elapsed().as_secs_f64()
    );

    let file = ParametersFile::from_toml(
        &std::fs::read_to_string(&cli.fitted).expect("the fitted parameters file reads"),
    )
    .expect("the fitted parameters file parses");
    let excess_by_name: BTreeMap<&str, f64> = file
        .inbreeding
        .by_sample
        .iter()
        .map(|row| (row.sample.as_str(), row.inbreeding_coefficient.value))
        .collect();
    let excess: Vec<f64> = names
        .iter()
        .map(|name| *excess_by_name.get(name.as_str()).unwrap_or(&0.0))
        .collect();

    // Each class the pooled fit has its own numbers for, with them held.
    let mut kept: Vec<StratumEvidence> = Vec::new();
    let mut held: Vec<HeldClass> = Vec::new();
    let tracts = &file.repeat_tracts;
    for stratum in strata {
        let key = (stratum.stratum.period, stratum.stratum.reference_repeats);
        let slippage = tracts
            .slippage_by_stratum_and_group
            .iter()
            .find(|row| (row.period, row.reference_repeats) == key && row.slippage_group == 0);
        let spectrum = tracts
            .length_spectrum_by_stratum
            .iter()
            .find(|row| (row.period, row.reference_repeats) == key);
        let (Some(slippage), Some(spectrum)) = (slippage, spectrum) else {
            continue;
        };
        held.push(HeldClass {
            slippage: vec![Slippage {
                level: slippage.share_of_reads_that_slip,
                shorter_share: slippage.shorter_share,
                fall_off: slippage.fall_off,
            }],
            length_spectrum: spectrum.shares_by_repeat_offset.clone(),
            concentration: spectrum.concentration,
        });
        kept.push(stratum);
    }
    let pooled_rates: Vec<f64> = held.iter().map(|class| class.slippage[0].level).collect();
    let config = SsrFitConfig::default();
    eprintln!("{} classes have the pooled fit's own numbers", kept.len());

    let fitting = Instant::now();
    let fit = fit_rates_and_multipliers(
        &mut kept,
        &mut held,
        &excess,
        &config,
        cli.pull,
        cli.rounds,
        |round, _, moved| {
            eprintln!(
                "round {round}: largest multiplier move {moved:.3} standard errors ({:.1} s)",
                fitting.elapsed().as_secs_f64()
            );
        },
    );

    // Library → sample, from the reads: every library belongs to one sample.
    let mut sample_of_library: BTreeMap<usize, usize> = BTreeMap::new();
    for stratum in &kept {
        for row in stratum.tracts.iter().flat_map(|tract| &tract.samples) {
            for (library, _) in &row.by_library {
                sample_of_library.insert(library.get(), row.sample as usize);
            }
        }
    }
    let truth: BTreeMap<String, (String, f64)> = read_tsv(&cli.truth_libraries)
        .into_iter()
        .map(|row| {
            (
                row[0].clone(),
                (row[1].clone(), row[2].parse().expect("a multiplier")),
            )
        })
        .collect();
    let libraries: Vec<usize> = (0..fit.multipliers.len())
        .filter(|library| sample_of_library.contains_key(library))
        .collect();
    let true_multiplier: Vec<f64> = libraries
        .iter()
        .map(|library| truth[&names[sample_of_library[library]]].1)
        .collect();
    let level_name = truth[&names[sample_of_library[&libraries[0]]]].0.clone();
    let weights: Vec<f64> = libraries
        .iter()
        .map(|library| fit.reads[*library].max(1e-9))
        .collect();
    let true_mean = weighted_geometric_mean(&true_multiplier, &weights);

    println!("library\tsample\ttrue\tfitted\tstandard_error\tz\tslipped_reads");
    let mut within_two = 0;
    let mut ratios = Vec::new();
    let mut pooled_ratios = Vec::new();
    for (index, library) in libraries.iter().enumerate() {
        let truth = true_multiplier[index] / true_mean;
        let fitted = fit.multipliers[*library];
        let error = fit.standard_errors[*library];
        let z = (fitted - truth) / error;
        within_two += usize::from(z.abs() <= 2.0);
        ratios.push(float::ln(fitted / truth).abs());
        pooled_ratios.push(float::ln(1.0 / truth).abs());
        println!(
            "{library}\t{}\t{truth:.3}\t{fitted:.3}\t{error:.3}\t{z:+.2}\t{:.0}",
            names[sample_of_library[library]], fit.slipped_reads[*library]
        );
    }

    let true_rates: BTreeMap<(u8, u64), f64> = read_tsv(&cli.truth_rates)
        .into_iter()
        .filter(|row| row[0] == level_name)
        .map(|row| {
            (
                (
                    row[1].parse().expect("a period"),
                    row[2].parse().expect("repeats"),
                ),
                row[3].parse().expect("a rate"),
            )
        })
        .collect();
    let mut rate_ratios = Vec::new();
    let mut pooled_rate_ratios = Vec::new();
    for ((stratum, fitted), pooled) in kept.iter().zip(&fit.class_rates).zip(&pooled_rates) {
        let key = (stratum.stratum.period, stratum.stratum.reference_repeats);
        if let Some(rate) = true_rates.get(&key) {
            let typical = rate * true_mean;
            rate_ratios.push(fitted / typical);
            pooled_rate_ratios.push(pooled / typical);
        }
    }

    println!();
    println!(
        "summary\tlibraries {}\tclasses {}\trounds {}\tlast move {:.3} errors\tfit {:.1} s",
        libraries.len(),
        kept.len(),
        fit.rounds,
        fit.last_move_in_errors,
        fitting.elapsed().as_secs_f64()
    );
    println!(
        "multipliers\twithin two errors {within_two} of {}\tmedian |ln(fitted/true)| {:.3}\tthe pooled fit's (every multiplier one) {:.3}",
        libraries.len(),
        median(ratios),
        median(pooled_ratios)
    );
    if !rate_ratios.is_empty() {
        println!(
            "class rates\t{} with a true rate\tmedian fitted/true {:.3}\tthe pooled fit's {:.3}",
            rate_ratios.len(),
            median(rate_ratios),
            median(pooled_rate_ratios)
        );
    }
}
