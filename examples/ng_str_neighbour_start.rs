//! **Does one walk from a neighbouring class's answer end where the best of the usual three ends?**
//! The owner's first idea for a faster repeat-tract fit (2026-10-10), measured.
//!
//! Every repeat class of a cohort is fitted as `estimate-parameters` fits it, from the usual
//! starting points, with the curves left undrawn so each answer is the class's own. Then each class
//! that has a fitted neighbour of the same motif length at a smaller repeat count is climbed once
//! more, from that neighbour's answer, and the two are compared: the rounds each took, and how much
//! lower the neighbour start's total log-likelihood ended than the best of the usual walks.
//!
//! ```text
//! cargo run --release --example ng_str_neighbour_start -- \
//!   --reference RUN/reference.fa --catalog RUN/reference.fa.repeats.parquet --psp PSPS \
//!   --output unused --fitted RUN/ng_fit.parameters.toml
//! ```

use std::collections::BTreeMap;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;
use rayon::prelude::*;

use pop_var_caller::calling::parameters_file::ParametersFile;
use pop_var_caller::cli::estimate_parameters::{
    EstimateParametersArgs, OpenedCensuses, contig_id_of, open_the_censuses,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::neighbour_start::{
    mean_log_likelihood_at, walk_from_a_neighbour,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::{
    Slippage, SsrFitConfig, StratumError, StratumFit, StratumOutcome, fit_strata,
};
use pop_var_caller::run::{
    CohortTractStrata, every_read_group_pooled, the_tract_strata_of_a_cohort,
};

#[derive(Parser)]
#[command(about = "Compare one walk from a neighbour's answer with the usual starting points")]
struct Cli {
    #[command(flatten)]
    args: EstimateParametersArgs,
    /// A parameters file whose homozygote excesses are used; every sample at zero without one.
    #[arg(long)]
    fitted: Option<PathBuf>,
    /// Points of the precise average both answers are scored again with.
    #[arg(long, default_value_t = 65_536)]
    precise_points: usize,
}

fn main() {
    let cli = Cli::parse();
    let OpenedCensuses {
        mut evidence, plan, ..
    } = open_the_censuses(&cli.args).unwrap_or_else(|error| panic!("{error:?}"));
    let names: Vec<String> = evidence.sample_names().map(str::to_owned).collect();
    let contig_of = contig_id_of(&plan);
    let pooled = every_read_group_pooled(&evidence);
    let CohortTractStrata { strata, .. } =
        the_tract_strata_of_a_cohort(&mut evidence, &plan.loci, &contig_of, &pooled)
            .unwrap_or_else(|error| panic!("{error:?}"));
    let excess: Vec<f64> = match &cli.fitted {
        Some(path) => {
            let file = ParametersFile::from_toml(
                &std::fs::read_to_string(path).expect("the parameters file reads"),
            )
            .expect("the parameters file parses");
            let by_name: BTreeMap<&str, f64> = file
                .inbreeding
                .by_sample
                .iter()
                .map(|row| (row.sample.as_str(), row.inbreeding_coefficient.value))
                .collect();
            names
                .iter()
                .map(|name| *by_name.get(name.as_str()).unwrap_or(&0.0))
                .collect()
        }
        None => vec![0.0; names.len()],
    };

    let mut config = SsrFitConfig::default();
    config.curve.draw_curves = false;
    config.strata_at_once = std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN);
    let started = Instant::now();
    let outcomes = fit_strata(&strata, &excess, &config);
    let usual_seconds = started.elapsed().as_secs_f64();
    let fitted: BTreeMap<(u8, u64), &StratumFit> = outcomes
        .iter()
        .filter_map(|outcome| match outcome {
            StratumOutcome::Fitted(fit) => Some((
                (fit.stratum.period, fit.stratum.reference_repeats),
                fit.as_ref(),
            )),
            _ => None,
        })
        .collect();

    // Each fitted class with a fitted neighbour below it in the same period.
    let pairs: Vec<(usize, &StratumFit, &StratumFit)> = strata
        .iter()
        .enumerate()
        .filter_map(|(index, evidence)| {
            let key = (evidence.stratum.period, evidence.stratum.reference_repeats);
            let own = fitted.get(&key)?;
            let neighbour = fitted
                .range((key.0, 0)..key)
                .next_back()
                .map(|(_, fit)| *fit)?;
            Some((index, *own, neighbour))
        })
        .collect();
    let started = Instant::now();
    let walks: Vec<_> = pairs
        .par_iter()
        .map(|(index, _, neighbour)| {
            walk_from_a_neighbour(&strata[*index], &excess, &config, neighbour)
        })
        .collect();
    let neighbour_seconds = started.elapsed().as_secs_f64();

    println!(
        "period\trepeats\tneighbour_repeats\ttracts\tusual_rounds\tneighbour_rounds\tneighbour_ending\tloss\tprecise_loss\tlevel_usual\tlevel_neighbour\tlevel_gap_in_errors"
    );
    let (mut usual_rounds, mut neighbour_rounds) = (0_u32, 0_u32);
    let mut losses: Vec<f64> = Vec::new();
    let mut precise_losses: Vec<f64> = Vec::new();
    let scoring = Instant::now();
    for ((index, own, neighbour), walk) in pairs.iter().zip(&walks) {
        let Some(walk) = walk else { continue };
        let evidence = &strata[*index];
        let tracts = evidence.tracts.len() as f64;
        let rounds: u32 = own.walks.iter().map(|record| record.rounds).sum();
        usual_rounds += rounds;
        neighbour_rounds += walk.record.rounds;
        let loss = (own.log_likelihood_a_tract - walk.log_likelihood_a_tract) * tracts;
        losses.push(loss);
        // Both answers scored again with the precise average, every slippage group the usual fit
        // left without numbers taking the neighbour walk's.
        let own_slippage: Vec<Slippage> = own
            .slippage
            .iter()
            .zip(&walk.slippage)
            .map(|(own, walked)| own.unwrap_or(*walked))
            .collect();
        let precise = |slippage: &[Slippage], spectrum: &[f64], concentration: f64| {
            mean_log_likelihood_at(
                evidence,
                slippage,
                spectrum,
                concentration,
                &excess,
                &config,
                cli.precise_points,
            )
        };
        let precise_loss = (precise(
            &own_slippage,
            &own.allele_length_variation,
            own.concentration
                .expect("the climb fills the concentration"),
        ) - precise(
            &walk.slippage,
            &walk.allele_length_variation,
            walk.concentration,
        )) * tracts;
        precise_losses.push(precise_loss);
        let usual_level = own.slippage[0].map_or(f64::NAN, |s| s.level);
        let error = own
            .standard_errors
            .as_ref()
            .and_then(|errors| errors.slippage[0].as_ref())
            .and_then(|errors| match errors.level {
                StratumError::Estimated(error) => Some(error),
                _ => None,
            })
            .unwrap_or(f64::NAN);
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:?}\t{:.3}\t{:.3}\t{:.5}\t{:.5}\t{:.2}",
            own.stratum.period,
            own.stratum.reference_repeats,
            neighbour.stratum.reference_repeats,
            tracts,
            rounds,
            walk.record.rounds,
            walk.record.ending,
            loss,
            precise_loss,
            usual_level,
            walk.slippage[0].level,
            (walk.slippage[0].level - usual_level) / error
        );
    }
    let summarise = |name: &str, values: &[f64]| {
        let at_most = |limit: f64| values.iter().filter(|value| **value <= limit).count();
        let at_least = |limit: f64| values.iter().filter(|value| **value >= limit).count();
        println!(
            "{name}\tneighbour worse by more than 1: {}\twithin 1: {}\tneighbour better by more than 1: {}\tlargest worse {:.2}\tlargest better {:.2}",
            values.len() - at_most(1.0),
            at_most(1.0) - (values.len() - at_least(-1.0)),
            values.len() - at_least(-1.0),
            values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            -values.iter().copied().fold(f64::INFINITY, f64::min)
        );
    };
    println!();
    println!(
        "summary\tclasses compared {}\tusual rounds {usual_rounds} ({usual_seconds:.1} s for every class)\tneighbour rounds {neighbour_rounds} ({neighbour_seconds:.1} s)\tprecise scoring {:.1} s at {} points",
        losses.len(),
        scoring.elapsed().as_secs_f64(),
        cli.precise_points
    );
    summarise("loss at 256 points", &losses);
    summarise("loss at the precise average", &precise_losses);
}
