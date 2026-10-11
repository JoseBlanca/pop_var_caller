//! **Every repeat class of a cohort fitted by counting** (`ssr_fit::composite_fit`): stutter,
//! allele-length shares and spread by expectation maximisation, with no average over tract
//! frequencies (the owner's routes 1 and 2, 2026-10-11). Writes one row a class in the columns
//! `ng_str_class_spectra` writes for the usual fit, so the two can be compared class by class; the
//! error columns are `NA`.
//!
//! ```text
//! cargo run --release --example ng_str_counting_fit -- --reference REF --catalog CATALOG \
//!   --psp PSPS --output unused --homozygote-excess 0.98 > counted.tsv
//! ```

use std::time::Instant;

use clap::Parser;

use pop_var_caller::cli::estimate_parameters::{
    EstimateParametersArgs, OpenedCensuses, contig_id_of, open_the_censuses,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::composite_fit::fit_class_by_counting;
use pop_var_caller::parameter_estimation::joint::ssr_fit::{
    DEFAULT_REFUSAL_FLOOR, Slippage, SsrFitConfig,
};
use pop_var_caller::run::{
    CohortTractStrata, every_read_group_pooled, the_tract_strata_of_a_cohort,
};

#[derive(Parser)]
#[command(about = "Fit every repeat class of a cohort by counting")]
struct Cli {
    #[command(flatten)]
    args: EstimateParametersArgs,
    /// One homozygote excess for every sample.
    #[arg(long, default_value_t = 0.0)]
    homozygote_excess: f64,
    /// Rounds of counting at most.
    #[arg(long, default_value_t = 2000)]
    rounds: usize,
}

fn main() {
    let cli = Cli::parse();
    let started = Instant::now();
    let OpenedCensuses {
        mut evidence, plan, ..
    } = open_the_censuses(&cli.args).unwrap_or_else(|error| panic!("{error:?}"));
    let samples = evidence.sample_names().count();
    let contig_of = contig_id_of(&plan);
    let pooled = every_read_group_pooled(&evidence);
    let CohortTractStrata { strata, .. } =
        the_tract_strata_of_a_cohort(&mut evidence, &plan.loci, &contig_of, &pooled)
            .unwrap_or_else(|error| panic!("{error:?}"));
    let config = SsrFitConfig::default();
    let excess = vec![cli.homozygote_excess; samples];
    let start = Slippage {
        level: 0.05,
        shorter_share: 0.5,
        fall_off: 0.3,
    };
    let span = config.allele_span;
    let mut header: Vec<String> = [
        "period",
        "repeats",
        "tracts",
        "reads",
        "level",
        "shorter_share",
        "fall_off",
        "concentration",
        "concentration_error",
    ]
    .map(str::to_owned)
    .to_vec();
    header.extend((-span..=span).map(|offset| format!("share_{offset}")));
    header.extend((-span..=span).map(|offset| format!("error_{offset}")));
    header.extend(["rounds".to_owned(), "seconds".to_owned()]);
    println!("{}", header.join("\t"));
    for evidence in &strata {
        if evidence.tracts_with_reads() < DEFAULT_REFUSAL_FLOOR {
            continue;
        }
        let fitting = Instant::now();
        let fit = fit_class_by_counting(evidence, start, &excess, &config, cli.rounds, 1e-7);
        let mut row = vec![
            evidence.stratum.period.to_string(),
            evidence.stratum.reference_repeats.to_string(),
            evidence.tracts_with_reads().to_string(),
            evidence.spanning_reads().to_string(),
            format!("{:.6}", fit.slippage.level),
            format!("{:.6}", fit.slippage.shorter_share),
            format!("{:.6}", fit.slippage.fall_off),
            format!("{:.6}", fit.concentration),
            "NA".to_owned(),
        ];
        row.extend(
            fit.length_spectrum
                .iter()
                .map(|share| format!("{share:.8}")),
        );
        row.extend((0..fit.length_spectrum.len()).map(|_| "NA".to_owned()));
        row.push(fit.rounds.to_string());
        row.push(format!("{:.2}", fitting.elapsed().as_secs_f64()));
        println!("{}", row.join("\t"));
    }
    eprintln!(
        "{} classes from {samples} samples in {:.1} s",
        strata.len(),
        started.elapsed().as_secs_f64()
    );
}
