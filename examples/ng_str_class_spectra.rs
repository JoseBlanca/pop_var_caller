//! **Every repeat class's own allele-length frequencies, from a real cohort**, for judging what shape
//! they have (the owner's "route 2", 2026-10-10).
//!
//! Gathers a cohort's repeat-tract evidence as `estimate-parameters` does and fits every class on
//! its own tracts, with the curves left undrawn, the SNP/indel half skipped and each sample's
//! homozygote excess read from a fitted parameters file. Writes one row a fitted class: its tracts,
//! its slippage numbers, its spread (concentration) and its share at each offset from the reference
//! length, each share with its standard error where the fit placed one.
//!
//! ```text
//! cargo run --release --example ng_str_class_spectra -- --reference REF --catalog CATALOG \
//!   --psp PSPS --output unused --fitted COHORT.parameters.toml > classes.tsv
//! ```

use std::collections::BTreeMap;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::time::Instant;

use clap::Parser;

use pop_var_caller::calling::parameters_file::ParametersFile;
use pop_var_caller::cli::estimate_parameters::{
    EstimateParametersArgs, OpenedCensuses, contig_id_of, open_the_censuses,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::{
    SsrFitConfig, StratumError, StratumOutcome, fit_strata,
};
use pop_var_caller::run::{
    CohortTractStrata, every_read_group_pooled, the_tract_strata_of_a_cohort,
};

#[derive(Parser)]
#[command(about = "Fit every repeat class of a cohort and write its allele-length frequencies")]
struct Cli {
    #[command(flatten)]
    args: EstimateParametersArgs,
    /// A parameters file whose homozygote excesses are used; every sample at zero without one.
    #[arg(long)]
    fitted: Option<PathBuf>,
    /// One homozygote excess for every sample, where no fitted file names them.
    #[arg(long, conflicts_with = "fitted")]
    homozygote_excess: Option<f64>,
    /// How far either side of the reference length alleles may sit (the fit's `ALLELE_SPAN` by
    /// default).
    #[arg(long)]
    allele_span: Option<i32>,
    /// Fit only these classes, written `period:repeats` and separated by commas, one at a time with
    /// every thread on each, and print each class's time; every class, several at once, without.
    #[arg(long, value_delimiter = ',')]
    classes: Vec<String>,
}

fn error_of(error: &StratumError) -> String {
    match error {
        StratumError::Estimated(value) => format!("{value:.6}"),
        _ => "NA".to_owned(),
    }
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
    let CohortTractStrata { mut strata, .. } =
        the_tract_strata_of_a_cohort(&mut evidence, &plan.loci, &contig_of, &pooled)
            .unwrap_or_else(|error| panic!("{error:?}"));
    if !cli.classes.is_empty() {
        let wanted: Vec<(u8, u64)> = cli
            .classes
            .iter()
            .map(|class| {
                let (period, repeats) = class.split_once(':').expect("a class is period:repeats");
                (
                    period.parse().expect("a whole period"),
                    repeats.parse().expect("a whole repeat count"),
                )
            })
            .collect();
        strata.retain(|stratum| {
            wanted.contains(&(stratum.stratum.period, stratum.stratum.reference_repeats))
        });
    }
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
            let missing = names
                .iter()
                .filter(|name| !by_name.contains_key(name.as_str()))
                .count();
            eprintln!(
                "homozygote excess from the file for {} of {} samples",
                names.len() - missing,
                names.len()
            );
            names
                .iter()
                .map(|name| *by_name.get(name.as_str()).unwrap_or(&0.0))
                .collect()
        }
        None => vec![cli.homozygote_excess.unwrap_or(0.0); names.len()],
    };

    let mut config = SsrFitConfig::default();
    config.curve.draw_curves = false;
    if let Some(span) = cli.allele_span {
        config.allele_span = span;
    }
    if cli.classes.is_empty() {
        config.strata_at_once = std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN);
    }
    eprintln!(
        "{} classes from {} samples gathered in {:.1} s; fitting",
        strata.len(),
        names.len(),
        started.elapsed().as_secs_f64()
    );
    let outcomes = if cli.classes.is_empty() {
        fit_strata(&strata, &excess, &config)
    } else {
        strata
            .iter()
            .flat_map(|stratum| {
                let class_started = Instant::now();
                let outcome = fit_strata(std::slice::from_ref(stratum), &excess, &config);
                eprintln!(
                    "class {}:{}, {} tracts, window ±{}: {:.1} s",
                    stratum.stratum.period,
                    stratum.stratum.reference_repeats,
                    stratum.tracts.len(),
                    config.allele_span,
                    class_started.elapsed().as_secs_f64()
                );
                outcome
            })
            .collect()
    };
    eprintln!("fitted in {:.1} s", started.elapsed().as_secs_f64());

    let span = config.allele_span;
    let mut header = vec![
        "period".to_owned(),
        "repeats".to_owned(),
        "tracts".to_owned(),
        "reads".to_owned(),
        "level".to_owned(),
        "shorter_share".to_owned(),
        "fall_off".to_owned(),
        "concentration".to_owned(),
        "concentration_error".to_owned(),
    ];
    header.extend((-span..=span).map(|offset| format!("share_{offset}")));
    header.extend((-span..=span).map(|offset| format!("error_{offset}")));
    println!("{}", header.join("\t"));
    for outcome in &outcomes {
        let StratumOutcome::Fitted(fit) = outcome else {
            continue;
        };
        let Some(slippage) = fit.slippage[0] else {
            continue;
        };
        let errors = fit.standard_errors.as_ref();
        let mut row = vec![
            fit.stratum.period.to_string(),
            fit.stratum.reference_repeats.to_string(),
            fit.tracts_of_its_own.to_string(),
            fit.reads_crossing.to_string(),
            format!("{:.6}", slippage.level),
            format!("{:.6}", slippage.shorter_share),
            format!("{:.6}", slippage.fall_off),
            format!("{:.6}", fit.concentration),
            errors.map_or("NA".to_owned(), |errors| error_of(&errors.concentration)),
        ];
        row.extend(
            fit.length_spectrum
                .iter()
                .map(|share| format!("{share:.8}")),
        );
        row.extend((0..fit.length_spectrum.len()).map(|class| {
            errors.map_or("NA".to_owned(), |errors| {
                error_of(&errors.length_spectrum[class])
            })
        }));
        println!("{}", row.join("\t"));
    }
}
