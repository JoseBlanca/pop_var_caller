//! **A large cohort's repeat-tract strata fitted on every sample and on the grown subset, side by
//! side** (`fit_precision.md` §4.5 item 3; plan step D3).
//!
//! The parameters fit reads each repeat-tract stratum of a cohort of more than 256 samples from a
//! subset of its samples, grown until the stratum's slippage level is measured to 2% of itself
//! (spec §4.4). This tool fits chosen strata of a real cohort both ways and prints, for each, the
//! three slippage numbers and the concentration with their standard errors, the samples the subset
//! reached, and how long each fit took — so the subset's answer can be judged against the whole
//! cohort's before it is trusted.
//!
//! ```text
//! cargo run --release --example ng_ssr_subset_against_every_sample -- \
//!     --reference ref.fa --psp cohort/ --output unused.toml [--inbreeding 0.9]
//! ```
//!
//! - **The flags are `estimate-parameters`'s**, and the cohort is opened and checked exactly as
//!   that command opens it. `--output` is required by the shared flags and not written to.
//! - **The homozygote excess** each sample's genotypes are drawn with is `--inbreeding` for every
//!   sample, or 0 when it is not given; the SNP/indel fit that measures it per sample is not run
//!   (on kimura it takes hours). Both fits of a stratum use the same value, so the comparison is
//!   fair whichever is given; tomato is mostly inbred, so give its usual value there.
//! - **Which strata**: `STRATA=1:8,2:10` (period in bases, reference repeats), or by default five
//!   spread across the cohort's sizes: the stratum with the most tracts with reads and those an
//!   eighth, a quarter, a half and three quarters of the way down, among the strata holding at
//!   least the 8 tracts with reads a fit needs.
//! - `FIRST_SUBSET`, `MIN_SAMPLES_A_GROUP` and `TARGET` set the first subset's size (256), the
//!   readers a slippage group is topped up to (8) and the level's relative error the subset grows to
//!   (0.02), so the tool can be tried on a small cohort.
//!
//! **What it costs.** Each stratum is fitted on the subset first, printed at once, then on every
//! sample, with a line before each fit; the run's state is never more than one fit old. On kimura
//! the fit on every sample is the long part: spec §1.1 records seven strata taking 11 h 20 min there
//! under the old five-round limit, and walks now get up to 40 rounds. Each fit on every sample says
//! the likelihood table it holds before it starts.
//!
//! **Reading "apart".** It is the subset's value minus the whole fit's, in the standard deviation
//! that difference has when the subset's samples are part of the whole, `√(subset error² − whole
//! error²)`. Drawing a subset moves a number by about one of these; over twenty numbers about one is
//! expected beyond two by chance, and several well beyond two say the subset is losing something.
//! It takes the errors at their word, so on an inbred cohort give `--inbreeding`. When the subset
//! grew to every sample both fits read the same samples, and the difference is given in the whole
//! fit's errors instead; "-" means one side has no error.
//!
//! Spec §4.5 item 3 asks that the two levels differ by less than their errors and that the time
//! saved be reported: the "apart" of the slippage level, and the last line's total, are those.

use std::time::Instant;

use clap::Parser;

use pop_var_caller::cli::estimate_parameters::{
    EstimateParametersArgs, OpenedCensuses, contig_id_of, open_the_censuses,
};
use pop_var_caller::parameter_estimation::joint::sample_order::{SAMPLE_ORDER_SEED, sample_order};
use pop_var_caller::parameter_estimation::joint::ssr_fit::{
    DEFAULT_REFUSAL_FLOOR, SsrFitConfig, StratumError, StratumEvidence, StratumFit, StratumOutcome,
    fit_strata_on_sample_subsets, fit_stratum,
};
use pop_var_caller::run::{
    CohortTractStrata, every_read_group_pooled, the_tract_strata_of_a_cohort,
};
use pop_var_caller::types::InbreedingF;

#[derive(Parser)]
#[command(about = "Fit chosen repeat-tract strata on every sample and on the grown subset")]
struct Cli {
    #[command(flatten)]
    args: EstimateParametersArgs,
}

/// The strata `STRATA` names, as (period, reference repeats); `None` when it is not set.
fn strata_asked_for() -> Option<Vec<(u8, u64)>> {
    let value = std::env::var("STRATA").ok()?;
    Some(
        value
            .split(',')
            .map(|item| {
                let (period, repeats) = item
                    .trim()
                    .split_once(':')
                    .unwrap_or_else(|| panic!("STRATA item {item:?} is not period:repeats"));
                (
                    period.parse().expect("a period in bases"),
                    repeats.parse().expect("a reference repeat count"),
                )
            })
            .collect(),
    )
}

/// A number and its error, as the fit gives them.
fn with_error(value: f64, error: StratumError) -> (f64, Option<f64>) {
    (value, error.value())
}

/// The three slippage numbers and the concentration of a fit, each with its error.
fn numbers_of(fit: &StratumFit) -> Vec<(&'static str, f64, Option<f64>)> {
    let slippage = fit.slippage[0].expect("a fitted stratum has its first group's numbers");
    let errors = fit.standard_errors.as_ref();
    let none = StratumError::NotIdentified;
    let group = errors.and_then(|errors| errors.slippage[0]);
    let (level, shorter, fall_off) = group.map_or((none, none, none), |group| {
        (group.level, group.shorter_share, group.fall_off)
    });
    let concentration = errors.map_or(none, |errors| errors.concentration);
    let rows = [
        ("slippage level", with_error(slippage.level, level)),
        ("shorter share", with_error(slippage.shorter_share, shorter)),
        ("fall-off", with_error(slippage.fall_off, fall_off)),
        (
            "concentration",
            with_error(fit.concentration, concentration),
        ),
    ];
    rows.into_iter()
        .map(|(name, (value, error))| (name, value, error))
        .collect()
}

/// An error and every cause beneath it, one per line.
fn with_its_causes(error: &dyn std::error::Error) -> String {
    let mut text = error.to_string();
    let mut source = error.source();
    while let Some(cause) = source {
        text.push_str(&format!("\n  because: {cause}"));
        source = cause.source();
    }
    text
}

fn format_error(error: Option<f64>) -> String {
    error.map_or("none".to_owned(), |error| format!("{error:.5}"))
}

/// **The strata compared when `STRATA` is not set: five spread across the cohort's sizes** — the
/// one holding the most tracts with reads, and those an eighth, a quarter, a half and three
/// quarters of the way down the ranking of those with at least the refusal floor's tracts with reads
/// ([`DEFAULT_REFUSAL_FLOOR`]): a thinner stratum is refused by both fits and shows nothing. The
/// largest strata cost the most to fit on every sample, and the ones that grow their subset past
/// the first are the ones `FIRST_SUBSET` is chosen from, so a spread shows both.
fn spread_across_sizes(strata: &[StratumEvidence]) -> Vec<&StratumEvidence> {
    let mut by_size: Vec<&StratumEvidence> = strata
        .iter()
        .filter(|evidence| evidence.tracts_with_reads() >= DEFAULT_REFUSAL_FLOOR)
        .collect();
    by_size.sort_by_key(|evidence| std::cmp::Reverse(evidence.tracts_with_reads()));
    let count = by_size.len();
    let mut ranks: Vec<usize> = [0, count / 8, count / 4, count / 2, 3 * count / 4]
        .into_iter()
        .filter(|rank| *rank < count)
        .collect();
    ranks.dedup();
    ranks.into_iter().map(|rank| by_size[rank]).collect()
}

/// The likelihood table one fit of `evidence` holds, as text: a row a (tract, sample with reads), of
/// 91 genotype pairs of 8 bytes at thirteen allele classes.
fn table_size(evidence: &StratumEvidence, classes: usize) -> String {
    let rows: usize = evidence
        .tracts
        .iter()
        .map(|tract| tract.samples.len())
        .sum();
    let genotypes = classes * (classes + 1) / 2;
    let bytes = (rows * genotypes * 8) as f64;
    let gib = bytes / f64::from(1_u32 << 30);
    if gib >= 1.0 {
        format!("{gib:.2} GiB")
    } else {
        format!("{:.1} MiB", bytes / f64::from(1_u32 << 20))
    }
}

fn main() {
    let Cli { args } = Cli::parse();
    if let Some(coefficient) = args.inbreeding {
        InbreedingF::try_new(coefficient).unwrap_or_else(|_| {
            panic!("--inbreeding {coefficient} is not a coefficient in 0 to 1")
        });
    }
    let OpenedCensuses {
        mut evidence, plan, ..
    } = open_the_censuses(&args).unwrap_or_else(|error| panic!("{}", with_its_causes(&error)));
    let names: Vec<String> = evidence.sample_names().map(str::to_owned).collect();
    let order = sample_order(&names, SAMPLE_ORDER_SEED);
    let excess = vec![args.inbreeding.unwrap_or(0.0); names.len()];
    let contig_of = contig_id_of(&plan);
    let pooled = every_read_group_pooled(&evidence);
    let CohortTractStrata { strata, .. } =
        the_tract_strata_of_a_cohort(&mut evidence, &plan.loci, &contig_of, &pooled)
            .unwrap_or_else(|error| panic!("{}", with_its_causes(&error)));

    let chosen: Vec<&StratumEvidence> = match strata_asked_for() {
        Some(asked) => asked
            .iter()
            .map(|&(period, repeats)| {
                strata
                    .iter()
                    .find(|evidence| {
                        evidence.stratum.period == period
                            && evidence.stratum.reference_repeats == repeats
                    })
                    .unwrap_or_else(|| panic!("no stratum {period}:{repeats} in this cohort"))
            })
            .collect(),
        None => spread_across_sizes(&strata),
    };

    let mut config = SsrFitConfig::default();
    config.curve.draw_curves = false;
    let setting = |name: &str| {
        std::env::var(name).ok().map(|value| {
            value
                .trim()
                .parse::<f64>()
                .unwrap_or_else(|_| panic!("{name}={value:?} is not a number"))
        })
    };
    if let Some(first) = setting("FIRST_SUBSET") {
        config.subsets.first = first as usize;
    }
    if let Some(least) = setting("MIN_SAMPLES_A_GROUP") {
        config.subsets.min_samples_a_group = least as usize;
    }
    if let Some(target) = setting("TARGET") {
        config.subsets.level_relative_error_target = target;
    }
    let classes = (2 * config.allele_span + 1) as usize;
    println!(
        "{} samples; homozygote excess {} for every sample; {} strata compared; first subset {}, \
         {} readers a slippage group at least, level target {}",
        names.len(),
        excess.first().copied().unwrap_or(0.0),
        chosen.len(),
        config.subsets.first,
        config.subsets.min_samples_a_group,
        config.subsets.level_relative_error_target,
    );
    let (mut subset_total, mut whole_total) = (0.0, 0.0);
    for evidence in chosen {
        let label = format!(
            "stratum {}:{} ({} tracts with reads)",
            evidence.stratum.period,
            evidence.stratum.reference_repeats,
            evidence.tracts_with_reads()
        );
        // **The subset first, printed at once**: it is the cheap fit, and on a large cohort the
        // fit on every sample can take hours.
        println!("{label}: fitting the grown subset");
        let started = Instant::now();
        let subset =
            fit_strata_on_sample_subsets(std::slice::from_ref(evidence), &excess, &order, &config);
        let subset_seconds = started.elapsed().as_secs_f64();
        let StratumOutcome::Fitted(subset) = &subset[0] else {
            println!(
                "{label}: the subset fit refused the stratum ({:?})",
                subset[0]
            );
            continue;
        };
        let reached = match subset.samples_fitted_on {
            None => "every sample, no subset drawn".to_owned(),
            Some(samples) if samples >= names.len() => format!("every sample ({samples})"),
            Some(samples) => format!("{samples} samples"),
        };
        println!("{label}: the subset reached {reached}, in {subset_seconds:.1} s");
        for (name, value, error) in numbers_of(subset) {
            println!("  {name}: subset {value:.5} ± {}", format_error(error));
        }
        println!(
            "{label}: fitting every sample (a likelihood table of about {})",
            table_size(evidence, classes)
        );
        let started = Instant::now();
        let Some(whole) = fit_stratum(evidence, &excess, &config) else {
            println!("{label}: the fit on every sample found no read to walk");
            continue;
        };
        let whole_seconds = started.elapsed().as_secs_f64();
        subset_total += subset_seconds;
        whole_total += whole_seconds;
        let every_sample_read = subset
            .samples_fitted_on
            .is_none_or(|samples| samples >= names.len());
        println!(
            "{label}: {subset_seconds:.1} s on the subset against {whole_seconds:.1} s on every \
             sample"
        );
        for ((name, whole_value, whole_error), (_, subset_value, subset_error)) in
            numbers_of(&whole).into_iter().zip(numbers_of(subset))
        {
            let apart = match (whole_error, subset_error) {
                // Both fits read every sample: the difference is the climb's, in the whole fit's
                // errors.
                (Some(whole_error), _) if every_sample_read => format!(
                    "{:+.2} of its error (both read every sample)",
                    (subset_value - whole_value) / whole_error
                ),
                (Some(whole_error), Some(subset_error)) if subset_error > whole_error => {
                    let spread = (subset_error * subset_error - whole_error * whole_error).sqrt();
                    format!("{:+.2}", (subset_value - whole_value) / spread)
                }
                _ => "-".to_owned(),
            };
            println!(
                "  {name}: every sample {whole_value:.5} ± {}, subset {subset_value:.5} ± {}, \
                 apart {apart}",
                format_error(whole_error),
                format_error(subset_error),
            );
        }
    }
    if whole_total > 0.0 {
        println!(
            "in all: {subset_total:.1} s on the subsets against {whole_total:.1} s on every \
             sample, {:.0}% of the time",
            100.0 * subset_total / whole_total
        );
    }
}
