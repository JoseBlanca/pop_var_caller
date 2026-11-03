//! **Each read group's stutter multiplier as the repeat-tract fit finds it**, for measuring the
//! fit against a cohort whose truth is known (`str_library_stutter.md` §8).
//!
//! ```text
//! cargo run --release --example ng_str_read_group_multipliers -- \
//!     --reference genome.fa --catalog catalog.parquet --psp psps/ --output multipliers.tsv \
//!     --fitted-parameters fit.toml
//! ```
//!
//! The tract evidence is gathered as `estimate-parameters` gathers it, and fitted as it fits it,
//! both rules on, with each sample's homozygote excess read back from `--fitted-parameters` (a
//! fit's inbreeding coefficients, which are that excess). Two tables:
//!
//! - `--output`: one line a motif length and read group: the sample whose read group it is, the
//!   multiplier, its standard error on the log scale (empty when it rests on the pull), and the
//!   read group's expected slipped reads;
//! - `<output>.classes.tsv`: one line a class: its slip rate, the rate of a read group of typical
//!   stutter.
//!
//! The options are `estimate-parameters`'s; only `--reference`, `--catalog` and `--psp` are read
//! from them, and `--output` names the first table.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use pop_var_caller::calling::parameters_file::ParametersFile;
use pop_var_caller::cli::estimate_parameters::{
    EstimateParametersArgs, contig_id_of, open_the_censuses,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::penalised_fit::{
    PenalisedFitConfig, fit_repeat_tracts,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::tract_exclusion::TractExclusion;
use pop_var_caller::run::census_fit::{every_read_group_pooled, the_tract_strata_of_a_cohort};

#[derive(Parser)]
struct Cli {
    #[command(flatten)]
    args: EstimateParametersArgs,
    /// A parameters file fitted from these psps, for each sample's inbreeding coefficient.
    #[arg(long)]
    fitted_parameters: PathBuf,
    /// Leave out, before the fit, every class with fewer loci with reads than this: for testing
    /// whether thin classes steer their motif length's curves. One keeps every class.
    #[arg(long, default_value_t = 1)]
    min_loci_a_class: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    let fitted = ParametersFile::from_toml(&std::fs::read_to_string(&cli.fitted_parameters)?)?;
    let coefficient: BTreeMap<&str, f64> = fitted
        .inbreeding
        .by_sample
        .iter()
        .map(|row| (row.sample.as_str(), row.inbreeding_coefficient.value))
        .collect();

    let opened = open_the_censuses(&cli.args)?;
    let mut evidence = opened.evidence;
    let samples: Vec<String> = evidence.sample_names().map(str::to_owned).collect();
    let excess: Vec<f64> = samples
        .iter()
        .map(|name| coefficient.get(name.as_str()).copied().unwrap_or(0.0))
        .collect();
    let contig_of = contig_id_of(&opened.plan);
    let pooled_groups = every_read_group_pooled(&evidence);
    let gathered =
        the_tract_strata_of_a_cohort(&mut evidence, &opened.plan.loci, &contig_of, &pooled_groups)?;
    let mut strata = gathered.strata;
    strata.retain(|class| class.tracts_with_reads() >= cli.min_loci_a_class);
    // Which sample each read group's reads came from.
    let mut sample_of: BTreeMap<usize, u32> = BTreeMap::new();
    for row in strata
        .iter()
        .flat_map(|class| &class.tracts)
        .flat_map(|tract| &tract.samples)
    {
        for (read_group, _) in &row.by_read_group {
            sample_of.insert(read_group.get(), row.sample);
        }
    }
    let fit = fit_repeat_tracts(
        &mut strata,
        &excess,
        &PenalisedFitConfig::production(),
        &TractExclusion::production(),
    );

    let mut out = BufWriter::new(File::create(&cli.args.output)?);
    writeln!(
        out,
        "period\tread_group\tsample\tmultiplier\tlog_error\tslipped_reads"
    )?;
    for (period, fitted) in &fit.read_group_multipliers {
        for (read_group, multiplier) in fitted.multipliers.iter().enumerate() {
            let sample = sample_of
                .get(&read_group)
                .map_or("", |sample| samples[*sample as usize].as_str());
            let error = fitted.log_errors[read_group].map_or(String::new(), |e| format!("{e}"));
            writeln!(
                out,
                "{period}\t{read_group}\t{sample}\t{multiplier}\t{error}\t{}",
                fitted.slipped_reads[read_group]
            )?;
        }
    }
    let mut classes = BufWriter::new(File::create(format!(
        "{}.classes.tsv",
        cli.args.output.display()
    ))?);
    writeln!(
        classes,
        "period\treference_repeats\tloci_with_reads\tslip_rate"
    )?;
    for (evidence, outcome) in strata.iter().zip(&fit.outcomes) {
        if let Some(Some(slippage)) = outcome.slippage().first() {
            writeln!(
                classes,
                "{}\t{}\t{}\t{}",
                evidence.stratum.period,
                evidence.stratum.reference_repeats,
                evidence.tracts_with_reads(),
                slippage.level
            )?;
        }
    }
    Ok(())
}
