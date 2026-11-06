//! **What do the repeat tracts the stutter fit leaves out look like, and do ordinary tracts'
//! off-length reads spread from tract to tract as a binomial does?** — measurements for
//! `parameter_prepass_ssr.md` §4.6.1 on a cohort.
//!
//! ```text
//! cargo run --release --example ng_str_tract_spread -- \
//!     --reference genome.fa --catalog catalog.parquet --psp psps/ --output tracts.tsv \
//!     --fitted-parameters fit.toml
//! ```
//!
//! The tract evidence is gathered exactly as `estimate-parameters` gathers it (the guard applied).
//! Then the two rules run as the fit runs them — the depth rule over every class, then the
//! penalised fit with the stutter-tail rule — with each sample's homozygote excess read back from
//! `--fitted-parameters` (a fit's inbreeding coefficients, which are that excess). **One line a
//! tract**, saying which rule left it out (`none`, `depth`, `tail`) and what its reads look like:
//!
//! - `samples`, `reads`: samples with reads, reads summed over them;
//! - `sample_lengths`: how many distinct lengths the samples' own most common lengths take, over
//!   the samples with a single most common length. Above one, the samples do not agree on the
//!   tract's length: in an inbred cohort that is a tract polymorphic between accessions;
//! - `off_shared`: reads at a length other than the cohort's most common one;
//! - `agree`: 1 when every sample has a single most common length and they are all the same, so
//!   in an inbred cohort every off-length read is stutter or an artefact — the tracts the binomial
//!   spread is measured on.
//!
//! **What the selection of `agree` costs.** A sample whose reads are mostly off-length has another
//! most common length and drops its tract, so the tracts that stutter most are under-represented.
//!
//! The options are `estimate-parameters`'s; only `--reference`, `--catalog` and `--psp` are read
//! from them, and `--output` names the table written.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use pop_var_caller::calling::parameters_file::ParametersFile;
use pop_var_caller::cli::estimate_parameters::{
    EstimateParametersArgs, contig_id_of, open_the_censuses,
};
use pop_var_caller::parameter_estimation::joint::census::RECORDED_OFFSET_RANGE;
use pop_var_caller::parameter_estimation::joint::ssr_fit::TractReads;
use pop_var_caller::parameter_estimation::joint::ssr_fit::penalised_fit::{
    PenalisedFitConfig, fit_repeat_tracts,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::tract_exclusion::{
    DepthRule, StutterTailRule, TractExclusion, leave_out_deep_tracts,
};
use pop_var_caller::run::census_fit::{every_read_group_pooled, the_tract_strata_of_a_cohort};

#[derive(Parser)]
struct Cli {
    #[command(flatten)]
    args: EstimateParametersArgs,
    /// A parameters file fitted from these psps, for each sample's inbreeding coefficient.
    #[arg(long)]
    fitted_parameters: PathBuf,
}

/// One sample's reads at a tract, its read groups pooled.
fn pooled(
    sample: &pop_var_caller::parameter_estimation::joint::ssr_fit::SampleTractReads,
) -> Vec<u64> {
    let mut row = vec![0_u64; (2 * RECORDED_OFFSET_RANGE + 1) as usize];
    for (_, reads) in &sample.by_read_group {
        for (bucket, count) in reads.iter().enumerate() {
            row[bucket] += u64::from(*count);
        }
    }
    row
}

/// The tract's line, after its period, repeats and which rule left it out.
fn describe(tract: &TractReads) -> String {
    let mut summed = vec![0_u64; (2 * RECORDED_OFFSET_RANGE + 1) as usize];
    let mut modes = std::collections::BTreeSet::new();
    let mut every_sample_single = true;
    for sample in &tract.samples {
        let row = pooled(sample);
        let most = *row.iter().max().expect("buckets");
        let at_most: Vec<usize> = (0..row.len()).filter(|b| row[*b] == most).collect();
        if at_most.len() == 1 {
            modes.insert(at_most[0]);
        } else {
            every_sample_single = false;
        }
        for (into, from) in summed.iter_mut().zip(&row) {
            *into += from;
        }
    }
    let reads: u64 = summed.iter().sum();
    let shared = *summed.iter().max().unwrap_or(&0);
    let agree = every_sample_single && modes.len() == 1 && tract.samples.len() >= 2;
    format!(
        "{}\t{}\t{}\t{}\t{}",
        tract.samples.len(),
        reads,
        modes.len(),
        reads - shared,
        u8::from(agree)
    )
}

/// Which of `before` are missing from `after`: `after` is `before` with some removed, in order.
fn removed(before: &[TractReads], after: &[TractReads]) -> Vec<bool> {
    let mut kept = after.iter().peekable();
    before
        .iter()
        .map(|tract| {
            if kept.peek() == Some(&tract) {
                kept.next();
                false
            } else {
                true
            }
        })
        .collect()
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
    let excess: Vec<f64> = evidence
        .sample_names()
        .map(|name| coefficient.get(name).copied().unwrap_or(0.0))
        .collect();
    let contig_of = contig_id_of(&opened.plan);
    let pooled_groups = every_read_group_pooled(&evidence);
    let gathered =
        the_tract_strata_of_a_cohort(&mut evidence, &opened.plan.loci, &contig_of, &pooled_groups)?;
    let original = gathered.strata.clone();

    let mut strata = gathered.strata;
    leave_out_deep_tracts(&mut strata, &DepthRule::production());
    let after_depth = strata.clone();
    let _ = fit_repeat_tracts(
        &mut strata,
        &excess,
        &PenalisedFitConfig::production(),
        &TractExclusion {
            depth: None,
            stutter_tail: Some(StutterTailRule::production()),
        },
    );

    let mut out = BufWriter::new(File::create(&cli.args.output)?);
    writeln!(
        out,
        "period\treference_repeats\tleft_out_by\tsamples\treads\tsample_lengths\toff_shared\tagree"
    )?;
    for ((first, deep), last) in original.iter().zip(&after_depth).zip(&strata) {
        let by_depth = removed(&first.tracts, &deep.tracts);
        let by_tail = removed(&deep.tracts, &last.tracts);
        let mut tail = by_tail.into_iter();
        for (tract, deep) in first.tracts.iter().zip(by_depth) {
            let rule = if deep {
                "depth"
            } else if tail
                .next()
                .expect("one verdict a tract the depth rule kept")
            {
                "tail"
            } else {
                "none"
            };
            writeln!(
                out,
                "{}\t{}\t{rule}\t{}",
                first.stratum.period,
                first.stratum.reference_repeats,
                describe(tract)
            )?;
        }
    }
    Ok(())
}
