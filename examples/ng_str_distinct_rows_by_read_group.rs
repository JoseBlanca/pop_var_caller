//! **How many distinct rows of reads does each repeat class hold when a row is keyed by its read
//! group as well as by its reads?** — the cost of fitting stutter per read group
//! (`str_library_stutter.md`) against the distinct-row speed-up of the penalised fit
//! (`str_penalised_fit.md` §2.8).
//!
//! ```text
//! cargo run --release --example ng_str_distinct_rows_by_read_group -- \
//!     --reference genome.fa --catalog catalog.parquet --psp psps/ --output rows.tsv
//! ```
//!
//! The evidence is gathered as `estimate-parameters` gathers it, the depth rule applied. **One line
//! a class**: its rows (one a (tract, sample) with reads), its distinct rows with a sample's read
//! groups pooled, as the fit holds them today, and its distinct rows keyed by the sample's read
//! groups and each read group's own reads, as a fit scoring each read group at its own stutter must.

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Write};

use clap::Parser;
use pop_var_caller::cli::estimate_parameters::{
    EstimateParametersArgs, contig_id_of, open_the_censuses,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::tract_exclusion::{
    DepthRule, leave_out_deep_tracts,
};
use pop_var_caller::run::census_fit::{every_read_group_pooled, the_tract_strata_of_a_cohort};

#[derive(Parser)]
struct Cli {
    #[command(flatten)]
    args: EstimateParametersArgs,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    let opened = open_the_censuses(&cli.args)?;
    let mut evidence = opened.evidence;
    let contig_of = contig_id_of(&opened.plan);
    let pooled_groups = every_read_group_pooled(&evidence);
    let gathered =
        the_tract_strata_of_a_cohort(&mut evidence, &opened.plan.loci, &contig_of, &pooled_groups)?;
    let mut strata = gathered.strata;
    leave_out_deep_tracts(&mut strata, &DepthRule::production());

    let mut out = BufWriter::new(File::create(&cli.args.output)?);
    writeln!(
        out,
        "period\treference_repeats\tread_groups\trows\tdistinct_pooled\tdistinct_by_read_group"
    )?;
    for stratum in &strata {
        let mut pooled: HashSet<Vec<u32>> = HashSet::new();
        let mut by_read_group: HashSet<Vec<(u32, Vec<u32>)>> = HashSet::new();
        let mut rows = 0_usize;
        for sample in stratum.tracts.iter().flat_map(|tract| &tract.samples) {
            rows += 1;
            let mut summed = vec![0_u32; (2 * stratum.read_span + 1) as usize];
            for (_, reads) in &sample.by_library {
                for (into, count) in summed.iter_mut().zip(reads) {
                    *into += count;
                }
            }
            pooled.insert(summed);
            by_read_group.insert(
                sample
                    .by_library
                    .iter()
                    .map(|(library, reads)| (library.0, reads.clone()))
                    .collect(),
            );
        }
        if rows == 0 {
            continue;
        }
        writeln!(
            out,
            "{}\t{}\t{}\t{rows}\t{}\t{}",
            stratum.stratum.period,
            stratum.stratum.reference_repeats,
            stratum.libraries.len(),
            pooled.len(),
            by_read_group.len()
        )?;
    }
    Ok(())
}
