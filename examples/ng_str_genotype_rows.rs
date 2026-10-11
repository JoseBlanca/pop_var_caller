//! **Writes each repeat class's genotype likelihoods, a row a (tract, sample), with the stutter held
//! at a class fit's values** — for estimating the class's allele-length frequencies, and fitting
//! shapes to them, outside the fit (the owner's "route 2" research, 2026-10-11).
//!
//! The stutter of each class is read from a table `ng_str_class_spectra` wrote (`--classes`). One
//! file a class, `<period>_<repeats>.bin`, little-endian: the row count and the row width as `u64`,
//! then each row's tract (`u32`), each row's homozygote excess (`f64`), then the rows (`f64`), each
//! divided by its largest value. `genotypes.tsv` gives the order of a row's genotype pairs, as
//! allele offsets from the reference length.
//!
//! ```text
//! cargo run --release --example ng_str_genotype_rows -- --reference REF --catalog CATALOG \
//!   --psp PSPS --output unused --classes classes.tsv --homozygote-excess 0.98 --out-dir DIR
//! ```

use std::collections::BTreeMap;
use std::io::Write;
use std::path::PathBuf;

use clap::Parser;

use pop_var_caller::cli::estimate_parameters::{
    EstimateParametersArgs, OpenedCensuses, contig_id_of, open_the_censuses,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::genotype_rows::genotype_rows;
use pop_var_caller::parameter_estimation::joint::ssr_fit::{Slippage, SsrFitConfig};
use pop_var_caller::run::{
    CohortTractStrata, every_read_group_pooled, the_tract_strata_of_a_cohort,
};

#[derive(Parser)]
#[command(about = "Write each repeat class's genotype likelihoods with its stutter held")]
struct Cli {
    #[command(flatten)]
    args: EstimateParametersArgs,
    /// The class table `ng_str_class_spectra` wrote, whose stutter is held.
    #[arg(long)]
    classes: PathBuf,
    /// One homozygote excess for every sample.
    #[arg(long, default_value_t = 0.0)]
    homozygote_excess: f64,
    /// How far either side of the reference length alleles may sit (the fit's `ALLELE_SPAN` by
    /// default).
    #[arg(long)]
    allele_span: Option<i32>,
    /// A parameters file whose per-sample homozygote excesses are used instead.
    #[arg(long)]
    fitted: Option<PathBuf>,
    /// Where the class files go.
    #[arg(long)]
    out_dir: PathBuf,
}

fn main() {
    let cli = Cli::parse();
    let OpenedCensuses {
        mut evidence, plan, ..
    } = open_the_censuses(&cli.args).unwrap_or_else(|error| panic!("{error:?}"));
    let names: Vec<String> = evidence.sample_names().map(str::to_owned).collect();
    let samples = names.len();
    let contig_of = contig_id_of(&plan);
    let pooled = every_read_group_pooled(&evidence);
    let CohortTractStrata { strata, .. } =
        the_tract_strata_of_a_cohort(&mut evidence, &plan.loci, &contig_of, &pooled)
            .unwrap_or_else(|error| panic!("{error:?}"));

    let text = std::fs::read_to_string(&cli.classes).expect("the class table reads");
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().expect("a header").split('\t').collect();
    let column = |name: &str| {
        header
            .iter()
            .position(|column| *column == name)
            .unwrap_or_else(|| panic!("no column {name}"))
    };
    let (period, repeats, level, shorter, fall_off) = (
        column("period"),
        column("repeats"),
        column("level"),
        column("shorter_share"),
        column("fall_off"),
    );
    let held: BTreeMap<(u8, u64), Slippage> = lines
        .map(|line| {
            let values: Vec<&str> = line.split('\t').collect();
            (
                (
                    values[period].parse().expect("a period"),
                    values[repeats].parse().expect("repeats"),
                ),
                Slippage {
                    level: values[level].parse().expect("a level"),
                    shorter_share: values[shorter].parse().expect("a share"),
                    fall_off: values[fall_off].parse().expect("a fall-off"),
                },
            )
        })
        .collect();

    std::fs::create_dir_all(&cli.out_dir).expect("the output directory");
    let mut config = SsrFitConfig::default();
    if let Some(span) = cli.allele_span {
        config.allele_span = span;
    }
    let excess: Vec<f64> = match &cli.fitted {
        Some(path) => {
            let file = pop_var_caller::calling::parameters_file::ParametersFile::from_toml(
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
                .map(|name| *by_name.get(name.as_str()).unwrap_or(&cli.homozygote_excess))
                .collect()
        }
        None => vec![cli.homozygote_excess; samples],
    };
    let mut wrote_genotypes = false;
    for evidence in &strata {
        let key = (evidence.stratum.period, evidence.stratum.reference_repeats);
        let Some(slippage) = held.get(&key) else {
            continue;
        };
        let rows = genotype_rows(evidence, &[*slippage], &excess, &config);
        if !wrote_genotypes {
            let span = (rows.classes / 2) as i64;
            let order: String = rows
                .genotypes
                .iter()
                .map(|(first, second)| {
                    format!("{}\t{}\n", *first as i64 - span, *second as i64 - span)
                })
                .collect();
            std::fs::write(
                cli.out_dir.join("genotypes.tsv"),
                format!("first\tsecond\n{order}"),
            )
            .expect("the genotype order writes");
            wrote_genotypes = true;
        }
        let path = cli.out_dir.join(format!("{}_{}.bin", key.0, key.1));
        let mut file = std::io::BufWriter::new(std::fs::File::create(&path).expect("a class file"));
        let count = rows.tract.len() as u64;
        file.write_all(&count.to_le_bytes()).expect("writes");
        file.write_all(&(rows.genotypes.len() as u64).to_le_bytes())
            .expect("writes");
        for tract in &rows.tract {
            file.write_all(&tract.to_le_bytes()).expect("writes");
        }
        for excess in &rows.homozygote_excess {
            file.write_all(&excess.to_le_bytes()).expect("writes");
        }
        for value in &rows.values {
            file.write_all(&value.to_le_bytes()).expect("writes");
        }
        eprintln!("{}:{} — {count} rows", key.0, key.1);
    }
}
