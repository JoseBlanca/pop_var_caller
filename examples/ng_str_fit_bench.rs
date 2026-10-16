//! **The repeat-tract fitting bench**: draw one motif length's classes from a known model, fit them
//! with each method asked for, and write how close each method came to the truth, how well it
//! predicts loci it was not fitted on, and what it cost (`doc/devel/ng/research/str_fit_research_plan.md`
//! §2, §3.1).
//!
//! One run is one **cell** of the plan's grid — a motif length, a number of samples, a read depth,
//! a truth and a seed. The drawing is the count-level simulator
//! (`ssr_fit::bench_fixtures::cohort_draw`); the methods and the held-out score are in
//! `ssr_fit::bench_fixtures::fit_bench`, whose documentation says how to add a method.
//!
//! ```text
//! cargo run --release --features bench-fixtures --example ng_str_fit_bench -- \
//!   --period 1 --samples 60 --depth 3 --seed 1 --truth pcr --methods B,D_tract_prior \
//!   --out rows.tsv --truth-out truth.tsv
//! ```
//!
//! # What it writes
//!
//! `--out`: one row per method, class and number, long format, tab-separated:
//!
//! | column | meaning |
//! |---|---|
//! | `cell`, `period`, `samples`, `depth`, `truth_variant`, `seed` | the cell |
//! | `method` | the method's name (`B`, `D_tract_prior`, …), with a fallback's suffix on held-out rows |
//! | `repeats` | the class's reference repeat count, or `all` |
//! | `parameter` | below |
//! | `truth`, `estimate`, `standard_error` | `NA` where there is none |
//!
//! Parameters: `level`, `shorter_share`, `fall_off`, `concentration`, `between_locus_share`
//! (`1 / (1 + concentration)`), `share_<offset>` for every allele class; the held-out score
//! `heldout_ln_per_read` and `heldout_ln_per_locus` (truth = the truth's own score on the same loci,
//! standard error = the Monte Carlo error) with `heldout_loci`, `heldout_reads` and
//! `heldout_loci_on_fallback` (loci of classes scored at calling's fallback rung, a flat shape at the
//! stated concentration, for want of their own shares or concentration); and the costs
//! `fit_seconds` (the fit on every locus), `heldout_fit_seconds`, `heldout_score_seconds`,
//! `draw_seconds` and `peak_rss_mib` (Linux only; elsewhere `NA`).
//!
//! `--fixture-out` writes none of that. It fits one penalised method on every locus and writes one
//! line of JSON holding everything the fit returned at full precision — each class's numbers and
//! standard errors, the curves and their strengths, and the round counts. That line is the unit of
//! the reference fixture which `scripts/research/str_fit_fixture.py` records and compares
//! against.
//!
//! **Peak memory is the process's**, so a run with several methods reports the largest of them. To
//! measure each method, run one method a process, as `scripts/research/str_fit_grid.sh` does; it
//! also records the process's peak through `/usr/bin/time`, which works on macOS too.
//!
//! # The psp mode: a real cohort, held out by whole locus
//!
//! With `--psp`, the bench scores the methods on a real cohort's reads instead of a drawn one
//! (spec `str_penalised_fit.md` §4.4; plan step D2a). There is no truth, so it computes only the
//! held-out score: how well each method's fit, made without a locus, predicts every sample's reads
//! at that locus.
//!
//! ```text
//! cargo run --release --features bench-fixtures --example ng_str_fit_bench -- \
//!   --psp PSPS --reference REF.fa --catalog CATALOG.parquet [--parameters FITTED.toml] \
//!   --methods B,D_tract_prior --out rows.tsv --locus-out loci.tsv
//! ```
//!
//! - **The evidence is `estimate-parameters`'s own**: the cohort is opened and checked by that
//!   command's `open_the_censuses` (samples in name order, read groups in declared-ID order), the
//!   selection checked by `refuse_another_selection`, and the classes gathered by
//!   `the_tract_strata_of_a_cohort` with every read group in one slippage group — so the
//!   repeat-tract guard and everything else that command applies apply here.
//! - **Each sample's homozygote excess** is the SNP/indel fit's, as `fit_a_cohort` hands it on
//!   (`homozygote_excess_in_cohort_order`); that fit takes about 11 minutes on four tomato psps and
//!   an hour on 59. `--parameters` reads the same numbers from a parameters file fitted on the same
//!   psps instead (its inbreeding coefficients, joined by sample name).
//! - **Each motif length is fitted on its own**, as `estimate-parameters` fits them; `--period`
//!   keeps one.
//! - **Folds are whole loci** within each class (`fold_plans_across_motif_lengths`): every locus is
//!   held out exactly once, `--seed` fixes the split, and at most `--heldout-loci` held-out loci a
//!   class are scored. Every method sees the same folds and the same Monte Carlo draws.
//! - **The rows** are the drawn mode's held-out rows, one set a motif length: `period` is the
//!   motif length, `depth` the mean spanning reads of a (locus, sample) pair with reads,
//!   `truth_variant` is `real`, `seed` the fold seed, and `truth` is `NA`. Each motif length adds
//!   `heldout_loci_to_score`, the loci picked for scoring; a locus whose class got no numbers from
//!   the other folds is scored by no method, and is the gap between that count and a variant's
//!   `heldout_loci`. Two methods' scores a locus are comparable only over the loci both scored:
//!   method `B` leaves a class with fewer than 8 loci with reads unfitted, and those loci
//!   unscored. `--locus-out` writes every scored locus, with its motif length in the last column,
//!   so the comparison can be paired by locus. The process's `peak_rss_mib` covers every motif
//!   length and is written with `period` `all`.
//!
//! # The proxy for calling
//!
//! `--calls-out` asks whether a fit's small errors change genotype calls (spec
//! `str_penalised_fit.md` §5.1, plan step D6), on cohorts too large to simulate reads for. Each
//! (locus, sample) with reads is called: the genotype its weights put highest, the weights being
//! the fit's own expectation step (`fit_bench::ClassCaller`) — a genotype's prior from the class's
//! allele-length shares, its spread and the sample's homozygote excess, times how well it explains
//! the sample's reads. The **reference** calls are made under the truth's numbers, and each of
//! these sets of numbers is compared with them:
//!
//! | `arm` | stutter | shares and spread |
//! |---|---|---|
//! | `level_0.995` | the truth's, slip rate 0.5% of its value low | the truth's |
//! | `fitted_level` | the truth's, slip rate the method's | the truth's |
//! | `fitted_stutter` | the method's three numbers | the truth's |
//! | `fitted_all` | the method's | the method's, at the held-out score's fallback where absent |
//!
//! The method is fitted once on every locus, as the rows' own estimates are. Columns: the cell's
//! six, then `method`, `arm`, `repeats` (or `all`), `true_level` and `arm_level` (`NA` on `all`
//! rows), `pairs` ((locus, sample) pairs with reads), `wrong_under_truth` (reference calls that
//! are not the genotype drawn, the error already there), `changed` (calls that differ from the
//! reference), and of those `right_to_wrong`, `wrong_to_right` and `wrong_to_wrong`. A genotype is
//! read over the fit's window, each copy past ±9 counted at ±9. A class a method gave no stutter
//! has no row for the arms that need it. **So an arm's `all` row sums only the classes that arm has
//! rows for**, and two arms' `all` rows can cover different classes: compare their `pairs` before
//! comparing their `changed`.

// The allocator `pvc` ships (`src/main.rs`), so that the bench's times and peak memory describe
// what `estimate-parameters` runs. A `#[global_allocator]` is per binary: without this line the
// bench ran the system allocator. Gated as `src/main.rs` gates it; this bench declares no dhat
// allocator, so the `dhat-heap` feature does not turn it off.
#[cfg(feature = "alloc-mimalloc")]
#[global_allocator]
static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

use clap::error::ErrorKind;
use clap::{CommandFactory, Parser};
use rayon::prelude::*;

use pop_var_caller::calling::parameters_file::{ParametersFile, Warrant};
use pop_var_caller::cli::estimate_parameters::{
    EstimateParametersArgs, OpenedCensuses, contig_id_of, open_the_censuses,
};
use pop_var_caller::parameter_estimation::joint::census::RECORDED_OFFSET_RANGE;
use pop_var_caller::parameter_estimation::joint::fit::{JointFitConfig, fit_jointly};
use pop_var_caller::parameter_estimation::joint::ssr_fit::bench_fixtures::cohort_draw::{
    ClassTruth, CohortModel, DrawnPeriod, HomozygoteExcessTruth, LociPerClass, ParameterCurve,
    StutterRegime, draw_period, write_truth_tsv,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::bench_fixtures::fit_bench::{
    CallChanges, ClassCaller, ClassEstimate, ClassNumbers, CountingFitErrorMethod, FoldPlan,
    HeldOutClasses, HeldOutRun, HeldOutSettings, HeldOutSums, PeriodFit, Sampler, ScoringFallback,
    TheCountingFit, ThePenalisedFit, call_changes, class_seed_across_motif_lengths,
    drawn_genotypes_in_window, fold_plans_across_motif_lengths, held_out_run, numbers_for_scoring,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::penalised_fit::{
    ParameterCurve as FittedCurve, PenalisedPeriodFit, stutter_prior_weight_problem,
};
use pop_var_caller::parameter_estimation::joint::ssr_fit::{
    ALLELE_SPAN, Slippage, StratumEvidence,
};
use pop_var_caller::run::{
    CohortTractStrata, every_read_group_pooled, homozygote_excess_in_cohort_order,
    refuse_another_selection, the_tract_strata_of_a_cohort,
};
use serde_json::{Value, json};

#[derive(Parser)]
#[command(
    about = "Score repeat-tract fitting methods on a cohort drawn from a known model, or \
                   on a real cohort's psps"
)]
struct Cli {
    /// The motif length; 1 by default. **With `--psp`, only this motif length is fitted**, and
    /// every one the cohort has when it is not given.
    #[arg(long)]
    period: Option<u8>,
    /// **The psp mode: a real cohort in place of a drawn one.** One psp per sample, or a directory
    /// holding them; repeat the flag. The censuses are read as `estimate-parameters` reads them,
    /// and only the held-out score is computed (see the module documentation). **The flags that
    /// shape a drawn cohort are refused beside it** (`--excess`, `--samples`, `--depth`, …), since
    /// a real cohort ignores them.
    #[arg(long = "psp", num_args = 1.., requires = "reference",
          conflicts_with_all = [
              "fixture_out", "truth_out", "excess", "samples", "depth", "truth", "loci",
              "loci_at_most", "repeats", "between_locus_share", "library_multipliers",
              "multipliers_told", "allele_reach",
          ])]
    psps: Vec<PathBuf>,
    /// The psp mode's reference FASTA, the one the psps were walked against.
    #[arg(long, requires = "psps")]
    reference: Option<PathBuf>,
    /// The psp mode's repeat catalog; `<reference>.repeats.parquet` by default.
    #[arg(long, requires = "psps")]
    catalog: Option<PathBuf>,
    /// The psp mode: each sample's homozygote excess read from this parameters file (its
    /// inbreeding coefficients, joined by sample name) instead of fitted from the cohort's
    /// ordinary positions. A file `estimate-parameters` fitted on the same psps holds the same
    /// numbers; reading it saves the SNP/indel fit's minutes. **Refused** when the file's samples
    /// are not the psps' or a coefficient is not `fitted_here` or `defaulted` (a run given
    /// `--inbreeding` writes `supplied`). **One case is not caught**: a sample whose excess came
    /// out at exactly 1.0 is written as 0.0 `defaulted` (a coefficient stops short of 1), and is
    /// read here as 0.0 where `estimate-parameters` fitted the tracts at 1.0. Check the file's
    /// `defaulted` rows, or leave this flag out to fit the excess here.
    #[arg(long, requires = "psps")]
    parameters: Option<PathBuf>,
    /// Samples in the cohort.
    #[arg(long)]
    samples: Option<usize>,
    /// Mean spanning reads per sample per locus.
    #[arg(long)]
    depth: Option<f64>,
    #[arg(long, default_value_t = 1)]
    seed: u64,
    /// `pcr` or `pcr_free` (the default truth, `CohortModel::tomato_like`), or a TOML model file.
    #[arg(long, default_value = "pcr")]
    truth: String,
    /// Loci per class: `tomato2` (as tomato2 has them) or a number.
    #[arg(long)]
    loci: Option<String>,
    /// At most this many loci per class.
    #[arg(long)]
    loci_at_most: Option<usize>,
    /// Only these reference repeat counts, e.g. `8-16` or `8,10,12`.
    #[arg(long)]
    repeats: Option<String>,
    /// Every sample's homozygote excess, in place of the truth's.
    #[arg(long)]
    excess: Option<f64>,
    /// One stutter multiplier per library, comma-separated.
    #[arg(long, value_delimiter = ',')]
    library_multipliers: Vec<f64>,
    /// Tell the methods the libraries' true multipliers.
    #[arg(long)]
    multipliers_told: bool,
    /// The methods, comma-separated: `B` (the counting fit, spread from pairs), `B_within` (the
    /// counting fit, spread from each sample's own two copies), `D` (penalised counting),
    /// `D_tract_prior` (the same as `estimate-parameters` runs it) and `D_working_F` (reporting at
    /// one sample the concentration it worked at); the `B_` variants in `method_named` are for
    /// measurements. `A`, the climb, went with it on 2026-10-16.
    #[arg(long, value_delimiter = ',', default_value = "B")]
    methods: Vec<String>,
    /// Folds for the held-out score; 0 switches it off.
    #[arg(long, default_value_t = 5)]
    folds: usize,
    /// At most this many held-out loci scored per class.
    #[arg(long, default_value_t = 100)]
    heldout_loci: usize,
    /// The Monte Carlo error a locus's held-out score is drawn down to.
    #[arg(long, default_value_t = 0.05)]
    target_se: f64,
    /// The most Monte Carlo draws a locus gets.
    #[arg(long, default_value_t = 65_536)]
    max_draws: usize,
    /// The held-out sampler beyond two samples: `automatic`, `dirichlet` or `sequential`.
    #[arg(long, default_value = "automatic")]
    sampler: String,
    /// Where the rows go.
    #[arg(long, required_unless_present = "fixture_out")]
    out: Option<PathBuf>,
    /// **Record the reference fixture instead of the rows**: fit the one penalised method named in
    /// `--methods` on every locus, write everything the fit returns to this file as one line of
    /// JSON — every number of every class at full precision, the curves, and the round counts —
    /// and stop. No held-out score is computed. `scripts/research/str_fit_fixture.py` runs the
    /// recorded cells this way and compares the result with the recorded fixture. **Not with
    /// `--out`**: this mode writes no rows. `D_working_F` records the same line as `D`, because
    /// the line holds both concentrations either way.
    #[arg(long, conflicts_with = "out")]
    fixture_out: Option<PathBuf>,
    /// Where the truth goes, one row a class.
    #[arg(long)]
    truth_out: Option<PathBuf>,
    /// The cell's name in the rows; built from the cell's settings by default.
    #[arg(long)]
    cell: Option<String>,
    /// Where every held-out locus's score goes, one row a method variant and locus, with the read
    /// offset most of the locus's reads report and the share the fit gave that offset; the last
    /// column is the class's motif length.
    #[arg(long)]
    locus_out: Option<PathBuf>,
    /// The shares' prior weight method B fits at, in tracts, in place of the counting fit's own
    /// (`SHARE_PRIOR_TRACTS`). For measuring what the prior does.
    #[arg(long)]
    share_prior_tracts: Option<f64>,
    /// Score a class with its own shares but no concentration at its own shares and the stated
    /// concentration, in place of calling's flat bottom rung. For measurement only.
    #[arg(long)]
    fallback_keeps_own_shares: bool,
    /// Score every method variant twice: at calling's fallback (as `<variant>`) and with the class's
    /// own shares kept (as `<variant>_own_shares`). A class whose numbers are the same both ways is
    /// scored once. For measurement only.
    #[arg(long, conflicts_with = "fallback_keeps_own_shares")]
    both_fallbacks: bool,
    /// The held-out score's stated concentration for a class without one, in place of the
    /// period's median (or one where no class has one). For measurement only.
    #[arg(long)]
    stated_concentration: Option<f64>,
    /// How far either side of the reference alleles are drawn, in place of the truth's (14). At
    /// most 127 repeats: the draw keeps each sample's genotype as offsets of one byte.
    #[arg(long)]
    allele_reach: Option<i32>,
    /// One between-locus share F for every class, in place of the truth's curve.
    #[arg(long)]
    between_locus_share: Option<f64>,
    /// The stutter prior's weight, in reads a class, for every penalised method (`D`,
    /// `D_working_F`, `D_tract_prior`), in place of the method's own weight (0 for `D` and
    /// `D_working_F`, one read for `D_tract_prior`). The method keeps its name. For measuring
    /// what the prior does: zero is the fit without it. A negative or non-finite weight is
    /// refused.
    #[arg(long, value_parser = a_stutter_prior_weight)]
    stutter_prior_reads: Option<f64>,
    /// **The proxy for calling** (spec `str_penalised_fit.md` §5.1; plan step D6): every sample
    /// at every locus called from its genotype weights under the truth's numbers and under each
    /// of a few others built from each method's full fit, and the calls compared. One row a
    /// method, set of numbers and class, then one over every class; see "The proxy for calling"
    /// in the module documentation. Drawn mode only, at one stutter multiplier.
    #[arg(long, conflicts_with_all = ["psps", "fixture_out"])]
    calls_out: Option<PathBuf>,
}

/// `--stutter-prior-reads` read as a number, refused unless it is a weight the fit can use
/// ([`stutter_prior_weight_problem`]).
fn a_stutter_prior_weight(text: &str) -> Result<f64, String> {
    let reads: f64 = text
        .parse()
        .map_err(|_| format!("{text:?} is not a number"))?;
    stutter_prior_weight_problem(reads).map_or(Ok(reads), Err)
}

/// The penalised fit a name on the command line stands for, or `None` when the name is another
/// method's; at `--stutter-prior-reads` when that is given.
fn penalised_method_named(name: &str, cli: &Cli) -> Option<ThePenalisedFit> {
    let mut fit = match name {
        "D" => ThePenalisedFit::default(),
        "D_working_F" => ThePenalisedFit {
            report_working_concentration: true,
            ..ThePenalisedFit::default()
        },
        "D_tract_prior" => ThePenalisedFit::production(),
        _ => return None,
    };
    if let Some(reads) = cli.stutter_prior_reads {
        fit.config.stutter_prior_reads = reads;
    }
    Some(fit)
}

/// The method a name on the command line stands for. **A new method adds a line here**, or in
/// [`penalised_method_named`] for a variant of the penalised fit.
fn method_named(name: &str, cli: &Cli) -> Box<dyn PeriodFit> {
    if let Some(fit) = penalised_method_named(name, cli) {
        return Box::new(fit);
    }
    match name {
        "B" => Box::new(TheCountingFit {
            share_prior_tracts: cli
                .share_prior_tracts
                .unwrap_or(TheCountingFit::default().share_prior_tracts),
            ..TheCountingFit::default()
        }),
        "B_within" => Box::new(TheCountingFit {
            pairs: false,
            ..TheCountingFit::default()
        }),
        "B_noerrors" => Box::new(TheCountingFit {
            errors: CountingFitErrorMethod::None,
            ..TheCountingFit::default()
        }),
        "B_tight" => Box::new(TheCountingFit {
            tolerance: 1e-10,
            max_rounds: 20_000,
            errors: CountingFitErrorMethod::None,
            ..TheCountingFit::default()
        }),
        "B_noprior" => Box::new(TheCountingFit {
            share_prior_tracts: 0.0,
            ..TheCountingFit::default()
        }),
        "B_bootstrap" => Box::new(TheCountingFit {
            errors: CountingFitErrorMethod::Bootstrap(50),
            ..TheCountingFit::default()
        }),
        other => panic!(
            "no method is called {other}; the methods are B, B_within, B_noerrors, \
             B_noprior, B_tight, B_bootstrap, D, D_working_F and D_tract_prior"
        ),
    }
}

fn repeats_wanted(text: &str) -> Vec<u64> {
    if let Some((low, high)) = text.split_once('-') {
        let low: u64 = low.parse().expect("a repeat count");
        let high: u64 = high.parse().expect("a repeat count");
        (low..=high).collect()
    } else {
        text.split(',')
            .map(|value| value.parse().expect("a repeat count"))
            .collect()
    }
}

fn model_of(cli: &Cli) -> CohortModel {
    let period = cli.period.unwrap_or(1);
    let mut model = match cli.truth.as_str() {
        "pcr" => CohortModel::tomato_like(period, StutterRegime::Pcr),
        "pcr_free" => CohortModel::tomato_like(period, StutterRegime::PcrFree),
        path => toml::from_str(&std::fs::read_to_string(path).expect("the model file reads"))
            .expect("the model file parses"),
    };
    if let Some(samples) = cli.samples {
        model.samples = samples;
    }
    if let Some(depth) = cli.depth {
        model.depth = depth;
    }
    model.seed = cli.seed;
    model.loci = match (cli.loci.as_deref(), model.loci.clone()) {
        (Some("tomato2"), _) => LociPerClass::Tomato2 {
            at_most: cli.loci_at_most,
        },
        (Some(number), _) => LociPerClass::Each(
            number
                .parse::<usize>()
                .expect("loci a class: tomato2 or a number")
                .min(cli.loci_at_most.unwrap_or(usize::MAX)),
        ),
        (None, LociPerClass::Tomato2 { at_most }) => LociPerClass::Tomato2 {
            at_most: cli.loci_at_most.or(at_most),
        },
        (None, LociPerClass::Each(loci)) => {
            LociPerClass::Each(loci.min(cli.loci_at_most.unwrap_or(usize::MAX)))
        }
        (None, LociPerClass::Listed(listed)) => LociPerClass::Listed(
            listed
                .into_iter()
                .map(|(repeats, loci)| (repeats, loci.min(cli.loci_at_most.unwrap_or(usize::MAX))))
                .collect(),
        ),
    };
    if let Some(repeats) = &cli.repeats {
        model.repeat_counts = Some(repeats_wanted(repeats));
    }
    if let Some(reach) = cli.allele_reach {
        model.true_allele_reach = reach;
    }
    if let Some(share) = cli.between_locus_share {
        model.between_locus_share = ParameterCurve::Constant(share);
    }
    if let Some(excess) = cli.excess {
        model.homozygote_excess = HomozygoteExcessTruth::Each(excess);
    }
    if !cli.library_multipliers.is_empty() {
        model
            .library_multipliers
            .clone_from(&cli.library_multipliers);
    }
    if cli.multipliers_told {
        model.multipliers_told = true;
    }
    model
}

/// The rows a run writes, with the cell's columns in front of each.
struct Rows {
    cell_columns: String,
    lines: Vec<String>,
}

fn number(value: Option<f64>) -> String {
    match value {
        Some(value) if value.is_finite() => format!("{value:.8}"),
        _ => "NA".to_owned(),
    }
}

impl Rows {
    fn push(
        &mut self,
        method: &str,
        repeats: &str,
        parameter: &str,
        truth: Option<f64>,
        estimate: Option<f64>,
        error: Option<f64>,
    ) {
        self.lines.push(format!(
            "{}\t{method}\t{repeats}\t{parameter}\t{}\t{}\t{}",
            self.cell_columns,
            number(truth),
            number(estimate),
            number(error)
        ));
    }
}

/// Every number of every class against its truth.
fn recovery_rows(rows: &mut Rows, method: &str, estimates: &[ClassEstimate], truth: &[ClassTruth]) {
    for (estimate, truth) in estimates.iter().zip(truth) {
        let repeats = truth.stratum.reference_repeats.to_string();
        let slippage = estimate.slippage;
        let errors = &estimate.errors;
        rows.push(
            method,
            &repeats,
            "level",
            Some(truth.slippage.level),
            slippage.map(|s| s.level),
            errors.level,
        );
        rows.push(
            method,
            &repeats,
            "shorter_share",
            Some(truth.slippage.shorter_share),
            slippage.map(|s| s.shorter_share),
            errors.shorter_share,
        );
        rows.push(
            method,
            &repeats,
            "fall_off",
            Some(truth.slippage.fall_off),
            slippage.map(|s| s.fall_off),
            errors.fall_off,
        );
        rows.push(
            method,
            &repeats,
            "concentration",
            Some(truth.concentration),
            estimate.concentration,
            errors.concentration,
        );
        rows.push(
            method,
            &repeats,
            "between_locus_share",
            Some(truth.between_locus_share()),
            estimate.concentration.map(|c| 1.0 / (1.0 + c)),
            errors
                .concentration
                .zip(estimate.concentration)
                .map(|(error, c)| error / ((1.0 + c) * (1.0 + c))),
        );
        let shares = truth.shares_in_window(ALLELE_SPAN);
        for (index, share) in shares.iter().enumerate() {
            let offset = index as i32 - ALLELE_SPAN;
            rows.push(
                method,
                &repeats,
                &format!("share_{offset}"),
                Some(*share),
                estimate
                    .allele_length_variation
                    .as_ref()
                    .and_then(|alv| alv.get(index).copied()),
                errors.allele_length_variation.get(index).copied().flatten(),
            );
        }
    }
}

/// The held-out rows of one method variant: one set a class, then one over every class (`all`).
/// Without a truth (the psp mode) the truth columns are `NA`.
fn heldout_rows(
    rows: &mut Rows,
    method: &str,
    per_class: &[HeldOutSums],
    classes: &[StratumEvidence],
) {
    let mut all = HeldOutSums::default();
    let write = |rows: &mut Rows, repeats: &str, sums: &HeldOutSums| {
        if sums.loci == 0 {
            return;
        }
        let reads = sums.reads.max(1) as f64;
        let loci = sums.loci as f64;
        let error = sums.variance.sqrt();
        rows.push(
            method,
            repeats,
            "heldout_ln_per_read",
            Some(sums.truth_ln / reads),
            Some(sums.ln / reads),
            Some(error / reads),
        );
        rows.push(
            method,
            repeats,
            "heldout_ln_per_locus",
            Some(sums.truth_ln / loci),
            Some(sums.ln / loci),
            Some(error / loci),
        );
        rows.push(method, repeats, "heldout_loci", None, Some(loci), None);
        rows.push(
            method,
            repeats,
            "heldout_loci_on_fallback",
            None,
            Some(sums.fallback_loci as f64),
            None,
        );
        rows.push(
            method,
            repeats,
            "heldout_reads",
            None,
            Some(sums.reads as f64),
            None,
        );
    };
    for (sums, evidence) in per_class.iter().zip(classes) {
        write(rows, &evidence.stratum.reference_repeats.to_string(), sums);
        all.ln += sums.ln;
        all.variance += sums.variance;
        all.truth_ln += sums.truth_ln;
        all.reads += sums.reads;
        all.loci += sums.loci;
        all.fallback_loci += sums.fallback_loci;
    }
    write(rows, "all", &all);
}

fn peak_rss_mib() -> Option<f64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    status.lines().find_map(|line| {
        line.strip_prefix("VmHWM:")?
            .split_whitespace()
            .next()?
            .parse::<f64>()
            .ok()
            .map(|kib| kib / 1024.0)
    })
}

/// A number for the fixture at full precision: JSON's shortest form that reads back to the same
/// `f64`. JSON has no NaN or infinity, so those are written as the strings `"NaN"`, `"inf"` and
/// `"-inf"`, which keeps them apart from an absent number (`null`).
fn exact(value: f64) -> Value {
    if value.is_finite() {
        Value::from(value)
    } else {
        Value::String(format!("{value}"))
    }
}

fn exact_or_null(value: Option<f64>) -> Value {
    value.map_or(Value::Null, exact)
}

fn exact_all(values: &[f64]) -> Value {
    Value::Array(values.iter().copied().map(exact).collect())
}

fn curve_json(curve: &FittedCurve) -> Value {
    json!({
        "first_repeats": curve.first_repeats,
        "values": exact_all(&curve.values),
        "smoothing": exact(curve.smoothing),
        "eased_at": curve.eased_at,
    })
}

/// **Everything one penalised fit of one motif length returned**, as the fixture's line for one
/// cell and seed. Times are kept apart under `seconds`, which a comparison ignores; `settings` is
/// the fit's configuration as Rust prints it, for the reader, and is not compared either.
///
/// **`input` records what the fit was given**: the sample-locus rows with reads over every class,
/// and each class's loci with reads. A comparison checks these first, so a change in the draw
/// (`cohort_draw.rs`, the `rand` crate) is reported as a moved input rather than as a moved fit.
fn fixture_line(
    cell_columns: Value,
    method: &ThePenalisedFit,
    classes: &[StratumEvidence],
    fitted: &PenalisedPeriodFit,
    draw_seconds: f64,
    fit_seconds: f64,
) -> Value {
    let class_lines: Vec<Value> = classes
        .iter()
        .zip(&fitted.classes)
        .map(|(evidence, class)| {
            let repeats = evidence.stratum.reference_repeats;
            match class {
                None => json!({ "repeats": repeats, "fitted": false }),
                Some(class) => json!({
                    "repeats": repeats,
                    "fitted": true,
                    "level": exact(class.slippage.level),
                    "shorter_share": exact(class.slippage.shorter_share),
                    "fall_off": exact(class.slippage.fall_off),
                    "allele_length_variation": exact_all(&class.allele_length_variation),
                    "concentration": exact_or_null(class.concentration),
                    "between_locus_share": exact_or_null(class.between_locus_share),
                    "working_concentration": exact(class.working_concentration),
                    "errors": class.errors.as_ref().map_or(Value::Null, |errors| json!({
                        "level": exact_or_null(errors.level),
                        "shorter_share": exact_or_null(errors.shorter_share),
                        "fall_off": exact_or_null(errors.fall_off),
                        "concentration": exact_or_null(errors.concentration),
                        "between_locus_share": exact_or_null(errors.between_locus_share),
                        "allele_length_variation": exact_all(&errors.allele_length_variation),
                    })),
                }),
            }
        })
        .collect();
    let curves = fitted.curves.as_ref().map_or(Value::Null, |curves| {
        json!({
            "level": curve_json(&curves.level),
            "shorter_share": curve_json(&curves.shorter_share),
            "fall_off": curve_json(&curves.fall_off),
            "between_locus_share": curves.between_locus_share.as_ref().map_or(Value::Null, curve_json),
            "allele_length_variation":
                Value::Array(curves.allele_length_variation.iter().map(curve_json).collect()),
        })
    });
    let input = json!({
        "rows_with_reads": classes.iter().map(StratumEvidence::rows_with_reads).sum::<usize>(),
        "classes": classes
            .iter()
            .map(|evidence| json!({
                "repeats": evidence.stratum.reference_repeats,
                "loci_with_reads":
                    evidence.tracts.iter().filter(|tract| !tract.samples.is_empty()).count(),
            }))
            .collect::<Vec<_>>(),
    });
    json!({
        "run": cell_columns,
        "input": input,
        "method": method.name(),
        "rounds": fitted.rounds,
        "rounds_per_counting": fitted.rounds_per_counting,
        "rounds_choosing_strength": fitted.rounds_choosing_strength,
        "pair_alternations": fitted.pair_alternations,
        "class_rounds_counted": fitted.class_rounds_counted,
        "last_change": exact(fitted.last_change),
        "composite_log_likelihood": exact(fitted.composite_log_likelihood),
        "classes": class_lines,
        "curves": curves,
        "seconds": {
            "draw": draw_seconds,
            "fit": fit_seconds,
            "expectation": fitted.expectation_seconds,
            "maximisation": fitted.maximisation_seconds,
            "strength_search": fitted.strength_search_seconds,
            "pair": fitted.pair_seconds,
            "errors": fitted.error_seconds,
        },
        "settings": format!("{:?}", method.config),
    })
}

/// How the held-out score is taken, from the command line. The fold seed is the bench's own
/// default here; the psp mode replaces it with `--seed`.
fn held_out_settings(cli: &Cli) -> HeldOutSettings {
    HeldOutSettings {
        folds: cli.folds,
        loci_scored_per_class: cli.heldout_loci,
        target_standard_error: cli.target_se,
        max_draws: cli.max_draws,
        sampler: match cli.sampler.as_str() {
            "automatic" => Sampler::Automatic,
            "dirichlet" => Sampler::DirichletOnly,
            "sequential" => Sampler::SequentialOnly,
            other => panic!("no sampler is called {other}"),
        },
        ..HeldOutSettings::default()
    }
}

/// The fallbacks every variant is scored at, each with the suffix its rows are named with.
fn fallbacks_asked(cli: &Cli) -> Vec<(&'static str, ScoringFallback)> {
    if cli.both_fallbacks {
        vec![
            (
                "",
                ScoringFallback {
                    keep_own_shares: false,
                    stated_concentration: cli.stated_concentration,
                },
            ),
            (
                "_own_shares",
                ScoringFallback {
                    keep_own_shares: true,
                    stated_concentration: cli.stated_concentration,
                },
            ),
        ]
    } else {
        vec![(
            "",
            ScoringFallback {
                keep_own_shares: cli.fallback_keeps_own_shares,
                stated_concentration: cli.stated_concentration,
            },
        )]
    }
}

/// A run's held-out rows for one method: its times, then every variant's scores.
fn held_out_rows(rows: &mut Rows, base: &str, run: &HeldOutRun, classes: &[StratumEvidence]) {
    rows.push(
        base,
        "all",
        "heldout_fit_seconds",
        None,
        Some(run.fit_seconds),
        None,
    );
    rows.push(
        base,
        "all",
        "heldout_score_seconds",
        None,
        Some(run.score_seconds),
        None,
    );
    for (variant, per_class) in run.variants.iter().zip(&run.sums) {
        heldout_rows(rows, variant, per_class, classes);
    }
}

/// A run's scored loci as `--locus-out` lines, one a locus and variant.
fn push_locus_lines(lines: &mut Vec<String>, run: &HeldOutRun, classes: &[StratumEvidence]) {
    for scored in &run.loci {
        let stratum = classes[scored.class].stratum;
        lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:e}\t{}",
            run.variants[scored.variant],
            stratum.reference_repeats,
            scored.locus,
            scored.ln_likelihood,
            scored.truth_ln,
            scored.reads,
            scored.most_reported_offset,
            scored.fitted_share,
            stratum.period
        ));
    }
}

/// The `--locus-out` file, if asked for.
fn write_locus_lines(cli: &Cli, locus_lines: &[String]) {
    if let Some(path) = &cli.locus_out {
        let mut out = BufWriter::new(File::create(path).expect("the locus file opens"));
        writeln!(
            out,
            "variant\trepeats\tlocus\theldout_ln\ttruth_ln\treads\tmost_reported_offset\t\
             fitted_share\tperiod"
        )
        .expect("the locus file writes");
        for line in locus_lines {
            writeln!(out, "{line}").expect("the locus file writes");
        }
    }
}

/// The rows file, written beside its final name and moved into place once whole.
fn write_rows(out_path: &std::path::Path, rows: &Rows) {
    let partial = out_path.with_extension("partial");
    {
        let mut out = BufWriter::new(File::create(&partial).expect("the rows file opens"));
        writeln!(
            out,
            "cell\tperiod\tsamples\tdepth\ttruth_variant\tseed\tmethod\trepeats\tparameter\ttruth\t\
             estimate\tstandard_error"
        )
        .expect("the rows file writes");
        for line in &rows.lines {
            writeln!(out, "{line}").expect("the rows file writes");
        }
        out.flush().expect("the rows file writes");
    }
    std::fs::rename(&partial, out_path).expect("the rows file moves into place");
}

/// An error and every cause under it, one a line.
fn with_its_causes(error: &dyn std::error::Error) -> String {
    let mut text = error.to_string();
    let mut cause = error.source();
    while let Some(next) = cause {
        text.push_str(&format!("\n  caused by: {next}"));
        cause = next.source();
    }
    text
}

/// **Each sample's homozygote excess from a parameters file**: its inbreeding coefficient, joined
/// by sample name, in the evidence's sample order.
///
/// The coefficient is the excess only where `estimate-parameters` wrote what its SNP/indel fit
/// found, so the file is refused when it cannot be that:
///
/// - **its samples are not the cohort's** — a file fitted on another cohort, a superset of these
///   psps included, holds that cohort's numbers;
/// - **a coefficient is not `fitted_here` or `defaulted`** — `supplied` means the run was given
///   `--inbreeding`, and the tract fit used the excess, not that number.
///
/// One case passes unseen: an excess of exactly 1.0 cannot be written as a coefficient (they stop
/// short of 1), so the file holds 0.0 `defaulted` where the tract fit used 1.0. The flag's help
/// says so.
fn homozygote_excess_from_file(path: &std::path::Path, samples: &[String]) -> Vec<f64> {
    let file = ParametersFile::from_toml(
        &std::fs::read_to_string(path).expect("the parameters file reads"),
    )
    .unwrap_or_else(|error| panic!("the parameters file parses: {}", with_its_causes(&error)));
    let instead = "fit the homozygote excess here by leaving out --parameters, or give a file \
                   `pvc estimate-parameters` fitted on exactly these psps without --inbreeding";
    let mut fitted_from: Vec<&str> = file
        .fitted_from
        .samples
        .iter()
        .map(String::as_str)
        .collect();
    fitted_from.sort_unstable();
    let mut cohort: Vec<&str> = samples.iter().map(String::as_str).collect();
    cohort.sort_unstable();
    assert!(
        fitted_from == cohort,
        "{} was fitted on other samples than these psps hold ({} against {}: in the file only: \
         {:?}; in the psps only: {:?}); {instead}",
        path.display(),
        fitted_from.len(),
        cohort.len(),
        fitted_from
            .iter()
            .filter(|name| cohort.binary_search(name).is_err())
            .collect::<Vec<_>>(),
        cohort
            .iter()
            .filter(|name| fitted_from.binary_search(name).is_err())
            .collect::<Vec<_>>(),
    );
    let not_fitted: Vec<&str> = file
        .inbreeding
        .by_sample
        .iter()
        .filter(|row| {
            !matches!(
                row.inbreeding_coefficient.warrant,
                Warrant::FittedHere | Warrant::Defaulted
            )
        })
        .map(|row| row.sample.as_str())
        .collect();
    assert!(
        not_fitted.is_empty(),
        "{} holds inbreeding coefficients that were not fitted (supplied with --inbreeding, or \
         borrowed) for sample(s) {not_fitted:?}, so they are not the homozygote excess the \
         repeat-tract fit used; {instead}",
        path.display()
    );
    let by_name: BTreeMap<&str, f64> = file
        .inbreeding
        .by_sample
        .iter()
        .map(|row| (row.sample.as_str(), row.inbreeding_coefficient.value))
        .collect();
    samples
        .iter()
        .map(|name| {
            *by_name.get(name.as_str()).unwrap_or_else(|| {
                panic!(
                    "{} gives no inbreeding coefficient for sample {name}",
                    path.display()
                )
            })
        })
        .collect()
}

/// **The psp mode**: a real cohort's repeat-tract evidence, read and arranged exactly as
/// `estimate-parameters` reads it, each motif length's classes folded by whole locus, each method
/// fitted on the kept folds and scored on the held-out loci. See the module documentation.
fn run_on_psps(cli: &Cli) {
    let out_path = cli.out.clone().expect("--out names where the rows go");
    let reading = Instant::now();
    // **`estimate-parameters`'s own arguments**, so the cohort is opened, checked and put in its
    // one order (samples by name, read groups by declared ID) by the command's own code. The
    // output is never written.
    let args = EstimateParametersArgs {
        reference: cli.reference.clone().expect("--psp requires --reference"),
        catalog: cli.catalog.clone(),
        psps: cli.psps.clone(),
        output: PathBuf::from("never-written.parameters.toml"),
        force: true,
        ploidy: 2,
        inbreeding: None,
        str_param_estimates_at_once: None,
        skip_contamination: false,
    };
    let OpenedCensuses {
        mut evidence, plan, ..
    } = open_the_censuses(&args).unwrap_or_else(|error| panic!("{}", with_its_causes(&error)));
    refuse_another_selection(&evidence, &plan.loci)
        .unwrap_or_else(|error| panic!("{}", with_its_causes(&error)));
    let samples: Vec<String> = evidence.sample_names().map(str::to_owned).collect();
    let excess_started = Instant::now();
    // The one number the repeat-tract fit takes from the SNP/indel fit, as `fit_a_cohort` takes
    // it — or from a file that fit wrote.
    let (homozygote_excess, excess_source) = match &cli.parameters {
        Some(path) => (
            homozygote_excess_from_file(path, &samples),
            "parameters_file",
        ),
        None => {
            let generic = fit_jointly(&mut evidence, &JointFitConfig::default())
                .unwrap_or_else(|error| panic!("{}", with_its_causes(&error)));
            (
                homozygote_excess_in_cohort_order(&evidence, &generic),
                "snp_indel_fit",
            )
        }
    };
    let excess_seconds = excess_started.elapsed().as_secs_f64();
    let contig_of = contig_id_of(&plan);
    let pooled = every_read_group_pooled(&evidence);
    let CohortTractStrata { strata, .. } =
        the_tract_strata_of_a_cohort(&mut evidence, &plan.loci, &contig_of, &pooled)
            .unwrap_or_else(|error| panic!("{}", with_its_causes(&error)));
    drop(evidence);
    let read_seconds = reading.elapsed().as_secs_f64();
    eprintln!(
        "{} psp sample(s): {} classes gathered in {read_seconds:.1} s, of which the homozygote \
         excess ({excess_source}) {excess_seconds:.1} s",
        samples.len(),
        strata.len(),
    );

    let mut by_period: BTreeMap<u8, Vec<StratumEvidence>> = BTreeMap::new();
    for evidence in strata {
        if cli
            .period
            .is_none_or(|period| period == evidence.stratum.period)
        {
            by_period
                .entry(evidence.stratum.period)
                .or_default()
                .push(evidence);
        }
    }
    assert!(
        !by_period.is_empty(),
        "the cohort has no repeat-tract class{}",
        cli.period
            .map_or_else(String::new, |period| format!(" of motif length {period}"))
    );
    let settings = HeldOutSettings {
        seed: cli.seed,
        ..held_out_settings(cli)
    };
    assert!(
        settings.folds > 0,
        "the psp mode computes the held-out score only; --folds 0 leaves it nothing to do"
    );
    let fallbacks = fallbacks_asked(cli);
    let cell = cli
        .cell
        .clone()
        .unwrap_or_else(|| format!("psps_n{}", samples.len()));
    let mut rows = Rows {
        cell_columns: String::new(),
        lines: Vec::new(),
    };
    let mut locus_lines: Vec<String> = Vec::new();
    for (period, classes) in &by_period {
        let rows_with_reads: usize = classes.iter().map(StratumEvidence::rows_with_reads).sum();
        let reads: u64 = classes.iter().map(StratumEvidence::spanning_reads).sum();
        // The depth column: the mean spanning reads of a (locus, sample) pair with reads.
        let depth = reads as f64 / rows_with_reads.max(1) as f64;
        rows.cell_columns = format!(
            "{cell}\t{period}\t{}\t{depth:.4}\treal\t{}",
            samples.len(),
            cli.seed
        );
        rows.push(
            "psps",
            "all",
            "rows_with_reads",
            None,
            Some(rows_with_reads as f64),
            None,
        );
        rows.push(
            "psps",
            "all",
            "read_seconds",
            None,
            Some(read_seconds),
            None,
        );
        let folded = HeldOutClasses::new(
            classes,
            &homozygote_excess,
            fold_plans_across_motif_lengths(classes, &settings),
            classes
                .iter()
                .map(|evidence| class_seed_across_motif_lengths(evidence.stratum))
                .collect(),
        );
        // A locus picked for scoring that no method scores (its class has no numbers from the
        // other folds) shows as the gap between this count and a variant's `heldout_loci`.
        eprintln!(
            "{cell}: motif length {period}: {} classes, {rows_with_reads} (locus, sample) pairs \
             with reads, {} held-out loci to score",
            classes.len(),
            folded.scored().len()
        );
        rows.push(
            "psps",
            "all",
            "heldout_loci_to_score",
            None,
            Some(folded.scored().len() as f64),
            None,
        );
        for name in cli.methods.iter().filter(|name| !name.is_empty()) {
            let method = method_named(name, cli);
            let base = method.name();
            let run = held_out_run(
                method.as_ref(),
                &folded,
                &settings,
                &fallbacks,
                cli.locus_out.is_some(),
            );
            eprintln!(
                "{cell}: motif length {period}: {base} held out {} folds: fits {:.1} s, scores \
                 {:.1} s",
                settings.folds, run.fit_seconds, run.score_seconds
            );
            held_out_rows(&mut rows, &base, &run, classes);
            push_locus_lines(&mut locus_lines, &run, classes);
        }
    }
    write_locus_lines(cli, &locus_lines);
    // The process's peak, over every motif length, under a `period` of `all`.
    rows.cell_columns = format!("{cell}\tall\t{}\tNA\treal\t{}", samples.len(), cli.seed);
    rows.push("process", "all", "peak_rss_mib", None, peak_rss_mib(), None);
    write_rows(&out_path, &rows);
}

/// The header of the `--calls-out` file.
const CALL_COLUMNS: &str = "cell\tperiod\tsamples\tdepth\ttruth_variant\tseed\tmethod\tarm\trepeats\t\
                            true_level\tarm_level\tpairs\twrong_under_truth\tchanged\tright_to_wrong\t\
                            wrong_to_right\twrong_to_wrong";

/// One `--calls-out` line.
fn call_line(
    cell_columns: &str,
    method: &str,
    arm: &str,
    repeats: &str,
    levels: Option<(f64, f64)>,
    changes: &CallChanges,
) -> String {
    format!(
        "{cell_columns}\t{method}\t{arm}\t{repeats}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        number(levels.map(|(truth, _)| truth)),
        number(levels.map(|(_, arm)| arm)),
        changes.pairs,
        changes.wrong_under_reference,
        changes.changed,
        changes.right_to_wrong,
        changes.wrong_to_right,
        changes.wrong_to_wrong()
    )
}

/// **The proxy for calling, for one method's full fit** (`estimates`, one a class): each class's
/// calls under the truth, against its calls under each arm of the module documentation's table,
/// as `--calls-out` lines: one line for each class and arm, then, for each arm, one line over
/// the classes that arm has lines for.
fn call_lines(
    lines: &mut Vec<String>,
    cell_columns: &str,
    method: &str,
    estimates: &[ClassEstimate],
    drawn: &DrawnPeriod,
) {
    const ARMS: [&str; 4] = [
        "level_0.995",
        "fitted_level",
        "fitted_stutter",
        "fitted_all",
    ];
    let scored = numbers_for_scoring(estimates, ScoringFallback::default());
    let mut over_all = [CallChanges::default(); ARMS.len()];
    for (class, evidence) in drawn.strata.iter().enumerate() {
        if evidence.rows_with_reads() == 0 {
            continue;
        }
        let truth = ClassNumbers::of_truth(&drawn.truth[class], ALLELE_SPAN);
        let fitted = estimates[class].slippage;
        let with_stutter = |slippage: Slippage| ClassNumbers {
            slippage,
            ..truth.clone()
        };
        let arms: [Option<ClassNumbers>; ARMS.len()] = [
            Some(with_stutter(Slippage {
                level: 0.995 * truth.slippage.level,
                ..truth.slippage
            })),
            fitted.map(|fitted| {
                with_stutter(Slippage {
                    level: fitted.level,
                    ..truth.slippage
                })
            }),
            fitted.map(with_stutter),
            scored[class].as_ref().map(|scored| scored.numbers.clone()),
        ];
        let caller = ClassCaller::of(evidence, &drawn.homozygote_excess);
        let genotypes = drawn_genotypes_in_window(&drawn.genotypes[class], ALLELE_SPAN);
        let under_truth = caller.under(&truth);
        let repeats = evidence.stratum.reference_repeats.to_string();
        for ((name, numbers), sum) in ARMS.iter().zip(&arms).zip(&mut over_all) {
            let Some(numbers) = numbers else { continue };
            let changes = call_changes(&under_truth, &caller.under(numbers), &genotypes);
            sum.add(&changes);
            lines.push(call_line(
                cell_columns,
                method,
                name,
                &repeats,
                Some((truth.slippage.level, numbers.slippage.level)),
                &changes,
            ));
        }
    }
    for (name, sum) in ARMS.iter().zip(&over_all) {
        lines.push(call_line(cell_columns, method, name, "all", None, sum));
    }
}

/// The `--calls-out` file, if asked for, written beside its final name and moved into place.
fn write_call_lines(cli: &Cli, lines: &[String]) {
    if let Some(path) = &cli.calls_out {
        let partial = path.with_extension("partial");
        {
            let mut out = BufWriter::new(File::create(&partial).expect("the calls file opens"));
            writeln!(out, "{CALL_COLUMNS}").expect("the calls file writes");
            for line in lines {
                writeln!(out, "{line}").expect("the calls file writes");
            }
            out.flush().expect("the calls file writes");
        }
        std::fs::rename(&partial, path).expect("the calls file moves into place");
    }
}

fn main() {
    let cli = Cli::parse();
    if !cli.psps.is_empty() {
        run_on_psps(&cli);
        return;
    }
    let model = model_of(&cli);
    // Checked on the model, not the flags, so that a truth file's numbers are refused too.
    if model.true_allele_reach > i32::from(i8::MAX) {
        Cli::command()
            .error(
                ErrorKind::ValueValidation,
                format!(
                    "the allele reach is {}, above the 127 repeats the draw can record \
                     (--allele-reach or the truth file's true_allele_reach)",
                    model.true_allele_reach
                ),
            )
            .exit();
    }
    // The reference calls are made at the truth's one stutter, which every sample's reads were
    // drawn through only when every library's multiplier is one.
    if cli.calls_out.is_some() && model.library_multipliers.iter().any(|m| *m != 1.0) {
        Cli::command()
            .error(
                ErrorKind::ArgumentConflict,
                "--calls-out compares calls at one stutter, so every library's multiplier must \
                 be 1 (--library-multipliers or the truth file's library_multipliers)",
            )
            .exit();
    }
    let mut truth_variant = match cli.truth.as_str() {
        "pcr" | "pcr_free" => cli.truth.clone(),
        path => PathBuf::from(path).file_stem().map_or_else(
            || "file".to_owned(),
            |stem| stem.to_string_lossy().into_owned(),
        ),
    };
    // A truth whose between-locus share was replaced is another truth, and its cells are named so.
    if let Some(share) = cli.between_locus_share {
        truth_variant = format!("{truth_variant}_F{share}");
    }
    let cell = cli.cell.clone().unwrap_or_else(|| {
        format!(
            "p{}_n{}_d{}_{}",
            model.period, model.samples, model.depth, truth_variant
        )
    });

    let drawing = Instant::now();
    let drawn = draw_period(&model, RECORDED_OFFSET_RANGE);
    let draw_seconds = drawing.elapsed().as_secs_f64();
    let rows_with_reads: usize = drawn
        .strata
        .iter()
        .map(StratumEvidence::rows_with_reads)
        .sum();
    eprintln!(
        "{cell} seed {}: drew {} classes, {} loci, {rows_with_reads} sample-locus rows with reads \
         in {draw_seconds:.1} s",
        model.seed,
        drawn.strata.len(),
        drawn.truth.iter().map(|t| t.loci).sum::<usize>(),
    );
    if let Some(path) = &cli.truth_out {
        let file = File::create(path).expect("the truth file opens");
        write_truth_tsv(&drawn.truth, ALLELE_SPAN, BufWriter::new(file))
            .expect("the truth file writes");
    }

    if let Some(path) = &cli.fixture_out {
        let [name] = cli.methods.as_slice() else {
            panic!("--fixture-out records one method; name exactly one in --methods");
        };
        let method = penalised_method_named(name, &cli)
            .unwrap_or_else(|| panic!("--fixture-out records a penalised fit; {name} is not one"));
        let started = Instant::now();
        let fitted =
            method.fit_period_in_full(&drawn.strata, &drawn.homozygote_excess, ALLELE_SPAN);
        let fit_seconds = started.elapsed().as_secs_f64();
        eprintln!(
            "{cell}: {} fitted every locus in {fit_seconds:.1} s, {} rounds",
            method.name(),
            fitted.rounds
        );
        let line = fixture_line(
            json!({
                "cell": cell,
                "period": model.period,
                "samples": model.samples,
                "depth": exact(model.depth),
                "truth_variant": truth_variant,
                "seed": model.seed,
            }),
            &method,
            &drawn.strata,
            &fitted,
            draw_seconds,
            fit_seconds,
        );
        let partial = path.with_extension("partial");
        std::fs::write(&partial, format!("{line}\n")).expect("the fixture file writes");
        std::fs::rename(&partial, path).expect("the fixture file moves into place");
        return;
    }
    let out_path = cli
        .out
        .clone()
        .expect("--out is required without --fixture-out");

    let mut rows = Rows {
        cell_columns: format!(
            "{cell}\t{}\t{}\t{}\t{truth_variant}\t{}",
            model.period, model.samples, model.depth, model.seed
        ),
        lines: Vec::new(),
    };
    rows.push(
        "draw",
        "all",
        "draw_seconds",
        None,
        Some(draw_seconds),
        None,
    );
    rows.push(
        "draw",
        "all",
        "rows_with_reads",
        None,
        Some(rows_with_reads as f64),
        None,
    );

    let settings = held_out_settings(&cli);
    // A drawn class is seeded by its repeat count alone (one motif length a run).
    let folded = HeldOutClasses::new(
        &drawn.strata,
        &drawn.homozygote_excess,
        drawn
            .strata
            .iter()
            .map(|evidence| FoldPlan::of(evidence, &settings, evidence.stratum.reference_repeats))
            .collect(),
        drawn
            .strata
            .iter()
            .map(|evidence| evidence.stratum.reference_repeats)
            .collect(),
    );
    // The truth's own score on every scored locus, once.
    let truth_scores: Vec<f64> = if settings.folds > 0 {
        let started = Instant::now();
        let scores: Vec<(f64, usize, f64)> = folded
            .scored()
            .par_iter()
            .map(|scored| {
                let score = folded.score(
                    scored.class,
                    scored.locus,
                    &ClassNumbers::of_truth(&drawn.truth[scored.class], ALLELE_SPAN),
                    &settings,
                );
                (score.ln_likelihood, score.draws, score.variance.sqrt())
            })
            .collect();
        let mut draws: Vec<usize> = scores.iter().map(|(_, draws, _)| *draws).collect();
        draws.sort_unstable();
        let short_of_target = scores
            .iter()
            .filter(|(_, _, error)| *error >= settings.target_standard_error)
            .count();
        let quantile = |q: f64| {
            draws
                .get(((draws.len() as f64 - 1.0) * q) as usize)
                .copied()
        };
        eprintln!(
            "{cell}: scored the truth on {} held-out loci in {:.1} s; Monte Carlo draws a locus: \
             median {:?}, 90th percentile {:?}, most {:?}; {short_of_target} loci short of the \
             error target",
            folded.scored().len(),
            started.elapsed().as_secs_f64(),
            quantile(0.5),
            quantile(0.9),
            draws.last(),
        );
        rows.push(
            "truth",
            "all",
            "heldout_median_draws",
            None,
            quantile(0.5).map(|d| d as f64),
            None,
        );
        rows.push(
            "truth",
            "all",
            "heldout_loci_short_of_target",
            None,
            Some(short_of_target as f64),
            None,
        );
        let scores: Vec<f64> = scores.into_iter().map(|(ln, _, _)| ln).collect();
        rows.push(
            "truth",
            "all",
            "heldout_score_seconds",
            None,
            Some(started.elapsed().as_secs_f64()),
            None,
        );
        scores
    } else {
        Vec::new()
    };

    let folded = if settings.folds > 0 {
        folded.with_truth_scores(&truth_scores)
    } else {
        folded
    };
    let fallbacks = fallbacks_asked(&cli);

    let mut locus_lines: Vec<String> = Vec::new();
    let mut calls: Vec<String> = Vec::new();
    // An empty name is no method: a run that only draws, to measure the draw's own memory.
    for name in cli.methods.iter().filter(|name| !name.is_empty()) {
        let method = method_named(name, &cli);
        let base = method.name();
        let started = Instant::now();
        let estimates = method.fit_period(&drawn.strata, &drawn.homozygote_excess, ALLELE_SPAN);
        let fit_seconds = started.elapsed().as_secs_f64();
        eprintln!("{cell}: {base} fitted every locus in {fit_seconds:.1} s");
        rows.push(&base, "all", "fit_seconds", None, Some(fit_seconds), None);
        recovery_rows(&mut rows, &base, &estimates, &drawn.truth);
        if cli.calls_out.is_some() {
            let started = Instant::now();
            call_lines(&mut calls, &rows.cell_columns, &base, &estimates, &drawn);
            eprintln!(
                "{cell}: {base} calls compared in {:.1} s",
                started.elapsed().as_secs_f64()
            );
        }

        if settings.folds > 0 {
            let run = held_out_run(
                method.as_ref(),
                &folded,
                &settings,
                &fallbacks,
                cli.locus_out.is_some(),
            );
            eprintln!(
                "{cell}: {base} held out {} folds: fits {:.1} s, scores {:.1} s",
                settings.folds, run.fit_seconds, run.score_seconds
            );
            held_out_rows(&mut rows, &base, &run, &drawn.strata);
            push_locus_lines(&mut locus_lines, &run, &drawn.strata);
        }
    }
    write_locus_lines(&cli, &locus_lines);
    write_call_lines(&cli, &calls);
    rows.push("process", "all", "peak_rss_mib", None, peak_rss_mib(), None);
    write_rows(&out_path, &rows);
}
