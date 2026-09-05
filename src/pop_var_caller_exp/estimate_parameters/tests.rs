//! `estimate-parameters` at the command line: what a person may type, what is refused before a
//! census is opened, and what the file it writes says about where its numbers came from.

use super::*;
use clap::Parser;
use std::path::Path;

use crate::pop_var_caller_exp::cli::{Cli, PopVarCallerExpCommand};
use crate::pop_var_caller_exp::generate_psps::{GeneratePspsArgs, run_generate_psps};
use crate::pop_var_caller_exp::test_fixtures::{AVaryingCohort, a_varying_cohort_on_disk};

/// Parse an argument vector into this subcommand's arguments, refusing any other subcommand.
fn args_of(argv: &[&str]) -> EstimateParametersArgs {
    match Cli::parse_from(argv).cmd {
        PopVarCallerExpCommand::EstimateParameters(args) => args,
        other => panic!("expected estimate-parameters, got {other:?}"),
    }
}

/// The shortest run a person can type.
fn a_shortest_run() -> Vec<&'static str> {
    vec![
        "pop_var_caller_exp",
        "estimate-parameters",
        "--reference",
        "ref.fa",
        "--census",
        "zeta.census",
        "--output",
        "cohort.parameters.toml",
    ]
}

/// **The varying fixture cohort**, walked, with its censuses beside its psps — which is the pair
/// this command is meant to be handed. The plain on-disk cohort has no repeat tract at all, so a
/// fit over it exercises only half of what this file records.
fn a_walked_cohort() -> (AVaryingCohort, PathBuf) {
    let cohort = a_varying_cohort_on_disk();
    let psps = cohort.directory.path().join("psps");
    run_generate_psps(&GeneratePspsArgs {
        reference: cohort.reference.clone(),
        catalog: Some(cohort.catalog.clone()),
        alignments: cohort.alignments.clone(),
        output_dir: psps.clone(),
        regions: None,
        force: false,
        build_index_if_missing: false,
        min_copies: MinCopies::default(),
        min_period: DEFAULT_MIN_PERIOD,
        max_period: DEFAULT_MAX_PERIOD,
        max_str_len: DEFAULT_MAX_STR_LEN,
        min_purity: DEFAULT_MIN_PURITY,
    })
    .expect("the cohort walks into psps");
    (cohort, psps)
}

fn args_over(cohort: &AVaryingCohort, psps: &Path, output: PathBuf) -> EstimateParametersArgs {
    EstimateParametersArgs {
        reference: cohort.reference.clone(),
        catalog: Some(cohort.catalog.clone()),
        censuses: vec![psps.to_path_buf()],
        output,
        force: false,
        ploidy: 2,
        inbreeding: 0.0,
        min_copies: MinCopies::default(),
        min_period: DEFAULT_MIN_PERIOD,
        max_period: DEFAULT_MAX_PERIOD,
        max_str_len: DEFAULT_MAX_STR_LEN,
        min_purity: DEFAULT_MIN_PURITY,
    }
}

#[test]
fn the_subcommand_is_spelled_estimate_parameters() {
    let args = args_of(&a_shortest_run());
    assert_eq!(args.reference, PathBuf::from("ref.fa"));
    assert_eq!(args.censuses, vec![PathBuf::from("zeta.census")]);
    assert_eq!(args.output, PathBuf::from("cohort.parameters.toml"));
    assert!(
        Cli::try_parse_from(["pop_var_caller_exp", SUBCOMMAND, "--help"])
            .expect_err("--help exits")
            .to_string()
            .contains(SUBCOMMAND),
        "the name this module records is the one clap answers to",
    );
}

/// **The command writes a parameters file from a cohort's own censuses.**
///
/// This is the first parameters file this tree produces from data. What it must say is where its
/// numbers came from — the reference, the samples, the read groups, and the census.
#[test]
fn it_writes_a_parameters_file_from_the_cohorts_censuses() {
    let (cohort, psps) = a_walked_cohort();
    let output = cohort.directory.path().join("cohort.parameters.toml");

    let (file, samples) =
        fit_and_assemble(&args_over(&cohort, &psps, output)).expect("the cohort fits");

    assert_eq!(samples, 2, "one entry a sample");
    assert_eq!(file.fitted_from.samples.len(), 2);
    assert!(
        !file.fitted_from.read_groups.is_empty(),
        "the file names the run's read groups",
    );
    assert!(
        !file.fitted_from.census.terms.is_empty(),
        "and names the census these numbers were fitted from, term by term",
    );
    assert_eq!(file.ploidy, 2);
}

/// **The same cohort, fitted twice, writes the same bytes.**
///
/// A parameters file is what a calling run scores with, so two runs over one cohort that
/// disagreed would make the calls depend on which one was kept.
#[test]
fn one_cohort_fitted_twice_writes_the_same_file() {
    let (cohort, psps) = a_walked_cohort();
    let output = cohort.directory.path().join("cohort.parameters.toml");

    let (first, _) =
        fit_and_assemble(&args_over(&cohort, &psps, output.clone())).expect("the cohort fits");
    let (again, _) = fit_and_assemble(&args_over(&cohort, &psps, output)).expect("and again");

    assert_eq!(
        first.to_toml(),
        again.to_toml(),
        "one cohort's parameters are one file",
    );
}

/// **The inbreeding coefficient is recorded as supplied, not fitted.**
///
/// It comes from a sample's own windowed genome histogram, which is the other pre-pass route. A
/// file that reported a declared value as fitted would make a run's own assumption look like a
/// measurement.
#[test]
fn the_inbreeding_coefficient_is_recorded_as_supplied() {
    let (cohort, psps) = a_walked_cohort();
    let output = cohort.directory.path().join("cohort.parameters.toml");
    let mut args = args_over(&cohort, &psps, output);
    args.inbreeding = 0.25;

    let (file, _) = fit_and_assemble(&args).expect("the cohort fits");

    assert!(!file.inbreeding.by_sample.is_empty(), "one row a sample");
    for row in &file.inbreeding.by_sample {
        assert_eq!(row.inbreeding_coefficient.value, 0.25);
        assert_eq!(
            row.inbreeding_coefficient.warrant,
            crate::ng::calling::parameters_file::Warrant::Supplied,
            "a declared coefficient is supplied, never fitted",
        );
    }
}

/// **An output already there is refused before a census is opened**, so a run that cannot write
/// its answer does not spend the fit finding out.
#[test]
fn an_output_already_there_is_refused_before_anything_is_fitted() {
    let (cohort, psps) = a_walked_cohort();
    let output = cohort.directory.path().join("cohort.parameters.toml");
    std::fs::write(&output, b"somebody's file").expect("the scratch dir is ours");

    let error = fit_and_assemble(&args_over(&cohort, &psps, output.clone()))
        .expect_err("the output is already there");

    assert!(
        matches!(&error, EstimateParametersCliError::OutputAlreadyThere { path } if path == &output),
        "{error:?}",
    );
    assert_eq!(
        std::fs::read(&output).expect("still there"),
        b"somebody's file",
        "and nothing was written over it",
    );
}

/// **A ploidy of zero is refused by name**, rather than reaching the fit as a genome with no
/// copies.
#[test]
fn a_ploidy_of_zero_is_refused() {
    let (cohort, psps) = a_walked_cohort();
    let output = cohort.directory.path().join("cohort.parameters.toml");
    let mut args = args_over(&cohort, &psps, output);
    args.ploidy = 0;

    let error = fit_and_assemble(&args).expect_err("zero copies is not a genome");

    assert!(
        matches!(&error, EstimateParametersCliError::NotAValue { what, .. } if *what == "--ploidy"),
        "{error:?}",
    );
}

/// **A directory of censuses is expanded in name order.**
#[test]
fn a_directory_contributes_every_census_inside_it() {
    let (cohort, psps) = a_walked_cohort();
    let args = args_over(&cohort, &psps, cohort.directory.path().join("out.toml"));

    let expanded = censuses_named_by(&args).expect("the directory lists");

    assert_eq!(expanded.len(), 2, "one census a sample, and nothing else");
    assert!(
        expanded.windows(2).all(|pair| pair[0] <= pair[1]),
        "sorted by name, so two runs read the same cohort in the same order: {expanded:?}",
    );
}

/// **A census whose psp is not beside it is refused**, and the refusal comes from the cohort's
/// own door rather than from the fit.
#[test]
fn a_census_without_its_psp_is_refused() {
    let (cohort, psps) = a_walked_cohort();
    let alone = cohort.directory.path().join("alone");
    std::fs::create_dir_all(&alone).expect("the scratch dir is ours");
    let orphan = alone.join("orphan.census");
    let first = censuses_named_by(&args_over(
        &cohort,
        &psps,
        cohort.directory.path().join("out.toml"),
    ))
    .expect("the directory lists")[0]
        .clone();
    std::fs::copy(first, &orphan).expect("a copy");

    let mut args = args_over(&cohort, &psps, cohort.directory.path().join("out.toml"));
    args.censuses = vec![orphan];

    let error = fit_and_assemble(&args).expect_err("nothing can check it");

    assert!(
        matches!(&error, EstimateParametersCliError::Cohort { .. }),
        "{error:?}",
    );
}
