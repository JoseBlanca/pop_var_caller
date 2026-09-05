//! `generate-census` at the command line: what a person may type, what is refused before a psp
//! is read, and what the run says about each file it read.

use super::*;
use clap::Parser;

use crate::pop_var_caller_exp::cli::{Cli, PopVarCallerExpCommand};
use crate::pop_var_caller_exp::generate_psps::{
    GeneratePspsArgs, census_path_for, psp_path_for, run_generate_psps,
};
use crate::pop_var_caller_exp::test_fixtures::{ACohortOnDisk, a_cohort_on_disk};

/// Parse an argument vector into this subcommand's arguments, refusing any other subcommand.
fn args_of(argv: &[&str]) -> GenerateCensusArgs {
    match Cli::parse_from(argv).cmd {
        PopVarCallerExpCommand::GenerateCensus(args) => args,
        other => panic!("expected generate-census, got {other:?}"),
    }
}

/// The shortest run a person can type.
fn a_shortest_run() -> Vec<&'static str> {
    vec![
        "pop_var_caller_exp",
        "generate-census",
        "--reference",
        "ref.fa",
        "--psp",
        "zeta.psp",
        "--output-dir",
        "out",
    ]
}

/// A walked cohort: two psps and their walk-time censuses, in a directory of their own.
fn a_walked_cohort() -> (ACohortOnDisk, PathBuf) {
    let cohort = a_cohort_on_disk();
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

/// This command's arguments over a walked cohort, writing into `output_dir`.
fn args_over(cohort: &ACohortOnDisk, psps: &Path, output_dir: PathBuf) -> GenerateCensusArgs {
    GenerateCensusArgs {
        reference: cohort.reference.clone(),
        catalog: Some(cohort.catalog.clone()),
        psps: vec![psps.to_path_buf()],
        output_dir,
        force: false,
        min_copies: MinCopies::default(),
        min_period: DEFAULT_MIN_PERIOD,
        max_period: DEFAULT_MAX_PERIOD,
        max_str_len: DEFAULT_MAX_STR_LEN,
        min_purity: DEFAULT_MIN_PURITY,
    }
}

#[test]
fn the_subcommand_is_spelled_generate_census() {
    let args = args_of(&a_shortest_run());
    assert_eq!(args.reference, PathBuf::from("ref.fa"));
    assert_eq!(args.psps, vec![PathBuf::from("zeta.psp")]);
    assert_eq!(args.output_dir, PathBuf::from("out"));
    assert!(
        Cli::try_parse_from(["pop_var_caller_exp", SUBCOMMAND, "--help"])
            .expect_err("--help exits")
            .to_string()
            .contains(SUBCOMMAND),
        "the name this module records is the one clap answers to",
    );
}

/// **There is no `--regions` flag, and that is a decision rather than an omission.**
///
/// The psps record the ground they were walked over, and the digest of that ground travels in
/// every census as one of its recording terms. A flag here could only let a person select over
/// ground the files were not walked over, producing censuses the cohort cannot be fitted from —
/// and the disagreement would surface hours later, at the fit.
#[test]
fn the_ground_cannot_be_narrowed_by_a_flag() {
    let mut argv = a_shortest_run();
    argv.extend(["--regions", "some.bed"]);
    let refused = Cli::try_parse_from(argv).expect_err("--regions is not a flag here");
    assert!(
        refused.to_string().contains("--regions"),
        "and got: {refused}",
    );
}

/// **The command's censuses are the walk's censuses, byte for byte.**
///
/// This is the end-to-end form of the producer's own agreement test: `generate-psps` writes a
/// census beside each psp as it walks, `generate-census` builds one from the psp afterwards, and
/// the two files are the same. It is what says the second route can stand in for the first.
#[test]
fn each_census_it_writes_equals_the_one_the_walk_wrote() {
    let (cohort, psps) = a_walked_cohort();
    let rebuilt = cohort.directory.path().join("rebuilt");

    let report =
        build_every_census(&args_over(&cohort, &psps, rebuilt.clone())).expect("the psps read");

    assert_eq!(report.samples.len(), 2, "one entry a sample");
    for sample in &report.samples {
        let walked = std::fs::read(census_path_for(&psps, &sample.sample))
            .expect("the walk wrote this sample's census");
        let built = std::fs::read(&sample.census).expect("this run wrote one too");
        assert_eq!(
            walked, built,
            "{}'s census differs between the walk and the rebuild",
            sample.sample,
        );
    }
}

/// **A census already at the path is refused, and nothing is written.**
///
/// The refusal comes before the first psp is read, so a cohort of sixty is never left with
/// forty replaced files and twenty originals.
#[test]
fn a_census_already_there_is_refused_and_nothing_is_replaced() {
    let (cohort, psps) = a_walked_cohort();
    // The psps' own directory already holds the walk's censuses, so writing there collides.
    let before = std::fs::read(census_path_for(&psps, "zeta")).expect("the walk wrote it");

    let error = build_every_census(&args_over(&cohort, &psps, psps.clone()))
        .expect_err("a census is already there");

    assert!(
        matches!(&error, GenerateCensusCliError::CensusAlreadyThere { .. }),
        "{error:?}",
    );
    let after = std::fs::read(census_path_for(&psps, "zeta")).expect("it is still there");
    assert_eq!(before, after, "the refused run replaced nothing");
}

/// **`--force` replaces them**, and what it writes is what was there — because both producers
/// agree.
#[test]
fn force_replaces_a_census_that_is_already_there() {
    let (cohort, psps) = a_walked_cohort();
    let before = std::fs::read(census_path_for(&psps, "zeta")).expect("the walk wrote it");

    let mut args = args_over(&cohort, &psps, psps.clone());
    args.force = true;
    let report = build_every_census(&args).expect("--force replaces them");

    assert_eq!(report.samples.len(), 2);
    let after = std::fs::read(census_path_for(&psps, "zeta")).expect("it is still there");
    assert_eq!(
        before, after,
        "the replacement is the same census, which is the point of the second producer",
    );
}

/// **A sample whose census holds no read is named, not omitted.**
///
/// The fixture's second sample carries no reads at all, so every kept position it has is a
/// zero — which is the denominator a fit needs and not an error. A run that left it out of its
/// report would leave somebody hunting for a file written exactly as asked.
#[test]
fn a_sample_with_no_reads_is_named_as_contributing_nothing() {
    let (cohort, psps) = a_walked_cohort();
    let rebuilt = cohort.directory.path().join("rebuilt");

    let report = build_every_census(&args_over(&cohort, &psps, rebuilt)).expect("the psps read");

    let alpha = report
        .samples
        .iter()
        .find(|sample| sample.sample == "alpha")
        .expect("the sample with no reads is in the report");
    assert!(
        alpha.tally.contributes_nothing(),
        "alpha has no reads, so no kept locus can have one: {:?}",
        alpha.tally,
    );
    assert!(
        alpha.line().contains("contributes nothing to a fit"),
        "and its line reads: {}",
        alpha.line(),
    );
    assert!(
        report
            .lines()
            .iter()
            .any(|line| line.contains("put nothing into the fit")),
        "the summary counts them too: {:?}",
        report.lines(),
    );
}

/// **Each sample's line names the psp it read and the census it wrote**, with both sizes, so a
/// person can pair the two files without knowing this command's naming rule.
#[test]
fn each_line_names_both_files_and_what_went_into_the_census() {
    let (cohort, psps) = a_walked_cohort();
    let rebuilt = cohort.directory.path().join("rebuilt");

    let report = build_every_census(&args_over(&cohort, &psps, rebuilt)).expect("the psps read");

    let zeta = report
        .samples
        .iter()
        .find(|sample| sample.sample == "zeta")
        .expect("the sample with reads is in the report");
    let line = zeta.line();
    assert!(line.contains("zeta.psp"), "and got: {line}");
    assert!(line.contains("zeta.census"), "and got: {line}");
    assert!(
        line.contains(&format!("{} stored loci read", zeta.records)),
        "and got: {line}",
    );
    assert!(
        line.contains(&format!("census {} bytes", zeta.census_bytes)),
        "and got: {line}",
    );
    assert_eq!(
        zeta.psp,
        psp_path_for(&psps, "zeta"),
        "the psp named is the one this sample's census came from",
    );
}

/// **A directory of psps is expanded in name order**, so two runs naming one directory read the
/// same cohort in the same order however the filesystem answers.
#[test]
fn a_directory_contributes_every_psp_inside_it() {
    let (cohort, psps) = a_walked_cohort();
    let rebuilt = cohort.directory.path().join("rebuilt");
    let args = args_over(&cohort, &psps, rebuilt);

    let expanded = psps_named_by(&args).expect("the directory lists");

    assert_eq!(
        expanded,
        vec![psp_path_for(&psps, "alpha"), psp_path_for(&psps, "zeta")],
        "sorted by name, and nothing but psps",
    );
}

/// **A directory holding no psp is refused by name**, rather than being read as an empty cohort.
#[test]
fn a_directory_with_no_psp_in_it_is_refused() {
    let (cohort, psps) = a_walked_cohort();
    let empty = cohort.directory.path().join("empty");
    std::fs::create_dir_all(&empty).expect("the scratch dir is ours");
    let mut args = args_over(&cohort, &psps, cohort.directory.path().join("rebuilt"));
    args.psps = vec![empty.clone()];

    let error = psps_named_by(&args).expect_err("an empty directory is not a cohort");

    assert!(
        matches!(&error, GenerateCensusCliError::NoPspsInDirectory { path } if path == &empty),
        "{error:?}",
    );
}

/// **A sample whose name could not be a file name is refused before anything is written.**
///
/// `@RG SM` is free header text and travels into the psp, so a name holding a separator would
/// otherwise put a census outside `--output-dir`.
#[test]
fn a_sample_name_that_is_not_a_file_name_is_refused() {
    for name in ["../escape", "a/b", "", ".", ".."] {
        let error = refuse_a_sample_name_that_is_not_a_file_name(name)
            .expect_err("{name:?} cannot be a file name");
        assert!(
            matches!(&error, GenerateCensusCliError::SampleNameNotAFileName { sample } if sample == name),
            "{name:?} gave {error:?}",
        );
    }
    refuse_a_sample_name_that_is_not_a_file_name("zeta").expect("an ordinary name is fine");
}
