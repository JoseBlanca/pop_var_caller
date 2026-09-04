//! The calling stage's command surface over stored files: what a person may type, what is
//! refused before a block is decoded, and what the run says about a cohort it read rather than
//! walked.

use super::*;
use clap::Parser;

use crate::pop_var_caller_exp::cli::{Cli, PopVarCallerExpCommand};
use crate::pop_var_caller_exp::generate_psps::{GeneratePspsArgs, psp_path_for, run_generate_psps};
use crate::pop_var_caller_exp::test_fixtures::{ACohortOnDisk, a_cohort_on_disk};

/// Parse an argument vector into this subcommand's arguments, refusing any other subcommand.
fn args_of(argv: &[&str]) -> CallFromPspsArgs {
    match Cli::parse_from(argv).cmd {
        PopVarCallerExpCommand::CallFromPsps(args) => args,
        other => panic!("expected call-from-psps, got {other:?}"),
    }
}

fn refusal_of(argv: &[&str]) -> clap::Error {
    Cli::try_parse_from(argv).expect_err("this argument vector must be refused")
}

/// The shortest run a person can type.
fn a_defaults_run() -> Vec<&'static str> {
    vec![
        "pop_var_caller_exp",
        "call-from-psps",
        "--reference",
        "ref.fa",
        "--psp",
        "zeta.psp",
        "--defaults",
        "--output",
        "calls.vcf",
    ]
}

#[test]
fn the_subcommand_is_spelled_call_from_psps() {
    let args = args_of(&a_defaults_run());
    assert_eq!(args.reference, PathBuf::from("ref.fa"));
    assert_eq!(args.psps, vec![PathBuf::from("zeta.psp")]);
    assert_eq!(args.output, PathBuf::from("calls.vcf"));
    assert!(
        Cli::try_parse_from(["pop_var_caller_exp", SUBCOMMAND, "--help"])
            .expect_err("--help exits")
            .to_string()
            .contains(SUBCOMMAND),
        "the name this module records is the one clap answers to",
    );
}

/// **A run naming neither a parameters file nor the defaults is refused, and told both
/// answers** — the same group direct mode uses, so the two commands cannot come to disagree
/// about how a run says where its numbers come from.
#[test]
fn a_run_that_names_neither_a_parameters_file_nor_the_defaults_is_refused() {
    let refused = refusal_of(&[
        "pop_var_caller_exp",
        "call-from-psps",
        "--reference",
        "ref.fa",
        "--psp",
        "zeta.psp",
        "--output",
        "calls.vcf",
    ]);
    let rendered = refused.to_string();
    assert!(rendered.contains("--parameters"), "and got: {rendered}");
    assert!(rendered.contains("--defaults"), "and got: {rendered}");
}

/// **There is no `--regions` flag, and that is a design decision rather than an omission**
/// (spec §5.3). A psp records the ground its walk covered and the cohort is refused unless
/// every file agrees about it, so a flag here could only let a person contradict the files.
#[test]
fn the_ground_cannot_be_narrowed_by_a_flag() {
    let mut argv = a_defaults_run();
    argv.extend(["--regions", "some.bed"]);

    let rendered = refusal_of(&argv).to_string();
    assert!(
        rendered.contains("--regions"),
        "the flag is refused by name, and got: {rendered}",
    );
}

/// **The psp flag repeats and keeps the order it was given**, which is the order of the VCF's
/// sample columns.
#[test]
fn the_psp_flag_repeats_and_keeps_the_order_it_was_given() {
    let mut argv = vec![
        "pop_var_caller_exp",
        "call-from-psps",
        "--reference",
        "ref.fa",
        "--defaults",
        "--output",
        "calls.vcf",
    ];
    for file in ["zeta.psp", "alpha.psp", "beta.psp"] {
        argv.extend(["--psp", file]);
    }

    assert_eq!(
        args_of(&argv).psps,
        vec![
            PathBuf::from("zeta.psp"),
            PathBuf::from("alpha.psp"),
            PathBuf::from("beta.psp"),
        ],
    );
}

/// **A directory contributes its psps in name order, and only its psps.**
///
/// Name order rather than the filesystem's, so two runs naming one directory open the same
/// cohort in the same order — the property §12.6's oracle is about, made true at the door
/// rather than relied on.
#[test]
fn a_directory_contributes_its_psps_in_name_order_and_nothing_else() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    for name in ["zeta.psp", "alpha.psp", "notes.txt", "alpha.psp.7.partial"] {
        std::fs::write(directory.path().join(name), b"x").expect("a file");
    }
    std::fs::create_dir(directory.path().join("inner.psp")).expect("a directory named like one");

    let mut args = args_of(&a_defaults_run());
    args.psps = vec![directory.path().to_path_buf()];

    let found = psps_named_by(&args).expect("the directory lists");

    assert_eq!(
        found,
        vec![
            directory.path().join("alpha.psp"),
            directory.path().join("zeta.psp"),
        ],
        "psps in name order; the text file, the half-written file and the directory are not \
         psps",
    );
}

/// **A file named directly is taken as typed**, whatever it is called — a psp that has been
/// renamed is still a psp, and the reader is what judges it.
#[test]
fn a_file_named_directly_is_taken_as_typed() {
    let mut args = args_of(&a_defaults_run());
    args.psps = vec![PathBuf::from("some/where/else.store")];

    assert_eq!(
        psps_named_by(&args).expect("a named file needs no listing"),
        vec![PathBuf::from("some/where/else.store")],
    );
}

/// **A directory with no psp in it is refused and told where to look.**
///
/// Refused rather than skipped: a directory typed by mistake and one whose walk has not
/// finished look the same from here, and calling a cohort short of a sample is not something to
/// do quietly.
#[test]
fn a_directory_with_no_psp_is_refused_and_says_what_to_do() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    let mut args = args_of(&a_defaults_run());
    args.psps = vec![directory.path().to_path_buf()];

    let refused = psps_named_by(&args).expect_err("an empty directory is not a cohort");

    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(
        rendered.contains("generate-psps"),
        "the message says which command makes them, and got: {rendered}",
    );
}

// ---------------------------------------------------------------------
// The command, driven over psps another command wrote
// ---------------------------------------------------------------------

/// The fixture cohort walked into psps, and the arguments that call them.
///
/// **The psps are named one by one rather than by their directory**, so the run's sample order
/// is the fixture's own — `zeta` then `alpha` — and a test can say which sample a line is about.
fn a_cohort_of_psps() -> (ACohortOnDisk, CallFromPspsArgs) {
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

    let args = CallFromPspsArgs {
        reference: cohort.reference.clone(),
        catalog: Some(cohort.catalog.clone()),
        psps: vec![psp_path_for(&psps, "zeta"), psp_path_for(&psps, "alpha")],
        output: cohort.directory.path().join("calls.vcf"),
        parameters: None,
        defaults: true,
        ploidy: None,
        max_cohort_locus_span: DEFAULT_MAX_COHORT_LOCUS_SPAN,
        max_candidate_alleles: DEFAULT_MAX_CANDIDATE_ALLELES.get(),
        cohort_locus_builder_regions_len: None,
        threads: 0,
        min_copies: MinCopies::default(),
        min_period: DEFAULT_MIN_PERIOD,
        max_period: DEFAULT_MAX_PERIOD,
        max_str_len: DEFAULT_MAX_STR_LEN,
        min_purity: DEFAULT_MIN_PURITY,
    };
    (cohort, args)
}

/// **The command writes a VCF and the parameters file beside it**, over psps another command
/// wrote — which is the whole of psp mode reaching a person's terminal.
///
/// This is the only test that drives `run_call_from_psps` itself, and it is what stops the
/// wiring from being proved only in pieces: direct mode's own review measured that gap and
/// found nine of fourteen mutations survived while every test called the helpers directly.
#[test]
fn the_command_writes_a_vcf_and_the_parameters_beside_it() {
    let (_cohort, args) = a_cohort_of_psps();

    run_call_from_psps(&args).expect("a cohort of stored samples calls");

    assert!(args.output.is_file(), "the VCF is at {:?}", args.output);
    let parameters = beside_the_vcf(&args.output);
    assert!(
        parameters.is_file(),
        "the parameters file is at {parameters:?}",
    );
    let vcf = std::fs::read_to_string(&args.output).expect("the VCF reads");
    let columns = vcf
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("the VCF has a column header");
    let samples: Vec<&str> = columns.split('\t').skip(9).collect();
    assert_eq!(
        samples,
        vec!["zeta", "alpha"],
        "the sample columns follow the order the psps were given",
    );
}

/// **What a person sees when the command finishes**, printed once so a reviewer of this step
/// reads the artefact rather than the code that builds it.
///
/// Ignored, because its value is the eye and not the assertion — the report's own tests
/// (`ng::run::report`) are what pin the lines. Run it deliberately:
///
/// ```text
/// cargo test --lib call_from_psps -- --ignored --nocapture the_report_a_person_sees
/// ```
#[test]
#[ignore = "prints the run report for a person to read; the assertions are in ng::run::report"]
fn the_report_a_person_sees() {
    let (_cohort, args) = a_cohort_of_psps();
    run_call_from_psps(&args).expect("a cohort of stored samples calls");
}

/// **A cohort walked under a different catalog is refused before a block is decoded**, naming
/// what differs — the refusal that makes stored evidence safe to call at all.
#[test]
fn a_run_whose_routing_is_not_the_walks_is_refused() {
    let (_cohort, mut args) = a_cohort_of_psps();
    // A tract that is a repeat to the walk and ordinary sequence to this run.
    args.min_purity = 0.99;

    let refused = run_call_from_psps(&args).expect_err("this run did not route as the walk did");

    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(
        rendered.contains("purity") || rendered.contains("criteria"),
        "the refusal names the field that differs, and got: {rendered}",
    );
}

// ---------------------------------------------------------------------
// The refusals the command must make, provoked THROUGH the command
// ---------------------------------------------------------------------

/// **Every refusal below is provoked through `run_call_from_psps` and not through the helper
/// that makes it**, which is the gap direct mode's own review measured: nine of its fourteen
/// mutations survived because every test called the helpers directly, and deleting the calls
/// from the command changed nothing. Measured here the same way — with these three helper calls
/// deleted from `run_call_from_psps`, every other test in this module stayed green.
///
/// **An `--output` naming a directory** is discovered by the writer only after the last locus
/// has been called — a whole cohort's work — and leaves the in-flight `<output>.tmp` beside the
/// directory the person named.
#[test]
fn an_output_that_names_a_directory_is_refused_by_the_command() {
    let (cohort, mut args) = a_cohort_of_psps();
    args.output = cohort.directory.path().to_path_buf();

    let refused = run_call_from_psps(&args).expect_err("a directory is not somewhere for a VCF");

    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(rendered.contains("is a directory"), "and got: {rendered}");
}

/// **An `--output` in a directory that does not exist** is a typing mistake, and a person
/// should not wait through a cohort's decoding to hear about it.
#[test]
fn an_output_in_a_missing_directory_is_refused_by_the_command() {
    let (cohort, mut args) = a_cohort_of_psps();
    args.output = cohort
        .directory
        .path()
        .join("no-such-directory")
        .join("calls.vcf");

    let refused = run_call_from_psps(&args).expect_err("there is no directory to write into");

    // **The command's own sentence, not the writer's.** The writer names the missing directory
    // too — after the whole cohort has been called — so a test that only asked for the name
    // would pass with this refusal deleted. Measured: it did.
    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(
        rendered.contains("cannot be written: there is no directory"),
        "and got: {rendered}",
    );
    assert!(
        rendered.contains("no-such-directory"),
        "the message names the directory that is missing, and got: {rendered}",
    );
}

/// **A run may not write its parameters file over the one it was handed** (spec §7 invites the
/// collision: copy the file your run wrote, change a line, re-run). What it would destroy is
/// the edit, and the numbers that came back would look ordinary.
#[test]
fn a_run_whose_output_would_overwrite_its_own_parameters_file_is_refused_by_the_command() {
    let (_cohort, mut args) = a_cohort_of_psps();
    run_call_from_psps(&args).expect("the first run writes both files");

    args.parameters = Some(beside_the_vcf(&args.output));
    args.defaults = false;

    let refused = run_call_from_psps(&args).expect_err("this run would write over its own input");

    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(
        rendered.contains("calls.parameters.toml"),
        "the message names the file that would be lost, and got: {rendered}",
    );
}

/// **A `--ploidy` that disagrees with the parameters file is refused** — a file's numbers were
/// fitted at one ploidy and mean nothing at another, so a run that wants the other number wants
/// a different fit. Provoked through the command, because the flag has to reach the check.
#[test]
fn a_ploidy_that_is_not_the_parameters_files_is_refused_by_the_command() {
    let (cohort, mut args) = a_cohort_of_psps();
    run_call_from_psps(&args).expect("a defaults run writes the parameters it scored with");

    // The file the first run wrote, moved out of the way so the second run may read it.
    let supplied = cohort.directory.path().join("supplied.parameters.toml");
    std::fs::rename(beside_the_vcf(&args.output), &supplied).expect("the file moves");
    args.parameters = Some(supplied);
    args.defaults = false;
    args.ploidy = Some(4);

    let refused = run_call_from_psps(&args).expect_err("the flag disagrees with the file");

    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(
        rendered.contains("--ploidy 4") && rendered.contains("ploidy 2"),
        "the message names both numbers, and got: {rendered}",
    );
}

/// **A run reads back the parameters a run over this cohort wrote** — psp mode's whole
/// supplied-file path, which the defaults runs above never touch: the by-name match of samples
/// and read groups (spec §6.2, §12.5), the reference binding, and the numbers reaching the
/// calling loop.
#[test]
fn the_parameters_a_previous_run_wrote_bind_to_this_cohort() {
    let (cohort, mut args) = a_cohort_of_psps();
    run_call_from_psps(&args).expect("a defaults run writes the parameters it scored with");
    let supplied = cohort.directory.path().join("supplied.parameters.toml");
    std::fs::rename(beside_the_vcf(&args.output), &supplied).expect("the file moves");

    args.parameters = Some(supplied);
    args.defaults = false;

    run_call_from_psps(&args).expect("its own parameters bind to the cohort that produced them");
}

/// **A parameters file naming another cohort's samples is refused, naming what is missing**
/// (spec §12.5). The join is by name, so a file whose sample list is not this run's cannot be
/// bound however its rows are ordered.
#[test]
fn parameters_that_name_another_cohorts_samples_are_refused() {
    let (cohort, mut args) = a_cohort_of_psps();
    run_call_from_psps(&args).expect("a defaults run writes the parameters it scored with");
    let supplied = cohort.directory.path().join("supplied.parameters.toml");
    std::fs::rename(beside_the_vcf(&args.output), &supplied).expect("the file moves");
    let text = std::fs::read_to_string(&supplied).expect("the file reads");
    std::fs::write(&supplied, text.replace("zeta", "someone-else")).expect("the file rewrites");

    args.parameters = Some(supplied);
    args.defaults = false;

    let refused = run_call_from_psps(&args).expect_err("these numbers are not this cohort's");

    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(
        rendered.contains("zeta") || rendered.contains("someone-else"),
        "the refusal names a sample the two sides disagree about, and got: {rendered}",
    );
}
