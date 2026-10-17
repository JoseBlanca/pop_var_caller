//! The calling stage's command surface over stored files: what a person may type, what is
//! refused before a block is decoded, and what the run says about a cohort it read rather than
//! walked.

use super::*;
use clap::Parser;

use crate::cli::command_line::{Cli, PopVarCallerCommand};
use crate::cli::generate_psps::{GeneratePspsArgs, psp_path_for, run_generate_psps};
use crate::cli::test_fixtures::{ACohortOnDisk, a_cohort_on_disk};
use crate::region_typing::segment_criteria::{
    DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
};
use crate::region_typing::{DEFAULT_MAX_STR_BUNDLE_LEN, DEFAULT_MAX_STR_LEN};

/// Parse an argument vector into this subcommand's arguments, refusing any other subcommand.
fn args_of(argv: &[&str]) -> CallFromPspsArgs {
    match Cli::parse_from(argv).cmd {
        PopVarCallerCommand::CallFromPsps(args) => args,
        other => panic!("expected call-from-psps, got {other:?}"),
    }
}

fn refusal_of(argv: &[&str]) -> clap::Error {
    Cli::try_parse_from(argv).expect_err("this argument vector must be refused")
}

/// The shortest run a person can type.
fn a_defaults_run() -> Vec<&'static str> {
    vec![
        "pvc",
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
        Cli::try_parse_from(["pvc", SUBCOMMAND, "--help"])
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
        "pvc",
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
        "pvc",
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
        max_str_bundle_len: DEFAULT_MAX_STR_BUNDLE_LEN,
        min_purity: DEFAULT_MIN_PURITY,
        max_depth: crate::locus_generation::pileup::DEFAULT_MAX_READS_PER_POSITION,
    })
    .expect("the cohort walks into psps");

    let args = CallFromPspsArgs {
        reference: cohort.reference.clone(),
        catalog: Some(cohort.catalog.clone()),
        psps: vec![psp_path_for(&psps, "zeta"), psp_path_for(&psps, "alpha")],
        output: cohort.directory.path().join("calls.vcf"),
        // Several tests call this cohort again into the same directory.
        force: true,
        log: None,
        tmp_dir: None,
        parameters: None,
        defaults: true,
        ploidy: None,
        max_locus_span: DEFAULT_MAX_COHORT_LOCUS_SPAN,
        max_generic_alleles: DEFAULT_MAX_CANDIDATE_ALLELES.get(),
        hidden_dup_fdr: 0.0,
        min_site_quality: crate::run::records::DEFAULT_MIN_SITE_QUALITY,
        hidden_dup_tag: false,
        locus_builder_regions_len: None,
        psp_prefetch_bytes: crate::run::psp_prefetch::DEFAULT_PSP_PREFETCH_BUDGET_BYTES,
        threads: 0,
        max_reads_per_position: None,
        max_depth: crate::run::depth_ceiling::DEFAULT_MAX_READ_GROUP_DEPTH,
        explain_loci: None,
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
fn the_command_writes_a_vcf_and_logs_its_parameters() {
    let (_cohort, args) = a_cohort_of_psps();

    run_call_from_psps(&args).expect("a cohort of stored samples calls");

    assert!(args.output.is_file(), "the VCF is at {:?}", args.output);
    let parameters =
        crate::cli::calling_run::parameters_the_run_logged(&args.output, args.log.as_deref());
    assert!(
        parameters.is_file(),
        "the run's log holds its parameters, copied out to {parameters:?}",
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
/// (`run::report`) are what pin the lines. Run it deliberately:
///
/// ```text
/// cargo test --lib call_from_psps -- --ignored --nocapture the_report_a_person_sees
/// ```
#[test]
#[ignore = "prints the run report for a person to read; the assertions are in run::report"]
fn the_report_a_person_sees() {
    let (_cohort, args) = a_cohort_of_psps();
    run_call_from_psps(&args).expect("a cohort of stored samples calls");
}

/// **A cohort whose psps disagree with each other about the repeat criteria is refused by the
/// command**, and not only by the library (`psp_census_pair.md` §6).
///
/// **The owner's ruling of 2026-09-09 is that this is a hard failure for every command that
/// opens a cohort**, calling included: two psps typed under different criteria hold records that
/// are not about the same loci, so scoring them together compares evidence gathered two ways.
/// The check lives in the cohort opener, which is why calling gets it without a line of its own —
/// and this test is what says calling actually goes through that opener.
///
/// The two psps come from two walks of the same reads under two purity floors. **On this
/// fixture's reference no tract is typed differently by the change** — it is one base repeated,
/// so every tract on it is perfectly pure and clears either floor. What differs is the criteria
/// record in the psp's header, which is what the opener compares, and that is the whole of what
/// this test is about.
#[test]
fn a_cohort_of_psps_that_disagree_about_the_criteria_is_refused_by_the_command() {
    let (cohort, mut args) = a_cohort_of_psps();
    let strictly_typed = cohort.directory.path().join("psps-strict");
    run_generate_psps(&GeneratePspsArgs {
        reference: cohort.reference.clone(),
        catalog: Some(cohort.catalog.clone()),
        // **`alpha` alone**, since that is the only psp of this walk the cohort below uses.
        alignments: vec![cohort.alignments[1].clone()],
        output_dir: strictly_typed.clone(),
        regions: None,
        force: false,
        build_index_if_missing: false,
        min_copies: MinCopies::default(),
        min_period: DEFAULT_MIN_PERIOD,
        max_period: DEFAULT_MAX_PERIOD,
        max_str_len: DEFAULT_MAX_STR_LEN,
        max_str_bundle_len: DEFAULT_MAX_STR_BUNDLE_LEN,
        min_purity: 0.99,
        max_depth: crate::locus_generation::pileup::DEFAULT_MAX_READS_PER_POSITION,
    })
    .expect("the cohort walks into psps a second time");
    // zeta from the first walk, alpha from the second: the disagreement is between the files.
    args.psps = vec![args.psps[0].clone(), psp_path_for(&strictly_typed, "alpha")];

    let refused = run_call_from_psps(&args).expect_err("two typings, one cohort");

    // **The whole sentence, not two words of it.** Without the opener's check this run reaches
    // the caller's own, which refuses the same cohort against *the run's* segmentation and names
    // one sample and the same field — so an assertion on the field alone cannot tell the two
    // refusals apart, and it is the pair of samples that says the check happened at the door.
    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(
        rendered.contains(
            "samples zeta and alpha do not agree on the set of repeat-tract \
                           criteria"
        ),
        "the refusal names the pair and the field, and got: {rendered}",
    );
}

/// **A cohort walked under its own repeat settings is called under them, with nothing to type**
/// (`doc/devel/ng/spec/command_line.md` §4.4). Before the settings were read from the psp
/// headers, this run — psps written at a purity of 0.99, called at the default 0.8 — was refused
/// for routing differently from its walk.
#[test]
fn a_cohort_is_called_under_the_repeat_settings_its_psps_were_written_with() {
    let cohort = a_cohort_on_disk();
    let psps = cohort.directory.path().join("psps-pure");
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
        max_str_bundle_len: DEFAULT_MAX_STR_BUNDLE_LEN,
        min_purity: 0.99,
        max_depth: crate::locus_generation::pileup::DEFAULT_MAX_READS_PER_POSITION,
    })
    .expect("the cohort walks into psps at a purity of 0.99");

    let (_defaults_cohort, args) = a_cohort_of_psps();
    let args = CallFromPspsArgs {
        catalog: Some(cohort.catalog.clone()),
        reference: cohort.reference.clone(),
        psps: vec![psp_path_for(&psps, "zeta"), psp_path_for(&psps, "alpha")],
        output: cohort.directory.path().join("pure.vcf"),
        log: Some(cohort.directory.path().join("pure.log")),
        tmp_dir: None,
        ..args
    };
    run_call_from_psps(&args).expect("the psps' own settings are the run's");
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

/// **A `--ploidy` that disagrees with the parameters file is refused** — a file's numbers were
/// fitted at one ploidy and mean nothing at another, so a run that wants the other number wants
/// a different fit. Provoked through the command, because the flag has to reach the check.
#[test]
fn a_ploidy_that_is_not_the_parameters_files_is_refused_by_the_command() {
    let (cohort, mut args) = a_cohort_of_psps();
    run_call_from_psps(&args).expect("a defaults run writes the parameters it scored with");

    // The file the first run wrote, moved out of the way so the second run may read it.
    let supplied = cohort.directory.path().join("supplied.parameters.toml");
    std::fs::rename(
        crate::cli::calling_run::parameters_the_run_logged(&args.output, args.log.as_deref()),
        &supplied,
    )
    .expect("the file moves");
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
    std::fs::rename(
        crate::cli::calling_run::parameters_the_run_logged(&args.output, args.log.as_deref()),
        &supplied,
    )
    .expect("the file moves");

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
    std::fs::rename(
        crate::cli::calling_run::parameters_the_run_logged(&args.output, args.log.as_deref()),
        &supplied,
    )
    .expect("the file moves");
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

/// **A run asking for the filter now runs it**, where it used to be refused.
///
/// The step that built the scoring and writing passes is the step that deletes the refusal, and
/// this is what says the deletion was not merely a deletion: the run finishes, writes its VCF, and
/// the header states what the filter was calibrated to — which is the one trace a dropped record
/// leaves.
#[test]
fn a_run_asking_for_the_filter_writes_a_calibrated_file() {
    let (cohort, mut args) = a_cohort_of_psps();
    args.hidden_dup_fdr = 0.01;

    run_call_from_psps(&args).expect("the filter finishes a run now");

    let written = std::fs::read_to_string(&args.output).expect("the calls are on disk");
    assert!(
        written.contains("##paralogFilter=target_fdr=0.0100;"),
        "the header must say what the filter was calibrated to: {}",
        written.lines().take(12).collect::<Vec<_>>().join("\n")
    );
    assert!(
        written.contains("##FILTER=<ID=hiddenDup,"),
        "and it must declare the filter it may have applied"
    );
    assert!(
        !std::path::Path::new(&format!("{}.paralog-spill.tmp", args.output.display())).exists(),
        "the spill must not outlive the run"
    );
    drop(cohort);
}

/// **A run that says nothing about the filter runs it, at spec §3.6's target.**
///
/// The default was `0` while the filter's scoring and writing passes did not exist. Every fixture
/// in this file sets the target explicitly, so nothing here would have noticed the constant going
/// back to `0.01` — or failing to. This reads it off the parsed command line, which is where an
/// operator's run gets it.
#[test]
fn a_run_that_says_nothing_about_the_filter_takes_the_specs_target() {
    let args = args_of(&a_defaults_run());

    assert!(
        (args.hidden_dup_fdr - 0.01).abs() < 1e-12,
        "the filter is on by default, at about one wrongly removed record in a hundred; got {}",
        args.hidden_dup_fdr
    );
    assert!(
        !args.hidden_dup_tag,
        "and it removes them rather than tagging them"
    );
}

/// **A target that is not a false-discovery rate is refused, whatever shape it takes.**
///
/// `clap` parses any `f64` into this flag, so a run can ask for `7`, for infinity, or for a
/// negative zero that compares equal to the off value while not being it. A false-discovery rate
/// is a probability; anything else is a mistake worth naming rather than clamping.
#[test]
fn a_paralog_target_that_is_not_a_fraction_is_refused() {
    let (cohort, mut args) = a_cohort_of_psps();

    for asked in [7.0, 1.0, -0.5, -0.0, f64::INFINITY, f64::NAN] {
        args.hidden_dup_fdr = asked;

        let refused = run_call_from_psps(&args).expect_err("not a false-discovery rate");

        assert!(
            matches!(
                refused,
                CallFromPspsCliError::ParalogTargetIsNotAFraction { .. }
            ),
            "expected {asked} to be refused as not a fraction, got {refused:?}"
        );
        assert!(
            !args.output.exists(),
            "the run was refused after opening the output"
        );
    }
    drop(cohort);
}

/// **A calling run writes back the SNP/indel fit's starts that the file it read carried** — the
/// whole command, from a supplied file on disk to the file it writes beside its VCF. A file read and
/// written again without them would drop the record of how the fit behind its numbers ended.
#[test]
fn a_run_writes_back_the_starts_its_file_carried() {
    use crate::calling::parameters_file::{ParametersFile, SnpIndelFitStart, StartOutcome};
    let (cohort, mut args) = a_cohort_of_psps();
    run_call_from_psps(&args).expect("a defaults run writes the parameters it scored with");
    let supplied = cohort.directory.path().join("supplied.parameters.toml");
    let mut file = ParametersFile::from_toml(
        &std::fs::read_to_string(crate::cli::calling_run::parameters_the_run_logged(
            &args.output,
            args.log.as_deref(),
        ))
        .expect("the run wrote its file"),
    )
    .expect("and it parses");
    assert_eq!(
        file.fitted_from.snp_indel_fit_starts, None,
        "no fit stands behind a defaults run"
    );
    let starts = vec![
        SnpIndelFitStart {
            start: 1,
            ended: StartOutcome::Converged,
            passes: 9,
            agreed_with_start: None,
        },
        SnpIndelFitStart {
            start: 2,
            ended: StartOutcome::AgreedWithAnEarlierStart,
            passes: 3,
            agreed_with_start: Some(1),
        },
    ];
    file.fitted_from.snp_indel_fit_starts = Some(starts.clone());
    std::fs::write(&supplied, file.to_toml()).expect("the supplied file is written");
    std::fs::remove_file(crate::cli::calling_run::parameters_the_run_logged(
        &args.output,
        args.log.as_deref(),
    ))
    .expect("the first run's file is removed");
    args.parameters = Some(supplied);
    args.defaults = false;
    run_call_from_psps(&args).expect("a run handed its own kind of file binds it");
    let written = ParametersFile::from_toml(
        &std::fs::read_to_string(crate::cli::calling_run::parameters_the_run_logged(
            &args.output,
            args.log.as_deref(),
        ))
        .expect("the run wrote its file"),
    )
    .expect("and it parses");
    assert_eq!(written.fitted_from.snp_indel_fit_starts, Some(starts));
}

/// **The parameters block a calling run writes into its log carries both informational tables the
/// file it read held** — `[ordinary_site_prior.fitted_density]` and `[duplicated_positions]`
/// (`parameter_prepass_joint_fit.md` §2.1.3, traps). Since a calling run records its parameters
/// in its log rather than beside its VCF, that block is the file a person copies out to re-run;
/// a block without the tables would silently strip them from every file passed on. The whole
/// command, from a supplied file on disk to the block read back out of the log.
#[test]
fn a_runs_logged_parameters_carry_the_informational_tables_its_file_held() {
    use crate::calling::parameters_file::{
        FittedDensity, FittedDuplicatedPositions, InformationalNumber, PairedNumber, ParametersFile,
    };
    let (cohort, mut args) = a_cohort_of_psps();
    run_call_from_psps(&args).expect("a defaults run writes the parameters it scored with");
    let mut file = ParametersFile::from_toml(
        &std::fs::read_to_string(crate::cli::calling_run::parameters_the_run_logged(
            &args.output,
            args.log.as_deref(),
        ))
        .expect("the run logged its parameters"),
    )
    .expect("and they parse");
    assert_eq!(file.ordinary_site_prior.fitted_density, None);
    assert_eq!(file.duplicated_positions, None);
    let number = |value: f64| InformationalNumber {
        value,
        standard_error: Some(value / 10.0),
    };
    let marked = |value: f64| PairedNumber {
        value,
        standard_error: None,
        separately_determined: false,
    };
    let density = FittedDensity {
        stated_frequency: 0.01,
        invariant_share: marked(0.84),
        fixed_non_reference_share: number(3e-5),
        segregating_shape_a: marked(0.02),
        segregating_shape_b: number(4.28),
        share_segregating_above_stated_frequency: number(8.4e-3),
        expected_heterozygosity: number(1.2e-3),
        expected_alternative_frequency: number(7.7e-4),
    };
    let duplicated = FittedDuplicatedPositions {
        stated_frequency: 0.01,
        duplicated_share: marked(0.015),
        carrier_shape_a: marked(0.02),
        carrier_shape_b: number(1.41),
        share_with_carriers_above_stated_frequency: number(1.2e-3),
    };
    file.ordinary_site_prior.fitted_density = Some(density.clone());
    file.duplicated_positions = Some(duplicated.clone());
    let supplied = cohort.directory.path().join("supplied.parameters.toml");
    std::fs::write(&supplied, file.to_toml()).expect("the supplied file is written");
    args.parameters = Some(supplied);
    args.defaults = false;
    args.output = cohort.directory.path().join("with_tables.vcf");
    args.log = None;
    run_call_from_psps(&args).expect("a run handed a file with both tables binds it");
    let logged = std::fs::read_to_string(crate::cli::calling_run::parameters_the_run_logged(
        &args.output,
        args.log.as_deref(),
    ))
    .expect("the run logged its parameters");
    let written = ParametersFile::from_toml(&logged).expect("and they parse");
    assert_eq!(written.ordinary_site_prior.fitted_density, Some(density));
    assert_eq!(written.duplicated_positions, Some(duplicated));
}

/// **`--explain-loci` explains and changes nothing that is called** (`spec/explain_loci.md`
/// §2): the VCF is the same file with the option as without it, and the explanation is a TSV
/// with the stated columns whose every row is a locus of the BED's.
#[test]
fn explaining_loci_writes_the_explanation_and_the_same_vcf() {
    let (cohort, plain) = a_cohort_of_psps();
    run_call_from_psps(&plain).expect("a cohort of stored samples calls");
    let plain_vcf = std::fs::read_to_string(&plain.output).expect("the VCF reads");

    let directory = cohort.directory.path();
    let bed = directory.join("explain.bed");
    let tsv = directory.join("explain.tsv");
    std::fs::write(&bed, "chr1\t0\t100\n").expect("the BED is written");
    let explained = CallFromPspsArgs {
        output: directory.join("explained.vcf"),
        explain_loci: Some(vec![bed, tsv.clone()]),
        ..plain.clone()
    };
    run_call_from_psps(&explained).expect("the same cohort calls while explaining");
    let explained_vcf = std::fs::read_to_string(&explained.output).expect("the VCF reads");

    let comparable = |vcf: &str| -> Vec<String> {
        vcf.lines()
            .filter(|line| !line.starts_with("##commandline") && !line.starts_with("##pvcLog"))
            .map(str::to_owned)
            .collect()
    };
    assert_eq!(comparable(&explained_vcf), comparable(&plain_vcf));

    let text = std::fs::read_to_string(&tsv).expect("the explanation reads");
    assert!(text.starts_with(&format!(
        "#explain_loci_columns_version={}\n",
        crate::run::explain::COLUMNS_VERSION
    )));
    let mut lines = text.lines().filter(|line| !line.starts_with('#'));
    assert_eq!(
        lines.next(),
        Some(crate::run::explain::COLUMNS.join("\t").as_str())
    );
    let rows: Vec<&str> = lines.collect();
    assert!(
        rows.iter()
            .any(|row| row.split('\t').nth(3) == Some("outcome")),
        "every explained locus ends in an outcome row: {rows:?}"
    );
    for row in &rows {
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(fields.len(), crate::run::explain::COLUMNS.len(), "{row}");
        assert_eq!(fields[0], "chr1", "{row}");
    }
}
