//! The command surface: what a person may type, what is refused before a byte is read, and
//! what the refusals tell them to do next.

use super::*;
use clap::Parser;

use crate::pop_var_caller_exp::cli::{Cli, PopVarCallerExpCommand};

/// Parse an argument vector into this subcommand's arguments, refusing any other subcommand.
fn args_of(argv: &[&str]) -> CallFromAlignmentsArgs {
    match Cli::parse_from(argv).cmd {
        PopVarCallerExpCommand::CallFromAlignments(args) => args,
        other => panic!("expected call-from-alignments, got {other:?}"),
    }
}

/// Try to parse, and hand back clap's own refusal.
fn refusal_of(argv: &[&str]) -> clap::Error {
    Cli::try_parse_from(argv).expect_err("these arguments are refused")
}

/// The shortest run a person can type: a reference, one alignment file, somewhere to write, and
/// the defaults.
fn a_defaults_run() -> Vec<&'static str> {
    vec![
        "pop_var_caller_exp",
        "call-from-alignments",
        "--reference",
        "ref.fa",
        "--alignment",
        "one.cram",
        "--output",
        "calls.vcf.gz",
        "--defaults",
    ]
}

/// **The subcommand is spelled `call-from-alignments`** — kebab-cased from its enum variant,
/// like the three beside it. The name is the command surface's, agreed for all four modes
/// (`doc/devel/ng/spec/typed_regions_cli.md`), so it is pinned rather than left to clap's
/// derive.
#[test]
fn the_subcommand_is_spelled_call_from_alignments() {
    let args = args_of(&a_defaults_run());

    assert_eq!(args.reference, PathBuf::from("ref.fa"));
    assert_eq!(args.output, PathBuf::from("calls.vcf.gz"));
    assert_eq!(args.alignments, vec![PathBuf::from("one.cram")]);
    assert!(args.defaults);
    assert!(args.parameters.is_none());
}

/// **A run has to say where its numbers come from**, and neither answer is the default: a run
/// that silently guessed its parameters would produce a file nothing on it says was guessed.
#[test]
fn a_run_that_names_neither_a_parameters_file_nor_the_defaults_is_refused() {
    let refusal = refusal_of(&[
        "pop_var_caller_exp",
        "call-from-alignments",
        "--reference",
        "ref.fa",
        "--alignment",
        "one.cram",
        "--output",
        "calls.vcf.gz",
    ]);

    assert!(
        refusal.to_string().contains("--defaults"),
        "the refusal names the flag that is missing, and got: {refusal}",
    );
}

/// **A run cannot both supply its numbers and ask for the defaults**, since the two are
/// different claims about what the genotypes rest on.
#[test]
fn a_run_that_names_both_a_parameters_file_and_the_defaults_is_refused() {
    let mut argv = a_defaults_run();
    argv.extend(["--parameters", "fitted.parameters.toml"]);

    let refusal = refusal_of(&argv);

    assert!(
        refusal.to_string().contains("--parameters"),
        "the refusal names the two that conflict, and got: {refusal}",
    );
}

/// **Every sample is one `--alignment`**, repeated, and a run with none is refused before a
/// reference is read.
#[test]
fn the_alignment_flag_repeats_and_keeps_the_order_it_was_given() {
    let args = args_of(&[
        "pop_var_caller_exp",
        "call-from-alignments",
        "--reference",
        "ref.fa",
        "--alignment",
        "zeta.cram",
        "--alignment",
        "alpha.cram",
        "--output",
        "calls.vcf",
        "--defaults",
    ]);

    assert_eq!(
        args.alignments,
        vec![PathBuf::from("zeta.cram"), PathBuf::from("alpha.cram")],
        "the run's sample order is the order the files were named in",
    );
    assert!(
        refusal_of(&[
            "pop_var_caller_exp",
            "call-from-alignments",
            "--reference",
            "ref.fa",
            "--output",
            "calls.vcf",
            "--defaults",
        ])
        .to_string()
        .contains("--alignment"),
    );
}

/// The ground defaults to every base of every contig, and the catalog to the file
/// `repeat-catalog` writes beside the reference.
#[test]
fn the_ground_and_the_catalog_have_defaults_a_person_need_not_type() {
    let args = args_of(&a_defaults_run());

    assert!(args.regions.is_none(), "no BED means the whole genome");
    assert!(args.catalog.is_none());
    assert_eq!(
        args.ploidy, None,
        "a ploidy nobody typed is absent, not two — which is what lets a supplied file's own \
         ploidy be taken without the flag's default contradicting it",
    );
    assert_eq!(
        ploidy_asked_for(&args).expect("two is a ploidy").get(),
        DEFAULT_PLOIDY,
    );
}

/// **A ploidy past what the read likelihood scores is refused before anything is read.**
///
/// It was a panic until 2026-09-01, and one reached after the whole cohort had been opened:
/// `Ploidy::try_new` turns down zero and nothing else, and the read likelihood's copy-share
/// table asserts at seventeen. A polyploid crop is an ordinary thing to call — sugarcane runs to
/// about twelve copies — so a person who types twenty gets a sentence naming the ceiling.
#[test]
fn a_ploidy_past_what_the_caller_scores_is_refused_and_names_the_ceiling() {
    let mut argv = a_defaults_run();
    argv.extend(["--ploidy", "20"]);
    let args = args_of(&argv);

    let refused = ploidy_asked_for(&args).expect_err("twenty copies cannot be scored");

    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(
        rendered.contains("16"),
        "the message names the ceiling, and got: {rendered}",
    );
    assert!(
        ploidy_asked_for(&args_of(&{
            let mut argv = a_defaults_run();
            argv.extend(["--ploidy", "16"]);
            argv
        }))
        .is_ok(),
        "sixteen is the ceiling and is scored, so the refusal is not off by one",
    );
}

/// A ploidy of zero is not a genome, and is refused with the type's own words.
#[test]
fn a_ploidy_of_zero_is_refused() {
    let mut argv = a_defaults_run();
    argv.extend(["--ploidy", "0"]);

    assert!(ploidy_asked_for(&args_of(&argv)).is_err());
}

/// **`--output` naming a directory is refused before the reference is read.**
///
/// Left to the writer it is discovered after the last locus has been called — a whole cohort's
/// work — and leaves the in-flight `<output>.tmp` beside the directory the person named.
#[test]
fn an_output_that_names_a_directory_is_refused_before_anything_is_read() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    let args = a_run_writing_to(directory.path().to_path_buf());

    let refused = refuse_an_output_that_cannot_be_written(&args)
        .expect_err("a directory is not somewhere to write a VCF");

    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(rendered.contains("is a directory"), "and got: {rendered}",);
}

/// **An `--output` in a directory that does not exist is refused before the reference is read**,
/// which is minutes of reading and opening a person does not have to wait through for a typo.
#[test]
fn an_output_in_a_missing_directory_is_refused_before_anything_is_read() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    let args = a_run_writing_to(directory.path().join("no-such-directory").join("calls.vcf"));

    let refused = refuse_an_output_that_cannot_be_written(&args)
        .expect_err("there is no directory to write into");

    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(
        rendered.contains("no-such-directory"),
        "the message names the directory that is missing, and got: {rendered}",
    );
}

/// **An output that can be written is admitted**, including a bare file name in the working
/// directory — whose parent is the empty path and must not be read as a missing directory.
#[test]
fn an_output_that_can_be_written_is_admitted() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    assert!(
        refuse_an_output_that_cannot_be_written(&a_run_writing_to(
            directory.path().join("calls.vcf.gz")
        ))
        .is_ok()
    );
    assert!(
        refuse_an_output_that_cannot_be_written(&a_run_writing_to(PathBuf::from("calls.vcf")))
            .is_ok(),
        "a bare file name is the working directory, not a missing one",
    );
}

/// A run's arguments with everything but the output fixed, so the output checks read as one
/// question.
fn a_run_writing_to(output: PathBuf) -> CallFromAlignmentsArgs {
    CallFromAlignmentsArgs {
        reference: PathBuf::from("ref.fa"),
        catalog: None,
        alignments: vec![PathBuf::from("one.cram")],
        output,
        regions: None,
        parameters: None,
        defaults: true,
        ploidy: None,
        build_index_if_missing: false,
    }
}

/// **A run with no catalog is told which file is missing and the command that builds it**, and
/// is refused before a single alignment file is opened.
///
/// The catalog is what says where the repeat tracts are; without it a run has no way to route
/// its ground, and a message naming only *not found* would leave an operator to work out that
/// `repeat-catalog` is what makes one.
#[test]
fn a_run_with_no_catalog_is_told_which_file_is_missing_and_how_to_build_it() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    let reference = directory.path().join("ref.fa");
    let args = CallFromAlignmentsArgs {
        reference: reference.clone(),
        catalog: None,
        alignments: vec![directory.path().join("one.cram")],
        output: directory.path().join("calls.vcf"),
        regions: None,
        parameters: None,
        defaults: true,
        ploidy: None,
        build_index_if_missing: false,
    };

    let refused = segments_over(
        &args,
        &GenomeRegions::whole_contigs(&[]),
        &ReferenceInfo {
            md5: None,
            contigs: Vec::new(),
            fasta_path: Some(reference),
        },
    )
    .expect_err("there is no catalog beside this reference");

    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(
        rendered.contains("ref.fa.repeats.parquet"),
        "the message names the file that is missing, and got: {rendered}",
    );
    assert!(
        rendered.contains("repeat-catalog"),
        "and the command that builds it, and got: {rendered}",
    );
}
