//! The command surface: what a person may type, what is refused before a byte is read, and
//! what the refusals tell them to do next.

use super::*;
use clap::Parser;
use std::path::Path;

use crate::ng::repeat_catalog::StrRepeatCriteria;
use crate::ng::types::InbreedingF;

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
        max_cohort_locus_span: DEFAULT_MAX_COHORT_LOCUS_SPAN,
        max_candidate_alleles: DEFAULT_MAX_CANDIDATE_ALLELES.get(),
        cohort_locus_builder_regions_len: None,
        threads: 0,
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
        max_cohort_locus_span: DEFAULT_MAX_COHORT_LOCUS_SPAN,
        max_candidate_alleles: DEFAULT_MAX_CANDIDATE_ALLELES.get(),
        cohort_locus_builder_regions_len: None,
        threads: 0,
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

// ---------------------------------------------------------------------
// The parameters the run wrote beside its VCF (spec §7)
// ---------------------------------------------------------------------

/// **The header names the file the run writes beside its output**, by name and not by path.
///
/// A reader holding the VCF has to be able to find what its genotypes rest on, and the two are
/// siblings by construction — so the name is enough, and a path would be wrong the moment
/// somebody moved the pair.
#[test]
fn the_headers_parameters_file_is_the_one_the_run_writes_beside_the_vcf() {
    for (output, expected) in [
        ("/data/run7/calls.vcf.gz", "calls.parameters.toml"),
        ("/data/run7/calls.vcf", "calls.parameters.toml"),
        ("tomato.bcf", "tomato.parameters.toml"),
    ] {
        assert_eq!(
            beside_the_vcf(Path::new(output))
                .file_name()
                .expect("a file name")
                .to_string_lossy(),
            expected,
            "from {output}",
        );
    }
}

/// **A run may not write its parameters over the file it was handed**, and the collision is the
/// one spec §7 invites: it tells a user to copy the file their run wrote and change a line, so
/// `--parameters calls.parameters.toml --output calls.vcf.gz` is the natural next command.
#[test]
fn a_run_whose_output_would_overwrite_its_own_parameters_file_is_refused() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    let supplied = directory.path().join("calls.parameters.toml");
    std::fs::write(&supplied, "# a file somebody edited\n").expect("a file to be handed");
    let mut args = a_run_writing_to(directory.path().join("calls.vcf.gz"));
    args.defaults = false;
    args.parameters = Some(supplied);

    let refused = refuse_an_output_whose_parameters_file_is_this_run_s_input(&args)
        .expect_err("the run would write over its own input");

    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(
        rendered.contains("calls.parameters.toml") && rendered.contains("calls.vcf.gz"),
        "the message names both files, and got: {rendered}",
    );
}

/// **Two spellings of one path are one file**, so a relative name does not slip past the
/// refusal.
#[test]
fn the_refusal_sees_through_two_spellings_of_one_path() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    std::fs::write(directory.path().join("calls.parameters.toml"), "# edited\n")
        .expect("a file to be handed");
    std::fs::create_dir(directory.path().join("sub")).expect("a subdirectory to walk back out of");
    let mut args = a_run_writing_to(directory.path().join("calls.vcf.gz"));
    args.defaults = false;
    args.parameters = Some(
        directory
            .path()
            .join(".")
            .join("sub")
            .join("..")
            .join("calls.parameters.toml"),
    );

    assert!(
        refuse_an_output_whose_parameters_file_is_this_run_s_input(&args).is_err(),
        "`<dir>/sub/../calls.parameters.toml` is `<dir>/calls.parameters.toml`",
    );
}

/// **A symlink pointing at the file the run would write is that file**, and following it would
/// have destroyed the target through the link — the same loss, one indirection away.
#[cfg(unix)]
#[test]
fn the_refusal_follows_a_symlink_to_the_file_the_run_would_write() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    let target = directory.path().join("calls.parameters.toml");
    std::fs::write(&target, "# edited by hand\n").expect("a file to be handed");
    let handy = directory.path().join("handy.toml");
    std::os::unix::fs::symlink(&target, &handy).expect("a symlink to it");
    let mut args = a_run_writing_to(directory.path().join("calls.vcf.gz"));
    args.defaults = false;
    args.parameters = Some(handy);

    assert!(
        refuse_an_output_whose_parameters_file_is_this_run_s_input(&args).is_err(),
        "a link to the destination is the destination",
    );
}

/// **A `--parameters` file that is not there is not this refusal's business.** The person
/// mistyped a name, and telling them to copy a file that does not exist sends them nowhere; the
/// message they need comes from the read that follows.
#[test]
fn a_parameters_file_that_does_not_exist_is_left_to_the_read_that_follows() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    let mut args = a_run_writing_to(directory.path().join("calls.vcf.gz"));
    args.defaults = false;
    args.parameters = Some(directory.path().join("calls.parameters.toml"));

    assert!(
        refuse_an_output_whose_parameters_file_is_this_run_s_input(&args).is_ok(),
        "nothing is there to be overwritten yet",
    );
}

/// A parameters file that is not the run's own output is admitted, and so is a defaults run,
/// which has no input to overwrite.
#[test]
fn a_parameters_file_that_is_not_the_runs_own_output_is_admitted() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    let elsewhere = directory.path().join("fitted.parameters.toml");
    std::fs::write(&elsewhere, "# a fit's own file\n").expect("a file to be handed");
    let mut args = a_run_writing_to(directory.path().join("calls.vcf.gz"));
    args.defaults = false;
    args.parameters = Some(elsewhere);
    assert!(refuse_an_output_whose_parameters_file_is_this_run_s_input(&args).is_ok());

    let defaults = a_run_writing_to(directory.path().join("calls.vcf.gz"));
    assert!(
        refuse_an_output_whose_parameters_file_is_this_run_s_input(&defaults).is_ok(),
        "a defaults run has no supplied file to overwrite",
    );
}

/// Two fixture alignment files, and the run's read-group table over them.
///
/// **A real table rather than a fabricated one**, because every axis of a parameters file is
/// keyed on it — one row a read group, one a sample, in the order this table fixes — and a
/// hand-built one would not exercise the joins `ParametersFile::of_run` holds in release.
fn a_cohorts_read_groups() -> (tempfile::TempDir, tempfile::TempDir, ReadGroups) {
    use crate::ng::read::input::test_fixtures::{header, indexed_named_bam, matching_contigs};

    let with_sample = |sample: &str, file: &str| {
        indexed_named_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some(sample))],
            ),
            &[],
            file,
        )
    };
    let (zeta_dir, zeta) = with_sample("zeta", "zeta.bam");
    let (alpha_dir, alpha) = with_sample("alpha", "alpha.bam");
    let read_groups = build_read_groups(&[zeta, alpha]).expect("the fixtures declare read groups");
    (zeta_dir, alpha_dir, read_groups)
}

/// **A run is reproducible from its own output** — spec §7's first purpose, as an assertion.
///
/// The file a defaults run writes beside its VCF reads back into the same numbers that run
/// scored with: the same ploidy, the same count of libraries and plants, and the same
/// coefficient for each of them. It is bound to this run at its own door on the way back in —
/// a file naming another reference or another cohort's read groups is refused there — so what
/// this shows is that a run's own file passes its own bindings and projects to its own numbers.
#[test]
fn what_a_defaults_run_writes_reads_back_into_the_numbers_it_scored_with() {
    let (_zeta_dir, _alpha_dir, read_groups) = a_cohorts_read_groups();
    let ploidy = Ploidy::try_new(2).expect("a diploid");
    let reference = ReferenceDigest([7; 16]);
    // **The two plants are given different coefficients, and that is what makes the round trip
    // mean something.** With both at the default the per-sample axis compares equal under any
    // permutation, so a writer or reader that reversed it would pass — measured, as a surviving
    // mutation. `sampleA` selfs and `sampleB` outcrosses, so the two rows are distinguishable
    // and their order is checked rather than assumed.
    let inbreeding = DeclaredInbreeding::nothing_said().and_this_sample(
        "zeta",
        InbreedingF::try_new(0.9).expect("a coefficient in [0, 1)"),
    );
    let scored_with = RunParameters::of_defaults(&read_groups, ploidy, &inbreeding);
    assert_ne!(
        scored_with.inbreeding_coefficient_by_sample()[0],
        scored_with.inbreeding_coefficient_by_sample()[1],
        "a fixture whose two samples score alike could not tell the axis's order",
    );

    let file = ParametersFile::of_run(
        &scored_with,
        &read_groups,
        &ReadsBehindEachCalibration::nothing_was_fitted(read_groups.len()),
        &inbreeding.of_each_sample(&read_groups),
        &reference,
        CensusIdentity::of_a_run_with_no_census(),
    );
    let directory = tempfile::tempdir().expect("a temporary directory");
    let vcf = directory.path().join("calls.vcf.gz");
    let at = file.write_beside_the_vcf(&vcf).expect("the write succeeds");

    assert_eq!(
        at,
        directory.path().join("calls.parameters.toml"),
        "beside the VCF, named after it",
    );
    let text = std::fs::read_to_string(&at).expect("the file is there");
    let read_back = ParametersFile::from_toml(&text)
        .expect("what a run wrote is what its reader reads")
        .to_run_parameters_for(&reference, &read_groups, None)
        .expect("and it binds to the run that wrote it")
        .from_file
        .parameters;

    assert_eq!(read_back.ploidy(), scored_with.ploidy());
    assert_eq!(read_back.read_group_count(), scored_with.read_group_count(),);
    assert_eq!(
        read_back.inbreeding_coefficient_by_sample(),
        scored_with.inbreeding_coefficient_by_sample(),
        "one coefficient a plant, in the run's own sample order — and the two differ, so this \
         compares the order and not only the multiset",
    );
    assert_eq!(
        read_back
            .calibration_by_read_group()
            .iter()
            .map(|calibration| (calibration.scale, calibration.provenance))
            .collect::<Vec<_>>(),
        scored_with
            .calibration_by_read_group()
            .iter()
            .map(|calibration| (calibration.scale, calibration.provenance))
            .collect::<Vec<_>>(),
        "one multiplier a library",
    );
}

/// **A defaults run's file says it fitted nothing**, which is the whole of why spec §7 makes
/// writing unconditional: a run that guessed its numbers is auditable in the same form as one
/// that measured them, and the file's own count is what tells a reader which they are holding.
#[test]
fn a_defaults_runs_file_says_it_fitted_nothing() {
    let (_zeta_dir, _alpha_dir, read_groups) = a_cohorts_read_groups();
    let inbreeding = DeclaredInbreeding::nothing_said();
    let scored_with = RunParameters::of_defaults(
        &read_groups,
        Ploidy::try_new(2).expect("a diploid"),
        &inbreeding,
    );

    let file = ParametersFile::of_run(
        &scored_with,
        &read_groups,
        &ReadsBehindEachCalibration::nothing_was_fitted(read_groups.len()),
        &inbreeding.of_each_sample(&read_groups),
        &ReferenceDigest([7; 16]),
        CensusIdentity::of_a_run_with_no_census(),
    );

    let fitted = file.what_the_run_fitted();
    assert!(fitted.nothing_was_fitted());
    assert_eq!(
        fitted.fitted().len(),
        0,
        "none of the {} groups rests on this cohort's reads",
        fitted.groups(),
    );
    assert!(
        fitted.groups() > 0,
        "a denominator of zero would make the count above say nothing",
    );
}

// ---------------------------------------------------------------------
// The whole command, driven
// ---------------------------------------------------------------------

/// **A reference, its catalog, and two samples' alignment files** — everything
/// [`run_call_from_alignments`] needs, built on disk.
///
/// **The reference is the shared fixture's**, a hundred `A`s on `chr1` and two hundred on
/// `chr2`, which is what the alignment fixtures declare in their `@SQ`. Every base of it is one
/// mononucleotide run, so the catalog routes the whole genome to the repeat-tract generator and
/// the run calls no locus at all. **That is what this fixture is for**: what it exercises is the
/// command's wiring — the files it writes and what they say about each other — and a run that
/// wrote no record still writes a header, a parameters file, and a summary.
fn a_cohort_on_disk() -> (
    tempfile::TempDir,
    tempfile::TempDir,
    tempfile::TempDir,
    CallFromAlignmentsArgs,
) {
    use crate::ng::read::input::test_fixtures::{header, indexed_named_bam, matching_contigs};
    use crate::ng::reference_info::{ReferenceSource, read_reference_info_observing};
    use crate::ng::repeat_catalog::RepeatCatalogBuilder;
    use crate::ng::tandem_repeat::ScanParams;
    use crate::pileup::per_sample::cram_files::{ContigSpec, build_fasta};

    let specs: Vec<ContigSpec> = crate::ng::read::input::test_fixtures::FIXTURE_CONTIGS
        .iter()
        .map(|(name, length)| ContigSpec {
            name: (*name).to_string(),
            length: *length as u64,
        })
        .collect();
    let (reference_dir, fasta) = build_fasta(&specs).expect("a reference on disk");

    let catalog_path = reference_dir.path().join("ref.fa.repeats.parquet");
    let criteria = StrRepeatCriteria::default();
    let mut builder = RepeatCatalogBuilder::create(
        &catalog_path,
        criteria.clone(),
        ScanParams {
            match_reward: 2,
            mismatch_penalty: 7,
            min_copies: 2,
        },
    )
    .expect("a catalog to build into");
    let reference = read_reference_info_observing(
        ReferenceSource::Fasta {
            fasta: fasta.clone(),
            fai: None,
        },
        &mut builder,
    )
    .expect("the reference reads");
    builder.finish(&reference).expect("the catalog is written");

    let with_sample = |sample: &str, file: &str| {
        indexed_named_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some(sample))],
            ),
            &[],
            file,
        )
    };
    let (zeta_dir, zeta) = with_sample("zeta", "zeta.bam");
    let (alpha_dir, alpha) = with_sample("alpha", "alpha.bam");

    let args = CallFromAlignmentsArgs {
        reference: fasta,
        catalog: Some(catalog_path),
        alignments: vec![zeta, alpha],
        output: reference_dir.path().join("calls.vcf"),
        regions: None,
        parameters: None,
        defaults: true,
        ploidy: None,
        build_index_if_missing: false,
        max_cohort_locus_span: DEFAULT_MAX_COHORT_LOCUS_SPAN,
        max_candidate_alleles: DEFAULT_MAX_CANDIDATE_ALLELES.get(),
        cohort_locus_builder_regions_len: None,
        threads: 0,
    };
    (reference_dir, zeta_dir, alpha_dir, args)
}

/// **The command writes both its files, and each says the right thing about the other.**
///
/// **This is the only test that drives `run_call_from_alignments` itself**, and the review that
/// asked for it measured what its absence cost: nine of fourteen mutations aimed at this step
/// survived, because every other test calls the helpers directly. Among them were a header
/// naming `calls.vcf` instead of `calls.parameters.toml` on every VCF a run writes, and two that
/// hand `ParametersFile::of_run` an axis of the wrong length — which panics at startup on any
/// real run and on nothing in the suite.
#[test]
fn the_command_writes_a_vcf_and_the_parameters_beside_it() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();

    run_call_from_alignments(&args).expect("the cohort runs");

    let parameters_at = beside_the_vcf(&args.output);
    assert!(args.output.is_file(), "the VCF is at {:?}", args.output);
    assert!(
        parameters_at.is_file(),
        "and the parameters beside it at {parameters_at:?}",
    );

    let vcf = std::fs::read_to_string(&args.output).expect("the VCF reads");
    assert!(
        vcf.contains("##parametersFile=calls.parameters.toml\n"),
        "the header names the file the run wrote beside it, and got:\n{}",
        vcf.lines()
            .take_while(|line| line.starts_with("##"))
            .collect::<Vec<_>>()
            .join("\n"),
    );

    let text = std::fs::read_to_string(&parameters_at).expect("the parameters read");
    let file = ParametersFile::from_toml(&text).expect("what the run wrote, its reader reads");
    assert_eq!(
        file.fitted_from.samples,
        vec!["zeta".to_string(), "alpha".to_string()],
        "the run's own samples, in the order its alignment files were named",
    );
    assert_eq!(
        file.fitted_from.read_groups.len(),
        2,
        "one row a library, which is what `of_run` joins every other axis against",
    );
    assert_eq!(file.inbreeding.by_sample.len(), 2);
    assert!(
        file.what_the_run_fitted().nothing_was_fitted(),
        "a `--defaults` run fitted nothing, and the file says so",
    );
}

/// **The second command a person types is refused by the driver**, not merely by the helper the
/// other tests call.
///
/// After a first run, `--parameters calls.parameters.toml --output calls.vcf` is the natural next
/// thing to type and would write over the file just read. Deleting the refusal's *call site* left
/// every other test green, which is what this closes.
#[test]
fn the_command_refuses_to_write_its_parameters_over_the_file_it_was_given() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();
    run_call_from_alignments(&args).expect("the first run writes both files");

    let again = CallFromAlignmentsArgs {
        parameters: Some(beside_the_vcf(&args.output)),
        defaults: false,
        ..args.clone()
    };
    let refused = run_call_from_alignments(&again).expect_err("the second run is refused");

    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(
        rendered.contains("calls.parameters.toml"),
        "the message names the file that would be lost, and got: {rendered}",
    );
}

/// **What a person sees when the command finishes**, printed once so a reviewer of this step
/// reads the artefact rather than the code that builds it.
///
/// Ignored, because its value is the eye and not the assertion — the report's own tests
/// (`ng::run::report`) are what pin the lines. Run it deliberately:
///
/// ```text
/// cargo test --lib call_from_alignments -- --ignored --nocapture the_report_a_person_sees
/// ```
#[test]
#[ignore = "prints the run report for a person to read; the assertions are in ng::run::report"]
fn the_report_a_person_sees() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();
    run_call_from_alignments(&args).expect("the cohort runs");
}

// ---------------------------------------------------------------------
// The round width chosen from the cohort's size
// ---------------------------------------------------------------------

/// **What a round costs is `width × samples`, so that is what the rule holds fixed.** A
/// round holds about one observation per covered base per sample, so a width that is right
/// at sixty-three samples holds sixteen times as much at a thousand. Pinning the product
/// rather than the width is what lets one rule serve both ends of the cohort range
/// (`design_principles.md` §0).
///
/// The two clamps are what the product cannot express: below the floor the merge's own
/// default takes over, and above the ceiling there is nothing left to buy — measured on four
/// accessions over 400 kb of SL4.0, 3.29 s at 8,000 bases, 3.18 s at 32,000 and 3.22 s at
/// 64,000, the last costing 407 MB of peak resident against 340.
#[test]
fn the_round_width_holds_one_rounds_observations_to_a_budget() {
    // Between the clamps, the product is the budget.
    for samples in [40_usize, 63, 100, 500] {
        let width = u64::from(round_width_for(samples).get());
        let held = width * samples as u64;
        assert!(
            held <= u64::from(ROUND_OBSERVATION_BUDGET),
            "{samples} samples got {width} bases, holding {held} observations a round",
        );
        assert!(
            held > u64::from(ROUND_OBSERVATION_BUDGET) / 2,
            "{samples} samples got {width} bases, which leaves most of the budget unspent",
        );
    }
}

/// **A cohort big enough to be memory-bound gets the number it has today.** The floor is the
/// merge's own compiled-in default, so nothing this rule does can make a thousand-sample run
/// hold more ground than it held before the rule existed.
#[test]
fn a_large_cohort_gets_the_merges_own_default() {
    assert_eq!(
        round_width_for(2_000).get(),
        DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN,
    );
    assert_eq!(
        round_width_for(1_000).get(),
        DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN,
    );
}

/// **A single sample gets the ceiling, not the budget.** One sample could hold half a million
/// bases of ground within the budget; the ceiling is there because the gain has saturated
/// long before that, and because a round that wide would make the run's memory jump on the
/// one input where it is least expected.
#[test]
fn a_small_cohort_gets_the_ceiling() {
    assert_eq!(round_width_for(1).get(), WIDEST_ROUND);
    assert_eq!(round_width_for(4).get(), WIDEST_ROUND);
    // Zero files cannot reach the command, but the rule must not divide by it.
    assert_eq!(round_width_for(0).get(), WIDEST_ROUND);
}

/// **The benchmark's own cohort gets the width the sweep found.** 63 accessions over the
/// whole 8 Mb of `benchmarks/tomato1/regions.bed` ran in 193.2 s at 500 bases and 115.3 s at
/// 8,000, writing the same VCF; the rule lands on 7,936.
#[test]
fn the_tomato_cohorts_width_is_the_one_that_was_measured() {
    assert_eq!(round_width_for(63).get(), 7_936);
}
