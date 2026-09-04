//! The command surface: what a person may type, what is refused before a byte is read, and
//! what the refusals tell them to do next.

use super::*;
use clap::Parser;
use std::path::Path;

use crate::ng::region_typing::{GenomeRegions, RegionKind, TypedRegion, TypedRegionConfig};
use crate::ng::repeat_catalog::{ReadScope, RepeatCatalog, StrRepeatCriteria};
use crate::ng::types::InbreedingF;
use crate::pop_var_caller_exp::run_ground::{GroundError, routing_criteria, segments_over};
use crate::regions::ContigBounds;

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
        min_copies: MinCopies::default(),
        min_period: DEFAULT_MIN_PERIOD,
        max_period: DEFAULT_MAX_PERIOD,
        max_str_len: DEFAULT_MAX_STR_LEN,
        min_purity: DEFAULT_MIN_PURITY,
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
        min_copies: MinCopies::default(),
        min_period: DEFAULT_MIN_PERIOD,
        max_period: DEFAULT_MAX_PERIOD,
        max_str_len: DEFAULT_MAX_STR_LEN,
        min_purity: DEFAULT_MIN_PURITY,
    };

    let refused = segments_over(
        &super::ground_request(&args),
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
    use crate::ng::read::input::test_fixtures::{
        header, indexed_named_bam, matching_contigs, read_group_for,
    };

    let with_sample = |sample: &str, file: &str| {
        indexed_named_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[(&read_group_for(sample), Some(sample))],
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
        &StrRepeatCriteria::default(),
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
        &StrRepeatCriteria::default(),
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
    use crate::ng::read::input::test_fixtures::{
        header, indexed_named_bam, matching_contigs, read_group_for,
    };
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
                &[(&read_group_for(sample), Some(sample))],
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
        min_copies: MinCopies::default(),
        min_period: DEFAULT_MIN_PERIOD,
        max_period: DEFAULT_MAX_PERIOD,
        max_str_len: DEFAULT_MAX_STR_LEN,
        min_purity: DEFAULT_MIN_PURITY,
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

// ---------------------------------------------------------------------
// What this run counts as a repeat
// ---------------------------------------------------------------------

/// **The defaults are ng's calling floors, not the catalog's storage floors.**
///
/// This is the whole of the routing change. The run asked the catalog with
/// `StrRepeatCriteria::default()`, which *is* what the file was built at, so every row the
/// file held became an STR locus of the run — on the human benchmark about seven times more
/// reference than ng's own floors would route. Parsed through clap rather than read off the
/// struct, so a `default_value` that drifted from the library constant fails here.
#[test]
fn the_default_routing_is_the_calling_floors_and_not_the_catalogs() {
    let asked = routing_criteria(&super::ground_request(&args_of(&a_defaults_run())).routing)
        .expect("the defaults make a range");

    assert_eq!(
        asked,
        StrRepeatCriteria::from(&TypedRegionConfig::default()),
        "the run's default policy is step 3's own, converted for a reader",
    );
    assert_ne!(
        asked,
        StrRepeatCriteria::default(),
        "and it is not the catalog's, which is what the run used to ask with",
    );

    // Stated as the gap itself: at every period the run now needs strictly more copies than
    // the file was written down at, which is the room the file exists to leave.
    let catalog = StrRepeatCriteria::default();
    for period in 1..=6 {
        assert!(
            asked.classification.min_copies.for_period(period)
                > catalog.classification.min_copies.for_period(period),
            "period {period}: the run asks for {} copies, the catalog holds from {}",
            asked.classification.min_copies.for_period(period),
            catalog.classification.min_copies.for_period(period),
        );
    }
    assert!(
        asked.max_str_len_bp.get() < catalog.max_str_len_bp.get(),
        "and the satellite cap is the calling one, {} bp against the file's {}",
        asked.max_str_len_bp.get(),
        catalog.max_str_len_bp.get(),
    );
}

/// **Every one of the five flags reaches the criteria the catalog is asked with.**
///
/// A flag that parsed and then went nowhere would route the whole genome on the defaults with
/// nothing crashing and every number in the report plausible, which is why each is moved to a
/// value nothing else in the run could produce and each is read back by itself.
#[test]
fn every_routing_flag_reaches_the_criteria_the_catalog_is_asked_with() {
    let mut argv = a_defaults_run();
    argv.extend([
        "--min-copies",
        "9,7,7,7,6,5",
        "--min-period",
        "2",
        "--max-period",
        "4",
        "--max-str-len",
        "60",
        "--min-purity",
        "0.95",
    ]);
    let asked = routing_criteria(&super::ground_request(&args_of(&argv)).routing)
        .expect("2..=4 is a range");

    let floors: Vec<u32> = (1..=6)
        .map(|period| asked.classification.min_copies.for_period(period))
        .collect();
    assert_eq!(floors, vec![9, 7, 7, 7, 6, 5]);
    assert_eq!(asked.classification.periods.min(), 2);
    assert_eq!(asked.classification.periods.max(), 4);
    assert_eq!(asked.max_str_len_bp.get(), 60);
    assert!((asked.classification.min_purity - 0.95).abs() < 1e-6);

    assert_eq!(
        asked.min_flank_bp.get(),
        StrRepeatCriteria::default().min_flank_bp.get(),
        "the flank floor is not a flag: the rows below the file's were never written",
    );
}

/// **A period range the wrong way round is a message, not a panic.** Clap bounds each end to
/// 1..=6 on its own, so this is the one way left to type a range that is not one.
#[test]
fn a_period_range_the_wrong_way_round_is_refused_before_the_catalog_is_opened() {
    let mut argv = a_defaults_run();
    argv.extend(["--min-period", "5", "--max-period", "3"]);

    let refused = routing_criteria(&super::ground_request(&args_of(&argv)).routing)
        .expect_err("5..=3 is not a range");
    assert!(
        matches!(refused, GroundError::PeriodRange { .. }),
        "got {refused:?}",
    );
}

/// A one-contig reference holding two homopolymers and nothing else a repeat scanner can
/// find, with the catalog a `repeat-catalog` run would have written beside it.
///
/// **The two tracts straddle the gap the catalog exists to leave.** The file is built at 5
/// copies for period 1 and ng calls from 8, so a run of **6** `A`s is in the file and below
/// the calling floor, and a run of **10** is above both. Everything between them and around
/// them is `CGTGCTG` repeated — period 7, outside the 1..=6 the scanner looks for, and
/// carrying no `A` at all, so neither tract can grow into its surroundings.
///
/// Positions, 1-based: filler 1–40, the six-base run 41–46, filler 47–86, the ten-base run
/// 87–96, filler 97–136. Forty bases of filler each side is more than the 15 bp of flank the
/// file requires and more than the 15 bp within which two tracts would be bundled together
/// instead of being loci.
fn the_fixture_references_bases() -> Vec<u8> {
    let filler: Vec<u8> = b"CGTGCTG".iter().copied().cycle().take(40).collect();
    let mut bases = Vec::new();
    bases.extend_from_slice(&filler);
    bases.extend(std::iter::repeat_n(b'A', 6));
    bases.extend_from_slice(&filler);
    bases.extend(std::iter::repeat_n(b'A', 10));
    bases.extend_from_slice(&filler);
    assert_eq!(bases.len(), 136, "the fixture's own geometry");
    bases
}

fn a_reference_with_a_tract_on_each_side_of_the_calling_floor()
-> (tempfile::TempDir, CallFromAlignmentsArgs, ReferenceInfo) {
    use crate::ng::reference_info::{ReferenceSource, read_reference_info_observing};
    use crate::ng::repeat_catalog::RepeatCatalogBuilder;
    use crate::ng::tandem_repeat::ScanParams;
    use std::io::Write;

    let bases = the_fixture_references_bases();
    let directory = tempfile::tempdir().expect("a temporary directory");
    let fasta = directory.path().join("ref.fa");
    let header = ">chr1\n";
    std::fs::write(
        &fasta,
        format!("{header}{}\n", String::from_utf8_lossy(&bases)),
    )
    .expect("the reference writes");
    let mut fai = std::fs::File::create(directory.path().join("ref.fa.fai")).expect("an index");
    writeln!(
        fai,
        "chr1\t{}\t{}\t{}\t{}",
        bases.len(),
        header.len(),
        bases.len(),
        bases.len() + 1
    )
    .expect("the index writes");

    // What `repeat-catalog` would have written: the file's own storage floors, which sit
    // below every floor a caller routes on.
    let catalog_path = directory.path().join("ref.fa.repeats.parquet");
    let mut builder = RepeatCatalogBuilder::create(
        &catalog_path,
        StrRepeatCriteria::default(),
        ScanParams::default(),
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

    let mut args = args_of(&a_defaults_run());
    args.reference = fasta;
    args.catalog = Some(catalog_path);
    (directory, args, reference)
}

/// **A repeat the catalog holds and this run's floors turn down becomes ordinary sequence, not
/// a hole** (`run_ssr_observations.md` §2.2) — and that is the whole of what the routing change
/// buys.
///
/// One switch, two settings, the same reference and the same catalog file. At the **catalog's**
/// floors — what every run asked with until now — both homopolymers are repeat tracts, and on
/// the human benchmark every truth variant on such ground was missed, because nothing calls a
/// tract yet. At **ng's calling floors**, the six-base run is generic ground and goes to the
/// SNP/indel caller, which finds variants there at 0.98 recall; the ten-base run is a tract
/// under either setting, so this is a test of the floor and not of the plumbing.
#[test]
fn a_tract_below_the_calling_floor_becomes_generic_ground_and_one_above_it_stays_a_tract() {
    let (_directory, calling_floors, reference) =
        a_reference_with_a_tract_on_each_side_of_the_calling_floor();

    // The same run asking with the file's own floors — which is what `StrRepeatCriteria::
    // default()` is, and what this command passed before it had flags.
    let mut catalog_floors = calling_floors.clone();
    catalog_floors.min_copies =
        crate::pop_var_caller_exp::cli::parsers::parse_min_copies("5,5,4,4,4,3")
            .expect("the catalog's own table");
    catalog_floors.max_str_len = StrRepeatCriteria::default().max_str_len_bp.get();

    let bounds = vec![ContigBounds {
        name: "chr1",
        length: 136,
    }];
    let analysed = GenomeRegions::whole_contigs(&bounds);
    let kind_covering = |args: &CallFromAlignmentsArgs, position: u64| {
        let segmentation = segments_over(&super::ground_request(args), &analysed, &reference)
            .expect("the catalog answers");
        segmentation
            .segments()
            .iter()
            .find(|segment| {
                segment.region.start.get() <= position && position <= segment.region.end.get()
            })
            .map(|segment| segment.kind.clone())
            .expect("the segments partition the contig, so every base is in one")
    };

    // The fixture is the claim's other half: if the six-base run were not a tract at the
    // file's floors there would be nothing for the calling floors to turn down, and the
    // assertion below would pass over a reference with no repeat in it at all.
    assert!(
        matches!(
            kind_covering(&catalog_floors, 43),
            RegionKind::SsrSegment(_)
        ),
        "the six-base run is in the catalog, which is what makes it a candidate",
    );
    assert!(
        matches!(
            kind_covering(&catalog_floors, 90),
            RegionKind::SsrSegment(_)
        ),
        "and so is the ten-base run",
    );

    assert_eq!(
        kind_covering(&calling_floors, 43),
        RegionKind::Generic,
        "six copies is below ng's period-1 floor of eight, so its bases go to the SNP/indel \
         caller rather than becoming ground nobody speaks for",
    );
    assert!(
        matches!(
            kind_covering(&calling_floors, 90),
            RegionKind::SsrSegment(_)
        ),
        "ten copies clears the floor under either setting — the fixture moves one tract, not \
         both",
    );
}

/// **The parameters file a run writes says what that run counted as a repeat**
/// (`parameters_file.md` §3.9), and it says what this run typed rather than what the binary
/// defaults to.
///
/// Two runs over the same reference and the same catalog with different floors here analyse
/// different ground, and nothing else in the file would show it. The flags are moved away from
/// every default before the run, so a writer that recorded `StrRepeatCriteria::default()`, or
/// the calling defaults, or the catalog's own header, fails.
#[test]
fn the_written_parameters_file_records_what_this_run_counted_as_a_repeat() {
    let (_reference_dir, _zeta_dir, _alpha_dir, mut args) = a_cohort_on_disk();
    args.min_copies = crate::pop_var_caller_exp::cli::parsers::parse_min_copies("9,7,7,7,6,5")
        .expect("six floors");
    args.min_period = 2;
    args.max_period = 5;
    args.max_str_len = 64;
    args.min_purity = 0.95;

    run_call_from_alignments(&args).expect("the cohort runs");

    let written = std::fs::read_to_string(beside_the_vcf(&args.output)).expect("it was written");
    let file = ParametersFile::from_toml(&written).expect("what this run wrote, it can read");
    let routing = file
        .repeat_routing
        .expect("every run this build makes records what it routed with");

    assert_eq!(routing.min_copies, [9, 7, 7, 7, 6, 5]);
    assert_eq!(routing.min_period, 2);
    assert_eq!(routing.max_period, 5);
    assert_eq!(routing.max_str_len, 64);
    assert!((routing.min_purity - 0.95).abs() < 1e-6);
    assert_eq!(
        routing.min_flank_bp,
        StrRepeatCriteria::default().min_flank_bp.get(),
        "not a flag, and pinned at the floor the catalog was built at",
    );

    assert!(
        written.contains("[repeat_routing]"),
        "and it is a section a person can find and edit, not only a field serde knows",
    );
}

// ---------------------------------------------------------------------
// Routing parity, and the no-regression pin
// ---------------------------------------------------------------------

/// **The run's ground partition is what the dump prints for the same reference at the same
/// floors** — spec §10's standing oracle, kept here as a fixture test.
///
/// `examples/ng_typed_region_dump.rs` is what a person runs to ask *which known variants fall
/// on ground the caller sends down the repeat path*, and every number quoted from it is a claim
/// about a run. That only holds while the two ask the catalog the same question. They resolve
/// their ground separately — the tool from a BED or from whole contigs, the run through
/// `analysed_regions` — and they build their criteria separately, so either could drift with
/// nothing failing.
#[test]
fn the_runs_ground_partition_is_the_dumps_at_the_same_floors() {
    let (_directory, args, reference) =
        a_reference_with_a_tract_on_each_side_of_the_calling_floor();
    let bounds = vec![ContigBounds {
        name: "chr1",
        length: 136,
    }];
    let analysed = GenomeRegions::whole_contigs(&bounds);

    for floors in [
        StrRepeatCriteria::from(&TypedRegionConfig::default()),
        StrRepeatCriteria::default(),
    ] {
        let mut asking_with = args.clone();
        asking_with.min_copies = the_min_copies_of(&floors);
        asking_with.min_period = floors.classification.periods.min();
        asking_with.max_period = floors.classification.periods.max();
        asking_with.max_str_len = floors.max_str_len_bp.get();
        asking_with.min_purity = floors.classification.min_purity;

        // The dump's own two lines, verbatim: open the catalog beside the reference, and walk
        // it at these criteria over the analysed ground.
        let catalog = RepeatCatalog::open_beside_reference(&args.reference, &reference)
            .expect("the catalog opens");
        let ground: Vec<_> = analysed.iter().collect();
        let dumped: Vec<TypedRegion> = catalog
            .genome_segments(&floors, ReadScope::Regions(&ground))
            .expect("the dump walks")
            .collect::<Result<_, _>>()
            .expect("every region reads");

        let run = segments_over(&super::ground_request(&asking_with), &analysed, &reference)
            .expect("the run walks");
        assert_eq!(
            run.segments(),
            dumped.as_slice(),
            "the run and the dump disagree about the ground at {:?}",
            floors.classification.min_copies,
        );

        // **And the partition is exact**, which the dump cannot say because it prints whatever
        // it is handed. Every base of the contig is in exactly one region, in order.
        let mut next = 1;
        for region in run.segments() {
            assert_eq!(
                region.region.start.get(),
                next,
                "a gap or an overlap before {region:?}",
            );
            next = region.region.end.get() + 1;
        }
        assert_eq!(next, 137, "the partition stops short of the contig's end");
    }
}

/// The six floors of `criteria`, in the shape the flag parses to.
fn the_min_copies_of(criteria: &StrRepeatCriteria) -> MinCopies {
    let floors: Vec<String> = (1..=6)
        .map(|period| {
            criteria
                .classification
                .min_copies
                .for_period(period)
                .to_string()
        })
        .collect();
    crate::pop_var_caller_exp::cli::parsers::parse_min_copies(&floors.join(","))
        .expect("six floors")
}

/// **Where the routing did not move, nothing about the run moved either** — spec §10's
/// no-regression claim, and the other half of the change B1 made.
///
/// The reference holds a six-base run of `A` that is a repeat tract under the catalog's floors
/// and ordinary sequence under ng's. **Every read here sits on the filler between the
/// homopolymers**, which is ordinary sequence under both — so the two runs analyse the same
/// ground *where there is any evidence*, and the VCF must be identical byte for byte. A run
/// whose reads covered the six-base run would legitimately differ, which is why they do not.
///
/// This is what makes B1 a change of *which* ground is generic rather than a change of what
/// happens on it — and a defect that made the criteria leak into the calling arithmetic, rather
/// than into the routing alone, is what it would catch.
#[test]
fn where_the_routing_did_not_move_the_vcf_is_byte_identical() {
    let (directory, mut args, _reference) =
        a_reference_with_a_tract_on_each_side_of_the_calling_floor();

    // One sample, reads over the filler at 5–34 and 101–130 only — never over either
    // homopolymer, whose routing is the thing that moves between the two runs. Thirty bases
    // because that is the shortest read the filters keep, and each is the reference with its
    // fifteenth base flipped, so there is one variant to call rather than thirty.
    let bases = the_fixture_references_bases();
    let header = crate::ng::read::input::test_fixtures::header(
        Some("coordinate"),
        &[("chr1", 136, None)],
        &[(
            &crate::ng::read::input::test_fixtures::read_group_for("one"),
            Some("one"),
        )],
    );
    let mut reads = Vec::new();
    for start in [5usize, 101] {
        for copy in 0..4 {
            let mut observed = bases[start - 1..start + 29].to_vec();
            observed[14] = if observed[14] == b'G' { b'T' } else { b'G' };
            reads.push(a_read_showing(
                &format!("r{start}-{copy}"),
                start,
                &observed,
            ));
        }
    }
    let (_bam_dir, bam) =
        crate::ng::read::input::test_fixtures::indexed_named_bam(&header, &reads, "one.bam");
    args.alignments = vec![bam];

    // **The same file name in two directories**, because the VCF header names the parameters
    // file beside it — two runs writing `calling.vcf` and `catalog.vcf` would differ on that
    // line and on nothing else, which is a difference in what the test asked for rather than in
    // what it is testing.
    let vcf_of = |args: &CallFromAlignmentsArgs, into: &str| {
        let mut args = args.clone();
        let directory = directory.path().join(into);
        std::fs::create_dir(&directory).expect("a directory to write into");
        args.output = directory.join("calls.vcf");
        run_call_from_alignments(&args).expect("the cohort runs");
        (
            std::fs::read(&args.output).expect("the VCF was written"),
            args,
        )
    };

    let (calling, at_the_calling_floors) = vcf_of(&args, "calling");

    let mut asking_the_catalogs_floors = args.clone();
    asking_the_catalogs_floors.min_copies = the_min_copies_of(&StrRepeatCriteria::default());
    asking_the_catalogs_floors.max_str_len = StrRepeatCriteria::default().max_str_len_bp.get();
    let (catalog, at_the_catalogs_floors) = vcf_of(&asking_the_catalogs_floors, "catalog");

    // The fixture's own half: the two runs really did route differently, or this compares two
    // identical runs and could not fail.
    let partition_at = |args: &CallFromAlignmentsArgs| {
        let bounds = vec![ContigBounds {
            name: "chr1",
            length: 136,
        }];
        segments_over(
            &super::ground_request(args),
            &GenomeRegions::whole_contigs(&bounds),
            &_reference,
        )
        .expect("the run walks")
        .segments()
        .iter()
        .filter(|region| matches!(region.kind, RegionKind::SsrSegment(_)))
        .count()
    };
    assert_eq!(
        partition_at(&at_the_calling_floors),
        1,
        "only the ten-base run clears ng's period-1 floor of eight",
    );
    assert_eq!(
        partition_at(&at_the_catalogs_floors),
        2,
        "both clear the catalog's floor of five, which is the difference under test",
    );

    // And there is something to compare: a header-only pair of files would be equal whatever
    // the caller did with the reads.
    let records = |vcf: &[u8]| {
        String::from_utf8_lossy(vcf)
            .lines()
            .filter(|line| !line.starts_with('#'))
            .count()
    };
    assert!(
        records(&calling) > 0,
        "the fixture's reads produced no record, so this compares two headers",
    );

    assert_eq!(
        String::from_utf8_lossy(&calling),
        String::from_utf8_lossy(&catalog),
        "the routing moved and the calls did not",
    );
}

/// One 12-base read starting at 1-based `start`, showing `observed`.
///
/// `read_named_with_length` writes a run of `A`, which on this fixture's reference is twelve
/// mismatches rather than one variant.
fn a_read_showing(name: &str, start: usize, observed: &[u8]) -> noodles_sam::alignment::RecordBuf {
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record::{Flags, MappingQuality};
    use noodles_sam::alignment::record_buf::{QualityScores, Sequence};

    RecordBuf::builder()
        .set_name(name.as_bytes())
        .set_reference_sequence_id(0usize)
        .set_flags(Flags::empty())
        .set_mapping_quality(MappingQuality::new(60).expect("mapq in range"))
        .set_alignment_start(noodles_core::Position::try_from(start).expect("a position"))
        .set_cigar([Op::new(Kind::Match, observed.len())].into_iter().collect())
        .set_sequence(Sequence::from(observed.to_vec()))
        .set_quality_scores(QualityScores::from(vec![30u8; observed.len()]))
        .build()
}
