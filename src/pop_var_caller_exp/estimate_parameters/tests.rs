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
    a_cohort_walked_asking_for(MinCopies::default())
}

/// The same, walked under a copy floor of the caller's choosing.
///
/// **What it is for: a cohort whose psps were not walked under this command's default flags.**
/// The fixture's repeat tract is ten copies of `GT`, so a walk asking for eleven leaves it
/// ordinary sequence, and a fit that rebuilt its selection under the default floors would keep
/// positions inside it that the censuses hold nothing at.
fn a_cohort_walked_asking_for(min_copies: MinCopies) -> (AVaryingCohort, PathBuf) {
    a_cohort_walked(min_copies, None)
}

/// The same again, with the ground the walk covered given as well.
///
/// **`regions` is what makes the cohort's ground different from the whole reference.** Every
/// other fixture here walks whole contigs, where *the ground the psps were walked over* and
/// *every base of the reference* are the same set of positions and no test can tell a command
/// that takes one from a command that takes the other.
fn a_cohort_walked(min_copies: MinCopies, regions: Option<&str>) -> (AVaryingCohort, PathBuf) {
    let cohort = a_varying_cohort_on_disk();
    let psps = cohort.directory.path().join("psps");
    let regions = regions.map(|bed| {
        let path = cohort.directory.path().join("walked.bed");
        std::fs::write(&path, bed).expect("the scratch dir is ours");
        path
    });
    let walk = GeneratePspsArgs {
        reference: cohort.reference.clone(),
        catalog: Some(cohort.catalog.clone()),
        alignments: cohort.alignments.clone(),
        output_dir: psps.clone(),
        regions,
        force: false,
        build_index_if_missing: false,
        min_copies,
        min_period: DEFAULT_MIN_PERIOD,
        max_period: DEFAULT_MAX_PERIOD,
        max_str_len: DEFAULT_MAX_STR_LEN,
        min_purity: DEFAULT_MIN_PURITY,
    };
    run_generate_psps(&walk).expect("the cohort walks into psps");
    // **The walk seals the census into each psp now** (`psp_census_pair.md` §3), and this
    // command still takes census *files* — so they are built here, by `generate-census`, which
    // is the producer that still writes them. Plan step C2 is what moves this command onto the
    // psps themselves and deletes the pair.
    crate::pop_var_caller_exp::test_fixtures::censuses_written_beside_the_psps(&walk, &psps);
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

/// **The repeat-tract criteria come out of the psps, and the flags that used to supply them are
/// not read** (spec §6).
///
/// One cohort, fitted twice with the copy floor set two different ways on the command line, and
/// the two files are the same bytes. Before this step the flags built the run's segmentation, so
/// the second fit rebuilt its selection over a different set of stretches from the one the
/// censuses were written against — eleven copies where the fixture's tract has ten leaves that
/// tract ordinary sequence, so the selection keeps positions inside it that the censuses hold
/// nothing at — and the fit refused the cohort as *built under another selection*.
///
/// **Why the flags cannot simply be believed.** They are values a person retypes from a walk that
/// already recorded them, and a mistyped one does not fail — it fits a plausible file over a
/// different set of loci from the one the censuses hold, and is caught, if at all, by the fit's
/// digest check reporting *another selection* twenty seconds later.
#[test]
fn the_criteria_flags_do_not_change_the_file_that_is_fitted() {
    let (cohort, psps) = a_walked_cohort();
    let output = cohort.directory.path().join("cohort.parameters.toml");

    let (as_walked, _) = fit_and_assemble(&args_over(&cohort, &psps, output.clone()))
        .expect("the cohort fits with the flags left at their defaults");
    let mut asking_for_eleven = args_over(&cohort, &psps, output);
    asking_for_eleven.min_copies = MinCopies::uniform(11);
    let (whatever_the_flags_said, _) =
        fit_and_assemble(&asking_for_eleven).expect("and fits with them set to anything else");

    assert_eq!(
        as_walked.to_toml(),
        whatever_the_flags_said.to_toml(),
        "the psps say what a repeat is; the flags are not read",
    );
}

/// **A cohort walked under criteria this command was never told is fitted all the same.**
///
/// The other half of the same rule: the psps here were walked asking for eleven copies, and the
/// fit is run with the flags at their defaults — the shape of every real run once the flags are
/// gone (plan step C2), and the one a person hits today by walking with `--min-copies` and
/// forgetting to repeat it here.
///
/// **It asserts the file differs from the default walk's**, because *it fitted something* is
/// satisfied by a fit that quietly used the defaults: the two walks type the fixture's tract two
/// different ways — a tract for one, ordinary sequence for the other — so two files that matched
/// would mean the walk's criteria never reached the arithmetic.
#[test]
fn a_cohort_walked_under_other_criteria_is_fitted_without_being_told_them() {
    let (cohort, psps) = a_cohort_walked_asking_for(MinCopies::uniform(11));
    let output = cohort.directory.path().join("cohort.parameters.toml");

    let (from_those_psps, samples) = fit_and_assemble(&args_over(&cohort, &psps, output))
        .expect("the criteria are read from the psps, so the default flags are no disagreement");

    assert_eq!(samples, 2, "one entry a sample");
    let (default_cohort, default_psps) = a_walked_cohort();
    let (from_a_default_walk, _) = fit_and_assemble(&args_over(
        &default_cohort,
        &default_psps,
        default_cohort
            .directory
            .path()
            .join("cohort.parameters.toml"),
    ))
    .expect("the default walk fits too");
    assert_ne!(
        from_those_psps.to_toml(),
        from_a_default_walk.to_toml(),
        "the two walks type the fixture's tract differently, so their fits cannot agree",
    );
}

/// **`--catalog` names the file that is read**, and not the one beside the reference.
///
/// The two are the same file in every other test here — the fixture's catalog is written to
/// `ref.fa.repeats.parquet`, which is exactly where a run looks when it is told nothing — so
/// nothing else in this file can tell a command that honours the flag from one that ignores it.
/// A person keeping catalogs in one directory and references in another (this project's own
/// tomato benchmark does, because the reference is on a read-only mount) is the ordinary case.
#[test]
fn the_catalog_flag_names_the_file_that_is_read() {
    let (cohort, psps) = a_walked_cohort();
    let elsewhere = cohort.directory.path().join("catalogs");
    std::fs::create_dir(&elsewhere).expect("the scratch dir is ours");
    let moved = elsewhere.join("the-only-catalog.parquet");
    std::fs::rename(&cohort.catalog, &moved).expect("the catalog moves");

    let output = cohort.directory.path().join("cohort.parameters.toml");
    let mut args = args_over(&cohort, &psps, output);
    args.catalog = Some(moved);

    let (from_where_it_was_moved_to, samples) =
        fit_and_assemble(&args).expect("the catalog it was pointed at is the one it reads");

    assert_eq!(samples, 2, "one entry a sample");
    let (beside_the_reference, psps_beside_it) = a_walked_cohort();
    let (from_the_sibling_path, _) = fit_and_assemble(&args_over(
        &beside_the_reference,
        &psps_beside_it,
        beside_the_reference
            .directory
            .path()
            .join("cohort.parameters.toml"),
    ))
    .expect("and the same cohort fits with its catalog left where a run looks by default");
    assert_eq!(
        from_where_it_was_moved_to.to_toml(),
        from_the_sibling_path.to_toml(),
        "one catalog read from two paths is one answer",
    );
}

/// **The selection is rebuilt over the ground the psps were walked on, not over the reference.**
///
/// This cohort is walked over the first 400 bases of a 600-base contig and fitted with no
/// `--regions` at all — which is every real run of this command, since it has no such flag: the
/// ground comes from the psp headers (spec §5.3, §6).
///
/// **Why it is worth its own fixture.** A fit that took every base of the reference does not fail
/// loudly: it rebuilds a plausible selection that keeps positions in the 200 bases nobody walked,
/// where the censuses hold nothing, and what the user sees is the cohort refused as *built under
/// another selection* — if the digest notices at all. Every other cohort in this file is walked
/// whole, so this is the only test here that can tell the two apart.
#[test]
fn the_selection_is_rebuilt_over_the_ground_the_psps_were_walked_on() {
    // `VARYING_CONTIG` is `("chrV", 600)`, and the fixture's repeat tract sits at 200..220, so
    // the walked part holds the tract and the unwalked part does not.
    let (cohort, psps) = a_cohort_walked(MinCopies::default(), Some("chrV\t0\t400\n"));
    let output = cohort.directory.path().join("cohort.parameters.toml");

    let (file, samples) = fit_and_assemble(&args_over(&cohort, &psps, output))
        .expect("the ground comes from the psps, so a partial walk fits");

    assert_eq!(samples, 2, "one entry a sample");
    assert!(
        !file.fitted_from.census.terms.is_empty(),
        "and the fit had evidence to fit from",
    );
}

/// **A catalog that cannot answer the psps' criteria is refused without naming a flag to move.**
///
/// The refusal for *the reader asks for tracts below what this file holds* names the flag that
/// asked — `--min-copies`, say — because on a walk that is what the person typed and what they
/// have to change (`run_ground::catalog_error_naming_the_flag`). On this command nobody typed it:
/// the criteria came out of the psps, and no flag here can move them. What the person can move is
/// which catalog they pointed at, which the general refusal names.
///
/// The catalog built here holds tracts of twenty copies and up; the walk's own calling floors ask
/// for six at period 2, so the file cannot serve them — the rows are not in it.
#[test]
fn a_catalog_that_cannot_serve_the_psps_criteria_is_refused_without_naming_a_flag() {
    use crate::ng::reference_info::{ReferenceSource, read_reference_info_observing};
    use crate::ng::region_typing::segment_criteria::SsrSegmentCriteria;
    use crate::ng::repeat_catalog::RepeatCatalogBuilder;
    use crate::ng::repeat_catalog::StrRepeatCriteria;
    use crate::ng::tandem_repeat::ScanParams;

    let (cohort, psps) = a_walked_cohort();
    let coarse = cohort
        .directory
        .path()
        .join("twenty-copies.repeats.parquet");
    let mut builder = RepeatCatalogBuilder::create(
        &coarse,
        StrRepeatCriteria {
            classification: SsrSegmentCriteria {
                min_copies: MinCopies::uniform(20),
                ..StrRepeatCriteria::default().classification
            },
            ..StrRepeatCriteria::default()
        },
        ScanParams {
            match_reward: 2,
            mismatch_penalty: 7,
            min_copies: 20,
        },
    )
    .expect("a catalog to build into");
    let reference_info = read_reference_info_observing(
        ReferenceSource::Fasta {
            fasta: cohort.reference.clone(),
            fai: None,
        },
        &mut builder,
    )
    .expect("the reference reads");
    builder.finish(&reference_info).expect("the catalog writes");

    let output = cohort.directory.path().join("cohort.parameters.toml");
    let mut args = args_over(&cohort, &psps, output);
    args.catalog = Some(coarse);

    let error = fit_and_assemble(&args).expect_err("that catalog does not hold what the psps ask");

    let rendered = crate::error_render::format_error_chain(&error);
    assert!(
        !rendered.contains("--min-copies"),
        "no flag on this command line can move the criteria, and got: {rendered}",
    );
    assert!(
        rendered.contains("twenty-copies.repeats.parquet"),
        "the catalog the run was pointed at is what the reader can change, and got: {rendered}",
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
