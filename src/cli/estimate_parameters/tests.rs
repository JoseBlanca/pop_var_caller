//! `estimate-parameters` at the command line: what a person may type, what is refused before a
//! census is opened, and what the file it writes says about where its numbers came from.

use super::*;
use clap::Parser;
use std::path::Path;

use crate::cli::command_line::{Cli, PopVarCallerCommand};
use crate::cli::generate_psps::{GeneratePspsArgs, run_generate_psps};
use crate::cli::test_fixtures::{AVaryingCohort, a_varying_cohort_on_disk};
use crate::region_typing::DEFAULT_MAX_STR_LEN;
use crate::region_typing::segment_criteria::{
    DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
};

/// Parse an argument vector into this subcommand's arguments, refusing any other subcommand.
fn args_of(argv: &[&str]) -> EstimateParametersArgs {
    match Cli::parse_from(argv).cmd {
        PopVarCallerCommand::EstimateParameters(args) => args,
        other => panic!("expected estimate-parameters, got {other:?}"),
    }
}

/// The shortest run a person can type.
fn a_shortest_run() -> Vec<&'static str> {
    vec![
        "pop_var_caller",
        "estimate-parameters",
        "--reference",
        "ref.fa",
        "--psp",
        "zeta.psp",
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
    (cohort, psps)
}

fn args_over(cohort: &AVaryingCohort, psps: &Path, output: PathBuf) -> EstimateParametersArgs {
    EstimateParametersArgs {
        reference: cohort.reference.clone(),
        catalog: Some(cohort.catalog.clone()),
        psps: vec![psps.to_path_buf()],
        output,
        force: false,
        ploidy: 2,
        inbreeding: None,
        str_param_estimates_at_once: std::num::NonZeroUsize::MIN,
    }
}

#[test]
fn the_subcommand_is_spelled_estimate_parameters() {
    let args = args_of(&a_shortest_run());
    assert_eq!(args.reference, PathBuf::from("ref.fa"));
    assert_eq!(args.psps, vec![PathBuf::from("zeta.psp")]);
    assert_eq!(args.output, PathBuf::from("cohort.parameters.toml"));
    assert!(
        Cli::try_parse_from(["pop_var_caller", SUBCOMMAND, "--help"])
            .expect_err("--help exits")
            .to_string()
            .contains(SUBCOMMAND),
        "the name this module records is the one clap answers to",
    );
}

/// **One stratum at a time unless the run says otherwise**, and the run says it as a count of at
/// least one.
///
/// The default is the smallest memory whatever the cohort, so a run told nothing cannot run out
/// of memory in this step for want of a flag; zero strata at once is not a schedule and is refused
/// when the command line is read.
#[test]
fn strata_are_fitted_one_at_a_time_unless_the_run_asks_for_more() {
    assert_eq!(
        args_of(&a_shortest_run()).str_param_estimates_at_once.get(),
        1
    );

    let mut asked = a_shortest_run();
    asked.extend(["--str-param-estimates-at-once", "4"]);
    assert_eq!(args_of(&asked).str_param_estimates_at_once.get(), 4);

    assert_eq!(
        repeat_tract_config(&args_of(&asked)).strata_at_once.get(),
        4,
        "the count the run asked for reaches the repeat-tract fit"
    );
    assert_eq!(
        repeat_tract_config(&args_of(&a_shortest_run())).strata_at_once,
        crate::parameter_estimation::joint::ssr_fit::DEFAULT_STRATA_AT_ONCE
    );

    let mut zero = a_shortest_run();
    zero.extend(["--str-param-estimates-at-once", "0"]);
    assert!(
        Cli::try_parse_from(zero).is_err(),
        "zero strata at once is refused, not read as one"
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

/// **A cohort walked under criteria this command cannot be told is fitted all the same** (spec
/// §6).
///
/// The psps here were walked asking for eleven motif copies where the default is six at period 2,
/// and this command has no flag that says so — there is nothing to type and nothing to get wrong.
/// Before plan step C1 the same run needed `--min-copies 11,11,11,11,11,11` and refused the
/// cohort without it.
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

/// **A catalog the psps were not walked with is refused as that, before a segment is cut** (plan
/// step C5, which makes spec §6's check exist).
///
/// The catalog built here is on the cohort's own reference and holds tracts of twenty copies and
/// up, where the psps' was built at the default floors — the same assembly, another file.
///
/// **How the order is shown.** Cut with the psps' criteria, which ask for six copies at period 2,
/// this catalog cannot serve them — the rows are not in it — so a run that cut segments before
/// comparing headers is refused as a catalog too coarse for its criteria: true, and no help to a
/// person whose fault is which file they named.
///
/// **And no flag to move is named, nor a census to regenerate.** The criteria came out of the
/// psps, so what this command line can change is which catalog it points at; and a census
/// regenerated against this catalog would be refused again by the next fit with the right one.
#[test]
fn a_catalog_the_psps_were_not_walked_with_is_refused_before_segments_are_cut() {
    use crate::reference_info::{ReferenceSource, read_reference_info_observing};
    use crate::region_typing::segment_criteria::SsrSegmentCriteria;
    use crate::repeat_catalog::RepeatCatalogBuilder;
    use crate::repeat_catalog::StrRepeatCriteria;
    use crate::tandem_repeat::ScanParams;

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
    args.catalog = Some(coarse.clone());

    let error = fit_and_assemble(&args).expect_err("the psps were walked with another catalog");

    assert!(
        matches!(
            &error,
            EstimateParametersCliError::WalkedWithAnotherCatalog { path, .. } if path == &coarse
        ),
        "the catalog is compared with the psps' before it is cut, and got: {error:?}",
    );
    let rendered = crate::error_render::format_error_chain(&error);
    assert!(
        rendered.contains("it was built under other repeat criteria"),
        "and what about it is not theirs: {rendered}",
    );
    assert!(
        rendered.contains("run this fit again with --catalog naming the one they were"),
        "and what to do: {rendered}",
    );
    assert!(
        !rendered.contains("--min-copies")
            && !rendered.contains(THE_COMMAND_THAT_REBUILDS_A_CENSUS),
        "no flag on this command line can move the criteria, and regenerating fixes nothing: \
         {rendered}",
    );
}

/// **A reference the psps were not walked against is refused as that, before the selection is
/// rebuilt** (plan step C5).
///
/// The reference here is the cohort's own with its first base changed, so its one contig has the
/// psps' name and length and only its bases differ — what a person meets when two builds of an
/// assembly share their contig names.
///
/// **How the order is shown: the catalog is the psps' own**, built on the right reference. A run
/// that opened the catalog before comparing the reference with the psps would be refused as a
/// catalog built on another reference, blaming the one file that is right; one that reached the
/// comparison of recorded settings would name every census as stale and tell the person to
/// regenerate them against the wrong reference.
#[test]
fn a_reference_the_psps_were_not_walked_against_is_refused_before_the_selection_is_rebuilt() {
    use crate::cli::test_fixtures::{VARYING_CONTIG, the_varying_cohorts_reference};

    let (cohort, psps) = a_walked_cohort();
    let mut bases = the_varying_cohorts_reference();
    bases[0] = match bases[0] {
        b'A' => b'C',
        _ => b'A',
    };
    let another = cohort.directory.path().join("another-build.fa");
    std::fs::write(
        &another,
        format!(
            ">{}\n{}\n",
            VARYING_CONTIG.0,
            std::str::from_utf8(&bases).expect("ACGT is text")
        ),
    )
    .expect("the scratch dir is ours");
    let mut args = args_over(&cohort, &psps, cohort.directory.path().join("out.toml"));
    args.reference = another.clone();

    let error =
        fit_and_assemble(&args).expect_err("the psps were walked against another reference");

    assert!(
        matches!(
            &error,
            EstimateParametersCliError::WalkedAgainstAnotherReference { path, .. } if path == &another
        ),
        "the reference is compared with the psps before the catalog is opened, and got: {error:?}",
    );
    let rendered = crate::error_render::format_error_chain(&error);
    assert!(
        rendered.contains(
            "they name ref.fa, so run this fit again with that one, which regenerates nothing"
        ),
        "what to do, and which file to do it with — the psps' headers name it: {rendered}",
    );
    assert!(
        rendered.contains("the psp for sample one was written against a different reference"),
        "and which psp says so, and how: {rendered}",
    );
    assert!(
        !rendered.contains(THE_COMMAND_THAT_REBUILDS_A_CENSUS),
        "regenerating against the wrong reference fixes nothing: {rendered}",
    );
}

/// **A stated coefficient overrides the fitted one and is recorded as supplied.**
///
/// The override is the owner's ruling of 2026-08-27 — a user who knows how their material was
/// bred knows it whatever the cohort size — and the warrant is what keeps it honest: a file that
/// reported a declared value as fitted would make a run's own assumption look like a measurement.
/// **The cohort here has two samples and its fit does produce coefficients**, so this pins the
/// precedence and not merely the absence of a fit.
#[test]
fn a_stated_inbreeding_coefficient_overrides_the_fit_and_is_recorded_as_supplied() {
    let (cohort, psps) = a_walked_cohort();
    let output = cohort.directory.path().join("cohort.parameters.toml");
    let mut args = args_over(&cohort, &psps, output);
    args.inbreeding = Some(0.25);

    let (file, _) = fit_and_assemble(&args).expect("the cohort fits");

    assert!(!file.inbreeding.by_sample.is_empty(), "one row a sample");
    for row in &file.inbreeding.by_sample {
        assert_eq!(row.inbreeding_coefficient.value, 0.25);
        assert_eq!(
            row.inbreeding_coefficient.warrant,
            crate::calling::parameters_file::Warrant::Supplied,
            "a declared coefficient is supplied, never fitted",
        );
    }
}

/// **With nothing stated, the file carries what the fit measured, marked `fitted_here`.**
///
/// This is the whole point of writing the coefficient into the parameters file: the calling
/// commands take no flag for it, so a coefficient the fit measured reaches a run only if it is
/// written here. Before 2026-09-11 the command declared zero on every row and threw the fit's own
/// number away.
///
/// **What pins that is the warrant and the observation count, not the value** — and it is worth
/// saying which, because this fixture's two samples are not inbred and their fitted homozygote
/// excess comes back at about 3 × 10⁻¹⁰. So an assertion on the number would pass against a
/// hard-coded zero and prove nothing.
///
/// The observation count is the second, independent discriminator: the writer leaves it out of the
/// file for a `supplied` or `defaulted` value on purpose — *"a zero count would claim a
/// measurement over no genome"* (`from_run_parameters`'s `warranted_value`) — so
/// `observations = { covered_positions = … }` can only have come from a fitted estimate.
#[test]
fn with_nothing_stated_the_fits_own_coefficient_is_written_and_marked_fitted() {
    let (cohort, psps) = a_walked_cohort();
    let output = cohort.directory.path().join("cohort.parameters.toml");
    let args = args_over(&cohort, &psps, output);
    assert!(
        args.inbreeding.is_none(),
        "the fixture must state nothing for this test to be about the fitted rung",
    );

    let (file, _) = fit_and_assemble(&args).expect("the cohort fits");

    assert!(!file.inbreeding.by_sample.is_empty(), "one row a sample");
    for row in &file.inbreeding.by_sample {
        assert_eq!(
            row.inbreeding_coefficient.warrant,
            crate::calling::parameters_file::Warrant::FittedHere,
            "two samples identify a homozygote excess, so it is fitted here: {:?}",
            row,
        );
        assert!(
            (0.0..1.0).contains(&row.inbreeding_coefficient.value),
            "a coefficient is a fraction below one: {:?}",
            row,
        );
        assert!(
            row.inbreeding_coefficient.observations.is_some(),
            "a fitted coefficient says how much data stood behind it, in its own unit: {:?}",
            row,
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

/// **A directory of psps is expanded in name order**, so two runs naming one directory read the
/// same cohort in the same order however the filesystem answers.
#[test]
fn a_directory_contributes_every_psp_inside_it() {
    let (cohort, psps) = a_walked_cohort();
    let args = args_over(&cohort, &psps, cohort.directory.path().join("out.toml"));

    let expanded = psps_named_by(&args).expect("the directory lists");

    assert_eq!(expanded.len(), 2, "one psp a sample, and nothing else");
    assert!(
        expanded
            .iter()
            .all(|path| path.extension().is_some_and(|it| it == "psp")),
        "and nothing that is not a psp: {expanded:?}",
    );
    assert!(
        expanded.windows(2).all(|pair| pair[0] <= pair[1]),
        "sorted by name: {expanded:?}",
    );
}

/// **A psp carrying no census stops the run, and the report names it and says what to run**
/// (spec §4.3).
///
/// This is every psp written before the census moved into the trailer, and every psp an `append`
/// has discarded the trailer of (spec §3.4).
#[test]
fn a_psp_carrying_no_census_is_refused_and_the_report_names_it() {
    let (cohort, psps) = a_walked_cohort();
    let paths = psps_named_by(&args_over(
        &cohort,
        &psps,
        cohort.directory.path().join("out.toml"),
    ))
    .expect("the directory lists");
    crate::psp::replace_trailer(&paths[1], b"").expect("the tail rewrites");

    let error = fit_and_assemble(&args_over(
        &cohort,
        &psps,
        cohort.directory.path().join("out.toml"),
    ))
    .expect_err("that psp carries no census");

    let EstimateParametersCliError::CohortCannotBeFitted { report } = &error else {
        panic!("a psp with no census is a regeneration report, not this: {error:?}");
    };
    let said = report.to_string();
    assert!(
        said.contains(&format!(
            "  two ({}) carries no census\n",
            paths[1].display()
        )),
        "one line naming the individual, its file and the cause: {said}",
    );
    assert!(
        said.contains(&format!(
            "{THE_COMMAND_THAT_REBUILDS_A_CENSUS} --reference {} --catalog {} --psp {}",
            cohort.reference.display(),
            cohort.catalog.display(),
            psps.display(),
        )),
        "and the command, with the arguments this run was given: {said}",
    );
    assert!(
        !said.contains("one ("),
        "the sample that is fine is counted, not listed: {said}",
    );
}

/// Replace the census in each of `paths` with one rebuilt from that psp under a selection keeping
/// `generic_target` positions — the census a build with another position budget would have
/// written, in this build's format.
fn recensus_under_a_budget_of(cohort: &AVaryingCohort, paths: &[PathBuf], generic_target: u64) {
    use crate::parameter_estimation::joint::census_file::write_census;
    use crate::run::census_from_psp;
    use crate::run::test_fixtures::a_census_plan_over_selecting;

    let (segmentation, plan) =
        a_census_plan_over_selecting(&cohort.reference, &cohort.catalog, generic_target);
    for path in paths {
        let rebuilt = census_from_psp(path, &plan, &segmentation).expect("the psp reads");
        let mut bytes = Vec::new();
        write_census(&rebuilt.evidence, &mut bytes).expect("the census encodes");
        crate::psp::replace_trailer(path, &bytes).expect("the tail rewrites");
    }
}

/// The line the report gives a sample whose census records a different position budget.
fn recorded_under_another_budget(sample: &str, psp: &Path) -> String {
    format!(
        "  {sample} ({}) carries a census recorded under settings this run does not use (the \
         first that differs: generic target position count)\n",
        psp.display()
    )
}

/// **A cohort whose censuses were recorded under another position budget is refused, naming every
/// sample, the setting and the command** (plan step C5; spec §4.2, third row).
///
/// **Half the shipped budget is the case that used to fit without a word.** This contig holds
/// fewer ordinary positions than either budget, so both selections keep every one of them, and the
/// fit — which compares the positions kept — had nothing to refuse, while every census's recorded
/// settings said the budget differed. **A budget of three changes the kept set as well**, and is
/// named the same way, because the budget comes before the kept positions in the order the
/// settings are compared.
///
/// **Every census here agrees with every other**, so a check that compared the samples only with
/// each other sees nothing to refuse; what is stale is the whole cohort against the run.
#[test]
fn a_cohort_recorded_under_another_position_budget_is_refused_naming_the_budget() {
    for budget in [CensusSelection::SHIPPED.generic_target / 2, 3] {
        let (cohort, psps) = a_walked_cohort();
        let args = args_over(&cohort, &psps, cohort.directory.path().join("out.toml"));
        let paths = psps_named_by(&args).expect("the directory lists");
        recensus_under_a_budget_of(&cohort, &paths, budget);

        let error =
            fit_and_assemble(&args).expect_err("those censuses were recorded under another budget");

        let EstimateParametersCliError::CohortCannotBeFitted { report } = &error else {
            panic!(
                "at a budget of {budget}, censuses recorded under other settings are a \
                 regeneration report, not this: {error:?}"
            );
        };
        let said = report.to_string();
        assert_eq!(report.stale_count(), 2, "at a budget of {budget}: {said}");
        for (sample, path) in [("one", &paths[0]), ("two", &paths[1])] {
            assert!(
                said.contains(&recorded_under_another_budget(sample, path)),
                "at a budget of {budget}, a line a sample naming the setting: {said}",
            );
        }
        assert!(
            said.ends_with(&the_command_that_regenerates(&args)),
            "at a budget of {budget}, and the command, last: {said}",
        );
    }
}

/// **A census that disagrees with the rest of its cohort is named alone, with the fix** (plan
/// step C5).
///
/// Before C5 this cohort was refused as *samples one and two disagree on generic target position
/// count; they did not record the same thing* — both samples named, and nothing to do.
///
/// **The stale census is the first sample's**, so a refusal naming the pair names the fresh
/// sample too, and one that took the first census as the standard names the fresh sample instead
/// of the stale one.
#[test]
fn a_census_that_disagrees_with_the_rest_of_its_cohort_is_named_alone() {
    let (cohort, psps) = a_walked_cohort();
    let args = args_over(&cohort, &psps, cohort.directory.path().join("out.toml"));
    let paths = psps_named_by(&args).expect("the directory lists");
    recensus_under_a_budget_of(&cohort, &paths[..1], 3);

    let error = fit_and_assemble(&args).expect_err("one census was recorded under another budget");

    let EstimateParametersCliError::CohortCannotBeFitted { report } = &error else {
        panic!("one stale census is a regeneration report, not this: {error:?}");
    };
    let said = report.to_string();
    assert_eq!(report.stale_count(), 1, "{said}");
    assert!(
        said.contains(&recorded_under_another_budget("one", &paths[0])),
        "the stale sample, its file and the setting: {said}",
    );
    assert!(
        !said.contains("two (")
            && said.contains("1 other sample carries a census this build reads."),
        "the fresh sample is counted, not listed: {said}",
    );
}

/// **The command the report names is one a person can run, and it is the repair.**
///
/// Until plan step D1 it named nothing that existed — deliberately, because the command it must
/// never name is the one that wrote a census *file* beside the psp, which no fit reads: following
/// that would cost the rebuild and leave the psp exactly as stale. Now the name is
/// `regenerate-census`'s own, and this parses the word the report prints, so a report that sent
/// somebody after a word the binary does not know fails here.
#[test]
fn the_report_names_the_command_that_rebuilds_a_psps_census() {
    // `regenerate_census::SUBCOMMAND` *is* this constant, so comparing the two is `X == X`;
    // what carries the claim is the parse below.
    let parsed = Cli::try_parse_from([
        "pop_var_caller",
        THE_COMMAND_THAT_REBUILDS_A_CENSUS,
        "--reference",
        "ref.fa",
        "--psp",
        "zeta.psp",
    ])
    .expect("the command this report prints is one clap answers to");

    assert!(
        matches!(parsed.cmd, PopVarCallerCommand::RegenerateCensus(_)),
        "and it is the repair rather than another subcommand: {:?}",
        parsed.cmd,
    );
}

/// **The command names no catalog when the run was not given one**, because the command it names
/// finds the same file beside the reference by the same rule.
#[test]
fn the_report_names_no_catalog_when_the_run_named_none() {
    let (cohort, psps) = a_walked_cohort();
    let paths = psps_named_by(&args_over(
        &cohort,
        &psps,
        cohort.directory.path().join("out.toml"),
    ))
    .expect("the directory lists");
    crate::psp::replace_trailer(&paths[1], b"").expect("the tail rewrites");
    let mut args = args_over(&cohort, &psps, cohort.directory.path().join("out.toml"));
    // The fixture's catalog is where a run looks when it is told nothing, so dropping the flag
    // changes what is printed and not what would be read.
    args.catalog = None;

    let error = fit_and_assemble(&args).expect_err("that psp carries no census");

    let EstimateParametersCliError::CohortCannotBeFitted { report } = &error else {
        panic!("a psp with no census is a regeneration report: {error:?}");
    };
    let said = report.to_string();
    assert!(!said.contains("--catalog"), "{said}");
    assert!(
        said.contains(&format!(
            "{THE_COMMAND_THAT_REBUILDS_A_CENSUS} --reference {} --psp {}",
            cohort.reference.display(),
            psps.display(),
        )),
        "and the rest of the line is as it was typed: {said}",
    );
}

/// **Two stale psps in a cohort are both named, in one message** — spec §4.1, and the whole point
/// of judging the cohort before refusing it.
///
/// **Rebuilding one census is a quarter of an hour** (spec §2, §4), so a refusal that named one
/// sample at a time would cost that wait once a stale sample, one after another, to learn a job
/// that fits in one message. The two here are stale for different reasons — one carries no census,
/// the other a census of a version this build does not read — because the report groups by cause
/// and a run that reported whatever the first stale psp was would pass a test with one.
#[test]
fn two_stale_psps_are_both_named_in_one_report() {
    let (cohort, psps) = a_walked_cohort();
    let paths = psps_named_by(&args_over(
        &cohort,
        &psps,
        cohort.directory.path().join("out.toml"),
    ))
    .expect("the directory lists");
    crate::psp::replace_trailer(&paths[0], b"").expect("the tail rewrites");
    let of_another_version = {
        use crate::parameter_estimation::joint::census_file::{
            BYTES_THAT_NAME_THE_VERSION, VERSION,
        };
        let mut census = {
            let mut psp = crate::psp::PspReader::open(&paths[1]).expect("the psp opens");
            psp.trailer().expect("its census reads")
        };
        let word_at = BYTES_THAT_NAME_THE_VERSION - size_of::<u16>();
        census[word_at..BYTES_THAT_NAME_THE_VERSION].copy_from_slice(&(VERSION - 1).to_le_bytes());
        census
    };
    crate::psp::replace_trailer(&paths[1], &of_another_version).expect("the tail rewrites");

    let error = fit_and_assemble(&args_over(
        &cohort,
        &psps,
        cohort.directory.path().join("out.toml"),
    ))
    .expect_err("neither psp carries a census this build reads");

    let EstimateParametersCliError::CohortCannotBeFitted { report } = &error else {
        panic!("two stale psps are a regeneration report: {error:?}");
    };
    let said = report.to_string();
    assert_eq!(report.stale_count(), 2, "both samples are in it: {said}");
    assert!(
        said.contains("one (") && said.contains("two ("),
        "and both are named: {said}",
    );
    assert!(
        said.contains("carries no census") && said.contains("older version"),
        "each with its own cause: {said}",
    );
}

/// **The refusal comes before the reference is opened** (spec §4.2), which is what makes a missing
/// census immediate rather than something a person waits for.
///
/// **How it is shown: the reference is a path that does not exist.** A run that read it first
/// would fail on that instead, and reading a human reference is minutes where judging a psp is one
/// seek and ten bytes.
#[test]
fn a_stale_cohort_is_refused_before_the_reference_is_read() {
    let (cohort, psps) = a_walked_cohort();
    let paths = psps_named_by(&args_over(
        &cohort,
        &psps,
        cohort.directory.path().join("out.toml"),
    ))
    .expect("the directory lists");
    crate::psp::replace_trailer(&paths[1], b"").expect("the tail rewrites");
    let mut args = args_over(&cohort, &psps, cohort.directory.path().join("out.toml"));
    args.reference = cohort.directory.path().join("no-such-reference.fa");

    let error = fit_and_assemble(&args).expect_err("that psp carries no census");

    assert!(
        matches!(
            &error,
            EstimateParametersCliError::CohortCannotBeFitted { .. }
        ),
        "the psps were judged before the reference was opened, and got: {error:?}",
    );
}
