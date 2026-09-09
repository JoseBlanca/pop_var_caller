//! `generate-census` at the command line: what a person may type, what is refused before a psp
//! is read, and what the run says about each file it read.

use super::*;
use clap::Parser;

use crate::ng::psp::PspReader;
use crate::pop_var_caller_exp::cli::{Cli, PopVarCallerExpCommand};
use crate::pop_var_caller_exp::generate_psps::{
    GeneratePspsArgs, census_path_for, psp_path_for, run_generate_psps,
};
use crate::pop_var_caller_exp::test_fixtures::{
    ACohortOnDisk, a_cohort_on_disk, censuses_written_beside_the_psps,
};

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

/// A walked cohort: two psps in a directory of their own, **and no census files** — the walk
/// writes the census into each psp's trailer now (`psp_census_pair.md` §3). Tests that need a
/// census beside a psp make one with `censuses_written_beside_the_psps`.
fn a_walked_cohort() -> (ACohortOnDisk, PathBuf) {
    let (cohort, psps, _) = a_walked_cohort_and_how_it_was_walked();
    (cohort, psps)
}

/// The same, and the arguments it was walked under — which is what a census built afterwards
/// has to be given, since a census is a function of the criteria the ground was cut with.
fn a_walked_cohort_and_how_it_was_walked() -> (ACohortOnDisk, PathBuf, GeneratePspsArgs) {
    let cohort = a_cohort_on_disk();
    let psps = cohort.directory.path().join("psps");
    let walk = GeneratePspsArgs {
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
    };
    run_generate_psps(&walk).expect("the cohort walks into psps");
    (cohort, psps, walk)
}

/// A walked cohort **with a census file beside each psp**, which is the pair the commands that
/// still take `--census` paths are meant to be handed.
fn a_walked_cohort_with_its_censuses() -> (ACohortOnDisk, PathBuf) {
    let (cohort, psps, walk) = a_walked_cohort_and_how_it_was_walked();
    censuses_written_beside_the_psps(&walk, &psps);
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

/// **The evidence this command writes is the evidence the walk wrote, byte for byte.**
///
/// This is the end-to-end form of the producer's own agreement test: `generate-psps` builds a
/// census as it walks and seals it into the psp's trailer, `generate-census` builds one from
/// that psp's records afterwards, and the two encode to the same bytes. It is what says the
/// second route can stand in for the first.
///
/// **What is compared is not the file this command wrote.** That file is decoded and re-encoded
/// with no pileup identity, because the trailer carries none — a census that *is* its psp's
/// trailer has nothing left to pair wrongly with (`psp_census_pair.md` §3) — while a census in a
/// file of its own does. So two things this comparison does *not* cover, and each is covered by
/// a test of its own: the file's own layout, which survives here only as far as `decode_census`
/// preserves it (`census_file.rs`'s `write_census_after_decode_census_returns_the_bytes_it_was_given`),
/// and the identity that was dropped (`the_census_it_writes_names_the_psp_it_read`, below).
#[test]
fn each_census_it_writes_equals_the_one_the_walk_wrote() {
    use crate::ng::parameter_estimation::joint::census_file::{decode_census, write_census};

    let (cohort, psps) = a_walked_cohort();
    let rebuilt = cohort.directory.path().join("rebuilt");

    let report =
        build_every_census(&args_over(&cohort, &psps, rebuilt.clone())).expect("the psps read");

    assert_eq!(report.samples.len(), 2, "one entry a sample");
    for sample in &report.samples {
        let walked = PspReader::open(&psp_path_for(&psps, &sample.sample))
            .expect("the psp opens")
            .trailer()
            .expect("its trailer reads");
        let built = decode_census(&std::fs::read(&sample.census).expect("this run wrote one"))
            .expect("this build's own census");
        let mut without_its_identity = Vec::new();
        write_census(&built.census, None, &mut without_its_identity)
            .expect("a vector accepts every write");
        assert_eq!(
            walked, without_its_identity,
            "{}'s census differs between the walk and the rebuild",
            sample.sample,
        );
    }
}

/// **A census file this command writes names the psp it read — its header and its record
/// count.**
///
/// The identity is what says a census file and the psp beside it are a pair; a census naming a
/// psp that is not there is refused for ever by every freshness check.
///
/// **The record count is the half nothing else checks.** The cohort opener compares only the
/// header digest (`freshness_by_header`), and the command-level test that counted a psp's
/// records was deleted when `generate-psps` stopped writing census files
/// (`psp_census_pair.md` §3) — so without this, a rebuild writing `records: 0` would pass the
/// whole suite.
#[test]
fn the_census_it_writes_names_the_psp_it_read() {
    use crate::ng::parameter_estimation::joint::census_file::{PileupIdentity, open_census};

    let (cohort, psps) = a_walked_cohort();
    let rebuilt = cohort.directory.path().join("rebuilt");

    let report = build_every_census(&args_over(&cohort, &psps, rebuilt)).expect("the psps read");

    for sample in &report.samples {
        let (_evidence, named) = open_census(&sample.census).expect("this build wrote it");
        let named = named.expect("a census file names the psp it was built from");
        let mut psp = PspReader::open(&psp_path_for(&psps, &sample.sample)).expect("the psp opens");
        let records = psp.records().expect("the walk starts").count() as u64;
        let expected = PileupIdentity::of_header(
            &psp.header().encode().expect("the header re-encodes"),
            records,
        );
        assert_eq!(
            named, expected,
            "{}'s census names its own psp — the header in the file, and the records the file \
             holds",
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
    // One run into the psps' own directory, so the run below collides with what it left.
    let (cohort, psps) = a_walked_cohort_with_its_censuses();
    let before = std::fs::read(census_path_for(&psps, "zeta")).expect("the first run wrote it");

    let error = build_every_census(&args_over(&cohort, &psps, psps.clone()))
        .expect_err("a census is already there");

    assert!(
        matches!(&error, GenerateCensusCliError::CensusAlreadyThere { .. }),
        "{error:?}",
    );
    let after = std::fs::read(census_path_for(&psps, "zeta")).expect("it is still there");
    assert_eq!(before, after, "the refused run replaced nothing");
}

/// **`--force` replaces them**, and what it writes is what was there — because this command is
/// a function of the psp it reads and nothing else.
///
/// **Both files here are this command's**, since the walk stopped writing census files
/// (`psp_census_pair.md` §3), so what this asserts is that two runs over one psp agree. The
/// cross-producer guarantee — that the walk and the rebuild agree — is
/// `each_census_it_writes_equals_the_one_the_walk_wrote` above, and only there.
#[test]
fn force_replaces_a_census_that_is_already_there() {
    let (cohort, psps) = a_walked_cohort_with_its_censuses();
    let before = std::fs::read(census_path_for(&psps, "zeta")).expect("the first run wrote it");

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

/// **The censuses this command writes assemble into a cohort a fit can read.**
///
/// This is what steps C1 and C2 exist for, end to end. Every census numbers its read groups from
/// zero, because a walk sees one sample — so before the names were recorded, two censuses of one
/// cohort both claimed read group 0 and were refused as libraries that would be fitted as one.
/// Now they are merged on the `@RG ID`s each declares and renumbered apart.
#[test]
fn the_censuses_it_writes_assemble_into_a_cohort() {
    use crate::ng::parameter_estimation::joint::census::CohortCensusEvidence;
    use crate::ng::parameter_estimation::joint::census_file::open_census;

    let (cohort, psps) = a_walked_cohort();
    let rebuilt = cohort.directory.path().join("rebuilt");
    let report = build_every_census(&args_over(&cohort, &psps, rebuilt)).expect("the psps read");

    let mut samples = Vec::new();
    for entry in &report.samples {
        let (evidence, _) = open_census(&entry.census).expect("this build wrote it");
        assert_eq!(
            evidence.read_groups(),
            vec![crate::ng::types::ReadGroupId(0)],
            "{}'s census numbers its own groups from zero, which is what makes the merge \
             necessary",
            entry.sample,
        );
        samples.push(evidence);
    }
    assert_eq!(samples.len(), 2, "the fixture has two samples");

    let built = CohortCensusEvidence::new(samples)
        .expect("two censuses of one cohort are a cohort, renumbered onto run-wide identifiers");

    assert_eq!(
        built.read_groups().len(),
        2,
        "the two samples' libraries end up under different identifiers, which is what keeps \
         their sequencing-error rates apart",
    );
    let declared: Vec<&str> = built
        .samples()
        .iter()
        .flat_map(|sample| {
            sample
                .declared_read_groups()
                .values()
                .map(|group| group.declared_id.as_str())
        })
        .collect();
    assert_eq!(
        declared.len(),
        2,
        "and each keeps the @RG ID it was declared under: {declared:?}",
    );
}
