//! `regenerate-census` at the command line: what a person may type, what is refused before a psp
//! is rewritten, and what the run says about each file it rebuilt.

use super::*;
use clap::Parser;

use crate::ng::psp::PspReader;
use crate::ng::region_typing::DEFAULT_MAX_STR_LEN;
use crate::ng::region_typing::segment_criteria::{
    DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
};
use crate::pop_var_caller_exp::cli::{Cli, PopVarCallerExpCommand};
use crate::pop_var_caller_exp::generate_psps::{GeneratePspsArgs, psp_path_for, run_generate_psps};
use crate::pop_var_caller_exp::test_fixtures::{
    ACohortOnDisk, a_cohort_on_disk, a_varying_cohort_on_disk,
};

/// Parse an argument vector into this subcommand's arguments, refusing any other subcommand.
fn args_of(argv: &[&str]) -> RegenerateCensusArgs {
    match Cli::parse_from(argv).cmd {
        PopVarCallerExpCommand::RegenerateCensus(args) => args,
        other => panic!("expected regenerate-census, got {other:?}"),
    }
}

/// The shortest run a person can type — **which is now every argument there is**, since the ground
/// and the criteria come from the psps.
fn a_shortest_run() -> Vec<&'static str> {
    vec![
        "pop_var_caller_exp",
        SUBCOMMAND,
        "--reference",
        "ref.fa",
        "--psp",
        "zeta.psp",
    ]
}

/// The walk that put the fixture cohorts into psps, under the default criteria.
fn a_walk_into(psps: &Path, reference: &Path, catalog: &Path, alignments: &[PathBuf]) {
    run_generate_psps(&GeneratePspsArgs {
        reference: reference.to_path_buf(),
        catalog: Some(catalog.to_path_buf()),
        alignments: alignments.to_vec(),
        output_dir: psps.to_path_buf(),
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
}

/// The plain cohort, walked: two psps in a directory of their own, each carrying the census its
/// walk wrote into its trailer. **Its second sample has no reads at all**, which is the shape a
/// producer that skipped empty sections would get wrong.
fn a_walked_cohort() -> (ACohortOnDisk, PathBuf) {
    let cohort = a_cohort_on_disk();
    let psps = cohort.directory.path().join("psps");
    a_walk_into(
        &psps,
        &cohort.reference,
        &cohort.catalog,
        &cohort.alignments,
    );
    (cohort, psps)
}

/// This command's arguments over a walked cohort.
fn args_over(reference: &Path, catalog: &Path, psps: &Path) -> RegenerateCensusArgs {
    RegenerateCensusArgs {
        reference: reference.to_path_buf(),
        catalog: Some(catalog.to_path_buf()),
        psps: vec![psps.to_path_buf()],
    }
}

/// Each sample's psp path beside the bytes of the census in its trailer.
fn the_trailers_in(psps: &Path, samples: &[&str]) -> Vec<(PathBuf, Vec<u8>)> {
    samples
        .iter()
        .map(|sample| {
            let path = psp_path_for(psps, sample);
            let trailer = PspReader::open(&path)
                .expect("the psp opens")
                .trailer()
                .expect("its trailer reads");
            (path, trailer)
        })
        .collect()
}

/// **The command is spelled by the constant the library's refusals name.**
///
/// Two messages in `ng::run` tell a person to run this — the report `estimate-parameters` refuses
/// a stale cohort with, and the fit's own backstop — and until this step they named a command that
/// did not exist. **What this pins is that the name they print is the one clap answers to**, in
/// both directions: the constant is this module's `SUBCOMMAND`, and parsing that word gives these
/// arguments.
#[test]
fn the_subcommand_is_spelled_regenerate_census() {
    let args = args_of(&a_shortest_run());
    assert_eq!(args.reference, PathBuf::from("ref.fa"));
    assert_eq!(args.psps, vec![PathBuf::from("zeta.psp")]);
    assert_eq!(
        SUBCOMMAND, THE_COMMAND_THAT_REBUILDS_A_CENSUS,
        "the name the library's refusals print is this command's own",
    );
    assert!(
        Cli::try_parse_from(["pop_var_caller_exp", SUBCOMMAND, "--help"])
            .expect_err("--help exits")
            .to_string()
            .contains(SUBCOMMAND),
        "and clap answers to it",
    );
}

/// **Nothing about the ground or about what counts as a repeat can be typed here** (spec §6, §8).
///
/// The psps record the ground they were walked over and the criteria they were cut with, and the
/// cohort is refused unless they agree, so every one of these flags could only say something the
/// files already say — and a person who typed one and was ignored would have a census rebuilt
/// under settings they thought they had changed. **`--output-dir` and `--force` go for a
/// different reason**: the census goes into the psp, and replacing it is the whole command.
#[test]
fn nothing_about_the_ground_or_the_criteria_can_be_typed_here() {
    for flag in [
        "--regions",
        "--min-copies",
        "--min-period",
        "--max-period",
        "--max-str-len",
        "--min-purity",
        "--output-dir",
    ] {
        let mut argv = a_shortest_run();
        argv.extend([flag, "whatever"]);
        let refused = Cli::try_parse_from(argv).expect_err("{flag} is not a flag here");
        assert!(
            refused.to_string().contains(flag),
            "{flag} was accepted, and got: {refused}",
        );
    }
    let mut argv = a_shortest_run();
    argv.push("--force");
    let refused = Cli::try_parse_from(argv).expect_err("--force is not a flag here");
    assert!(
        refused.to_string().contains("--force"),
        "and got: {refused}"
    );
}

/// **The census this command writes into a psp is the census the walk wrote there, byte for
/// byte** — on the cohort that carries a repeat tract and three read groups.
///
/// This is the parity oracle in the shape the trailer gives it (`psp_census_pair.md` §11):
/// `generate-psps` builds a census as it walks and seals it into the psp, this command builds one
/// from that psp's records afterwards and replaces the trailer with it, and the bytes are the same
/// ones. It is what says the repair can stand in for the walk.
///
/// **What this covers that `census_from_psp`'s own parity tests cannot: the selection is derived
/// twice.** Those hand one plan to both producers, so a disagreement about how a plan is *built*
/// is invisible there. Here the walk and this command each assemble their own segmentation and
/// their own `CensusPlan` — the walk from its flags, this command from the psp headers — and a
/// divergence changes which loci are kept and therefore the bytes.
///
/// **The fixture is the one with a tract in it**, because the half of a census keyed by stratum as
/// well as by read group is empty on both sides in a cohort whose selection keeps no tract, and
/// then a producer that dropped every tract read would pass. The assertions below say what the
/// comparison was actually over.
///
/// **What no parity oracle can catch**, here or in the sibling: a defect both producers share.
/// They build their writer through one `CensusPlan::writer_for`, so a fault in `CensusWriter` or
/// in `write_census` corrupts both sides identically. "Byte for byte" is a statement about the psp
/// being complete, not about the census being right.
#[test]
fn the_census_it_writes_is_the_one_the_walk_wrote_on_a_cohort_with_a_repeat_tract() {
    let cohort = a_varying_cohort_on_disk();
    let psps = cohort.directory.path().join("psps");
    a_walk_into(
        &psps,
        &cohort.reference,
        &cohort.catalog,
        &cohort.alignments,
    );
    let walked = the_trailers_in(&psps, &["one", "two"]);

    let report = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect("the psps read");

    assert_eq!(report.samples.len(), 2, "one entry a sample");
    for (path, before) in &walked {
        let after = PspReader::open(path)
            .expect("the psp still opens")
            .trailer()
            .expect("its new trailer reads");
        assert_eq!(
            before,
            &after,
            "{} carries a different census after the rebuild",
            path.display(),
        );
    }

    // **What the comparison was actually over.** Both halves of a census can be present and
    // empty, and two empty halves agree — so a fixture that stopped putting reads on kept tracts
    // would leave the comparison above passing while covering nothing this step is for.
    assert!(
        report
            .samples
            .iter()
            .any(|sample| sample.tally.tracts_with_reads > 0),
        "no sample has a read at a kept repeat tract, so the tract half of the comparison is two \
         never-walked section sets agreeing with each other: {:?}",
        report.samples.iter().map(|it| it.tally).collect::<Vec<_>>(),
    );
    let declared: Vec<usize> = ["one", "two"]
        .iter()
        .map(|sample| {
            PspReader::open(&psp_path_for(&psps, sample))
                .expect("the psp opens")
                .header()
                .read_groups
                .len()
        })
        .collect();
    assert_eq!(
        declared,
        vec![2, 1],
        "this fixture's two samples declare two read groups and one, and a census keys its \
         sections by group — so the comparison spans more than one key a sample",
    );
}

/// **And on the plain cohort, whose second sample carries no reads at all** — the case where every
/// kept position of a census is a zero, which is a denominator the fit needs and the shape a
/// producer that skipped empty sections would get wrong.
#[test]
fn the_census_it_writes_is_the_one_the_walk_wrote_for_a_sample_with_no_reads() {
    let (cohort, psps) = a_walked_cohort();
    let walked = the_trailers_in(&psps, &["alpha", "zeta"]);

    let report = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect("the psps read");

    for (path, before) in &walked {
        let after = PspReader::open(path)
            .expect("the psp still opens")
            .trailer()
            .expect("its new trailer reads");
        assert_eq!(before, &after, "{} changed", path.display());
    }

    // **Both shapes have to be present, or this is not the test it is named for.** One sample with
    // nothing is what it covers; one sample with something is what stops it passing on two
    // producers that both read no records at all.
    let tallies: Vec<_> = report.samples.iter().map(|it| it.tally).collect();
    assert!(
        tallies.iter().any(|tally| tally.contributes_nothing()),
        "the fixture's empty sample stopped being empty: {tallies:?}",
    );
    assert!(
        tallies.iter().any(|tally| !tally.contributes_nothing()),
        "every sample is empty, so the two producers agreed about nothing: {tallies:?}",
    );
}

/// **The census written into a trailer names no psp**, where a census in a file of its own names
/// the one it was built from.
///
/// A census that *is* its psp's trailer has nothing left to pair wrongly with
/// (`psp_census_pair.md` §3.1), and writing an identity into it would make this command's output
/// differ from the walk's for a field neither needs — which the byte comparison above would catch,
/// but as *the censuses differ* rather than as the one field that does.
#[test]
fn the_census_it_writes_into_the_trailer_names_no_psp() {
    use crate::ng::parameter_estimation::joint::census_file::decode_census;

    let (cohort, psps) = a_walked_cohort();

    let report = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect("the psps read");

    for sample in &report.samples {
        let trailer = PspReader::open(&sample.psp)
            .expect("the psp opens")
            .trailer()
            .expect("its trailer reads");
        let census = decode_census(&trailer).expect("this build's own census");
        assert!(
            census.pileup.is_none(),
            "{}'s census is its psp's trailer, so it names no psp",
            sample.sample,
        );
    }
}

/// **Only the trailer is rewritten**: the header the psp declares and the records it holds are
/// what they were.
///
/// The whole-file oracle is plan step D4 — a copied psp regenerated and `cmp`-identical, header,
/// blocks, index and all. This is the half that can be asserted through the reader: a rebuild that
/// re-encoded the header, or that dropped a block, would still pass a comparison of trailers.
#[test]
fn the_psps_header_and_records_are_untouched() {
    let (cohort, psps) = a_walked_cohort();
    let before: Vec<(String, u64)> = ["alpha", "zeta"]
        .iter()
        .map(|sample| {
            let mut psp = PspReader::open(&psp_path_for(&psps, sample)).expect("the psp opens");
            let records = psp.records().expect("the walk starts").count() as u64;
            let header = format!("{:?}", psp.header());
            (header, records)
        })
        .collect();

    regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect("the psps read");

    let after: Vec<(String, u64)> = ["alpha", "zeta"]
        .iter()
        .map(|sample| {
            let mut psp = PspReader::open(&psp_path_for(&psps, sample)).expect("the psp opens");
            let records = psp.records().expect("the walk starts").count() as u64;
            let header = format!("{:?}", psp.header());
            (header, records)
        })
        .collect();

    assert_eq!(before, after, "the rebuild rewrote more than the trailer");
}

/// **A sample whose census holds no read is named, not omitted.**
///
/// The plain fixture's second sample carries no reads at all, so every kept position it has is a
/// zero — which is the denominator a fit needs and not an error. A run that left it out of its
/// report would leave somebody hunting for a psp rewritten exactly as asked.
#[test]
fn a_sample_with_no_reads_is_named_as_contributing_nothing() {
    let (cohort, psps) = a_walked_cohort();

    let report = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect("the psps read");

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

/// **Each sample's line names the psp it read and rewrote, and how much of it the census is** —
/// one file now, where the old command wrote a second one beside it.
#[test]
fn each_line_names_the_psp_and_what_went_into_its_census() {
    let (cohort, psps) = a_walked_cohort();

    let report = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect("the psps read");

    let zeta = report
        .samples
        .iter()
        .find(|sample| sample.sample == "zeta")
        .expect("the sample with reads is in the report");
    let line = zeta.line();
    assert!(line.contains("zeta.psp"), "and got: {line}");
    assert!(
        line.contains(&format!("{} stored loci read", zeta.records)),
        "and got: {line}",
    );
    assert!(
        line.contains(&format!(
            "census {} bytes written into its trailer",
            zeta.census_bytes
        )),
        "and got: {line}",
    );
    assert!(
        !line.contains(".census"),
        "no census file is written any more, so no line may name one: {line}",
    );
    assert_eq!(
        zeta.psp,
        psp_path_for(&psps, "zeta"),
        "the psp named is the one this sample's census came from",
    );
}

/// **A directory of psps is expanded in name order**, so two runs naming one directory rebuild the
/// same cohort in the same order however the filesystem answers.
#[test]
fn a_directory_contributes_every_psp_inside_it() {
    let (cohort, psps) = a_walked_cohort();
    let args = args_over(&cohort.reference, &cohort.catalog, &psps);

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
    let mut args = args_over(&cohort.reference, &cohort.catalog, &psps);
    args.psps = vec![empty.clone()];

    let error = psps_named_by(&args).expect_err("an empty directory is not a cohort");

    assert!(
        matches!(&error, RegenerateCensusCliError::NoPspsInDirectory { path } if path == &empty),
        "{error:?}",
    );
}

/// **What this command writes is what a fit reads**: the rebuilt trailers assemble into a cohort.
///
/// Every census numbers its read groups from zero, because a walk sees one sample — so two
/// censuses of one cohort both claim read group 0 and are merged on the `@RG ID`s they declare.
/// This is the end of the repair: the psps a fit refused are read back by the fit's own reader.
#[test]
fn the_censuses_it_writes_assemble_into_a_cohort() {
    use crate::ng::run::every_census_in_the_cohorts_psps;

    let (cohort, psps) = a_walked_cohort();
    let args = args_over(&cohort.reference, &cohort.catalog, &psps);
    let paths = psps_named_by(&args).expect("the directory lists");
    regenerate_every_census(&args).expect("the psps read");

    let open = OpenPspCohort::open(&paths).expect("the rebuilt psps are one cohort");
    let evidence = every_census_in_the_cohorts_psps(&open)
        .expect("each rebuilt trailer is a census, and the two are one cohort");

    assert_eq!(evidence.len(), 2, "the fixture has two samples");
    assert_eq!(
        evidence.read_groups().len(),
        2,
        "this fixture's two samples declare one library each, and the two end up under different \
         identifiers — which is what keeps their sequencing-error rates apart, and what every \
         census numbering its own groups from zero would otherwise lose",
    );
}

/// **A reference the psps were not walked against is refused, and no psp is rewritten.**
///
/// This is the failure this command must not commit: rebuilding a census against another reference
/// keeps other positions and records that reference's digest, so every later fit refuses the
/// cohort — after this run has spent a quarter of an hour a sample producing it. Before this step
/// the command it replaces made no such check.
///
/// The reference here is the cohort's own with its first base changed, so its contig table is the
/// psps' and only the bases differ.
#[test]
fn a_reference_the_psps_were_not_walked_against_is_refused_before_anything_is_rewritten() {
    use crate::pop_var_caller_exp::test_fixtures::{VARYING_CONTIG, the_varying_cohorts_reference};

    let cohort = a_varying_cohort_on_disk();
    let psps = cohort.directory.path().join("psps");
    a_walk_into(
        &psps,
        &cohort.reference,
        &cohort.catalog,
        &cohort.alignments,
    );
    let walked = the_trailers_in(&psps, &["one", "two"]);
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

    let error = regenerate_every_census(&args_over(&another, &cohort.catalog, &psps))
        .expect_err("the psps were walked against another reference");

    assert!(
        matches!(
            &error,
            RegenerateCensusCliError::WalkedAgainstAnotherReference { path, .. } if path == &another
        ),
        "{error:?}",
    );
    let said = crate::error_render::format_error_chain(&error);
    assert!(
        said.contains("they name ref.fa, so run this again with that one"),
        "and it names the reference the psps do: {said}",
    );
    for (path, before) in &walked {
        let after = PspReader::open(path)
            .expect("the psp opens")
            .trailer()
            .expect("its trailer reads");
        assert_eq!(
            before,
            &after,
            "{} was rewritten against a reference it was not walked against",
            path.display(),
        );
    }
}

/// **A catalog the psps were not walked with is refused, and no psp is rewritten**, for the reason
/// a wrong reference is: the census records what the catalog was built at.
///
/// The catalog built here is on the cohort's own reference and holds tracts of twenty copies and
/// up, where the psps' was built at the default floors — the same assembly, another file.
#[test]
fn a_catalog_the_psps_were_not_walked_with_is_refused_before_anything_is_rewritten() {
    use crate::ng::reference_info::{ReferenceSource, read_reference_info_observing};
    use crate::ng::region_typing::segment_criteria::SsrSegmentCriteria;
    use crate::ng::repeat_catalog::{RepeatCatalogBuilder, StrRepeatCriteria};
    use crate::ng::tandem_repeat::ScanParams;

    let (cohort, psps) = a_walked_cohort();
    let walked = the_trailers_in(&psps, &["alpha", "zeta"]);
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

    let error = regenerate_every_census(&args_over(&cohort.reference, &coarse, &psps))
        .expect_err("the psps were walked with another catalog");

    assert!(
        matches!(
            &error,
            RegenerateCensusCliError::WalkedWithAnotherCatalog { path, .. } if path == &coarse
        ),
        "{error:?}",
    );
    let said = crate::error_render::format_error_chain(&error);
    assert!(
        said.contains("it was built under other repeat criteria")
            && said.contains("run this again with --catalog naming the one they were"),
        "it says what differs and what to do: {said}",
    );
    for (path, before) in &walked {
        let after = PspReader::open(path)
            .expect("the psp opens")
            .trailer()
            .expect("its trailer reads");
        assert_eq!(
            before,
            &after,
            "{} was rewritten against a catalog it was not walked with",
            path.display(),
        );
    }
}

/// **A psp that will not take its new census says which sample, which file, and what state that
/// file is in** (spec §8).
///
/// The state is what the person acts on: a psp left byte for byte what it was needs this run
/// again, and one already cut back to its trailer has to be rewritten before anything can read
/// it. **The message is pinned rather than the filesystem failure that produces it**, because
/// this suite runs as root inside the container, where a read-only file is not read-only — so a
/// test that made one and expected a refusal would pass for the wrong reason on one machine and
/// fail on another.
#[test]
fn a_psp_that_will_not_take_its_census_names_the_file_and_its_state() {
    use crate::ng::psp::{FileAfterAFailedReplacement, PspWriteError};

    for (state, expected) in [
        (
            FileAfterAFailedReplacement::Unchanged,
            "the file is exactly as it was, and the call can be made again",
        ),
        (
            FileAfterAFailedReplacement::Torn,
            "the file is cut short and has no footer; the sample has to be written again",
        ),
    ] {
        let error = RegenerateCensusCliError::TrailerNotReplaced {
            sample: "zeta".to_string(),
            psp: PathBuf::from("/psps/zeta.psp"),
            source: Box::new(TrailerReplacementFailure {
                file: state,
                source: PspWriteError::Io {
                    path: PathBuf::from("/psps/zeta.psp"),
                    while_doing: "writing the trailer",
                    source: std::io::Error::from(std::io::ErrorKind::PermissionDenied),
                },
            }),
        };

        let said = crate::error_render::format_error_chain(&error);
        assert!(
            said.contains("replacing zeta's census in /psps/zeta.psp"),
            "which sample and which file: {said}",
        );
        assert!(said.contains(expected), "and what state it is in: {said}");
    }
}
