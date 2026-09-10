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
/// walk wrote into its trailer. **The sample `alpha` has no reads at all**, which is the shape a
/// producer that skipped empty sections would get wrong — and it is the *first* psp this command
/// reads, since psps arrive in name order.
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

/// Empty every psp's trailer, after `walked` has captured what the walk had written there.
///
/// **What it is for: making the bytes a comparison asserts come from the run under test.** Since
/// plan step D2 a psp whose census is the one this run would write is skipped and its records are
/// never read, so a test that walks a cohort and regenerates it compares the walk's own bytes with
/// themselves and passes whatever this command does — which is the defect D1's review found in
/// every test in this file. Emptying the trailers first makes each psp owed a rebuild.
fn empty_every_trailer(walked: &[(PathBuf, Vec<u8>)]) {
    for (path, _) in walked {
        crate::ng::psp::replace_trailer(path, b"").expect("the tail rewrites");
    }
}

/// **Overwrite 32 bytes just below where the psp's index begins, which is inside its last block.**
///
/// **What it is for: a psp whose records cannot be read, while everything the open pass reads is
/// intact.** `PspReader::open` reads the header, the index and the footer and checks the index's
/// own checksum, so corrupting any of those refuses the file instead — and the trailer is above the
/// index, so a census in it still decodes. What is left is a psp that opens, judges and skips
/// cleanly, and whose rebuild fails the moment it reads a block.
///
/// **Permissions would be the plainer way to make a psp unusable and cannot be used here**: this
/// suite runs as root inside the dev container, where a read-only file is not read-only.
fn corrupt_a_block_of(psp: &Path) {
    use std::io::{Seek, SeekFrom, Write};

    let index_begins_at = PspReader::open(psp)
        .expect("the psp opens")
        .footer()
        .index_offset;
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .open(psp)
        .expect("the scratch dir is ours");
    file.seek(SeekFrom::Start(index_begins_at - 32))
        .expect("the blocks end where the index begins");
    file.write_all(&[0xFF; 32]).expect("32 bytes of nonsense");
    file.sync_all().expect("the bytes land");
}

/// **The command is spelled by the constant the library's refusals name.**
///
/// The refusals in `ng::run` tell a person to run this command, and until this step they named one
/// that did not exist. **What this pins is that the word they print is the one clap answers to**:
/// the argument list parsed here is built from that constant, and clap answers `--help` for it
/// rather than reporting a subcommand it does not know.
#[test]
fn the_subcommand_is_spelled_regenerate_census() {
    let args = args_of(&a_shortest_run());
    assert_eq!(args.reference, PathBuf::from("ref.fa"));
    assert_eq!(args.psps, vec![PathBuf::from("zeta.psp")]);
    // **`SUBCOMMAND` is *defined as* that constant, so comparing them is `X == X`.** What
    // carries the claim is the parse above, whose argument list is built from the constant the
    // library's refusals print: if clap did not answer to that word, `args_of` would not return.
    let refused = Cli::try_parse_from([
        "pop_var_caller_exp",
        THE_COMMAND_THAT_REBUILDS_A_CENSUS,
        "--help",
    ])
    .expect_err("--help exits rather than running");
    assert_eq!(
        refused.kind(),
        clap::error::ErrorKind::DisplayHelp,
        "clap printed help for this name, rather than failing to recognise it — which its \
         unrecognised-subcommand error would also echo back",
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
        let refused = Cli::try_parse_from(argv)
            .expect_err("this command takes no such flag")
            .to_string();
        assert!(
            refused.contains(flag),
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
    empty_every_trailer(&walked);

    let report = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect("the psps read");

    assert_eq!(report.samples.len(), 2, "one entry a sample");
    assert!(
        report.skipped.is_empty(),
        "both psps were emptied, so neither can have needed nothing: {:?}",
        report.skipped,
    );
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
    empty_every_trailer(&walked);

    let report = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect("the psps read");

    assert_eq!(
        report.samples.len(),
        2,
        "both psps were emptied, so both are owed a rebuild"
    );
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

/// **One psp owed a rebuild, one that needs nothing: the first is written and the second is
/// skipped** — plan step D2's own case, and the test that says the write happens at all.
///
/// **Why it is both at once.** Emptying one trailer makes that psp's new bytes provably this run's,
/// which is what D1's review found nothing asserted: every comparison read a trailer the walk had
/// already put there, so all fourteen tests passed with the write removed. Seven of them empty
/// their trailers now and would catch that too; **what only this one covers is the pair** — one psp
/// rebuilt and one skipped in the same run, which is the case spec §8 asks for.
#[test]
fn a_psp_with_no_census_is_rebuilt_and_the_fresh_one_is_skipped() {
    let (cohort, psps) = a_walked_cohort();
    let walked = the_trailers_in(&psps, &["alpha", "zeta"]);
    let emptied = psp_path_for(&psps, "alpha");
    crate::ng::psp::replace_trailer(&emptied, b"").expect("the tail rewrites");
    assert_eq!(
        PspReader::open(&emptied)
            .expect("the psp opens")
            .footer()
            .trailer_bytes,
        0,
        "this fixture starts with one psp carrying no census at all",
    );

    let report = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect("the psps read");

    assert_eq!(
        report
            .samples
            .iter()
            .map(|sample| sample.sample.as_str())
            .collect::<Vec<_>>(),
        vec!["alpha"],
        "the emptied psp is the one rebuilt",
    );
    assert_eq!(
        report
            .skipped
            .iter()
            .map(|it| it.sample.as_str())
            .collect::<Vec<_>>(),
        vec!["zeta"],
        "and the other needed nothing",
    );
    for (path, before) in &walked {
        let after = PspReader::open(path)
            .expect("the psp opens")
            .trailer()
            .expect("its trailer reads");
        assert_eq!(
            before,
            &after,
            "{} does not carry the census its walk wrote",
            path.display(),
        );
    }
    // **The one run that renders every singular the first line has**, since it is the only shape
    // with exactly one of each — and each of those three pieces is a `match` a test can otherwise
    // leave to whichever arm the fixtures happen to take.
    let said = report.lines();
    assert!(
        said[0].contains("regenerated 1 census over")
            && said[0].contains("and skipped 1 psp that needed nothing"),
        "one of each, in the singular: {said:?}",
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
    empty_every_trailer(&the_trailers_in(&psps, &["alpha", "zeta"]));

    let report = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect("the psps read");

    assert_eq!(
        report.samples.len(),
        2,
        "the censuses asserted below are this run's"
    );
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
/// The whole-file comparison is `a_regenerated_psp_is_the_walked_one_byte_for_byte`, below. This is
/// the half that can be asserted through the reader, and it is what catches a rebuild that dropped
/// a block: the record count comes from reading them.
#[test]
fn the_psps_header_and_records_are_untouched() {
    let (cohort, psps) = a_walked_cohort();
    empty_every_trailer(&the_trailers_in(&psps, &["alpha", "zeta"]));
    let before: Vec<(String, u64)> = ["alpha", "zeta"]
        .iter()
        .map(|sample| {
            let mut psp = PspReader::open(&psp_path_for(&psps, sample)).expect("the psp opens");
            let records = psp.records().expect("the walk starts").count() as u64;
            let header = format!("{:?}", psp.header());
            (header, records)
        })
        .collect();

    let report = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect("the psps read");
    assert_eq!(
        report.samples.len(),
        2,
        "a psp this run skipped would be untouched for a reason this test is not about",
    );

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
    empty_every_trailer(&the_trailers_in(&psps, &["alpha", "zeta"]));

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
    empty_every_trailer(&the_trailers_in(&psps, &["alpha", "zeta"]));

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
    assert!(
        !report.lines()[0].contains("skipped"),
        "every psp here was owed a rebuild, so the first line has nothing to say about \
         skipping: {:?}",
        report.lines(),
    );
    assert_eq!(
        zeta.psp,
        psp_path_for(&psps, "zeta"),
        "the psp named is the one this sample's census came from",
    );
    // **The count is the psp's own, and nothing else here checks it.** The line is built from the
    // same field it is compared against, so a rebuild reporting no records at all would pass every
    // other test in this file — which the command this replaced said in as many words.
    let counted = {
        let mut psp = PspReader::open(&psp_path_for(&psps, "zeta")).expect("the psp opens");
        psp.records().expect("the walk starts").count() as u64
    };
    assert_eq!(
        zeta.records, counted,
        "the line reports {} stored loci where the psp holds {counted}",
        zeta.records,
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
    empty_every_trailer(&the_trailers_in(&psps, &["alpha", "zeta"]));
    let report = regenerate_every_census(&args).expect("the psps read");
    assert_eq!(
        report.samples.len(),
        2,
        "the censuses assembled below are the ones this run wrote",
    );

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

/// **A cohort that needs nothing is skipped whole, and its records are never read** (spec §8).
///
/// Freshness here is all three of spec §4.2's causes: the two the psps' heads answer, and — since
/// this command rebuilds the selection anyway — whether each census recorded the settings this run
/// records under. A freshly walked cohort passes all three, so there is nothing to do.
///
/// **What is shown, and it is narrower than "its records are never read".** One psp's *blocks* are
/// corrupted before the run: 32 bytes overwritten just below where the index begins, which is
/// inside the last block ([`corrupt_a_block_of`]). A rebuild reads every record and would fail on
/// them; this run succeeds, so **no rebuild pass happened for that psp** — which is the pass that
/// costs a quarter of an hour a sample at 50× human.
///
/// **What it does not pin: that nothing read those blocks at all.** Measured by mutation — code
/// added to the skip arm that reads the psp's records and discards the result passes every test in
/// this file, because a discarded failure is invisible and nothing counts block bytes the way
/// [`trailer_bytes_read`](crate::ng::psp::trailer_bytes_read) counts trailer bytes.
///
/// **The control is the second half of the test**: with that psp's trailer emptied, the same
/// corrupted file is owed a rebuild, and then the run does fail on it — which is what says the
/// corruption was fatal to a record pass rather than harmless.
#[test]
fn a_cohort_that_needs_nothing_is_skipped_whole_and_its_records_are_not_read() {
    let (cohort, psps) = a_walked_cohort();
    let corrupted = psp_path_for(&psps, "zeta");
    corrupt_a_block_of(&corrupted);

    let report = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect("every census is the one this run would write, so no psp is read past its census");

    assert!(
        report.samples.is_empty(),
        "nothing was rebuilt: {:?}",
        report
            .samples
            .iter()
            .map(|it| &it.sample)
            .collect::<Vec<_>>(),
    );
    assert_eq!(
        report
            .skipped
            .iter()
            .map(|it| it.sample.as_str())
            .collect::<Vec<_>>(),
        vec!["alpha", "zeta"],
        "both psps needed nothing, in the order they were read",
    );
    let said = report.lines();
    assert!(
        said[0].contains("regenerated 0 censuses")
            && said[0].contains("skipped 2 psps that needed nothing"),
        "and the run says so: {said:?}",
    );
    assert!(
        said.iter()
            .any(|line| line.contains("zeta: skipped, its census is the one this run would write")),
        "a line a sample, either way: {said:?}",
    );

    // **The control.** With its census gone, the same corrupted psp is owed a rebuild — and now
    // the run fails on it, which is what makes the pass above a statement about records not read
    // rather than about a psp that would have read cleanly anyway.
    crate::ng::psp::replace_trailer(&corrupted, b"").expect("the tail rewrites");

    let error = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect_err("this psp's records are nonsense, and now they have to be read");

    assert!(
        matches!(
            &error,
            RegenerateCensusCliError::Build { sample, psp, .. } if sample == "zeta" && psp == &corrupted
        ),
        "and it names the sample and the file: {error:?}",
    );
}

/// **A cohort whose censuses were recorded under another selection is rebuilt whole, without being
/// told to** (spec §8).
///
/// This is the case a changed build makes: the selection's constants move, every census on disk
/// records the old ones, and no flag says so. The fixture reaches it the way C5's tests do — each
/// census rebuilt from its own psp under a selection keeping a handful of positions instead of the
/// shipped budget — and what the command must do is rebuild every one and leave the trailers as
/// the walk had them.
/// **Both budgets are run, and the first is the one that matters.** At half the shipped budget
/// this fixture's short contig keeps every ordinary position either way, so the set of positions is
/// the same and only the recorded settings differ — which is the case a check comparing the kept
/// positions alone would skip, and the case C5 measured on the fit. At three the set differs too.
#[test]
fn a_cohort_recorded_under_another_selection_is_rebuilt_whole() {
    use crate::ng::parameter_estimation::joint::census_file::write_census;
    use crate::ng::run::census_from_psp;
    use crate::ng::run::test_fixtures::a_census_plan_over_selecting;

    for budget in [CensusSelection::SHIPPED.generic_target / 2, 3] {
        let (cohort, psps) = a_walked_cohort();
        let walked = the_trailers_in(&psps, &["alpha", "zeta"]);
        let (segmentation, under_another_budget) =
            a_census_plan_over_selecting(&cohort.reference, &cohort.catalog, budget);
        for (path, _) in &walked {
            let rebuilt =
                census_from_psp(path, &under_another_budget, &segmentation).expect("the psp reads");
            let mut bytes = Vec::new();
            write_census(&rebuilt.evidence, None, &mut bytes).expect("the census encodes");
            crate::ng::psp::replace_trailer(path, &bytes).expect("the tail rewrites");
        }

        let report = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
            .expect("the psps read");

        assert_eq!(
            report.samples.len(),
            2,
            "at a budget of {budget}, every census recorded another selection, so every one is \
             owed a rebuild",
        );
        assert!(
            report.skipped.is_empty(),
            "at a budget of {budget}, none needed nothing: {:?}",
            report.skipped,
        );
        for (path, before) in &walked {
            let after = PspReader::open(path)
                .expect("the psp opens")
                .trailer()
                .expect("its trailer reads");
            assert_eq!(
                before,
                &after,
                "at a budget of {budget}, {} did not come back to the census its walk wrote",
                path.display(),
            );
        }
    }
}

/// **A census damaged past its version word is rebuilt, not skipped and not refused.**
///
/// The head's two reads cannot tell it from a whole census — the magic and the version word are
/// this build's — and no cheap read can. What settles it is the census reader, and regenerating
/// rewrites exactly the bytes that are damaged, so this command's answer is to rebuild rather than
/// to stop (`CensusVerdict`'s own note on what is not a verdict).
///
/// **The premise is asserted rather than assumed**: the head still says fresh, and the census still
/// will not read. Without both, a test that saw the sample rebuilt would not know which of the two
/// judgements had done it.
#[test]
fn a_census_damaged_past_its_version_word_is_rebuilt() {
    use std::io::{Seek, SeekFrom, Write};

    use crate::ng::parameter_estimation::joint::census_file::BYTES_THAT_NAME_THE_VERSION;
    use crate::ng::run::{
        the_census_in_a_psp, what_the_footer_and_the_trailers_head_say_about_a_census,
    };

    let (cohort, psps) = a_walked_cohort();
    let walked = the_trailers_in(&psps, &["alpha", "zeta"]);
    let damaged = psp_path_for(&psps, "zeta");
    let census_begins_at = PspReader::open(&damaged)
        .expect("the psp opens")
        .footer()
        .trailer_offset;
    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .open(&damaged)
        .expect("the scratch dir is ours");
    file.seek(SeekFrom::Start(
        census_begins_at + BYTES_THAT_NAME_THE_VERSION as u64,
    ))
    .expect("past the magic and the version word");
    file.write_all(&[0xFF; 16]).expect("16 bytes of nonsense");
    file.sync_all().expect("the bytes land");
    drop(file);

    let mut psp = PspReader::open(&damaged).expect("the psp opens");
    assert_eq!(
        what_the_footer_and_the_trailers_head_say_about_a_census(&mut psp).expect("its head reads"),
        crate::ng::run::CensusVerdict::Fresh,
        "the damage is past the version word, so the head cannot see it",
    );
    assert!(
        the_census_in_a_psp(&damaged, &psp).is_err(),
        "and the census reader can",
    );
    drop(psp);

    let report = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect("a damaged census is rebuilt rather than refused");

    assert_eq!(
        report
            .samples
            .iter()
            .map(|sample| sample.sample.as_str())
            .collect::<Vec<_>>(),
        vec!["zeta"],
        "the damaged psp is the one rebuilt",
    );
    for (path, before) in &walked {
        let after = PspReader::open(path)
            .expect("the psp opens")
            .trailer()
            .expect("its trailer reads");
        assert_eq!(
            before,
            &after,
            "{} does not carry the census its walk wrote",
            path.display(),
        );
    }
}

/// **The cohort is opened before the reference is read** — spec §8's own order, and the finding
/// carried out of D1's review.
///
/// **How it is shown: both are wrong at once.** The `--psp` names a file that is not a psp and the
/// `--reference` names a file that is not there, so whichever is read first is the one that
/// refuses. Before this step it was the reference, and a person who mistyped a psp path paid a
/// full reference read — minutes on a human FASTA — to be told about it.
#[test]
fn a_cohort_that_will_not_open_is_refused_before_the_reference_is_read() {
    let (cohort, _psps) = a_walked_cohort();
    let not_a_psp = cohort.directory.path().join("notes.psp");
    std::fs::write(&not_a_psp, b"not a psp").expect("the scratch dir is ours");

    let error = regenerate_every_census(&RegenerateCensusArgs {
        reference: cohort.directory.path().join("no-such-reference.fa"),
        catalog: Some(cohort.catalog.clone()),
        psps: vec![not_a_psp],
    })
    .expect_err("neither the psp nor the reference is usable");

    assert!(
        matches!(&error, RegenerateCensusCliError::Cohort { .. }),
        "the cohort is opened first, so it is the cohort that refuses: {error:?}",
    );
}

/// **A run stopped part-way does only what is left when it runs again** (spec §8, §10).
///
/// This is what the skip rule is for. The first run rebuilds one psp and then fails on the second,
/// leaving the cohort half repaired; the second run skips the psp already done and rebuilds the one
/// still owed, so nobody has to remember which those were.
///
/// **The property nothing else here covers: the second run skips a census *this command* wrote.**
/// Every other skip test skips the bytes the walk wrote. If what `write_census` records and what
/// `CensusPlan::recording_terms` is compared against came apart, this command would rebuild its own
/// output for ever and a stopped run would never converge, however many times it was run.
///
/// **How the failure is made, and why it lands after the first psp.** The psps are read in name
/// order, so `alpha` is rebuilt before `zeta` is reached, and `zeta`'s blocks are corrupted, so its
/// rebuild — the only pass that reads records — fails. **Permissions would have been the plainer
/// way and cannot be used**: this suite runs as root inside the container, where a read-only file
/// is not read-only. What the corruption stands in for is a read that failed and then stopped
/// failing, which is why the bytes are put back between the runs.
///
/// **Two psps rather than the plan's three**, because the fixture cohorts have two samples and the
/// property needs one of each: a psp skipped and a psp rebuilt in the second run. **What two cannot
/// see** is a run that pressed on past the failure and rebuilt *later* samples before returning the
/// first error — with the failing psp last there is nothing after it to have been touched. That is
/// the case the plan's three covered, and closing it needs a three-sample fixture.
#[test]
fn a_run_stopped_part_way_does_only_what_is_left() {
    let (cohort, psps) = a_walked_cohort();
    let walked = the_trailers_in(&psps, &["alpha", "zeta"]);
    empty_every_trailer(&walked);
    let breaks_on_its_records = psp_path_for(&psps, "zeta");
    let with_its_blocks_intact = std::fs::read(&breaks_on_its_records).expect("the psp reads");
    corrupt_a_block_of(&breaks_on_its_records);

    let error = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect_err("zeta's records are nonsense, and a rebuild has to read them");

    assert!(
        matches!(
            &error,
            RegenerateCensusCliError::Build { sample, psp, .. }
                if sample == "zeta" && psp == &breaks_on_its_records
        ),
        "the run stops on the sample it could not rebuild, and names its file: {error:?}",
    );
    let alphas_census = the_trailers_in(&psps, &["alpha"])
        .pop()
        .expect("one entry")
        .1;
    assert_eq!(
        alphas_census,
        walked
            .iter()
            .find(|(path, _)| path == &psp_path_for(&psps, "alpha"))
            .expect("alpha was walked")
            .1,
        "the sample before the failure was rebuilt and kept",
    );

    // **The bytes are put back**, which is the part a person does by mending whatever made the read
    // fail — and this psp's trailer is still empty, so it is still owed.
    std::fs::write(&breaks_on_its_records, &with_its_blocks_intact)
        .expect("the scratch dir is ours");

    let report = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &psps))
        .expect("the psps read");

    assert_eq!(
        report
            .samples
            .iter()
            .map(|sample| sample.sample.as_str())
            .collect::<Vec<_>>(),
        vec!["zeta"],
        "the second run rebuilds the one still owed",
    );
    assert_eq!(
        report
            .skipped
            .iter()
            .map(|it| it.sample.as_str())
            .collect::<Vec<_>>(),
        vec!["alpha"],
        "and skips the one the first run finished — whose census this command wrote",
    );
    for (path, before) in &walked {
        let after = PspReader::open(path)
            .expect("the psp opens")
            .trailer()
            .expect("its trailer reads");
        assert_eq!(
            before,
            &after,
            "{} does not carry the census its walk wrote",
            path.display(),
        );
    }
}

/// **A regenerated psp is the walked one, byte for byte, whole** — header, blocks, index, trailer
/// and footer (plan step D4 of `psp_census_pair.md`; spec §11's parity oracle).
///
/// **What it covers that the trailer comparisons in this file do not**, and it is two things rather
/// than the four a first draft claimed: the header's exact **bytes**, where the sibling test
/// compares its `Debug` rendering, and every record's **content**, where that test compares only
/// how many there are.
///
/// **And the honest limit.** Against this command as it stands there is no small defect this
/// catches first: the only write is `replace_trailer`, which touches nothing below the trailer's
/// offset, and a psp whose index or footer was disturbed is refused by the reader before any
/// comparison here. What this is for is the change that would stop that being true — a rebuild that
/// wrote the psp through `PspWriter` again, re-compressed its blocks, or re-encoded its header.
///
/// **The copies' trailers are emptied first** (the plan's own note on this step). A psp whose
/// census this run would write is skipped since plan step D2, so copies handed straight to the
/// command come back untouched and a byte comparison passes without a rebuild having happened.
/// Emptying them also makes the comparison stronger than the plan asks: the footer's trailer length
/// goes to zero and has to come back.
///
/// **Both psps of the fixture with a repeat tract**, which is what gets all three of that cohort's
/// read groups — its samples declare two and one — and a cohort of more than one, where a run that
/// rebuilt only its first psp would pass over a cohort of one.
#[test]
fn a_regenerated_psp_is_the_walked_one_byte_for_byte() {
    let cohort = a_varying_cohort_on_disk();
    let psps = cohort.directory.path().join("psps");
    a_walk_into(
        &psps,
        &cohort.reference,
        &cohort.catalog,
        &cohort.alignments,
    );

    // The copies live in a directory of their own, so the command is handed them and not the
    // originals.
    let copies = cohort.directory.path().join("copies");
    std::fs::create_dir(&copies).expect("the scratch dir is ours");
    let mut as_the_walk_sealed_them = Vec::new();
    for sample in ["one", "two"] {
        let walked = psp_path_for(&psps, sample);
        let copy = copies.join(format!("{sample}.psp"));
        std::fs::copy(&walked, &copy).expect("the scratch dir is ours");
        as_the_walk_sealed_them
            .push((copy.clone(), std::fs::read(&walked).expect("the psp reads")));
        crate::ng::psp::replace_trailer(&copy, b"").expect("the tail rewrites");
        assert_ne!(
            std::fs::read(&copy).expect("the copy reads"),
            as_the_walk_sealed_them.last().expect("just pushed").1,
            "{sample}'s copy differs from the walk's file before the rebuild, or this comparison \
             is a file against itself",
        );
    }

    let report = regenerate_every_census(&args_over(&cohort.reference, &cohort.catalog, &copies))
        .expect("the copies read");

    assert_eq!(
        report
            .samples
            .iter()
            .map(|sample| sample.sample.as_str())
            .collect::<Vec<_>>(),
        vec!["one", "two"],
        "both copies were rebuilt",
    );
    // **What the comparison was over.** Two psps with no read at a kept tract would compare equal
    // while covering none of the census's stratum-keyed half.
    assert!(
        report
            .samples
            .iter()
            .any(|sample| sample.tally.tracts_with_reads > 0),
        "no sample has a read at a kept repeat tract, so this fixture stopped being the one with a \
         tract in it: {:?}",
        report.samples.iter().map(|it| it.tally).collect::<Vec<_>>(),
    );
    for (copy, as_walked) in &as_the_walk_sealed_them {
        assert_eq!(
            &std::fs::read(copy).expect("the copy reads"),
            as_walked,
            "{} is not the walked psp byte for byte",
            copy.display(),
        );
    }
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
