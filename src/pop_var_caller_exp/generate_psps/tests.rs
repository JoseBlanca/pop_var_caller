//! The walk stage's command surface: what a person may type, what is refused before a read
//! is decoded, and what the files it writes say about themselves.

use super::*;
use clap::Parser;

use crate::ng::psp::PspReader;
use crate::ng::repeat_catalog::StrRepeatCriteria;
use crate::pop_var_caller_exp::cli::{Cli, PopVarCallerExpCommand};

/// Parse an argument vector into this subcommand's arguments, refusing any other subcommand.
fn args_of(argv: &[&str]) -> GeneratePspsArgs {
    match Cli::parse_from(argv).cmd {
        PopVarCallerExpCommand::GeneratePsps(args) => args,
        other => panic!("expected generate-psps, got {other:?}"),
    }
}

fn refusal_of(argv: &[&str]) -> clap::Error {
    Cli::try_parse_from(argv).expect_err("this argument vector must be refused")
}

/// The shortest walk a person can type.
fn a_walk() -> Vec<&'static str> {
    vec![
        "pop_var_caller_exp",
        "generate-psps",
        "--reference",
        "ref.fa",
        "--alignment",
        "a.bam",
        "--output-dir",
        "psps",
    ]
}

#[test]
fn the_subcommand_is_spelled_generate_psps() {
    let args = args_of(&a_walk());
    assert_eq!(args.reference, PathBuf::from("ref.fa"));
    assert_eq!(args.output_dir, PathBuf::from("psps"));
}

/// **The alignment flag repeats and keeps the order it was given**, because the order the
/// samples are first seen in is the order their psps are written.
#[test]
fn the_alignment_flag_repeats_and_keeps_the_order_it_was_given() {
    let mut argv = vec![
        "pop_var_caller_exp",
        "generate-psps",
        "--reference",
        "ref.fa",
        "--output-dir",
        "psps",
    ];
    for file in ["zeta.bam", "alpha.bam", "beta.bam"] {
        argv.extend(["--alignment", file]);
    }
    let args = args_of(&argv);

    assert_eq!(
        args.alignments,
        vec![
            PathBuf::from("zeta.bam"),
            PathBuf::from("alpha.bam"),
            PathBuf::from("beta.bam"),
        ],
    );
}

/// A walk with no alignment file has nothing to do and is refused by clap.
#[test]
fn a_walk_with_no_alignment_is_refused() {
    let refused = refusal_of(&[
        "pop_var_caller_exp",
        "generate-psps",
        "--reference",
        "ref.fa",
        "--output-dir",
        "psps",
    ]);
    assert_eq!(
        refused.kind(),
        clap::error::ErrorKind::MissingRequiredArgument
    );
}

/// **The catalog defaults to the file `repeat-catalog` writes beside the reference**, the
/// same default direct mode takes, so the two modes read one file without being told twice.
#[test]
fn the_catalog_and_the_ground_have_defaults_a_person_need_not_type() {
    let args = args_of(&a_walk());

    assert!(args.catalog.is_none());
    assert!(args.regions.is_none());
    assert_eq!(
        ground_request(&args).catalog_path(),
        PathBuf::from("ref.fa.repeats.parquet"),
    );
}

/// **The routing defaults are ng's calling floors, not the catalog's storage floors** — the
/// same values direct mode routes with. A drift between the two commands' `default_value`s
/// would give one mode more repeat ground than the other over one catalog, and every psp
/// would then be refused by a calling run for disagreeing about the criteria (spec §6.2).
#[test]
fn the_routing_defaults_are_the_same_as_direct_modes() {
    let args = args_of(&a_walk());
    let asked = run_ground::routing_criteria(&ground_request(&args).routing)
        .expect("the defaults make a range");

    assert_ne!(
        asked,
        StrRepeatCriteria::default(),
        "the run's floors are its own, not the file's storage floors",
    );

    let direct = match Cli::parse_from([
        "pop_var_caller_exp",
        "call-from-alignments",
        "--reference",
        "ref.fa",
        "--alignment",
        "a.bam",
        "--output",
        "calls.vcf",
        "--defaults",
    ])
    .cmd
    {
        PopVarCallerExpCommand::CallFromAlignments(args) => args,
        other => panic!("expected call-from-alignments, got {other:?}"),
    };
    // **The whole criteria struct, not a tuple of the axes someone remembered.** An earlier
    // draft compared four of the five and left `min_purity` out, and a mutation moving that
    // one default alone passed the suite — so the comparison is against the value direct
    // mode's own defaults produce, which cannot go stale as axes are added.
    let direct_routing = run_ground::RepeatRouting {
        min_copies: direct.min_copies,
        min_period: direct.min_period,
        max_period: direct.max_period,
        max_str_len: direct.max_str_len,
        min_purity: direct.min_purity,
    };
    assert_eq!(
        asked,
        run_ground::routing_criteria(&direct_routing).expect("direct mode's defaults too"),
        "the two modes must route one catalog the same way",
    );
}

/// **A psp records the subcommand a person can actually type.** The name is written into
/// every header from one constant; clap derives the command's own spelling from the enum
/// variant. Nothing but this ties them, and a rename of the variant would otherwise leave
/// every psp naming a subcommand that no longer exists.
#[test]
fn the_recorded_subcommand_is_the_one_a_person_types() {
    let parsed = Cli::try_parse_from([
        "pop_var_caller_exp",
        SUBCOMMAND,
        "--reference",
        "ref.fa",
        "--alignment",
        "a.bam",
        "--output-dir",
        "psps",
    ])
    .expect("the recorded subcommand name is the one clap answers to");
    assert!(matches!(
        parsed.cmd,
        PopVarCallerExpCommand::GeneratePsps(_)
    ));
}

/// **An alignment file that is not there names itself**, and it is refused at the door —
/// before the reference is read, because reading a reference is minutes and a path that was
/// never going to open should cost none of them.
#[test]
fn an_alignment_file_that_is_not_there_names_itself_in_the_refusal() {
    let (_reference_dir, _zeta_dir, _alpha_dir, mut args) = a_cohort_on_disk();
    args.alignments.push(args.output_dir.join("absent.bam"));

    let refused = run_generate_psps(&args).expect_err("that file is not on disk");
    assert!(
        matches!(refused, GeneratePspsCliError::ReadGroups { .. }),
        "{refused:?}",
    );
    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(
        rendered.contains("absent.bam"),
        "the chain names the file that could not be read, and got: {rendered}",
    );
}

/// **A walk with no catalog is told which file is missing and the command that builds it** —
/// the shared refusal, so this reads exactly as it does from `call-from-alignments`.
#[test]
fn a_walk_with_no_catalog_is_told_which_file_is_missing_and_how_to_build_it() {
    let (_reference_dir, _zeta_dir, _alpha_dir, mut args) = a_cohort_on_disk();
    std::fs::remove_file(args.catalog.as_ref().expect("the fixture built one"))
        .expect("the catalog goes away");
    args.catalog = None;

    let refused = run_generate_psps(&args).expect_err("there is no catalog beside this reference");

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

/// **A catalog built on another reference of the same shape is refused.** Contig names and
/// lengths match, only the bases differ — the one case the digest comparison exists for, and
/// the reason the run checks the catalog against the reference it *read* rather than against
/// the `.fai`'s digest-free view of it.
///
/// A psp written against the wrong catalog puts every repeat tract at the wrong coordinates
/// for the whole genome, in a file that opens cleanly and says nothing is amiss.
#[test]
fn a_catalog_built_on_another_reference_of_the_same_shape_is_refused() {
    use crate::ng::read::input::test_fixtures::FIXTURE_CONTIGS;
    use crate::ng::reference_info::{ReferenceSource, read_reference_info_observing};
    use crate::ng::repeat_catalog::RepeatCatalogBuilder;
    use crate::ng::tandem_repeat::ScanParams;
    use crate::pileup::per_sample::cram_files::{ContigSpec, build_fasta};

    let (_reference_dir, _zeta_dir, _alpha_dir, mut args) = a_cohort_on_disk();

    // A second reference with the same contig names and lengths and different bases: every
    // `A` becomes a `C`, which leaves each line's width — and so the `.fai` — untouched.
    let specs: Vec<ContigSpec> = FIXTURE_CONTIGS
        .iter()
        .map(|(name, length)| ContigSpec {
            name: (*name).to_string(),
            length: *length as u64,
        })
        .collect();
    let (other_dir, other_fasta) = build_fasta(&specs).expect("a second reference on disk");
    let bases = std::fs::read(&other_fasta).expect("the second reference reads");
    let recoded: Vec<u8> = bases
        .into_iter()
        .map(|base| if base == b'A' { b'C' } else { base })
        .collect();
    std::fs::write(&other_fasta, recoded).expect("the second reference rewrites");

    let other_catalog = other_dir.path().join("other.repeats.parquet");
    let mut builder = RepeatCatalogBuilder::create(
        &other_catalog,
        StrRepeatCriteria::default(),
        ScanParams {
            match_reward: 2,
            mismatch_penalty: 7,
            min_copies: 2,
        },
    )
    .expect("a catalog to build into");
    let other = read_reference_info_observing(
        ReferenceSource::Fasta {
            fasta: other_fasta,
            fai: None,
        },
        &mut builder,
    )
    .expect("the second reference reads");
    builder.finish(&other).expect("the catalog is written");

    args.catalog = Some(other_catalog);
    let refused =
        run_generate_psps(&args).expect_err("this catalog was built on another reference");
    let rendered = crate::error_render::format_error_chain(&refused);
    assert!(
        rendered.contains("catalog does not describe this reference"),
        "the refusal is about the reference's digests, and got: {rendered}",
    );
}

/// **A walk that stops names its sample, and the samples before it keep their psps.** That
/// is the whole point of walking one sample at a time: the repair is to re-run one sample.
///
/// It also pins the order — beta is named last, so under any other order the run would stop
/// before zeta and alpha were written.
#[test]
fn a_walk_that_stops_names_its_sample_and_leaves_the_earlier_samples_psps_written() {
    use crate::ng::read::input::test_fixtures::{header, matching_contigs, named_bam};

    let (_reference_dir, _zeta_dir, _alpha_dir, mut args) = a_cohort_on_disk();
    // beta's file carries no index and the run was not told to build one, so its walk stops
    // at the open. zeta and alpha are named before it.
    let (_beta_dir, beta) = named_bam(
        &header(
            Some("coordinate"),
            &matching_contigs(),
            &[("rg9", Some("beta"))],
        ),
        &[],
        "beta.bam",
    );
    args.alignments.push(beta);

    let stopped = run_generate_psps(&args).expect_err("beta's file has no index");
    match &stopped {
        GeneratePspsCliError::Walk { sample, path, .. } => {
            assert_eq!(sample, "beta");
            assert_eq!(path, &psp_path_for(&args.output_dir, "beta"));
        }
        other => panic!("expected a stopped walk, got {other:?}"),
    }
    assert!(
        psp_path_for(&args.output_dir, "zeta").is_file(),
        "zeta is named first, so its psp was finished before beta stopped the run",
    );
    assert!(
        psp_path_for(&args.output_dir, "alpha").is_file(),
        "and alpha's too",
    );
    assert!(
        !psp_path_for(&args.output_dir, "beta").exists(),
        "and the sample that stopped left nothing at its own path",
    );
}

/// **A stopped re-walk leaves the psp it was replacing intact.** The command's advertised
/// repair is *re-run the one sample that failed*; if the walk wrote straight to the final
/// path, `PspWriter::create` would truncate a good psp at the first byte and a second failure
/// would leave a stump every reader refuses.
#[test]
fn a_stopped_rewalk_does_not_destroy_the_psp_it_was_replacing() {
    use crate::ng::read::input::test_fixtures::{
        header, matching_contigs, named_bam, read_group_for,
    };

    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();
    run_generate_psps(&args).expect("the cohort walks");
    let zeta_psp = psp_path_for(&args.output_dir, "zeta");
    let written = std::fs::read(&zeta_psp).expect("zeta's psp is on disk");

    // Re-walk zeta alone, from a file that cannot be opened: same sample, same output path.
    let (_broken_dir, broken) = named_bam(
        &header(
            Some("coordinate"),
            &matching_contigs(),
            &[(&read_group_for("zeta"), Some("zeta"))],
        ),
        &[],
        "zeta.bam",
    );
    let mut rewalk = a_cohort_on_disk().3;
    rewalk.output_dir = args.output_dir.clone();
    rewalk.alignments = vec![broken];
    rewalk.reference = args.reference.clone();
    rewalk.catalog = args.catalog.clone();
    // **`--force`, or this never reaches the walk at all**: zeta's psp is already there, so
    // without it the run stops at the overwrite check and this test would prove nothing about
    // what a stopped *walk* leaves. (It did exactly that until the check landed.)
    rewalk.force = true;

    let stopped = run_generate_psps(&rewalk).expect_err("the re-walk's file has no index");
    assert!(
        matches!(stopped, GeneratePspsCliError::Walk { .. }),
        "the re-walk must fail in the walk, not at the door: {stopped:?}",
    );

    assert_eq!(
        std::fs::read(&zeta_psp).expect("zeta's psp is still on disk"),
        written,
        "the psp the re-walk was replacing is untouched, byte for byte",
    );
}

/// **`--regions` narrows the ground the walk analyses**, and the psp records the narrowed
/// ground — which is what a later calling run compares across the cohort.
#[test]
fn the_regions_flag_narrows_the_ground_the_walk_analyses() {
    let (reference_dir, _zeta_dir, _alpha_dir, mut args) = a_cohort_on_disk();
    let bed = reference_dir.path().join("ground.bed");
    std::fs::write(&bed, "chr1\t0\t10\n").expect("a bed on disk");
    args.regions = Some(bed);

    run_generate_psps(&args).expect("the cohort walks");

    let reader =
        PspReader::open(&psp_path_for(&args.output_dir, "zeta")).expect("a finished psp opens");
    let spans: Vec<_> = reader
        .header()
        .segmentation_inputs
        .analysed_regions
        .iter()
        .collect();
    assert_eq!(
        spans.len(),
        1,
        "one BED interval was asked for, not both contigs: {spans:?}",
    );
}

/// **`--catalog` is read, rather than the file beside the reference.** The fixture writes its
/// catalog exactly where the default looks, so only moving it can tell the two apart.
#[test]
fn the_catalog_flag_is_read_rather_than_the_file_beside_the_reference() {
    let (reference_dir, _zeta_dir, _alpha_dir, mut args) = a_cohort_on_disk();
    let moved = reference_dir.path().join("elsewhere.parquet");
    std::fs::rename(args.catalog.as_ref().expect("a catalog"), &moved).expect("the catalog moves");
    args.catalog = Some(moved);

    run_generate_psps(&args).expect("the cohort walks from the catalog it was handed");
}

/// **`--min-purity` reaches the criteria the catalog is asked with**, and the psp records
/// what was asked. The four other routing flags are pinned at their defaults against direct
/// mode's; this is the one whose value is read back from a written file.
#[test]
fn the_min_purity_flag_reaches_the_criteria_the_catalog_is_asked_with() {
    let (_reference_dir, _zeta_dir, _alpha_dir, mut args) = a_cohort_on_disk();
    args.min_purity = 0.95;

    run_generate_psps(&args).expect("the cohort walks");

    let reader =
        PspReader::open(&psp_path_for(&args.output_dir, "zeta")).expect("a finished psp opens");
    assert_eq!(
        reader
            .header()
            .segmentation_inputs
            .repeat_tract_criteria
            .classification
            .min_purity,
        0.95_f32,
    );
}

/// **`--build-index-if-missing` reaches the open that would use it.** Without the flag this
/// sample's walk stops at the open, which is what the stopped-walk test above relies on.
#[test]
fn the_build_index_flag_reaches_the_open_that_would_use_it() {
    use crate::ng::read::input::test_fixtures::{header, matching_contigs, named_bam};

    let (_reference_dir, _zeta_dir, _alpha_dir, mut args) = a_cohort_on_disk();
    let (_beta_dir, beta) = named_bam(
        &header(
            Some("coordinate"),
            &matching_contigs(),
            &[("rg9", Some("beta"))],
        ),
        &[],
        "beta.bam",
    );
    args.alignments.push(beta);
    args.build_index_if_missing = true;

    run_generate_psps(&args).expect("beta's index is built rather than demanded");
    assert!(psp_path_for(&args.output_dir, "beta").is_file());
}

/// **The names on disk, read from the directory rather than recomputed by the code under
/// test.** Every other test locates a psp by calling `psp_path_for` again, so none of them
/// can see the name change.
#[test]
fn each_psp_is_named_for_its_sample_with_the_psp_extension() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();

    run_generate_psps(&args).expect("the cohort walks");

    let mut names: Vec<String> = std::fs::read_dir(&args.output_dir)
        .expect("the output directory")
        .map(|entry| {
            entry
                .expect("an entry")
                .file_name()
                .to_string_lossy()
                .into_owned()
        })
        .collect();
    names.sort();
    assert_eq!(names, vec!["alpha.psp".to_string(), "zeta.psp".to_string()]);
}

/// **A sample whose `@RG SM` cannot be a file name is refused before anything is read.**
/// The tag is free header text, so a sample called `../elsewhere` would otherwise write
/// outside `--output-dir`, and one called `lane/1` would fail at the write with a whole
/// walk's decoding already spent.
#[test]
fn a_sample_name_that_cannot_be_a_file_name_is_refused_at_the_door() {
    for name in ["../elsewhere", "lane/1", "", ".", "..", "/absolute"] {
        let refused = refuse_a_sample_name_that_is_not_a_file_name(name)
            .expect_err("this name cannot be a psp's file name");
        assert!(
            matches!(refused, GeneratePspsCliError::SampleNameNotAFileName { .. }),
            "{name:?} got {refused:?}",
        );
    }
    for name in ["zeta", "NA12878", "sample.1", "a-b_c"] {
        refuse_a_sample_name_that_is_not_a_file_name(name)
            .unwrap_or_else(|error| panic!("{name:?} is a usable file name, got {error:?}"));
    }
}

fn a_cohort_on_disk() -> (
    tempfile::TempDir,
    tempfile::TempDir,
    tempfile::TempDir,
    GeneratePspsArgs,
) {
    let cohort = crate::pop_var_caller_exp::test_fixtures::a_cohort_on_disk();
    let args = GeneratePspsArgs {
        reference: cohort.reference,
        catalog: Some(cohort.catalog),
        alignments: cohort.alignments,
        output_dir: cohort.directory.path().join("psps"),
        regions: None,
        force: false,
        build_index_if_missing: false,
        min_copies: MinCopies::default(),
        min_period: DEFAULT_MIN_PERIOD,
        max_period: DEFAULT_MAX_PERIOD,
        max_str_len: DEFAULT_MAX_STR_LEN,
        min_purity: DEFAULT_MIN_PURITY,
    };
    (cohort.directory, cohort.zeta, cohort.alpha, args)
}

/// **The command writes one psp per sample, named for the sample, and each one opens.**
///
/// This is the only test that drives `run_generate_psps` itself, and it is what stops the
/// wiring from being proved only in pieces: the direct-mode command's own review measured
/// that gap and found nine of fourteen mutations survived while every test called the
/// helpers directly.
#[test]
fn the_command_writes_one_psp_per_sample_named_for_the_sample() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();

    run_generate_psps(&args).expect("the cohort walks");

    for sample in ["zeta", "alpha"] {
        let psp = psp_path_for(&args.output_dir, sample);
        assert!(psp.is_file(), "{sample}'s psp is at {psp:?}");
        let reader = PspReader::open(&psp).expect("a finished psp opens");
        assert_eq!(
            reader.header().sample,
            sample,
            "and it names the sample it holds",
        );
    }
}

/// **A sample's reads reach its psp.** Without this the suite could not tell a psp holding a
/// sample's evidence from a well-formed empty one — every other command-level test here is
/// satisfied by a file with no records at all.
#[test]
fn a_samples_reads_reach_its_psp() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();

    run_generate_psps(&args).expect("the cohort walks");

    let mut reader =
        PspReader::open(&psp_path_for(&args.output_dir, "zeta")).expect("a finished psp opens");
    let records = reader.records().expect("the walk starts").count();
    assert!(
        records > 0,
        "zeta's three reads must put records in its psp, and got {records}",
    );

    // And the sample with no reads still gets a whole, openable, empty psp — §12.9's
    // analysed-but-empty case reaching the command.
    let mut empty =
        PspReader::open(&psp_path_for(&args.output_dir, "alpha")).expect("alpha's psp opens");
    assert_eq!(empty.records().expect("the walk starts").count(), 0);
}

/// **The header records what a calling run will check** (spec §6.1, §6.2): the ground this
/// walk analysed and the criteria it routed with. Two psps written by one invocation must
/// agree on both, or no cohort could ever be assembled from them.
#[test]
fn every_psp_records_the_ground_and_the_criteria_the_walk_used() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();

    run_generate_psps(&args).expect("the cohort walks");

    let inputs_of = |sample: &str| {
        PspReader::open(&psp_path_for(&args.output_dir, sample))
            .expect("a finished psp opens")
            .header()
            .segmentation_inputs
            .clone()
    };
    let zeta = inputs_of("zeta");
    assert_eq!(
        zeta,
        inputs_of("alpha"),
        "one invocation's psps must be joinable into one cohort",
    );
    assert_eq!(
        zeta.repeat_tract_criteria,
        run_ground::routing_criteria(&ground_request(&args).routing).expect("a range"),
        "and the criteria recorded are the ones the walk actually routed with",
    );
    assert!(
        !zeta.analysed_regions.iter().collect::<Vec<_>>().is_empty(),
        "the whole reference is the analysed ground when no BED is given",
    );
}

/// **The provenance says which command wrote the file**, so a psp found on disk months later
/// can be traced to the invocation that made it (spec §6.1).
#[test]
fn every_psp_records_the_command_that_wrote_it() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();

    run_generate_psps(&args).expect("the cohort walks");

    let reader =
        PspReader::open(&psp_path_for(&args.output_dir, "zeta")).expect("a finished psp opens");
    let writer = &reader.header().writer;
    assert_eq!(writer.tool, "ng");
    assert_eq!(writer.subcommand, "generate-psps");
    assert_eq!(
        writer.input_alignments,
        vec!["zeta.bam"],
        "the gatherer overwrote the placeholder with the file it opened",
    );
    assert!(
        writer.parameters.contains_key("read-filter-min-mapq"),
        "and the read filters the walk applied are recorded: {:?}",
        writer.parameters,
    );
}

/// **The output directory is created rather than demanded.** A walk is often the first thing
/// run on a fresh machine, and refusing because a directory does not exist yet would be a
/// refusal with nothing behind it.
#[test]
fn the_output_directory_is_created_when_it_does_not_exist() {
    let (_reference_dir, _zeta_dir, _alpha_dir, mut args) = a_cohort_on_disk();
    args.output_dir = args.output_dir.join("not").join("there").join("yet");
    assert!(!args.output_dir.exists());

    run_generate_psps(&args).expect("the cohort walks");

    assert!(psp_path_for(&args.output_dir, "zeta").is_file());
}

/// **Two files naming one sample are one walk into one psp**, which is what `@RG SM` means
/// (`read_groups.md` §4) — and the file the walk was given twice is opened once.
#[test]
fn two_files_naming_one_sample_become_one_psp() {
    use crate::ng::read::input::test_fixtures::{header, indexed_named_bam, matching_contigs};

    let (_reference_dir, _zeta_dir, _alpha_dir, mut args) = a_cohort_on_disk();
    let with_second_read_group = |file: &str| {
        indexed_named_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg2", Some("zeta"))],
            ),
            &[],
            file,
        )
    };
    let (_second_dir, second) = with_second_read_group("zeta-lane2.bam");
    // zeta now has two files; alpha still has one.
    args.alignments.push(second);

    run_generate_psps(&args).expect("the cohort walks");

    let reader =
        PspReader::open(&psp_path_for(&args.output_dir, "zeta")).expect("a finished psp opens");
    assert_eq!(
        reader.header().writer.input_alignments,
        vec!["zeta.bam", "zeta-lane2.bam"],
        "both of the sample's files went into its one psp",
    );
    assert_eq!(
        reader.header().read_groups.len(),
        2,
        "and its walk-local table holds both read groups: {:?}",
        reader.header().read_groups,
    );
    assert_eq!(
        std::fs::read_dir(&args.output_dir)
            .expect("the output directory")
            .count(),
        2,
        "two samples, two psps",
    );
}

/// **A file declaring several of one sample's read groups is opened once**, which is the
/// property [`files_of`] exists for: the read-group table holds one row per group, so a
/// three-lane file appears three times in it and handing all three to one walk would open
/// the same file three times.
///
/// The fixture table mints every read group against one synthetic path, which is exactly
/// this shape — without the deduplication the list would come back three long.
#[test]
fn a_file_holding_several_of_a_samples_read_groups_is_listed_once() {
    use crate::ng::read::input::read_groups::ReadGroups;

    let read_groups = ReadGroups::of_libraries(&[
        ("rg1", "zeta"),
        ("rg2", "zeta"),
        ("rg3", "alpha"),
        ("rg4", "zeta"),
    ]);
    assert_eq!(
        read_groups
            .iter()
            .filter(|(_, group)| group.sample.as_ref() == "zeta")
            .count(),
        3,
        "the fixture must give zeta three read groups, or this proves no deduplication",
    );

    let zeta = read_groups
        .read_groups_per_sample()
        .iter()
        .find(|sample| sample.sample.as_ref() == "zeta")
        .expect("zeta is in the table");
    let files = alignment_files_of(zeta, &read_groups);
    assert_eq!(
        files.len(),
        1,
        "zeta's three read groups live in one file, so the walk opens one: {files:?}",
    );
}

// ---------------------------------------------------------------------
// C2 — what the run says about itself
// ---------------------------------------------------------------------

/// **The report names every sample, what it stored, and how much of the ground it could
/// speak for.** The lines are a value rather than something only a terminal sees, so this
/// can hold them.
#[test]
fn the_report_says_what_each_sample_stored_and_how_much_ground_it_covered() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();

    let report = walk_every_sample(&args).expect("the cohort walks");
    let lines = report.lines();

    assert_eq!(
        report.samples.len(),
        2,
        "one outcome per sample: {:?}",
        report.samples.iter().map(|o| &o.sample).collect::<Vec<_>>(),
    );
    let whole = lines.join("\n");
    for sample in ["zeta", "alpha"] {
        assert!(
            whole.contains(sample),
            "the report names {sample}, and got:\n{whole}",
        );
    }
    assert!(
        whole.contains("bases asked for"),
        "and says over how much ground, and got:\n{whole}",
    );

    // The sample with reads stored observations; the empty one stored none — and the report
    // says so rather than reporting one number for the cohort.
    let zeta = report
        .samples
        .iter()
        .find(|outcome| outcome.sample == "zeta")
        .expect("zeta walked");
    let alpha = report
        .samples
        .iter()
        .find(|outcome| outcome.sample == "alpha")
        .expect("alpha walked");
    assert!(zeta.stats.records > 0, "zeta has reads");
    assert_eq!(alpha.stats.records, 0, "alpha has none");
    assert!(
        zeta.counts.regions_handled > 0,
        "and the ground it spoke for is counted: {:?}",
        zeta.counts,
    );
}

/// **The report's per-sample line carries the numbers, not just the name.** A line that named
/// a sample and said nothing measurable about it would pass a `contains` check while telling
/// a person nothing.
#[test]
fn the_reports_sample_line_carries_what_that_sample_stored() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();

    let report = walk_every_sample(&args).expect("the cohort walks");
    let zeta = report
        .samples
        .iter()
        .find(|outcome| outcome.sample == "zeta")
        .expect("zeta walked");
    let line = report
        .lines()
        .into_iter()
        .find(|line| line.contains("zeta"))
        .expect("zeta has a line");

    assert!(
        line.contains(&zeta.stats.records.to_string()),
        "its observation count is in its line: {line}",
    );
    assert!(
        line.contains(&zeta.stats.bytes.to_string()),
        "and how big its psp is: {line}",
    );
    assert!(
        line.contains(&zeta.counts.regions_handled.to_string()),
        "and how many segments it spoke for: {line}",
    );
}

// ---------------------------------------------------------------------
// C3 — a psp already there, and what a stopped run leaves
// ---------------------------------------------------------------------

/// **A second run refuses rather than replacing a psp**, naming the sample, the file and the
/// flag that means it.
#[test]
fn a_rerun_refuses_to_replace_a_psp_that_is_already_there() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();
    run_generate_psps(&args).expect("the cohort walks");
    let zeta_psp = psp_path_for(&args.output_dir, "zeta");
    let first = std::fs::read(&zeta_psp).expect("zeta's psp is on disk");

    let refused = run_generate_psps(&args).expect_err("zeta's psp is already there");
    match &refused {
        GeneratePspsCliError::PspAlreadyThere { sample, path } => {
            assert_eq!(sample, "zeta");
            assert_eq!(path, &zeta_psp);
        }
        other => panic!("expected a refusal to replace, got {other:?}"),
    }
    assert!(
        crate::error_render::format_error_chain(&refused).contains("--force"),
        "and the message says how to mean it",
    );
    assert_eq!(
        std::fs::read(&zeta_psp).expect("still there"),
        first,
        "and the psp it refused to replace is untouched",
    );
}

/// **The refusal comes before any sample is walked.** A cohort whose *second* psp is already
/// there must not spend the first sample's walk before saying so.
#[test]
fn the_refusal_comes_before_the_first_sample_is_walked() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();
    std::fs::create_dir_all(&args.output_dir).expect("the output directory");
    // Only alpha's psp is there, and alpha is walked second.
    std::fs::write(psp_path_for(&args.output_dir, "alpha"), b"not really a psp")
        .expect("a file in the way");

    let refused = run_generate_psps(&args).expect_err("alpha's psp is already there");
    assert!(
        matches!(refused, GeneratePspsCliError::PspAlreadyThere { .. }),
        "{refused:?}",
    );
    assert!(
        !psp_path_for(&args.output_dir, "zeta").exists(),
        "zeta is named first and must not have been walked",
    );
}

/// **`--force` is how a person says they mean it**, and then the psp is replaced.
#[test]
fn force_replaces_a_psp_that_is_already_there() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();
    std::fs::create_dir_all(&args.output_dir).expect("the output directory");
    let zeta_psp = psp_path_for(&args.output_dir, "zeta");
    std::fs::write(&zeta_psp, b"not really a psp").expect("a file in the way");

    let mut forced = a_cohort_on_disk().3;
    forced.reference = args.reference.clone();
    forced.catalog = args.catalog.clone();
    forced.alignments = args.alignments.clone();
    forced.output_dir = args.output_dir.clone();
    forced.force = true;

    run_generate_psps(&forced).expect("--force replaces it");
    PspReader::open(&zeta_psp).expect("what is there now is a whole psp");
}

/// **A file that is not a whole psp is refused by a reader, whatever it is called.** The
/// format guarantees it; what this pins is that the command's own output directory is not a
/// place where a half-written file can masquerade as a finished one.
#[test]
fn a_truncated_psp_is_refused_as_interrupted() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();
    run_generate_psps(&args).expect("the cohort walks");

    let zeta_psp = psp_path_for(&args.output_dir, "zeta");
    let whole = std::fs::read(&zeta_psp).expect("zeta's psp is on disk");
    assert!(whole.len() > 64, "the fixture psp has a body to cut");
    std::fs::write(&zeta_psp, &whole[..whole.len() / 2]).expect("cut it in half");

    PspReader::open(&zeta_psp).expect_err("half a psp is not a psp");
}

/// **A stopped walk leaves nothing in the output directory that a later run could open** —
/// no psp at the sample's own path, and no scratch file beside it either.
#[test]
fn a_stopped_walk_leaves_no_file_at_the_samples_path_or_beside_it() {
    use crate::ng::read::input::test_fixtures::{header, matching_contigs, named_bam};

    let (_reference_dir, _zeta_dir, _alpha_dir, mut args) = a_cohort_on_disk();
    let (_beta_dir, beta) = named_bam(
        &header(
            Some("coordinate"),
            &matching_contigs(),
            &[("rg9", Some("beta"))],
        ),
        &[],
        "beta.bam",
    );
    args.alignments.push(beta);

    // A scratch file from an earlier stopped run of this same process, so the cleanup has
    // something to remove: beta stops inside `open`, before a walk of its own creates one.
    std::fs::create_dir_all(&args.output_dir).expect("the output directory");
    let stale = psp_path_for(&args.output_dir, "beta")
        .with_extension(format!("psp.{}.partial", std::process::id()));
    std::fs::write(&stale, b"a stopped walk's leavings").expect("a stale partial");

    run_generate_psps(&args).expect_err("beta's file has no index");

    assert!(
        !stale.exists(),
        "the stopped sample's scratch file was cleared, and got {stale:?}",
    );
    let left: Vec<String> = std::fs::read_dir(&args.output_dir)
        .expect("the output directory")
        .map(|entry| {
            entry
                .expect("an entry")
                .file_name()
                .to_string_lossy()
                .into_owned()
        })
        .filter(|name| name.starts_with("beta"))
        .collect();
    assert!(
        left.is_empty(),
        "the stopped sample left nothing behind, and got {left:?}",
    );
}

/// **Every number the per-sample line prints is the number it claims to be**, checked on a
/// walk outcome built by hand so that no two of them are equal — the command's own fixture
/// covers all of its ground, which makes a swap of the two region counts invisible.
///
/// Three of these eight were asserted before; mutations that swapped the region counts,
/// blanked the denominator or deleted the uncovered-ground clause all passed.
#[test]
fn the_per_sample_line_prints_the_numbers_it_names() {
    let outcome = SampleWalkOutcome {
        sample: "zeta".to_string(),
        psp: PathBuf::from("psps/zeta.psp"),
        stats: crate::ng::psp::WriteStats {
            records: 41,
            blocks: 3,
            bytes: 6007,
        },
        counts: crate::ng::locus_generation::LocusCounts {
            regions_in: 17,
            regions_handled: 11,
            regions_handled_bp: 500,
            loci_emitted: 41,
            unhandled_not_implemented: 4,
            unhandled_not_implemented_bp: 300,
            unhandled_out_of_scope: 2,
            unhandled_out_of_scope_bp: 200,
        },
    };
    let line = outcome.line();

    assert_eq!(
        outcome.bases_walked(),
        1_000,
        "the whole is the three parts' sum, not the ground asked for",
    );
    assert!(line.contains("41 loci stored"), "{line}");
    assert!(line.contains("6007 bytes"), "{line}");
    assert!(line.contains("psps/zeta.psp"), "{line}");
    assert!(
        line.contains("spoke for 11 of 17 typed regions"),
        "handled before dispatched, in that order: {line}",
    );
    assert!(
        line.contains("500 of 1000 bases walked, 50.0%"),
        "the share's denominator is the ground the walk was handed: {line}",
    );
    assert!(
        line.contains("clean flanks: 300 bases (30.0%)"),
        "the ground no generator is built for yet is named and sized: {line}",
    );
    assert!(
        line.contains("callable: 200 bases (20.0%)"),
        "and so is the ground no caller will ever speak for: {line}",
    );
}

/// **A walk that covered everything says so by saying nothing about what it did not.** The
/// two clauses are conditional, so a clean run's line is not two zeros a reader must skip.
#[test]
fn a_walk_that_covered_its_whole_ground_carries_no_uncovered_clause() {
    let outcome = SampleWalkOutcome {
        sample: "zeta".to_string(),
        psp: PathBuf::from("psps/zeta.psp"),
        stats: crate::ng::psp::WriteStats {
            records: 41,
            blocks: 3,
            bytes: 6007,
        },
        counts: crate::ng::locus_generation::LocusCounts {
            regions_in: 11,
            regions_handled: 11,
            regions_handled_bp: 1_000,
            loci_emitted: 41,
            ..Default::default()
        },
    };
    let line = outcome.line();

    assert!(line.contains("1000 of 1000 bases walked, 100.0%"), "{line}");
    assert!(!line.contains("not stored"), "nothing to report: {line}");
}

/// **The shares are of the ground the walk was handed, not of the ground asked for.** A
/// typed region is walked whole even where a BED cuts one, so the two totals differ — and
/// dividing by the ask is how the sibling report once printed 200.0%.
#[test]
fn the_shares_are_of_the_ground_the_walk_was_handed() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();
    let report = walk_every_sample(&args).expect("the cohort walks");

    for outcome in &report.samples {
        let counts = &outcome.counts;
        assert_eq!(
            outcome.bases_walked(),
            counts.regions_handled_bp
                + counts.unhandled_not_implemented_bp
                + counts.unhandled_out_of_scope_bp,
            "the whole is the three parts' own sum",
        );
        assert!(
            counts.regions_handled_bp <= outcome.bases_walked(),
            "so no share can exceed one hundred per cent",
        );
    }
}

/// **The report names the ground it walked**, so a psp found later can be matched to the run
/// that made it. The sibling report records this same omission being found once already.
#[test]
fn the_report_names_the_ground_it_walked() {
    let (reference_dir, _zeta_dir, _alpha_dir, mut args) = a_cohort_on_disk();
    let bed = reference_dir.path().join("ground.bed");
    std::fs::write(&bed, "chr1\t0\t10\n").expect("a bed on disk");
    args.regions = Some(bed);

    let report = walk_every_sample(&args).expect("the cohort walks");
    let opening = report.lines().first().cloned().expect("a first line");

    assert!(
        opening.contains("chr1"),
        "the ground is named by its chromosome, and got: {opening}",
    );
}

/// **The refusal names the sample that is actually blocked**, not the first one in the run.
/// With only the second sample's psp in the way, naming the first would send a person to the
/// wrong file.
#[test]
fn the_refusal_names_the_sample_whose_psp_is_in_the_way() {
    let (_reference_dir, _zeta_dir, _alpha_dir, args) = a_cohort_on_disk();
    std::fs::create_dir_all(&args.output_dir).expect("the output directory");
    std::fs::write(psp_path_for(&args.output_dir, "alpha"), b"not really a psp")
        .expect("a file in the way");

    let refused = run_generate_psps(&args).expect_err("alpha's psp is already there");
    match &refused {
        GeneratePspsCliError::PspAlreadyThere { sample, path } => {
            assert_eq!(
                sample, "alpha",
                "zeta is walked first but is not the blocked one"
            );
            assert_eq!(path, &psp_path_for(&args.output_dir, "alpha"));
        }
        other => panic!("expected a refusal to replace, got {other:?}"),
    }
}

/// **The flag the refusal tells a person to type is the flag that exists.** Nothing else ties
/// the message's `--force` to the argument's own spelling.
#[test]
fn the_flag_the_refusal_names_is_the_flag_clap_answers_to() {
    let mut argv = a_walk();
    argv.push("--force");
    let args = args_of(&argv);

    assert!(args.force, "clap answers to the flag the refusal names");
}
