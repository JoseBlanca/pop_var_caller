//! What a run says about itself, as text a test can hold.
//!
//! **The summary was the one part of this command a mutation could change with the whole suite
//! green**, because it went straight to `println!`. It is lines now, and these are what pin
//! them: the arithmetic a reader would otherwise do by hand, the three kinds of nothing kept
//! apart, and every span named by its chromosome rather than by an index.

use super::*;
use crate::fasta::{ContigEntry, ContigList};
use crate::ng::calling::parameters_file::{
    CensusIdentity, DeclaredInbreeding, ReadsBehindEachCalibration,
};
use crate::ng::calling::run_parameters::RunParameters;
use crate::ng::locus_generation::LocusCounts;
use crate::ng::parameter_estimation::joint::loci::ReferenceDigest;
use crate::ng::read::input::read_groups::{ReadGroups, build_read_groups};
use crate::ng::region_typing::GenomeRegions;
use crate::ng::repeat_catalog::StrRepeatCriteria;
use crate::ng::run::AssemblyCheckOutcome;
use crate::ng::run::callers::TractOutcomes;
use crate::ng::run::callers::{CohortCallingTallies, CohortWalkTallies, SampleWalkTallies};
use crate::ng::types::{ContigId, Ploidy, Position, ReadGroupId};

/// The two contigs every fixture here names its spans on.
fn contigs() -> ContigList {
    ContigList {
        entries: vec![
            ContigEntry {
                name: "chr1".to_owned(),
                length: 1_000,
                md5: None,
            },
            ContigEntry {
                name: "chr2".to_owned(),
                length: 2_000,
                md5: None,
            },
        ],
    }
}

/// A span on one of those contigs, 1-based and inclusive.
fn region(contig: u32, start: u64, end: u64) -> GenomeRegion {
    GenomeRegion {
        contig: ContigId(contig),
        start: Position(start),
        end: Position(end),
    }
}

/// One sample's walk, with the ground it covered and what its filters did.
fn walked(
    sample_name: &str,
    regions: LocusCounts,
    read_filters: Vec<(Option<ReadGroupId>, ReadFilterCounts)>,
) -> SampleWalkTallies {
    SampleWalkTallies {
        sample_name: sample_name.to_owned(),
        regions,
        read_filters,
        // **`None` is the shape of a sample whose generator counted nothing**, which is what a
        // sample with no reads in the analysed ground reports — and the case the "with none"
        // line is about.
        snp_indel: None,
    }
}

/// The ground counts of a run that called most of what it was given.
fn ground_mostly_called() -> LocusCounts {
    LocusCounts {
        regions_in: 10,
        regions_handled: 7,
        regions_handled_bp: 800,
        loci_emitted: 40,
        unhandled_not_implemented: 2,
        unhandled_not_implemented_bp: 150,
        unhandled_out_of_scope: 1,
        unhandled_out_of_scope_bp: 50,
    }
}

/// What a run produced, with everything the report reads named.
fn a_run(
    records_written: u64,
    loci_called_but_not_written: u64,
    too_wide: Vec<GenomeRegion>,
    nobody: Vec<GenomeRegion>,
    per_sample: Vec<SampleWalkTallies>,
) -> WrittenCohort {
    WrittenCohort {
        calling: CohortCallingTallies {
            records_written,
            loci_called_but_not_written,
            loci_too_wide_to_assemble: too_wide,
            loci_with_nobody_to_call: nobody,
            tracts: TractOutcomes::default(),
        },
        walk: CohortWalkTallies {
            per_sample,
            assembly_check: AssemblyCheckOutcome::NothingCouldBeChecked {
                because: crate::ng::run::callers::NoChecksums::TheReferenceCarriesNone,
            },
        },
    }
}

/// The parameters a defaults run writes, over a real read-group table.
fn a_defaults_runs_parameters(read_groups: &ReadGroups) -> ParametersFile {
    let inbreeding = DeclaredInbreeding::nothing_said();
    let scored_with = RunParameters::of_defaults(
        read_groups,
        Ploidy::try_new(2).expect("a diploid"),
        &inbreeding,
    );
    ParametersFile::of_run(
        &scored_with,
        read_groups,
        &ReadsBehindEachCalibration::nothing_was_fitted(read_groups.len()),
        &inbreeding.of_each_sample(read_groups),
        &ReferenceDigest([7; 16]),
        CensusIdentity::of_a_run_with_no_census(),
        &StrRepeatCriteria::default(),
    )
}

/// Two fixture alignment files, and the run's read-group table over them.
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

/// The ground a report is rendered over: one interval of `bases` on `chr1`.
fn ground_of(bases: u64) -> GenomeRegions {
    GenomeRegions::whole_contigs(&[crate::regions::ContigBounds {
        name: "chr1",
        length: bases as u32,
    }])
}

/// The bounds the fixtures' runs called under — the shipped ones, named so a report's advice can
/// quote them.
fn shipped_bounds() -> BoundsTheRunCalledUnder {
    BoundsTheRunCalledUnder {
        max_cohort_locus_span: 50,
        max_candidate_alleles: 6,
    }
}

/// The whole report as one string, for the tests that ask what it does and does not say.
fn rendered(
    written: &WrittenCohort,
    read_groups: &ReadGroups,
    parameters: &ParametersFile,
    analysed_bases: u64,
) -> String {
    lines_of(written, read_groups, parameters, analysed_bases).join("\n")
}

/// The report's own lines, for the tests that ask about one.
fn lines_of(
    written: &WrittenCohort,
    read_groups: &ReadGroups,
    parameters: &ParametersFile,
    analysed_bases: u64,
) -> Vec<String> {
    let contigs = contigs();
    let ground = ground_of(analysed_bases);
    RunReport::of(
        written,
        &contigs,
        read_groups,
        parameters,
        &ground,
        shipped_bounds(),
    )
    .lines()
}

/// **Called and written are different counts, and the report states both and their difference.**
///
/// A locus no written genotype carries an alternative at establishes no variant and is left out
/// (`vcf_output.md` §9), so a reader holding only the record count cannot tell how much ground
/// was called. A run whose two counts are far apart called a great deal that came back matching
/// the reference, which is ordinary at low depth and worth being able to see.
#[test]
fn the_report_states_what_was_called_and_what_of_it_was_written() {
    let (_zeta, _alpha, read_groups) = a_cohorts_read_groups();
    let written = a_run(
        120,
        45,
        Vec::new(),
        Vec::new(),
        vec![walked("zeta", ground_mostly_called(), Vec::new())],
    );

    let text = rendered(
        &written,
        &read_groups,
        &a_defaults_runs_parameters(&read_groups),
        1_000,
    );

    assert!(text.contains("records written: 120"), "{text}");
    assert!(
        text.contains("loci called: 165 — 120 written, 45 establishing no variant"),
        "the two parts and their sum, so a reader need not add them: {text}",
    );
}

/// **The analysed ground partitions into what was called and the two kinds of what was not**,
/// each with the share of the whole that a bare base count does not give.
///
/// The two kinds are kept apart because a reader acts on each differently: ground this caller
/// has not built a generator for **yet** is a gap that will close, and a satellite never will.
#[test]
fn the_report_partitions_the_analysed_ground_and_gives_each_part_its_share() {
    let (_zeta, _alpha, read_groups) = a_cohorts_read_groups();
    let written = a_run(
        10,
        0,
        Vec::new(),
        Vec::new(),
        vec![walked("zeta", ground_mostly_called(), Vec::new())],
    );

    let text = rendered(
        &written,
        &read_groups,
        &a_defaults_runs_parameters(&read_groups),
        1_000,
    );

    assert!(
        text.contains("analysed ground: chr1:1-1000 — 1000 bases asked for, in 10 typed regions"),
        "the ground is named, not only counted: {text}",
    );
    assert!(text.contains("called: 800 bases (80.0%)"), "{text}");
    assert!(
        text.contains(
            "clusters of repeats too close together to have clean flanks: 150 bases (15.0%)"
        ),
        "{text}",
    );
    assert!(
        text.contains("tandem arrays longer than this run types as callable: 50 bases (5.0%)"),
        "not \"will never call\", which claims a permanent refusal from a tunable threshold: \
         {text}",
    );
    assert_eq!(
        800 + 150 + 50,
        1_000,
        "the fixture's three parts are the whole, which is what the shares above add to",
    );
}

/// **Every span a person reads names its chromosome.**
///
/// `GenomeRegion`'s own `Display` writes `contig 0:120-190`, because a region carries no
/// reference and cannot spend a contig table it does not have. A run has one. Somebody whose
/// genome is `SL4.0ch01`…`SL4.0ch13` should not have to count from zero in a `.fai` to use a
/// report.
#[test]
fn every_span_the_report_shows_is_named_by_its_chromosome() {
    let (_zeta, _alpha, read_groups) = a_cohorts_read_groups();
    let written = a_run(
        3,
        0,
        vec![region(0, 120, 190)],
        vec![region(1, 40, 40)],
        vec![walked("zeta", ground_mostly_called(), Vec::new())],
    );

    let text = rendered(
        &written,
        &read_groups,
        &a_defaults_runs_parameters(&read_groups),
        1_000,
    );

    assert!(text.contains("chr1:120-190"), "{text}");
    assert!(text.contains("chr2:40-40"), "{text}");
    assert!(
        !text.contains("contig 0") && !text.contains("contig 1"),
        "no index reaches the page: {text}",
    );
}

/// **A refusal that did not happen gets a count and no advice**, because a line telling somebody
/// how to fix a thing that did not occur is a line they read and discard.
#[test]
fn a_refusal_that_did_not_happen_is_a_count_and_nothing_else() {
    let (_zeta, _alpha, read_groups) = a_cohorts_read_groups();
    let written = a_run(
        3,
        0,
        Vec::new(),
        Vec::new(),
        vec![walked("zeta", ground_mostly_called(), Vec::new())],
    );

    let text = rendered(
        &written,
        &read_groups,
        &a_defaults_runs_parameters(&read_groups),
        1_000,
    );

    assert!(
        text.contains("loci the merge declined to assemble for being too wide: 0"),
        "{text}",
    );
    assert!(
        !text.contains("--max-cohort-locus-span"),
        "nothing was refused, so nothing tells the reader to raise the bound: {text}",
    );
}

/// **A refusal that did happen says what to do about it**, which is what
/// `cohort_merge.md` §3.3 asks a non-zero count to lead a reader to — and shows a handful of
/// spans rather than all of them, since a badly-set bound refuses thousands.
#[test]
fn a_refusal_that_happened_shows_a_handful_of_spans_and_what_to_do() {
    let (_zeta, _alpha, read_groups) = a_cohorts_read_groups();
    let too_wide: Vec<GenomeRegion> = (0..8)
        .map(|n| region(0, 100 + n * 10, 160 + n * 10))
        .collect();
    let written = a_run(
        3,
        0,
        too_wide,
        Vec::new(),
        vec![walked("zeta", ground_mostly_called(), Vec::new())],
    );

    let text = rendered(
        &written,
        &read_groups,
        &a_defaults_runs_parameters(&read_groups),
        1_000,
    );

    assert!(
        text.contains("loci the merge declined to assemble for being too wide: 8"),
        "{text}",
    );
    assert_eq!(
        // **The indented lines alone.** The ground line above names `chr1:1-1000` too, so
        // counting every mention would count that one and read six where five were shown.
        text.lines()
            .filter(|line| line.starts_with("  chr1:"))
            .count(),
        SPANS_A_REPORT_SHOWS,
        "five spans shown, not eight: {text}",
    );
    assert!(text.contains("… and 3 more"), "{text}");
    assert!(
        text.contains("(61 bases)"),
        "each span carries its length, since the advice asks it to be compared with a bound: \
         {text}",
    );
    assert!(
        text.contains("the bound is --max-cohort-locus-span 50"),
        "and the bound is named as the flag a person types, with its value: {text}",
    );
}

/// **A sample with no reads in the analysed ground is named**, because in the VCF it looks
/// exactly like a sample that matched the reference everywhere.
///
/// It is not an error and the run does not stop for one: such a sample is still a sample of the
/// cohort and still carries a genotype, which the prior alone produced.
#[test]
fn a_sample_with_no_reads_is_named_and_the_run_does_not_treat_it_as_a_failure() {
    let (_zeta, _alpha, read_groups) = a_cohorts_read_groups();
    let written = a_run(
        3,
        0,
        Vec::new(),
        Vec::new(),
        vec![
            walked("zeta", ground_mostly_called(), Vec::new()),
            walked("alpha", ground_mostly_called(), Vec::new()),
        ],
    );

    let text = rendered(
        &written,
        &read_groups,
        &a_defaults_runs_parameters(&read_groups),
        1_000,
    );

    assert!(
        text.contains(
            "samples: 2 — 0 whose reads the caller used, 0 whose reads the filters took all of, \
             2 that contributed none"
        ),
        "{text}",
    );
    assert!(
        text.contains("no read reached the caller: zeta, alpha"),
        "{text}"
    );
    assert!(
        text.contains("written ./."),
        "and what the file says about them, which is a no-call and not a genotype: {text}",
    );
    assert!(
        !text.contains("from the prior alone"),
        "the loop calls such a sample; the file does not write that call (spec §7.1): {text}",
    );
}

/// **The read filters say why they dropped reads, per read group, and only the reasons that
/// fired.**
///
/// A read filter has nine reasons and a run trips two or three; printing the other six as zeros
/// is eight numbers a reader scans past to find the one that matters.
#[test]
fn the_report_names_only_the_filter_reasons_that_fired() {
    let (_zeta, _alpha, read_groups) = a_cohorts_read_groups();
    let filters = ReadFilterCounts {
        kept: 900,
        duplicate: 60,
        low_mapq: 40,
        ..ReadFilterCounts::default()
    };
    let written = a_run(
        3,
        0,
        Vec::new(),
        Vec::new(),
        vec![walked(
            "zeta",
            ground_mostly_called(),
            vec![(Some(ReadGroupId(0)), filters)],
        )],
    );

    let text = rendered(
        &written,
        &read_groups,
        &a_defaults_runs_parameters(&read_groups),
        1_000,
    );

    assert!(
        text.contains("zeta: 900 reads kept, 100 dropped by the read filters"),
        "the two totals, and 100 is 60 + 40: {text}",
    );
    assert!(
        text.contains("library rg-zeta of zeta: 60 duplicate, 40 mapping quality too low"),
        "the library as its file named it, not an index: {text}",
    );
    assert!(
        !text.contains("supplementary") && !text.contains("unmapped"),
        "the reasons that did not fire are not printed: {text}",
    );
}

/// **A read group that dropped nothing is not printed at all**, which is the same rule one step
/// up: a line saying a filter did nothing is a line to discard.
#[test]
fn a_read_group_that_dropped_nothing_gets_no_line_of_its_own() {
    let (_zeta, _alpha, read_groups) = a_cohorts_read_groups();
    let written = a_run(
        3,
        0,
        Vec::new(),
        Vec::new(),
        vec![walked(
            "zeta",
            ground_mostly_called(),
            vec![(
                Some(ReadGroupId(0)),
                ReadFilterCounts {
                    kept: 900,
                    ..ReadFilterCounts::default()
                },
            )],
        )],
    );

    let text = rendered(
        &written,
        &read_groups,
        &a_defaults_runs_parameters(&read_groups),
        1_000,
    );

    assert!(
        text.contains("zeta: 900 reads kept, 0 dropped by the read filters"),
        "{text}",
    );
    assert!(!text.contains("library rg1 of zeta:"), "{text}");
}

/// **A sample that read nothing gets no line about its read filters**, because the line naming
/// it as having no reads already said everything two zeroes would.
#[test]
fn a_sample_that_read_nothing_gets_no_filter_line() {
    let (_zeta, _alpha, read_groups) = a_cohorts_read_groups();
    let written = a_run(
        3,
        0,
        Vec::new(),
        Vec::new(),
        vec![walked("zeta", ground_mostly_called(), Vec::new())],
    );

    let text = rendered(
        &written,
        &read_groups,
        &a_defaults_runs_parameters(&read_groups),
        1_000,
    );

    assert!(text.contains("no read reached the caller: zeta"), "{text}");
    assert!(
        !text.contains("zeta: 0 reads kept"),
        "and nothing beyond it: {text}",
    );
}

/// **The line naming the parameters file and the line counting what was fitted are different
/// lines about different things**, and neither opens on the other's word — a reader scanning
/// for one must not land on the other.
#[test]
fn the_parameters_path_and_the_fitted_count_do_not_share_an_opening_word() {
    let (_zeta, _alpha, read_groups) = a_cohorts_read_groups();
    let written = a_run(
        3,
        0,
        Vec::new(),
        Vec::new(),
        vec![walked("zeta", ground_mostly_called(), Vec::new())],
    );

    let lines = lines_of(
        &written,
        &read_groups,
        &a_defaults_runs_parameters(&read_groups),
        1_000,
    );

    assert!(
        !lines.iter().any(|line| line.starts_with("parameters:")),
        "the report's own lines leave `parameters:` to the driver, which prints the path: {lines:#?}",
    );
}

/// **Which of the run's numbers rest on a measurement, and the ones that do not are named** —
/// `parameters_file.md` §8's question.
///
/// Named rather than counted, because a run whose contamination and slippage are compiled-in
/// constants is a different claim from one whose base-quality calibration is, and a count of
/// five says neither.
#[test]
fn the_report_names_the_groups_of_numbers_that_were_not_measured_here() {
    let (_zeta, _alpha, read_groups) = a_cohorts_read_groups();
    let written = a_run(
        3,
        0,
        Vec::new(),
        Vec::new(),
        vec![walked("zeta", ground_mostly_called(), Vec::new())],
    );

    let text = rendered(
        &written,
        &read_groups,
        &a_defaults_runs_parameters(&read_groups),
        1_000,
    );

    assert!(
        text.contains("numbers behind the calls: 0 of 7 groups the file says were fitted"),
        "a defaults run fitted none of the seven: {text}",
    );
    assert!(
        text.contains("the base-quality calibration")
            && text.contains("contamination")
            && text.contains("the inbreeding coefficients"),
        "and the ones it did not are named, in a reader's words: {text}",
    );
}

/// **A run given no ground to call over does not divide by it.** Every share is against the
/// analysed bases, and a run whose BED named nothing would otherwise report each part as a
/// not-a-number.
#[test]
fn a_run_over_no_ground_reports_no_share_rather_than_a_not_a_number() {
    assert_eq!(share_of(0, 0), "—");
    assert_eq!(share_of(150, 1_000), "15.0%");
}

/// **A BED that cuts a typed region does not make a share exceed a hundred per cent.**
///
/// A repeat tract is typed and walked **whole** even where a BED asks for part of it (spec §4.2:
/// findings whole, generic clipped), so the bases the walk was handed can exceed the bases asked
/// for. Measured in review: a BED of 120 bases inside two tracts charged 240 to *not built yet*,
/// and dividing by the 120 printed **200.0%**. The parts were right and the denominator was
/// wrong; the shares are now of the parts' own sum, and the two totals are printed together
/// whenever they differ.
#[test]
fn a_bed_that_cuts_a_typed_region_does_not_produce_a_share_above_a_hundred() {
    let (_zeta, _alpha, read_groups) = a_cohorts_read_groups();
    let cut_across_two_tracts = LocusCounts {
        regions_in: 4,
        regions_handled: 2,
        regions_handled_bp: 48,
        loci_emitted: 0,
        unhandled_not_implemented: 2,
        unhandled_not_implemented_bp: 192,
        unhandled_out_of_scope: 0,
        unhandled_out_of_scope_bp: 0,
    };
    let written = a_run(
        0,
        0,
        Vec::new(),
        Vec::new(),
        vec![walked("zeta", cut_across_two_tracts, Vec::new())],
    );

    let text = rendered(
        &written,
        &read_groups,
        &a_defaults_runs_parameters(&read_groups),
        120,
    );

    assert!(
        text.contains("120 bases asked for"),
        "what was asked for is still stated: {text}",
    );
    assert!(
        text.contains("it spoke for 240 bases; the shares below are of that"),
        "and so is the difference, and why: {text}",
    );
    assert!(text.contains("called: 48 bases (20.0%)"), "{text}");
    assert!(
        text.contains(
            "clusters of repeats too close together to have clean flanks: 192 bases (80.0%)"
        ),
        "not 160.0%: {text}",
    );
}

/// **Where the two totals agree, the report does not explain a difference there is none of.**
#[test]
fn a_whole_contig_run_says_nothing_about_two_totals() {
    let (_zeta, _alpha, read_groups) = a_cohorts_read_groups();
    let written = a_run(
        3,
        0,
        Vec::new(),
        Vec::new(),
        vec![walked("zeta", ground_mostly_called(), Vec::new())],
    );

    let text = rendered(
        &written,
        &read_groups,
        &a_defaults_runs_parameters(&read_groups),
        1_000,
    );

    assert!(
        !text.contains("it spoke for"),
        "800 + 150 + 50 is the 1,000 asked for, so there is nothing to explain: {text}",
    );
}

/// **A sample the filters emptied is not a sample with no reads**, and the report says which.
///
/// Three situations were written with one sentence and two of them were wrong: a sample whose
/// reads were all duplicates was called *no reads* four lines above the line saying it had 720 of
/// them. A geneticist checks the duplicate marking in one case and the sample sheet in the other.
#[test]
fn a_sample_the_filters_emptied_is_told_apart_from_one_that_had_no_reads() {
    let (_zeta, _alpha, read_groups) = a_cohorts_read_groups();
    let all_duplicates = ReadFilterCounts {
        kept: 0,
        duplicate: 720,
        ..ReadFilterCounts::default()
    };
    let written = a_run(
        3,
        0,
        Vec::new(),
        Vec::new(),
        vec![
            walked(
                "zeta",
                ground_mostly_called(),
                vec![(Some(ReadGroupId(0)), all_duplicates)],
            ),
            walked("alpha", ground_mostly_called(), Vec::new()),
        ],
    );

    let text = rendered(
        &written,
        &read_groups,
        &a_defaults_runs_parameters(&read_groups),
        1_000,
    );

    assert!(
        text.contains(
            "samples: 2 — 0 whose reads the caller used, 1 whose reads the filters took all of, \
             1 that contributed none"
        ),
        "{text}",
    );
    assert!(text.contains("every read filtered out: zeta"), "{text}");
    assert!(text.contains("no read reached the caller: alpha"), "{text}");
    assert!(
        !text.contains("no read reached the caller: zeta"),
        "zeta had 720 reads and the filters took them; saying it had none is the defect this \
         closes: {text}",
    );
}

/// **A run that built repeat-tract loci and could not score them says so, in loci**, and a run
/// that built none says nothing rather than printing a zero.
///
/// The two facts are easy to run together and a reader acts on each differently. The base lines
/// above say what ground *no generator looked at* — clusters of repeats with no clean flanks,
/// and tandem arrays too long to call. This one says what was looked at, built, merged across
/// the cohort, and then not scored, because nothing in the run scores a repeat tract yet
/// (`run_ssr_observations.md` §5). Since the tract slot was filled, a tract's bases count as
/// *called* on the base lines — so a run that printed only those would say a tract's ground was
/// called when nothing was called there.
#[test]
fn tract_loci_the_run_could_not_score_are_a_line_of_their_own_and_only_when_there_are_some() {
    let (_zeta, _alpha, read_groups) = a_cohorts_read_groups();
    let parameters = a_defaults_runs_parameters(&read_groups);
    let ground = vec![walked(
        "zeta",
        LocusCounts {
            regions_in: 10,
            regions_handled: 9,
            regions_handled_bp: 950,
            unhandled_not_implemented: 1,
            unhandled_not_implemented_bp: 50,
            ..LocusCounts::default()
        },
        Vec::new(),
    )];

    let none_built = a_run(3, 0, Vec::new(), Vec::new(), ground.clone());
    let text = rendered(&none_built, &read_groups, &parameters, 1_000);
    assert!(
        !text.contains("repeat tracts:"),
        "a run that built no tract prints no line for them, rather than a row of zeros: {text}",
    );

    let mut some_built = a_run(3, 0, Vec::new(), Vec::new(), ground.clone());
    some_built.calling.tracts = TractOutcomes {
        called: 40,
        not_periodic: 5,
        too_many_alleles: 3,
        without_whole_repeats: 1,
        bundles_set_aside: 7,
    };
    let text = rendered(&some_built, &read_groups, &parameters, 1_000);

    // **The headline is the sum and the share of it that was called**, so a reader who stops
    // after one line still knows how much of the tract ground the run spoke for.
    assert!(
        text.contains("repeat tracts: 56 built, of which 40 called"),
        "the five outcomes sum to the tracts built: {text}",
    );
    for (line, what) in [
        (
            "(notPeriodic): 5",
            "the reads do not vary in whole motif units",
        ),
        (
            "(tooManyAlleles): 3",
            "more sequences segregate than the cap admits",
        ),
        (
            "shorter than one copy of the motif: 1",
            "no rung on the stutter ladder",
        ),
        (
            "no clean flanks, which nothing builds a caller for yet: 7",
            "a bundle",
        ),
    ] {
        assert!(
            text.contains(line),
            "the report says {what} and how many, and got: {text}",
        );
    }
    assert!(
        text.contains("called: 950 bases (95.0%)"),
        "the tract's bases are called ground — which is exactly why the lines above have to \
         exist: {text}",
    );

    // **Each refusal's line appears only when it happened.** Four lines of zeros under a
    // headline is a report a reader stops reading.
    let mut only_called = a_run(3, 0, Vec::new(), Vec::new(), ground);
    only_called.calling.tracts = TractOutcomes {
        called: 40,
        ..TractOutcomes::default()
    };
    let text = rendered(&only_called, &read_groups, &parameters, 1_000);
    assert!(text.contains("repeat tracts: 40 built, of which 40 called"));
    for absent in [
        "notPeriodic",
        "tooManyAlleles",
        "shorter than one copy",
        "clean flanks, which",
    ] {
        assert!(
            !text.contains(absent),
            "a run where {absent} never happened does not print a zero for it: {text}",
        );
    }
}

/// **The five outcomes are a partition** — every tract-kind locus the merge built is in exactly
/// one of them, which is what makes the headline's sum a fact rather than an addition.
#[test]
fn the_tract_outcomes_sum_to_the_tracts_built() {
    let outcomes = TractOutcomes {
        called: 40,
        not_periodic: 5,
        too_many_alleles: 3,
        without_whole_repeats: 1,
        bundles_set_aside: 7,
    };
    assert_eq!(outcomes.built(), 56);
    assert_eq!(outcomes.refused_by_a_filter(), 8);
    assert_eq!(TractOutcomes::default().built(), 0);
}
