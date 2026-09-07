//! **What the run's report says the filter did, and what its header states.**
//!
//! The passes themselves are tested next door; this is about the words and the numbers a person
//! reads afterwards. Two things matter here. **A rejected sample is named with its reason**, not
//! counted — "the coverage model rested on 2 of 3" and "on 2 of 3, and the one that dropped out
//! was never covered" are different statements, and only the second tells an operator whether to
//! look at that sample's reads. And **the header states what the run was calibrated to**, which is
//! the only trace a dropped record leaves anywhere.

use std::fs;
use std::path::PathBuf;

use super::{WhatTheOperatorAskedFor, fit_score_and_write_the_calls, what_to_tell_the_operator};
use crate::ng::paralog::CoverageFitConfig;
use crate::ng::run::paralog_filter::{
    GenericLocusSample, SpillEntry, SpillFile, SpilledSamples, TargetFdr, WindowCoverage,
};
use crate::ng::types::{ContigId, InbreedingF, Ploidy, Position};
use crate::ng::vcf::{HeaderContig, VcfHeaderMetadata};
use crate::ng::window_coverage::{CoverageByGcHistogram, SampleHistogram};

/// A histogram a sample's depth model can be fitted from.
fn a_fittable_histogram() -> SampleHistogram {
    let depth_bins = 40u32;
    let mut counts = vec![0u32; depth_bins as usize + 1];
    for (bin, count) in [(8usize, 200u32), (9, 900), (10, 2000), (11, 900), (12, 200)] {
        counts[bin] = count;
    }
    SampleHistogram::Fitted(CoverageByGcHistogram {
        window_bp: 500,
        gc_bins: 1,
        depth_bin_width: 0.5,
        depth_bins,
        windows_folded: 4200,
        windows_under_the_floor: 0,
        counts,
    })
}

const ONE_COPY_DEPTH: f32 = 5.25;

fn a_record(position: u64, samples: usize, depth: f32) -> SpillEntry {
    let columns: String = (0..samples).map(|_| "\t0/1:5,5".to_string()).collect();
    SpillEntry {
        contig: ContigId(0),
        position: Position(position),
        is_repeat_tract: false,
        line: format!("chr0\t{position}\t.\tA\tG\t42.5\tPASS\tAN=2;DP=10\tGT:AD{columns}")
            .into_bytes(),
        samples: SpilledSamples::GenericLocus(
            (0..samples)
                .map(|_| GenericLocusSample {
                    window: WindowCoverage {
                        gc_fraction: 0.4,
                        mean_depth: depth,
                    },
                    ref_reads: 20,
                    alt_reads: 10,
                })
                .collect(),
        ),
    }
}

/// A record no sample can speak for: every window absent, so nothing is scored and its ratio is
/// not a number.
fn a_record_no_sample_covered(position: u64, samples: usize) -> SpillEntry {
    let mut entry = a_record(position, samples, ONE_COPY_DEPTH);
    entry.samples = SpilledSamples::GenericLocus(
        (0..samples)
            .map(|_| GenericLocusSample {
                window: WindowCoverage {
                    gc_fraction: f32::NAN,
                    mean_depth: f32::NAN,
                },
                ref_reads: 20,
                alt_reads: 10,
            })
            .collect(),
    );
    entry
}

fn scratch(name: &str) -> PathBuf {
    let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tmp")
        .join("paralog_finish_tests")
        .join(name);
    fs::create_dir_all(&directory).expect("a scratch directory");
    directory
}

fn a_spill_and_an_output(test_name: &str, entries: Vec<SpillEntry>) -> (SpillFile, PathBuf) {
    let output = scratch(test_name).join("cohort.vcf");
    let _ = fs::remove_file(&output);
    let _ = fs::remove_file(SpillFile::beside(&output).path());
    let mut spill = SpillFile::beside(&output);
    for entry in &entries {
        spill.append(entry).expect("the entry is appended");
    }
    spill.finish_writing().expect("the spill is flushed");
    (spill, output)
}

fn metadata_over(samples: usize) -> VcfHeaderMetadata {
    VcfHeaderMetadata::try_new(
        vec![HeaderContig {
            name: "chr0".to_string(),
            length: 1_000_000,
            md5: None,
        }],
        // **Names an operator would recognise**, because the report prints the name and not the
        // index — `sample0` would read as "sample sample0" and hide whether the name is used.
        (0..samples).map(|i| format!("SRR72794{i:02}")).collect(),
        String::new(),
        String::new(),
        String::new(),
    )
    .expect("a header")
}

fn outbred(n: usize) -> Vec<InbreedingF> {
    vec![InbreedingF::try_new(0.0).expect("zero is a coefficient"); n]
}

fn diploid() -> Ploidy {
    Ploidy::try_new(2).expect("two copies")
}

/// What an operator taking the defaults asks for: remove the records, at spec §3.6's target.
fn an_ordinary_request() -> WhatTheOperatorAskedFor {
    WhatTheOperatorAskedFor {
        target_fdr: TargetFdr::try_new(0.01).expect("a target"),
        tag_instead_of_dropping: false,
    }
}

/// **The header says what the filter was calibrated to, and declares what it may have applied.**
///
/// A dropped record leaves no trace in the file, so without these lines a reader cannot tell a run
/// that removed nothing from one that was never asked to remove anything.
#[test]
fn the_header_states_the_run_s_calibration_and_declares_the_filter() {
    let (spill, output) = a_spill_and_an_output(
        "the_header_states_the_run_s_calibration_and_declares_the_filter",
        vec![
            a_record(100, 3, ONE_COPY_DEPTH),
            a_record(200, 3, ONE_COPY_DEPTH * 2.0),
        ],
    );

    let run = fit_score_and_write_the_calls(
        &spill,
        (0..3).map(|_| a_fittable_histogram()).collect(),
        &outbred(3),
        an_ordinary_request(),
        &output,
        metadata_over(3),
        diploid(),
    )
    .expect("the run finishes");

    let written = fs::read_to_string(&output).expect("the calls are on disk");
    assert!(
        written.contains("##paralogFilter=target_fdr=0.0100;pi="),
        "the header must carry the calibration: {}",
        written.lines().take(10).collect::<Vec<_>>().join("\n")
    );
    assert!(
        written.contains("samples_with_coverage_model=3/3"),
        "including how much of the cohort the coverage evidence rested on"
    );
    assert!(written.contains("##INFO=<ID=PARALOG_LR,"));
    assert!(written.contains("##INFO=<ID=PARALOG_POST,"));
    assert!(written.contains("##FILTER=<ID=hiddenParalog,"));
    assert_eq!(run.did.written + run.did.dropped, 2);
}

/// **Zero means no filter, and it is the type that says so.**
///
/// **This is not decoration.** Handed a target of zero the scoring would not remove nothing: a
/// strongly duplicated record's tail false-discovery value underflows to exactly zero, and zero
/// is not above zero, so the most extreme records would go — measured next door in
/// `a_looser_target_removes_more_records`, where the strictest target a run can ask for still
/// removes four of seven. Until this rule lived in one place, all that stood between that and an
/// operator who asked for no filtering was a `> 0.0` comparison copied into two call sites.
#[test]
fn a_target_of_zero_is_no_filter_rather_than_a_filter_that_removes_nothing() {
    assert_eq!(
        WhatTheOperatorAskedFor::from_the_flags(0.0, false).expect("zero is a legal flag value"),
        None,
        "zero asks for no filter at all"
    );
    assert_eq!(
        WhatTheOperatorAskedFor::from_the_flags(0.0, true).expect("zero is a legal flag value"),
        None,
        "and the tag flag has nothing to act on, so it changes nothing"
    );

    let asked = WhatTheOperatorAskedFor::from_the_flags(0.01, true)
        .expect("a hundredth is a target")
        .expect("and it asks for the filter");
    assert!((asked.target_fdr.get() - 0.01).abs() < 1e-12);
    assert!(asked.tag_instead_of_dropping);

    // **A negative zero is a mistake, not a way of saying "off".** It compares equal to zero, so
    // a sign-blind rule would answer a wrong target with "the filter did not run".
    for refused in [-0.0, -1.0, 1.0, 5.0, f64::NAN, f64::INFINITY] {
        assert!(
            WhatTheOperatorAskedFor::from_the_flags(refused, false).is_err(),
            "{refused} is neither a target nor a way of turning the filter off"
        );
    }
}

/// **A run with the filter off declares none of it** — which is what keeps its header the header
/// it wrote before the filter existed, and spec §10's first oracle a property of the code.
#[test]
fn a_header_without_the_filter_declares_none_of_it() {
    let text = crate::ng::vcf::header_text(&metadata_over(2));

    assert!(!text.contains("##paralogFilter="));
    assert!(!text.contains("PARALOG_LR"));
    assert!(!text.contains("PARALOG_POST"));
    assert!(!text.contains("hiddenParalog"));
}

/// **A sample with no coverage model is named in the report, with the reason.**
///
/// The two samples that fail here fail differently — one was never covered, one had every window
/// refused — and a report that only counted them would say the same thing about both.
#[test]
fn the_report_names_each_sample_that_has_no_coverage_model_and_why() {
    let (spill, output) = a_spill_and_an_output(
        "the_report_names_each_sample_that_has_no_coverage_model_and_why",
        vec![a_record(100, 3, ONE_COPY_DEPTH)],
    );

    let run = fit_score_and_write_the_calls(
        &spill,
        vec![
            a_fittable_histogram(),
            SampleHistogram::NoWindowFinalised,
            SampleHistogram::EveryWindowUnderTheFloor,
        ],
        &outbred(3),
        an_ordinary_request(),
        &output,
        metadata_over(3),
        diploid(),
    )
    .expect("the run finishes");

    let lines = what_to_tell_the_operator(&run).join("\n");
    assert!(
        lines.contains("1 of 3 sample(s) with a coverage model"),
        "the report must say how much of the cohort was fitted: {lines}"
    );
    // **By name, not by index.** An index is a fact about the order the run's arguments were
    // typed in, and an operator would have to count their own command line to use it.
    assert!(
        lines.contains("sample SRR7279401 has no coverage model: no window was finalised"),
        "and name the first that was not, with its reason: {lines}"
    );
    assert!(
        lines.contains(
            "sample SRR7279402 has no coverage model: every window this sample finalised held"
        ),
        "and the second, whose reason is a different one: {lines}"
    );
    assert!(
        !lines.contains("SRR7279400 has no coverage model"),
        "the sample that was fitted is not named: {lines}"
    );
}

/// **A run whose rate could not be fitted says so in words**, rather than reporting the fallback
/// as though it had been measured.
#[test]
fn the_report_says_when_the_rate_is_the_fallback_and_not_this_run_s() {
    let (spill, output) = a_spill_and_an_output(
        "the_report_says_when_the_rate_is_the_fallback_and_not_this_run_s",
        vec![a_record(100, 2, ONE_COPY_DEPTH)],
    );

    // Every sample's coverage model is refused, so nothing can be scored and there is no rate to
    // fit — the state an operator most needs told about, because the file looks ordinary.
    let run = fit_score_and_write_the_calls(
        &spill,
        vec![SampleHistogram::NoWindowFinalised; 2],
        &outbred(2),
        an_ordinary_request(),
        &output,
        metadata_over(2),
        diploid(),
    )
    .expect("the run finishes");

    assert_eq!(run.verdicts.records_in_the_fit, 0);
    assert_eq!(run.did.unscored, 1);
    let lines = what_to_tell_the_operator(&run).join("\n");
    assert!(
        lines.contains("the documented fallback"),
        "the report must not present the fallback as a measurement: {lines}"
    );
    // **The warning C3 wrote, and the half the rate line cannot carry**: that the records this
    // run removed were calibrated against a constant rather than against this cohort.
    assert!(
        lines.contains("warning: the hidden-duplication rate could not be fitted"),
        "the warning must be printed, not merely available: {lines}"
    );
    assert!(
        lines.contains("calibrated against that rate and not against this run"),
        "and it must say what that means for the file: {lines}"
    );
    assert!(
        lines.contains("0 of 2 sample(s) with a coverage model"),
        "and the report must say the evidence was absent: {lines}"
    );
}

/// **The report's four counts are four different numbers, and each is the one it says it is.**
///
/// The fixture is deliberately asymmetric — one record removed, one written, one of them
/// unscorable — because on a fixture where two of the counts happen to be equal, swapping them
/// changes nothing. Two such swaps survived the first version of this suite.
#[test]
fn the_report_counts_what_went_and_what_stayed() {
    let entries = vec![
        // Duplication-shaped: over-covered, reads split at a third. This is the one that goes.
        a_record(100, 3, ONE_COPY_DEPTH * 2.0),
        // One-copy depth: kept, and scored.
        a_record(200, 3, ONE_COPY_DEPTH),
        // No usable window anywhere: kept, and unscored.
        a_record_no_sample_covered(300, 3),
    ];

    let (dropping, dropping_output) = a_spill_and_an_output(
        "the_report_counts_what_went_and_what_stayed__drop",
        entries.clone(),
    );
    let dropped = fit_score_and_write_the_calls(
        &dropping,
        (0..3).map(|_| a_fittable_histogram()).collect(),
        &outbred(3),
        an_ordinary_request(),
        &dropping_output,
        metadata_over(3),
        diploid(),
    )
    .expect("the run finishes");

    assert_eq!(dropped.did.dropped, 1, "the duplication-shaped record goes");
    assert_eq!(dropped.did.tagged, 0, "and nothing is tagged when dropping");
    assert_eq!(dropped.did.written, 2, "the other two are written");
    assert_eq!(
        dropped.did.unscored, 1,
        "one of which no sample could speak for"
    );

    let lines = what_to_tell_the_operator(&dropped).join("\n");
    assert!(
        lines.contains("1 record(s) dropped, 0 tagged, 2 written; 2 were scored, 1 could not be"),
        "each count must be the one it is named for: {lines}"
    );

    // **The same records with the tag flag**: nothing is dropped, the same one is tagged, and the
    // file holds all three. A pass three called with the flag hardcoded would give the drop counts.
    let (tagging, tagging_output) =
        a_spill_and_an_output("the_report_counts_what_went_and_what_stayed__tag", entries);
    let tagged = fit_score_and_write_the_calls(
        &tagging,
        (0..3).map(|_| a_fittable_histogram()).collect(),
        &outbred(3),
        WhatTheOperatorAskedFor {
            tag_instead_of_dropping: true,
            ..an_ordinary_request()
        },
        &tagging_output,
        metadata_over(3),
        diploid(),
    )
    .expect("the run finishes");

    assert_eq!(tagged.did.dropped, 0, "tagging drops nothing");
    assert_eq!(
        tagged.did.tagged, 1,
        "and tags what dropping would have dropped"
    );
    assert_eq!(tagged.did.written, 3, "so every record is in the file");
}

/// **The header says how much of the cohort the coverage evidence rested on, in that order.**
///
/// The fixture fits one sample of three, so the two numbers differ and a swap reads as `3/1`.
/// It also fits a rate, so `em_converged` is `true` here where every other fixture in this file
/// has it `false` — a negation would otherwise be invisible.
#[test]
fn the_header_says_how_much_of_the_cohort_was_fitted_and_whether_the_rate_was() {
    let (spill, output) = a_spill_and_an_output(
        "the_header_says_how_much_of_the_cohort_was_fitted_and_whether_the_rate_was",
        vec![
            a_record(100, 3, ONE_COPY_DEPTH),
            a_record(200, 3, ONE_COPY_DEPTH * 2.0),
        ],
    );

    let run = fit_score_and_write_the_calls(
        &spill,
        vec![
            a_fittable_histogram(),
            SampleHistogram::NoWindowFinalised,
            SampleHistogram::NoWindowFinalised,
        ],
        &outbred(3),
        an_ordinary_request(),
        &output,
        metadata_over(3),
        diploid(),
    )
    .expect("the run finishes");

    let written = fs::read_to_string(&output).expect("the calls are on disk");
    assert!(
        written.contains("samples_with_coverage_model=1/3"),
        "one of three was fitted, in that order: {}",
        written.lines().take(10).collect::<Vec<_>>().join("\n")
    );
    assert!(
        written.contains("em_converged=true"),
        "and this run's rate was fitted, not fallen back to"
    );
    assert!(run.verdicts.calibration.prior.converged);
}

/// **The report writes the rate and the cut to the precisions the header uses**, so the two can
/// be read against each other.
#[test]
fn the_report_writes_the_rate_and_the_cut_to_their_precisions() {
    let (spill, output) = a_spill_and_an_output(
        "the_report_writes_the_rate_and_the_cut_to_their_precisions",
        vec![a_record(100, 2, ONE_COPY_DEPTH)],
    );

    let run = fit_score_and_write_the_calls(
        &spill,
        vec![SampleHistogram::NoWindowFinalised; 2],
        &outbred(2),
        an_ordinary_request(),
        &output,
        metadata_over(2),
        diploid(),
    )
    .expect("the run finishes");

    let lines = what_to_tell_the_operator(&run).join("\n");
    assert!(
        lines.contains("duplication rate 0.030000 "),
        "the rate carries six decimals, as the header's does: {lines}"
    );
    let cut = lines
        .split("cut at ")
        .nth(1)
        .and_then(|rest| rest.split(',').next())
        .expect("the report states a cut");
    assert_eq!(
        cut.split_once('.').expect("a decimal point").1.len(),
        4,
        "and the cut carries four, as a record's own ratio does: {cut}"
    );
}

/// **A run whose ratios ran past the ends of the histogram says how many did.**
///
/// The score is a sum over the samples weighed, so it grows with the cohort while the range does
/// not; a run where every record saturates has one bin for the cut to work with. Nothing acts on
/// the count yet, and the plan's D2 and D3 are what report it on real data — but a line that is
/// never printed cannot be read there either.
#[test]
fn the_report_says_when_ratios_sat_past_the_ends_of_the_range() {
    let (spill, output) = a_spill_and_an_output(
        "the_report_says_when_ratios_sat_past_the_ends_of_the_range",
        vec![a_record(100, 6, ONE_COPY_DEPTH * 2.0)],
    );

    let run = fit_score_and_write_the_calls(
        &spill,
        (0..6).map(|_| a_fittable_histogram()).collect(),
        &outbred(6),
        an_ordinary_request(),
        &output,
        metadata_over(6),
        diploid(),
    )
    .expect("the run finishes");

    assert_eq!(
        run.verdicts.ratios_outside_the_histogram, 1,
        "six samples all at twice one copy score past the range's top; the ratio was {:?}",
        run.verdicts.ratios
    );
    let lines = what_to_tell_the_operator(&run).join("\n");
    assert!(
        lines.contains("1 scored record(s) sat past the ends of the ratio range"),
        "and the report must say so: {lines}"
    );
}

/// **A fit configured with knobs that are not a fit fails the run**, rather than reporting that
/// every sample's coverage was bad.
#[test]
fn a_fit_that_is_not_configured_fails_the_run_rather_than_the_cohort() {
    let (spill, output) = a_spill_and_an_output(
        "a_fit_that_is_not_configured_fails_the_run_rather_than_the_cohort",
        vec![a_record(100, 2, ONE_COPY_DEPTH)],
    );

    // A single-copy band that is not a band.
    let refused = crate::ng::run::paralog_filter::ParalogScoringContext::new(
        vec![a_fittable_histogram(); 2],
        &outbred(2),
        &crate::ng::paralog::ParalogModelParams::default(),
        &CoverageFitConfig {
            single_copy_lo: 1.6,
            single_copy_hi: 0.4,
            ..CoverageFitConfig::default()
        },
    );
    assert!(
        refused.is_err(),
        "a band whose ends are the wrong way round is not a configuration"
    );
    drop((spill, output));
}
