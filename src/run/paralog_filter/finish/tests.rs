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
use crate::paralog::CoverageFitConfig;
use crate::run::paralog_filter::{
    GenericLocusSample, SpillEntry, SpillFile, SpilledSamples, TargetFdr, WindowCoverage,
};
use crate::types::{ContigId, InbreedingF, Ploidy, Position};
use crate::vcf::{HeaderContig, VcfHeaderMetadata};
use crate::window_coverage::{CoverageByGcHistogram, SampleHistogram};

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
    // **The calling pass is given something to measure.** `calling` is the interval from the
    // spill being named to the fit, and in a fixture that interval is microseconds — small
    // enough that a clock started and read on one line looks the same as the real one. It did:
    // replacing the span with `Instant::now().elapsed()` gave 83 ns and passed every test.
    // `Instant` is monotonic, so a sleep here is a hard lower bound the mutation cannot reach.
    std::thread::sleep(THE_CALLING_PASS_IS_AT_LEAST);
    (spill, output)
}

/// How long [`a_spill_and_an_output`] holds the spill open before handing it over, so that the
/// calling pass has a duration a test can assert a floor on.
const THE_CALLING_PASS_IS_AT_LEAST: std::time::Duration = std::time::Duration::from_millis(20);

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
    let text = crate::vcf::header_text(&metadata_over(2));

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
    let refused = crate::run::paralog_filter::ParalogScoringContext::new(
        vec![a_fittable_histogram(); 2],
        &outbred(2),
        &crate::paralog::ParalogModelParams::default(),
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

/// A histogram whose one-copy depth peak sits at `peak_bin`, so a cohort built from several of
/// them has a **different** fitted depth per sample.
///
/// **The shape is symmetric about the peak on purpose**: the fit refuses a sample whose mode and
/// median disagree by more than half (spec §3.1's guard), and a symmetric histogram puts the two
/// together whichever bin the peak is moved to.
fn a_fittable_histogram_peaking_at(peak_bin: usize) -> SampleHistogram {
    let depth_bins = 40u32;
    let mut counts = vec![0u32; depth_bins as usize + 1];
    for (offset, count) in [(-2i64, 200u32), (-1, 900), (0, 2000), (1, 900), (2, 200)] {
        let bin = usize::try_from(i64::try_from(peak_bin).expect("a small bin") + offset)
            .expect("the peak is at least two bins from the bottom");
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

/// **The report says what each sample's coverage fit came to** — plan step D1 asks for "the fit's
/// outcome", and before this the report said only how many fits were accepted.
///
/// Without it the only route to the depth a record was compared against is inverting the run's own
/// ratios, and the obvious substitute — the median depth of the records the run wrote — is a
/// variant-site subset that runs well above one copy.
#[test]
fn the_report_says_what_each_sample_s_coverage_fit_came_to() {
    let (spill, output) = a_spill_and_an_output(
        "the_report_says_what_each_sample_s_coverage_fit_came_to",
        vec![a_record(100, 3, ONE_COPY_DEPTH)],
    );

    // **Two samples fitted at different depths, and one refused.** Different depths, because a
    // cohort whose fits were all equal could not tell a line that prints each sample's own number
    // from one that prints the first sample's number three times.
    let run = fit_score_and_write_the_calls(
        &spill,
        vec![
            a_fittable_histogram_peaking_at(10),
            a_fittable_histogram_peaking_at(16),
            SampleHistogram::NoWindowFinalised,
        ],
        &outbred(3),
        an_ordinary_request(),
        &output,
        metadata_over(3),
        diploid(),
    )
    .expect("the run finishes");

    let fits = run.scoring.what_each_fit_came_to();
    let first = fits[0].expect("the first sample was fitted");
    let second = fits[1].expect("the second sample was fitted");
    assert!(fits[2].is_none(), "the third sample's fit was refused");
    assert!(
        second.one_copy_depth > first.one_copy_depth + 2.0,
        "the fixture must give the two samples different depths, or this test cannot tell a \
         per-sample line from a repeated one: {} against {}",
        first.one_copy_depth,
        second.one_copy_depth
    );

    let lines = what_to_tell_the_operator(&run).join("\n");
    assert!(
        lines.contains(&format!(
            "sample SRR7279400 fitted one copy at {:.2} reads a window, scatter {:.3}",
            first.one_copy_depth, first.single_copy_depth_sd
        )),
        "the first sample's own fitted numbers must be in the report: {lines}"
    );
    assert!(
        lines.contains(&format!(
            "sample SRR7279401 fitted one copy at {:.2} reads a window, scatter {:.3}",
            second.one_copy_depth, second.single_copy_depth_sd
        )),
        "and the second sample's, which are different numbers: {lines}"
    );
    // A sample with no model has no fit to report, and says why instead.
    assert!(
        !lines.contains("SRR7279402 fitted one copy"),
        "a refused sample has no fitted numbers to print: {lines}"
    );
    assert!(
        lines.contains("sample SRR7279402 has no coverage model:"),
        "and it still says why it has none: {lines}"
    );
}

/// **A cohort too large to name sample by sample is summarised, not truncated silently.**
///
/// A run is up to several thousand samples (spec §4) and a report is read by a person, so past a
/// cap the line gives the spread instead. What it must not do is print the first ten and leave a
/// reader believing that is the cohort.
#[test]
fn a_cohort_larger_than_the_cap_is_reported_as_a_spread() {
    const SAMPLES: usize = 12;
    let (spill, output) = a_spill_and_an_output(
        "a_cohort_larger_than_the_cap_is_reported_as_a_spread",
        vec![a_record(100, SAMPLES, ONE_COPY_DEPTH)],
    );

    // Twelve samples, each fitted at its own depth, so the lowest, the median and the highest are
    // three different numbers and a line that confused them would fail.
    let run = fit_score_and_write_the_calls(
        &spill,
        (0..SAMPLES)
            .map(|sample| a_fittable_histogram_peaking_at(6 + sample))
            .collect(),
        &outbred(SAMPLES),
        an_ordinary_request(),
        &output,
        metadata_over(SAMPLES),
        diploid(),
    )
    .expect("the run finishes");

    let fits = run.scoring.what_each_fit_came_to();
    let mut depths: Vec<f64> = fits
        .iter()
        .map(|fit| fit.expect("every sample was fitted").one_copy_depth)
        .collect();
    assert_eq!(depths.len(), SAMPLES);
    depths.sort_by(f64::total_cmp);
    let (lowest, median, highest) = (depths[0], depths[(SAMPLES - 1) / 2], depths[SAMPLES - 1]);
    assert!(
        lowest < median && median < highest,
        "the fixture must spread the depths, or the spread line cannot be checked: {depths:?}"
    );

    let lines = what_to_tell_the_operator(&run).join("\n");
    assert!(
        lines.contains(&format!(
            "across the {SAMPLES} fitted sample(s): one copy is {median:.2} reads a window at \
             the median ({lowest:.2} to {highest:.2})"
        )),
        "past the cap the report gives the spread, with the median between the ends: {lines}"
    );
    assert!(
        !lines.contains("fitted one copy at"),
        "and names no sample individually, so a reader cannot mistake ten for the cohort: {lines}"
    );
}

/// Pull the duration and the share the report printed for one named pass, so a test can check
/// that a pass's name sits beside **its own** number.
///
/// **Why this exists.** An earlier version checked only that the four pass names appeared and that
/// the four shares summed to about a hundred. Both survive swapping two passes' figures, which is
/// the likeliest way this line goes wrong.
fn what_the_line_says_about(lines: &str, pass: &str) -> (String, f64) {
    // **Anchored on the separator before the pass name**, because the line opens with "time from
    // the start of the calling pass" — a bare search for "calling " finds that instead, and the
    // test then compares the total against the calling pass and fails for the wrong reason.
    let after = lines
        .split_once(pass)
        .unwrap_or_else(|| panic!("the time line names no pass {pass:?}: {lines}"))
        .1;
    let (duration, rest) = after
        .split_once(" (")
        .unwrap_or_else(|| panic!("no duration after {pass:?}: {lines}"));
    let share = rest
        .split_once("%)")
        .unwrap_or_else(|| panic!("no share after {pass:?}: {lines}"))
        .0;
    (
        duration.to_string(),
        share
            .parse()
            .unwrap_or_else(|_| panic!("share {share:?} is not a number: {lines}")),
    )
}

/// **The report says where the run's time went, and each pass's figure sits beside its own name** —
/// plan step D2 asks for the wall per pass, and spec §8 defers parallel scoring until a run says
/// whether it is worth threads.
#[test]
fn the_report_says_where_the_time_went_and_what_share_each_pass_took() {
    let (spill, output) = a_spill_and_an_output(
        "the_report_says_where_the_time_went_and_what_share_each_pass_took",
        vec![
            a_record(100, 1, ONE_COPY_DEPTH),
            a_record(200, 1, ONE_COPY_DEPTH),
        ],
    );

    let run = fit_score_and_write_the_calls(
        &spill,
        vec![a_fittable_histogram()],
        &outbred(1),
        an_ordinary_request(),
        &output,
        metadata_over(1),
        diploid(),
    )
    .expect("the run finishes");

    // **The calling pass really was clocked from the spill's birth.** The fixture holds the spill
    // open for a known interval first, so a clock started inside `fit_score_and_write_the_calls`
    // — which is what the span looked like to every earlier test — cannot reach this floor.
    let spent = run.spent;
    assert!(
        spent.calling >= THE_CALLING_PASS_IS_AT_LEAST,
        "the calling pass is timed from the spill being named, so it covers the {:?} the fixture \
         held it open; got {:?}",
        THE_CALLING_PASS_IS_AT_LEAST,
        spent.calling
    );
    // **Nothing here asserts an order between the four passes**, though the sleep above usually
    // makes the calling pass the largest. Under the full suite the writing pass contends for the
    // same filesystem as six thousand other tests and has been seen at 77 ms against calling's
    // 23 ms, so an assertion on their order fails on a loaded machine and passes on an idle one.
    // The per-pass checks below compare each printed figure against *that pass's own* measured
    // duration, which tells the four apart without needing an order between them.

    let lines = what_to_tell_the_operator(&run).join("\n");
    assert!(
        lines.contains("time from the start of the calling pass"),
        "the line must say where its clock starts: {lines}"
    );
    assert!(
        lines.contains("startup before the calling pass"),
        "and what it leaves out: {lines}"
    );

    // **Each pass's number, checked against that pass's own duration.** Swapping two passes'
    // figures leaves the shares summing to a hundred and every name present, so a check on the
    // sum alone cannot see it.
    let mut shares = Vec::new();
    for (pass, measured) in [
        (": calling ", spent.calling),
        (", fitting the coverage models ", spent.fitting),
        (", scoring ", spent.scoring),
        (", writing ", spent.writing),
    ] {
        let (printed, share) = what_the_line_says_about(&lines, pass);
        assert_eq!(
            printed,
            super::how_long(measured),
            "the duration printed for {pass} is not {pass}'s own: {lines}"
        );
        let expected = 100.0 * measured.as_secs_f64()
            / (spent.calling + spent.fitting + spent.scoring + spent.writing).as_secs_f64();
        assert!(
            (share - expected).abs() <= 1.0,
            "the share printed for {pass} is {share} and its own share is {expected:.1}: {lines}"
        );
        shares.push(share);
    }
    // The four are shares of one total: a line dividing by anything else would not sum to a
    // hundred, whichever pass it divided by.
    let summed: f64 = shares.iter().sum();
    assert!(
        (summed - 100.0).abs() <= 2.0,
        "the four shares must sum to about 100, got {summed}: {lines}"
    );
}

/// **The report says how much disk the parked records took**, which nothing else in the run does —
/// the spill holds every record's line uncompressed plus about ten bytes a sample, so an operator
/// cannot read it off the compressed output (B2 review M6).
#[test]
fn the_report_says_how_much_disk_the_parked_records_took() {
    // **Enough records that the printed figure has a digit that can be wrong.** Two records give
    // a few hundred bytes, and at that size every plausible unit and divisor renders the same
    // thing — so a test on two records pins the words and not the number.
    let records: Vec<SpillEntry> = (0u64..1_200)
        .map(|at| a_record(100 + at, 2, ONE_COPY_DEPTH))
        .collect();
    let (spill, output) = a_spill_and_an_output(
        "the_report_says_how_much_disk_the_parked_records_took",
        records,
    );

    let run = fit_score_and_write_the_calls(
        &spill,
        vec![a_fittable_histogram(), a_fittable_histogram()],
        &outbred(2),
        an_ordinary_request(),
        &output,
        metadata_over(2),
        diploid(),
    )
    .expect("the run finishes");

    // **The size the file really had, measured here rather than taken from the field.** Comparing
    // the printed figure against the field it was formatted from is self-consistent by
    // construction; the spill is gone by now, so the check is that the field is a real size and
    // that the line renders it in the unit it claims.
    let bytes = run
        .spill_bytes_on_disk
        .expect("pass one finished, so the size was read");
    assert!(
        bytes > 64 * 1024,
        "the fixture must exceed the 64 KiB write buffer, or a size read before the flush would \
         pass this test too; got {bytes} bytes"
    );
    let lines = what_to_tell_the_operator(&run).join("\n");
    let printed = lines
        .split_once("the parked records took ")
        .expect("the report gives the spill's size")
        .1
        .split_once(" on disk")
        .expect("the size line ends as it began")
        .0;
    // Rendered as kibibytes, and the digits are the real ones: a divisor of 1000 rather than 1024
    // moves the first decimal place at this size, and calling it MB would move the unit.
    assert_eq!(
        printed,
        format!("{:.1} kiB", bytes as f64 / 1024.0),
        "the size must be the file's own, in the unit it names: {lines}"
    );
    assert!(
        lines.contains("and are gone"),
        "and say the file does not outlive the run: {lines}"
    );
}

/// **A spill whose size could not be read says so**, rather than leaving the line out — the figure
/// is there so an operator can size a filesystem, and silence reads as an empty spill.
#[test]
fn a_spill_whose_size_is_unknown_says_so_rather_than_printing_nothing() {
    let (spill, output) = a_spill_and_an_output(
        "a_spill_whose_size_is_unknown_says_so_rather_than_printing_nothing",
        vec![a_record(100, 1, ONE_COPY_DEPTH)],
    );
    let mut run = fit_score_and_write_the_calls(
        &spill,
        vec![a_fittable_histogram()],
        &outbred(1),
        an_ordinary_request(),
        &output,
        metadata_over(1),
        diploid(),
    )
    .expect("the run finishes");

    run.spill_bytes_on_disk = None;
    let lines = what_to_tell_the_operator(&run).join("\n");
    assert!(
        lines.contains("size on disk could not be read"),
        "an unreadable size is said, not omitted: {lines}"
    );
    assert!(
        !lines.contains("the parked records took"),
        "and it does not also claim a size: {lines}"
    );
}

/// **Finishing pass one twice keeps the size it measured the first time.** The early return for an
/// already-finished spill sits above the `stat`; a refactor hoisting the `stat` out of the match
/// would clobber the size with a second read, and nothing else would notice.
#[test]
fn finishing_pass_one_twice_keeps_the_size_it_measured() {
    let (mut spill, _output) = a_spill_and_an_output(
        "finishing_pass_one_twice_keeps_the_size_it_measured",
        vec![a_record(100, 2, ONE_COPY_DEPTH)],
    );
    let first = spill.bytes_on_disk().expect("pass one finished");
    assert!(first > 0, "a parked record occupies bytes");
    spill.finish_writing().expect("finishing twice is allowed");
    assert_eq!(
        spill.bytes_on_disk(),
        Some(first),
        "the size is the one pass one measured, not a second reading"
    );
}
