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
        (0..samples).map(|i| format!("sample{i}")).collect(),
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
    assert!(
        lines.contains("sample 1 has no coverage model: no window was finalised"),
        "and name the first that was not, with its reason: {lines}"
    );
    assert!(
        lines.contains("sample 2 has no coverage model: every window this sample finalised held"),
        "and the second, whose reason is a different one: {lines}"
    );
    assert!(
        !lines.contains("sample 0 has no coverage model"),
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
        lines.contains("the documented fallback — this run's could not be fitted"),
        "the report must not present the fallback as a measurement: {lines}"
    );
    assert!(
        lines.contains("0 of 2 sample(s) with a coverage model"),
        "and must say the evidence was absent: {lines}"
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
