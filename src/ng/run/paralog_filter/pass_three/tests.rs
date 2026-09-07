//! **What pass three has to get right about a line it did not write.**
//!
//! Two things carry this step. **A record the verdict does not touch comes back byte for byte** —
//! spec §10's oracle rests on it, and `a_record_the_verdict_does_not_touch_is_written_unchanged`
//! is where it is checked against the input's own bytes rather than against a re-encoding. And
//! **a verdict reaches a record only by position**: the ratios carry no contig and no position, so
//! `the_verdict_follows_the_file_s_order_and_not_the_ratios_sorted_order` parks records whose
//! evidence does not rise with the file, and checks that the *right* one goes.
//!
//! **Which dimensions these fixtures vary.** Whether a record is removed, kept or tagged; whether
//! it was scored at all; whether the operator asked to drop or to tag; the contig, so a place
//! built from the wrong entry shows up as an ordering refusal; the `FILTER` column's prior
//! contents, since a record already carrying one is joined and not replaced; and the number of
//! records, at zero, one and several.

use std::fs;
use std::path::PathBuf;

use super::{PassThreeError, WhatTheFilterDid, write_the_records_the_filter_kept};
use crate::ng::paralog::{
    CalibrationConfig, ParalogCalibration, ParalogFdrCurve, ParalogPrior,
    calibrate_from_the_ratio_histogram,
};
use crate::ng::run::paralog_filter::{
    GenericLocusSample, LrHistogramShape, ParalogVerdicts, SpillEntry, SpillFile, SpilledSamples,
    WindowCoverage,
};
use crate::ng::types::{ContigId, Ploidy, Position};
use crate::ng::vcf::{HeaderContig, VcfHeaderMetadata, VcfWriter};

// ---------------------------------------------------------------- the records

/// One sample's row; the numbers never reach the scorer here, because pass three does not score.
fn a_row() -> GenericLocusSample {
    GenericLocusSample {
        window: WindowCoverage {
            gc_fraction: 0.41,
            mean_depth: 6.25,
        },
        ref_reads: 5,
        alt_reads: 5,
    }
}

/// A parked record whose line is a real nine-column VCF line with one sample column.
fn a_record(contig: u32, position: u64, filter: &str) -> SpillEntry {
    SpillEntry {
        contig: ContigId(contig),
        position: Position(position),
        is_repeat_tract: false,
        line: format!(
            "chr{contig}\t{position}\t.\tA\tG\t42.5\t{filter}\tAN=2;DP=10\tGT:AD\t0/1:5,5"
        )
        .into_bytes(),
        samples: SpilledSamples::GenericLocus(vec![a_row()]),
    }
}

// ---------------------------------------------------------------- the verdicts

/// **A calibration built from the ratios the fixture will carry**, so the cut is a real cut over
/// this run's own numbers rather than a hand-written threshold.
fn verdicts_over(ratios: Vec<f64>, target_fdr: f64) -> ParalogVerdicts {
    let lr_histogram = LrHistogramShape::shipped();
    let mut histogram = lr_histogram.histogram();
    for &ratio in &ratios {
        histogram.push(ratio);
    }
    let records_in_the_fit = histogram.total();
    ParalogVerdicts {
        calibration: calibrate_from_the_ratio_histogram(
            &histogram,
            target_fdr,
            &CalibrationConfig::default(),
        ),
        ratios,
        records_in_the_fit,
        ratios_outside_the_histogram: 0,
        config: CalibrationConfig::default(),
        lr_histogram,
    }
}

/// **A calibration that removes every finite ratio**, whatever it is: the prior is a certainty, so
/// the tail false-discovery value is zero everywhere and the target reaches it.
fn verdicts_that_remove_everything(ratios: Vec<f64>) -> ParalogVerdicts {
    let lr_histogram = LrHistogramShape::shipped();
    let mut histogram = lr_histogram.histogram();
    for &ratio in &ratios {
        histogram.push(ratio);
    }
    let records_in_the_fit = histogram.total();
    let prior = ParalogPrior {
        prior_probability: 1.0 - f64::EPSILON,
        converged: true,
    };
    let curve = ParalogFdrCurve::from_histogram(&histogram, &prior);
    ParalogVerdicts {
        calibration: ParalogCalibration {
            prior,
            curve,
            lr_threshold: Some(-100.0),
            target_fdr: 0.99,
        },
        ratios,
        records_in_the_fit,
        ratios_outside_the_histogram: 0,
        config: CalibrationConfig::default(),
        lr_histogram,
    }
}

/// A calibration that removes nothing: an unreachable target.
fn verdicts_that_remove_nothing(ratios: Vec<f64>) -> ParalogVerdicts {
    verdicts_over(ratios, 0.0)
}

// ---------------------------------------------------------------- the files

fn scratch(name: &str) -> PathBuf {
    let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tmp")
        .join("paralog_pass_three_tests")
        .join(name);
    fs::create_dir_all(&directory).expect("a scratch directory");
    directory
}

/// A spill holding `entries`, and the output path the VCF will be written to.
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

fn a_writer(output: &std::path::Path) -> VcfWriter {
    let metadata = VcfHeaderMetadata::try_new(
        vec![
            HeaderContig {
                name: "chr0".to_string(),
                length: 1_000_000,
                md5: None,
            },
            HeaderContig {
                name: "chr1".to_string(),
                length: 900_000,
                md5: None,
            },
        ],
        vec!["one".to_string()],
        String::new(),
        String::new(),
        String::new(),
    )
    .expect("a header");
    VcfWriter::create(output, metadata, Ploidy::try_new(2).expect("two copies")).expect("a writer")
}

/// The record lines of a written VCF, without the header.
fn record_lines(output: &std::path::Path) -> Vec<String> {
    fs::read_to_string(output)
        .expect("the calls are on disk")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(ToString::to_string)
        .collect()
}

/// Run pass three over a spill and hand back what it did and what it wrote.
fn write_and_read_back(
    spill: &SpillFile,
    output: &std::path::Path,
    verdicts: &ParalogVerdicts,
    tag: bool,
) -> (WhatTheFilterDid, Vec<String>) {
    let mut writer = a_writer(output);
    let did = write_the_records_the_filter_kept(spill, verdicts, tag, &mut writer)
        .expect("the spill is readable and the lines patch");
    writer.finish().expect("the calls are finished");
    (did, record_lines(output))
}

// ---------------------------------------------------------------- the oracle

/// **A record the verdict does not touch is written byte for byte as it was parked.**
///
/// This is spec §10's oracle at the level of one line: with the two `INFO` keys stripped, a
/// filtered run's records must be the unfiltered run's exactly. The comparison is against the
/// entry's own bytes, not against a re-encoding, because a re-encoding would agree with a writer
/// that rebuilt the line wrongly in the same way.
#[test]
fn a_record_the_verdict_does_not_touch_is_written_unchanged() {
    let parked = a_record(0, 100, "PASS");
    let untouched = String::from_utf8(parked.line.clone()).expect("the line is text");
    let (spill, output) = a_spill_and_an_output(
        "a_record_the_verdict_does_not_touch_is_written_unchanged",
        vec![parked],
    );

    // An unscored record: no verdict, no INFO fields, nothing to add.
    let verdicts = verdicts_that_remove_nothing(vec![f64::NAN]);
    let (did, lines) = write_and_read_back(&spill, &output, &verdicts, false);

    assert_eq!(lines, vec![untouched], "the line must come back unchanged");
    assert_eq!(
        did,
        WhatTheFilterDid {
            written: 1,
            dropped: 0,
            tagged: 0,
            unscored: 1,
        }
    );
}

/// **A scored record the cut keeps gains both `INFO` fields and nothing else.**
#[test]
fn a_scored_record_the_cut_keeps_carries_its_ratio_and_its_probability() {
    let parked = a_record(0, 100, "PASS");
    let (spill, output) = a_spill_and_an_output(
        "a_scored_record_the_cut_keeps_carries_its_ratio_and_its_probability",
        vec![parked],
    );

    let verdicts = verdicts_that_remove_nothing(vec![-4.25]);
    let (did, lines) = write_and_read_back(&spill, &output, &verdicts, false);

    assert_eq!(did.written, 1);
    assert_eq!(did.dropped, 0);
    assert_eq!(did.unscored, 0, "the record was scored");
    let line = &lines[0];
    assert!(
        line.contains("PARALOG_LR=-4.2500"),
        "the ratio goes in as written, to the decimals the header's cut uses: {line}"
    );
    assert!(
        line.contains("PARALOG_POST="),
        "and the probability it implies: {line}"
    );
    assert!(
        line.contains("AN=2;DP=10;PARALOG_LR="),
        "appended to the INFO already there, not replacing it: {line}"
    );
    assert!(
        line.split('\t').nth(6) == Some("PASS"),
        "and the FILTER column is untouched where nothing was flagged: {line}"
    );
}

// ---------------------------------------------------------------- drop and tag

/// **Removed means the line is not in the file**, and the count says so.
#[test]
fn a_removed_record_is_not_written() {
    let (spill, output) = a_spill_and_an_output(
        "a_removed_record_is_not_written",
        vec![a_record(0, 100, "PASS"), a_record(0, 200, "PASS")],
    );

    let verdicts = verdicts_that_remove_everything(vec![12.0, f64::NAN]);
    let (did, lines) = write_and_read_back(&spill, &output, &verdicts, false);

    assert_eq!(
        did,
        WhatTheFilterDid {
            written: 1,
            dropped: 1,
            tagged: 0,
            unscored: 1,
        },
        "the scored record goes, the unscored one stays"
    );
    assert_eq!(lines.len(), 1);
    assert!(
        lines[0].contains("\t200\t"),
        "the record that stayed is the one that was not scored: {}",
        lines[0]
    );
}

/// **The tag-mode file is the drop-mode file plus the removed lines, on the filter's id** —
/// spec §10's third relation, at the level of these two records.
#[test]
fn tagging_writes_what_dropping_leaves_out_on_the_filter_s_id() {
    let entries = vec![a_record(0, 100, "PASS"), a_record(0, 200, "PASS")];
    let ratios = vec![12.0, f64::NAN];

    let (dropping_spill, dropping_output) = a_spill_and_an_output(
        "tagging_writes_what_dropping_leaves_out__drop",
        entries.clone(),
    );
    let (_, dropped_lines) = write_and_read_back(
        &dropping_spill,
        &dropping_output,
        &verdicts_that_remove_everything(ratios.clone()),
        false,
    );

    let (tagging_spill, tagging_output) =
        a_spill_and_an_output("tagging_writes_what_dropping_leaves_out__tag", entries);
    let (did, tagged_lines) = write_and_read_back(
        &tagging_spill,
        &tagging_output,
        &verdicts_that_remove_everything(ratios),
        true,
    );

    assert_eq!(
        did,
        WhatTheFilterDid {
            written: 2,
            dropped: 0,
            tagged: 1,
            unscored: 1,
        },
        "nothing is dropped when the operator asked to tag"
    );
    let kept: Vec<&String> = tagged_lines
        .iter()
        .filter(|line| !line.contains("hiddenParalog"))
        .collect();
    assert_eq!(
        kept,
        dropped_lines.iter().collect::<Vec<_>>(),
        "the tag-mode file with its tagged lines removed is the drop-mode file"
    );
    assert_eq!(
        tagged_lines[0].split('\t').nth(6),
        Some("hiddenParalog"),
        "and the tagged line carries the filter's id: {}",
        tagged_lines[0]
    );
}

/// **A record already carrying a filter is joined, not replaced.**
///
/// The calling loop marks a record `EMNoConv` when its pass cap ran out; a record that is also a
/// hidden duplication is both, and replacing the first would lose a fact the run established.
#[test]
fn a_record_that_already_carries_a_filter_is_joined_not_replaced() {
    let (spill, output) = a_spill_and_an_output(
        "a_record_that_already_carries_a_filter_is_joined_not_replaced",
        vec![a_record(0, 100, "EMNoConv")],
    );

    let verdicts = verdicts_that_remove_everything(vec![12.0]);
    let (did, lines) = write_and_read_back(&spill, &output, &verdicts, true);

    assert_eq!(did.tagged, 1);
    assert_eq!(
        lines[0].split('\t').nth(6),
        Some("EMNoConv;hiddenParalog"),
        "both filters, in the order they were applied: {}",
        lines[0]
    );
}

// ---------------------------------------------------------------- the pairing

/// **A verdict reaches a record by position in the file, and by nothing else.**
///
/// The ratios carry no contig and no position, so the only thing pairing them with records is the
/// order of the walk. Here the removed record is neither first nor last and its ratio is not the
/// largest, so a pass that sorted, reversed or shifted the ratios removes a different record —
/// and the file says which.
#[test]
fn the_verdict_follows_the_file_s_order_and_not_the_ratios_sorted_order() {
    let (spill, output) = a_spill_and_an_output(
        "the_verdict_follows_the_file_s_order_and_not_the_ratios_sorted_order",
        vec![
            a_record(0, 100, "PASS"),
            a_record(0, 200, "PASS"),
            a_record(0, 300, "PASS"),
        ],
    );

    // **Measured on this fixture**, the three records' tail false-discovery values are
    // `0.667`, `1.17e-5` and `0.500`; at a target of three in ten only the middle one is
    // reached. The third is deliberately close — at a target of one in two it goes too — so the
    // fixture is not passing because the other two are far away.
    let verdicts = verdicts_over(vec![-30.0, 12.0, -8.0], 0.3);

    let (did, lines) = write_and_read_back(&spill, &output, &verdicts, true);

    assert_eq!(did.written, 3, "tagging writes them all");
    let tagged: Vec<&String> = lines
        .iter()
        .filter(|line| line.contains("hiddenParalog"))
        .collect();
    assert_eq!(tagged.len(), 1, "one record is removed by this cut");
    assert!(
        tagged[0].contains("\t200\t"),
        "and it is the middle one, whose ratio it is: {}",
        tagged[0]
    );
}

/// **A spill and a ratio vector that do not describe the same run are refused.**
///
/// Nothing in the types says the two came from one run, and pairing them by position would give
/// every record after the difference its neighbour's verdict — a file that looks entirely normal.
#[test]
fn a_ratio_vector_that_is_not_the_spill_s_length_is_refused() {
    let (spill, output) = a_spill_and_an_output(
        "a_ratio_vector_that_is_not_the_spill_s_length_is_refused",
        vec![a_record(0, 100, "PASS"), a_record(0, 200, "PASS")],
    );

    let verdicts = verdicts_that_remove_nothing(vec![-1.0]);
    let mut writer = a_writer(&output);
    let refused = write_the_records_the_filter_kept(&spill, &verdicts, false, &mut writer)
        .expect_err("one ratio cannot answer for two records");

    let PassThreeError::RatiosDoNotMatchTheSpill { records, ratios } = refused else {
        panic!("expected a pairing failure, got {refused:?}")
    };
    assert_eq!(records, 2);
    assert_eq!(ratios, 1);
}

// ---------------------------------------------------------------- the ends

/// **A run that called nothing writes a header and no records**, rather than failing.
#[test]
fn a_spill_with_no_records_writes_an_empty_file() {
    let (spill, output) =
        a_spill_and_an_output("a_spill_with_no_records_writes_an_empty_file", vec![]);

    let verdicts = verdicts_that_remove_nothing(Vec::new());
    let (did, lines) = write_and_read_back(&spill, &output, &verdicts, false);

    assert_eq!(did, WhatTheFilterDid::default());
    assert!(lines.is_empty());
}

/// **Records on a second contig are written in the order the spill holds them**, and the place
/// each is written at comes from its own entry.
///
/// Every other fixture here sits on one contig, where a place built from the wrong entry would
/// still be non-decreasing and pass the writer's ordering check.
#[test]
fn records_on_a_second_contig_are_written_in_the_spill_s_order() {
    let (spill, output) = a_spill_and_an_output(
        "records_on_a_second_contig_are_written_in_the_spill_s_order",
        vec![
            a_record(0, 900, "PASS"),
            a_record(1, 100, "PASS"),
            a_record(1, 200, "PASS"),
        ],
    );

    let verdicts = verdicts_that_remove_nothing(vec![-1.0, -2.0, -3.0]);
    let (did, lines) = write_and_read_back(&spill, &output, &verdicts, false);

    assert_eq!(did.written, 3);
    let positions: Vec<&str> = lines
        .iter()
        .map(|line| line.split('\t').next().expect("a contig"))
        .collect();
    assert_eq!(positions, vec!["chr0", "chr1", "chr1"]);
}
