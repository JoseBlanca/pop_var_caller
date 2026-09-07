//! **The one thing this pass must not do: turn "we could not score this" into "this is not a
//! duplication".**
//!
//! A record no sample can speak for is unscored — kept, never flagged, and out of the fit. It
//! carries `NaN`, and the histogram the duplication rate is fitted from refuses to fold a `NaN`.
//! A zero in its place is a record voting against duplication, so a run where nothing could be
//! scored would come back with a confident rate fitted from evidence that does not exist, and a
//! VCF that looks right (spec §6 trap 4). Four tests here are about that single value:
//! `a_record_no_sample_can_speak_for_is_unscored`,
//! `a_record_whose_samples_have_models_but_no_usable_spread_is_unscored`,
//! `the_histogram_folds_exactly_the_finite_ratios`, and
//! `an_empty_spill_falls_back_rather_than_fitting_a_rate_from_nothing`.
//!
//! **Which dimensions these fixtures vary, and which they do not.** They vary whether a record
//! is scorable, where the unscorable one sits in the file, whether a record is a repeat tract or
//! a generic locus, how many samples the run has (one and six), and whether the rate estimate
//! settles. They hold GC content constant — every histogram here has one GC bin — because GC
//! enters through the coverage model and nothing in this pass reads it; the tests that vary GC
//! are in `scoring_context/tests.rs`, where a one-bin fixture was the defect that hid a wrong
//! argument.

use std::fs;
use std::path::PathBuf;

use super::{
    PassTwoError, calibrate_from_the_ratio_histogram, score_the_parked_records_and_resolve_the_cut,
};
use crate::ng::paralog::{
    CalibrationConfig, CoverageFitConfig, DEFAULT_FALLBACK_PARALOG_PRIOR, EmConfig,
    ParalogLrHistogram, ParalogModelParams,
};
use crate::ng::run::paralog_filter::{
    GenericLocusSample, ParalogScoringContext, RepeatTractSample, SpillEntry, SpillFile,
    SpilledSamples,
};
use crate::ng::types::{ContigId, InbreedingF, Position};
use crate::ng::window_coverage::{CoverageByGcHistogram, SampleHistogram, WindowCoverage};

// ---------------------------------------------------------------- the run's samples

/// A histogram a sample's depth model can be fitted from: a single-copy peak well inside the
/// regular depth bins, at one GC bin — see the module note on what this fixture cannot test.
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

/// The depth this histogram's model calls one copy — the peak bin's centre.
const ONE_COPY_DEPTH: f32 = 5.25;

fn outbred(n: usize) -> Vec<InbreedingF> {
    vec![InbreedingF::try_new(0.0).expect("zero is a coefficient"); n]
}

fn a_context_over(n: usize) -> ParalogScoringContext {
    a_context_over_with(n, &CoverageFitConfig::default())
}

fn a_context_over_with(n: usize, fit: &CoverageFitConfig) -> ParalogScoringContext {
    ParalogScoringContext::new(
        (0..n).map(|_| a_fittable_histogram()).collect(),
        &outbred(n),
        &ParalogModelParams::default(),
        fit,
    )
    .expect("the fit configuration is a usable one")
}

// ---------------------------------------------------------------- the parked records

fn a_window(mean_depth: f32) -> WindowCoverage {
    WindowCoverage {
        gc_fraction: 0.4,
        mean_depth,
    }
}

fn no_window() -> WindowCoverage {
    WindowCoverage {
        gc_fraction: f32::NAN,
        mean_depth: f32::NAN,
    }
}

/// A generic locus at `position`, one row per `(window, ref_reads, alt_reads)`.
fn a_generic_locus(position: u64, samples: Vec<(WindowCoverage, u32, u32)>) -> SpillEntry {
    SpillEntry {
        contig: ContigId(0),
        position: Position(position),
        is_repeat_tract: false,
        line: format!("chr1\t{position}\t.\tA\tG\t42\tPASS\t.\tGT\t0/1").into_bytes(),
        samples: SpilledSamples::GenericLocus(
            samples
                .into_iter()
                .map(|(window, ref_reads, alt_reads)| GenericLocusSample {
                    window,
                    ref_reads,
                    alt_reads,
                })
                .collect(),
        ),
    }
}

/// A repeat tract at `position`, whose rows carry windows and — by construction — no reads.
fn a_repeat_tract(position: u64, windows: Vec<WindowCoverage>) -> SpillEntry {
    SpillEntry {
        contig: ContigId(0),
        position: Position(position),
        is_repeat_tract: true,
        line: format!("chr1\t{position}\t.\tAT\tATAT\t42\tPASS\t.\tGT\t0/1").into_bytes(),
        samples: SpilledSamples::RepeatTract(
            windows
                .into_iter()
                .map(|window| RepeatTractSample { window })
                .collect(),
        ),
    }
}

/// **A record shaped like a hidden duplication**: every sample over-covered, and the reads split
/// at a third rather than a half — the share two of three collapsed copies would give.
fn a_record_shaped_like_a_duplication(position: u64, n: usize) -> SpillEntry {
    a_generic_locus(
        position,
        (0..n)
            .map(|_| (a_window(ONE_COPY_DEPTH * 2.0), 20, 10))
            .collect(),
    )
}

/// **A record shaped like a real heterozygote**: one copy's depth, reads split near a half.
fn a_record_shaped_like_a_variant(position: u64, n: usize) -> SpillEntry {
    a_generic_locus(
        position,
        (0..n).map(|_| (a_window(ONE_COPY_DEPTH), 8, 8)).collect(),
    )
}

/// A record where `carriers` of the `n` samples look duplicated and the rest look like a plain
/// one-copy heterozygote — a dial from "no duplication signal" to "every sample says so".
fn a_record_with_carriers(position: u64, n: usize, carriers: usize) -> SpillEntry {
    a_generic_locus(
        position,
        (0..n)
            .map(|i| {
                if i < carriers {
                    (a_window(ONE_COPY_DEPTH * 2.0), 20, 10)
                } else {
                    (a_window(ONE_COPY_DEPTH), 8, 8)
                }
            })
            .collect(),
    )
}

/// A record nothing can be said about: no sample has a usable window.
fn a_record_no_sample_covered(position: u64, n: usize) -> SpillEntry {
    a_generic_locus(position, (0..n).map(|_| (no_window(), 8, 8)).collect())
}

// ---------------------------------------------------------------- the spill on disk

/// A scratch directory inside the project, per `CLAUDE.md` — never the system temp.
fn scratch(name: &str) -> PathBuf {
    let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tmp")
        .join("paralog_pass_two_tests")
        .join(name);
    fs::create_dir_all(&directory).expect("a scratch directory");
    directory
}

/// A spill holding `entries`, ready to read. **Anything an earlier run of this test left at the
/// path is removed first**, because the file is created with `create_new` and a leftover would
/// fail the run rather than the assertion.
fn a_spill_holding(test_name: &str, entries: Vec<SpillEntry>) -> SpillFile {
    let output = scratch(test_name).join("cohort.vcf");
    let _ = fs::remove_file(SpillFile::beside(&output).path());
    let mut spill = SpillFile::beside(&output);
    for entry in &entries {
        spill.append(entry).expect("the entry is appended");
    }
    spill.finish_writing().expect("the spill is flushed");
    spill
}

fn default_config() -> CalibrationConfig {
    CalibrationConfig::default()
}

/// A fit that cannot settle: one iteration at a tolerance nothing reaches.
fn a_config_whose_estimate_cannot_settle() -> CalibrationConfig {
    CalibrationConfig {
        em: EmConfig {
            max_iter: 1,
            tol: 1e-300,
            ..EmConfig::default()
        },
        ..CalibrationConfig::default()
    }
}

// ---------------------------------------------------------------- the trap

/// **A record no sample can speak for carries `NaN`, and nothing about it reaches the fit.**
///
/// The scorer's own answer here is a ratio of `0.0` — the same number it returns when the two
/// stories are exactly balanced. Folded, that record would be one more vote for "not a
/// duplication" (spec §6 trap 4).
#[test]
fn a_record_no_sample_can_speak_for_is_unscored() {
    let context = a_context_over(6);
    let spill = a_spill_holding(
        "a_record_no_sample_can_speak_for_is_unscored",
        vec![a_record_no_sample_covered(100, 6)],
    );

    let verdicts =
        score_the_parked_records_and_resolve_the_cut(&spill, &context, 0.01, &default_config())
            .expect("the spill is readable");

    assert_eq!(verdicts.ratios.len(), 1, "the record still gets a slot");
    assert!(
        verdicts.ratios[0].is_nan(),
        "an unscored record must carry NaN, not the scorer's neutral 0.0; it carried {}",
        verdicts.ratios[0]
    );
    assert_eq!(
        verdicts.records_scored, 0,
        "nothing was scored, so nothing was folded"
    );
}

/// **A record every sample has a model for, and still nothing to weigh.**
///
/// The samples here are not absent: each has a fitted coverage model, a usable window and reads,
/// so every observation the scorer is handed is `Some`. What they have not got is a usable spread
/// for one copy's depth — σ₀ is zero — so the scorer drops all six and returns its neutral `0.0`.
/// Screening on "were any observations built?" would fold that zero; screening on "did the scorer
/// weigh anything?" does not, and this fixture is the difference between the two rules.
#[test]
fn a_record_whose_samples_have_models_but_no_usable_spread_is_unscored() {
    let context = a_context_over_with(
        6,
        &CoverageFitConfig {
            single_copy_depth_sd_override: Some(0.0),
            ..CoverageFitConfig::default()
        },
    );
    assert_eq!(
        context.how_many_samples_have_a_coverage_model(),
        6,
        "the fixture is only about σ₀: every sample must still have a model"
    );
    let spill = a_spill_holding(
        "a_record_whose_samples_have_models_but_no_usable_spread_is_unscored",
        vec![a_record_shaped_like_a_duplication(100, 6)],
    );

    let verdicts =
        score_the_parked_records_and_resolve_the_cut(&spill, &context, 0.01, &default_config())
            .expect("the spill is readable");

    assert!(
        verdicts.ratios[0].is_nan(),
        "the scorer weighed no sample, so the record is unscored; it carried {}",
        verdicts.ratios[0]
    );
    assert_eq!(verdicts.records_scored, 0);
}

/// **The plan's test: the histogram holds exactly the finite ratios, no more and no fewer.**
///
/// The unscorable record sits in the *middle* of the file on purpose — at either end, a pass that
/// stopped early or started late would give the same counts.
#[test]
fn the_histogram_folds_exactly_the_finite_ratios() {
    let context = a_context_over(6);
    let spill = a_spill_holding(
        "the_histogram_folds_exactly_the_finite_ratios",
        vec![
            a_record_shaped_like_a_variant(100, 6),
            a_record_no_sample_covered(200, 6),
            a_record_shaped_like_a_duplication(300, 6),
            a_repeat_tract(400, vec![a_window(ONE_COPY_DEPTH); 6]),
        ],
    );

    let verdicts =
        score_the_parked_records_and_resolve_the_cut(&spill, &context, 0.01, &default_config())
            .expect("the spill is readable");

    let finite = verdicts.ratios.iter().filter(|r| r.is_finite()).count();
    assert_eq!(
        verdicts.ratios.len(),
        4,
        "every parked record keeps a slot, scored or not"
    );
    assert_eq!(finite, 3, "one of the four could not be scored");
    assert_eq!(
        verdicts.records_scored, finite as u64,
        "the fit rests on exactly the finite ratios"
    );
    assert!(
        verdicts.ratios[1].is_nan(),
        "the middle record is the unscored one"
    );
}

// ---------------------------------------------------------------- order, kind and range

/// **Each record's ratio is its own, and they are kept in the order the file holds them.**
///
/// The seven records differ only in how many of the six samples look duplicated — none, then
/// one, up to all six — so the evidence for a hidden duplication rises across the file and the
/// ratios must rise with it. **Measured on this fixture** they run
/// `-31.15, -3.43, 24.29, 52.01, 79.74, 115.19, 156.26`. A pass that scored the wrong entry, or
/// kept the ratios out of step with the records, breaks the ordering rather than merely moving a
/// number; the assertion is on the order, so it does not go stale if the scorer's arithmetic
/// changes.
#[test]
fn the_ratios_are_in_spill_order_and_each_belongs_to_its_own_record() {
    let context = a_context_over(6);
    let spill = a_spill_holding(
        "the_ratios_are_in_spill_order_and_each_belongs_to_its_own_record",
        (0..=6)
            .map(|carriers| a_record_with_carriers(100 + carriers as u64 * 10, 6, carriers))
            .collect(),
    );

    let verdicts =
        score_the_parked_records_and_resolve_the_cut(&spill, &context, 0.01, &default_config())
            .expect("the spill is readable");

    assert_eq!(verdicts.ratios.len(), 7);
    assert!(
        verdicts.ratios.windows(2).all(|pair| pair[1] > pair[0]),
        "each record has one more sample that looks duplicated than the one before it, so its \
         ratio must be larger; got {:?}",
        verdicts.ratios
    );
}

/// **A repeat tract is scored like any other record, on its coverage alone.**
///
/// Every fixture but this one holds at least one generic locus, so a pass that quietly dropped
/// tracts would still fold something and still look calibrated. Here there is nothing else to
/// fold: if tracts do not score, the fit rests on nothing.
#[test]
fn a_spill_of_nothing_but_repeat_tracts_is_still_scored() {
    let context = a_context_over(6);
    let spill = a_spill_holding(
        "a_spill_of_nothing_but_repeat_tracts_is_still_scored",
        vec![
            a_repeat_tract(100, vec![a_window(ONE_COPY_DEPTH); 6]),
            a_repeat_tract(200, vec![a_window(ONE_COPY_DEPTH * 2.0); 6]),
        ],
    );

    let verdicts =
        score_the_parked_records_and_resolve_the_cut(&spill, &context, 0.01, &default_config())
            .expect("the spill is readable");

    assert_eq!(
        verdicts.records_scored, 2,
        "both tracts must be scored: their ratios are {:?}",
        verdicts.ratios
    );
    assert!(verdicts.ratios.iter().all(|r| r.is_finite()));
}

/// **One sample is a run, not a special case** (spec §4). The score self-gates rather than
/// refusing, so the ratio is a number and the pass completes.
#[test]
fn a_run_of_one_sample_scores_finitely() {
    let context = a_context_over(1);
    let spill = a_spill_holding(
        "a_run_of_one_sample_scores_finitely",
        vec![a_record_shaped_like_a_duplication(100, 1)],
    );

    let verdicts =
        score_the_parked_records_and_resolve_the_cut(&spill, &context, 0.01, &default_config())
            .expect("the spill is readable");

    assert_eq!(verdicts.records_scored, 1);
    assert!(
        verdicts.ratios[0].is_finite(),
        "a one-sample run must produce a number, not NaN; got {}",
        verdicts.ratios[0]
    );
}

// ---------------------------------------------------------------- the rate and the cut

/// **A run with nothing to fit uses the documented rate and says it did not fit one.**
#[test]
fn an_empty_spill_falls_back_rather_than_fitting_a_rate_from_nothing() {
    let context = a_context_over(6);
    let spill = a_spill_holding(
        "an_empty_spill_falls_back_rather_than_fitting_a_rate_from_nothing",
        vec![],
    );

    let verdicts =
        score_the_parked_records_and_resolve_the_cut(&spill, &context, 0.01, &default_config())
            .expect("an empty spill is a spill");

    assert!(verdicts.ratios.is_empty());
    assert_eq!(verdicts.records_scored, 0);
    assert!(!verdicts.calibration.prior.converged);
    assert!(
        (verdicts.calibration.prior.prior_probability - DEFAULT_FALLBACK_PARALOG_PRIOR).abs()
            < 1e-12,
        "an unfitted rate must be the fallback, not the last iterate; got {}",
        verdicts.calibration.prior.prior_probability
    );
    assert!(verdicts.why_the_paralog_rate_is_not_fitted().is_some());
}

/// **An estimate that runs out of iterations is replaced, not used** — and the words the run
/// report prints say how many records it had.
#[test]
fn an_estimate_that_cannot_settle_is_replaced_by_the_fallback() {
    let context = a_context_over(6);
    let spill = a_spill_holding(
        "an_estimate_that_cannot_settle_is_replaced_by_the_fallback",
        vec![
            a_record_shaped_like_a_variant(100, 6),
            a_record_shaped_like_a_duplication(200, 6),
        ],
    );

    let verdicts = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        0.01,
        &a_config_whose_estimate_cannot_settle(),
    )
    .expect("the spill is readable");

    assert!(!verdicts.calibration.prior.converged);
    assert!(
        (verdicts.calibration.prior.prior_probability - DEFAULT_FALLBACK_PARALOG_PRIOR).abs()
            < 1e-12
    );
    let warning = verdicts
        .why_the_paralog_rate_is_not_fitted()
        .expect("an unfitted rate is worth a warning");
    assert!(
        warning.contains("2 scored record"),
        "the warning must say how much evidence there was: {warning}"
    );
}

/// **An estimate that settles is used, and there is nothing to warn about.**
#[test]
fn a_settled_estimate_is_used_and_warns_about_nothing() {
    let context = a_context_over(6);
    let spill = a_spill_holding(
        "a_settled_estimate_is_used_and_warns_about_nothing",
        vec![
            a_record_shaped_like_a_variant(100, 6),
            a_record_shaped_like_a_duplication(200, 6),
        ],
    );

    let verdicts =
        score_the_parked_records_and_resolve_the_cut(&spill, &context, 0.01, &default_config())
            .expect("the spill is readable");

    assert!(
        verdicts.calibration.prior.converged,
        "two records are enough for the estimate to settle"
    );
    assert!(verdicts.why_the_paralog_rate_is_not_fitted().is_none());
}

/// **The operator's target is what decides how many records go**, and a looser one takes more.
///
/// The seven records span the evidence from "nothing looks duplicated" to "every sample does",
/// so the cut has somewhere to move. **Measured on this fixture**: 4 records removed at a target
/// of zero, 5 at one in a hundred, 7 — all of them — at one in two. The assertion is that the
/// counts rise, not that they are those three numbers, so it survives a change to the scorer's
/// arithmetic; the numbers are here because a test asserting only "some are flagged" would pass
/// on a fixture where the target did nothing at all, which is what the first draft of this test
/// was written on.
#[test]
fn a_looser_target_removes_more_records() {
    let context = a_context_over(6);
    let spill = a_spill_holding(
        "a_looser_target_removes_more_records",
        (0..=6)
            .map(|carriers| a_record_with_carriers(100 + carriers as u64 * 10, 6, carriers))
            .collect(),
    );

    let removed_at = |target_fdr: f64| {
        let verdicts = score_the_parked_records_and_resolve_the_cut(
            &spill,
            &context,
            target_fdr,
            &default_config(),
        )
        .expect("the spill is readable");
        assert_eq!(
            verdicts.calibration.target_fdr, target_fdr,
            "the calibration must carry the target it was asked for"
        );
        assert_eq!(
            verdicts.ratios.len(),
            7,
            "the target does not change which records were scored"
        );
        verdicts
            .ratios
            .iter()
            .filter(|ratio| verdicts.calibration.flags(**ratio))
            .count()
    };

    let strict = removed_at(0.0);
    let ordinary = removed_at(0.01);
    let loose = removed_at(0.5);
    assert!(
        strict < ordinary && ordinary < loose,
        "a looser false-discovery target must remove more records; got {strict}, {ordinary}, \
         {loose}"
    );
    assert!(
        strict > 0,
        "the fixture must reach the cut even at the strictest target"
    );
}

// ---------------------------------------------------------------- what stops the pass

/// **A record that does not describe the run's cohort stops the pass, naming the record.**
///
/// The scorer answers a cohort-length mismatch with a *neutral score* rather than an error, so a
/// pass that let it through would fold zeros for every short record and finish looking calibrated.
#[test]
fn a_record_that_is_not_the_runs_cohort_stops_the_pass() {
    let context = a_context_over(3);
    let spill = a_spill_holding(
        "a_record_that_is_not_the_runs_cohort_stops_the_pass",
        vec![a_record_shaped_like_a_variant(700, 2)],
    );

    let failure =
        score_the_parked_records_and_resolve_the_cut(&spill, &context, 0.01, &default_config())
            .expect_err("a two-sample record in a three-sample run is a wiring error");

    match failure {
        PassTwoError::CohortSize(mismatch) => {
            assert_eq!(mismatch.record, 2);
            assert_eq!(mismatch.cohort, 3);
            assert_eq!(mismatch.position, 700, "the failure names the record");
        }
        other => panic!("expected a cohort-size failure, got {other:?}"),
    }
}

/// **A spill that lost its tail is refused, not calibrated on what survived.**
///
/// Pass one counted the records it wrote, so a file holding fewer is a shortfall the reader can
/// see. Without that, a truncated spill reads as a complete, shorter run: the rate would be
/// fitted from a prefix and pass three would write a prefix of the VCF, with nothing anywhere
/// saying so.
#[test]
fn a_spill_that_lost_its_tail_is_refused() {
    let context = a_context_over(6);
    let spill = a_spill_holding(
        "a_spill_that_lost_its_tail_is_refused",
        vec![
            a_record_shaped_like_a_variant(100, 6),
            a_record_shaped_like_a_duplication(200, 6),
        ],
    );
    let whole = fs::metadata(spill.path()).expect("the spill exists").len();
    fs::OpenOptions::new()
        .write(true)
        .open(spill.path())
        .expect("the spill opens")
        .set_len(whole / 2)
        .expect("the spill is truncated");

    let failure =
        score_the_parked_records_and_resolve_the_cut(&spill, &context, 0.01, &default_config())
            .expect_err("half a spill is not a spill");

    let PassTwoError::Spill(why) = &failure else {
        panic!("a short spill is a read failure, got {failure:?}")
    };
    let said = why.to_string();
    assert!(
        said.contains("paralog-spill"),
        "the failure must name the file it happened on; it said {said}"
    );
}

// ---------------------------------------------------------------- against production

/// **The fallback, the curve and the cut agree with production's, bit for bit.**
///
/// The three pieces underneath are copies checked against production elsewhere
/// (`src/ng/paralog/production_parity.rs`); what is ng's own is the four lines that put them
/// together and substitute the documented rate for an estimate that did not settle. Production's
/// `calibrate_from_histogram` does the same four, so it is the oracle. Both settled and unsettled
/// estimates are compared, because the fallback only shows up in one of them.
#[test]
fn the_fallback_and_the_cut_agree_with_productions_bit_for_bit() {
    use crate::paralog::{EmConfig as TheirEmConfig, ParalogLrHistogram as TheirHistogram};
    use crate::var_calling::paralog_filter::calibrate::{
        CalibrationConfig as TheirConfig, calibrate_from_histogram,
    };

    // Two runs of the estimate: one that settles, and one given a single iteration at a
    // tolerance nothing meets. Only the second reaches the fallback.
    let settling = (
        CalibrationConfig::default(),
        TheirConfig {
            em: TheirEmConfig::default(),
            fallback_prior: DEFAULT_FALLBACK_PARALOG_PRIOR,
        },
    );
    let cramped_em = TheirEmConfig {
        max_iter: 1,
        tol: 1e-300,
        ..TheirEmConfig::default()
    };
    let cramped = (
        a_config_whose_estimate_cannot_settle(),
        TheirConfig {
            em: cramped_em,
            fallback_prior: DEFAULT_FALLBACK_PARALOG_PRIOR,
        },
    );

    let mut compared = 0usize;
    for (ours_config, theirs_config) in [settling, cramped] {
        for ratios in [
            Vec::new(),
            vec![-4.0, -2.0, -1.0, 0.5, 1.0],
            (0..500)
                .map(|i| f64::from(i % 50) - 25.0 + f64::from(i % 7) * 0.1)
                .collect(),
        ] {
            let mut ours = ParalogLrHistogram::with_defaults();
            let mut theirs = TheirHistogram::with_defaults();
            for &lr in &ratios {
                ours.push(lr);
                theirs.push(lr);
            }
            for target_fdr in [0.0, 0.01, 0.05, 0.2, 0.5] {
                let ours = calibrate_from_the_ratio_histogram(&ours, target_fdr, &ours_config);
                let theirs = calibrate_from_histogram(&theirs, target_fdr, &theirs_config);
                let case = format!("{} ratios, target {target_fdr}", ratios.len());
                assert_eq!(
                    ours.prior.prior_probability.to_bits(),
                    theirs.prior.prior_probability.to_bits(),
                    "{case}: the two trees fitted different rates — ng {}, src/var_calling/ {}",
                    ours.prior.prior_probability,
                    theirs.prior.prior_probability,
                );
                assert_eq!(ours.prior.converged, theirs.prior.converged, "{case}");
                assert_eq!(
                    ours.lr_threshold.map(f64::to_bits),
                    theirs.lr_threshold.map(f64::to_bits),
                    "{case}: the two trees cut at different ratios — ng {:?}, \
                     src/var_calling/ {:?}",
                    ours.lr_threshold,
                    theirs.lr_threshold,
                );
                for lr in [-100.0, -8.0, -1.0, 0.0, 1.0, 8.0, 30.0, 100.0, f64::NAN] {
                    assert_eq!(
                        ours.flags(lr),
                        theirs.flags(lr),
                        "{case}: the two trees disagree about dropping a record at ratio {lr}"
                    );
                    compared += 1;
                }
            }
        }
    }
    assert_eq!(
        compared,
        2 * 3 * 5 * 9,
        "every combination is compared: {compared}"
    );
}
