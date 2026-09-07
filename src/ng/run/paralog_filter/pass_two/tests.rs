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
//! **Which dimensions these fixtures vary.** Whether a record can be scored; where the
//! unscorable one sits in the file; whether the file's order and the evidence's order agree;
//! repeat tract against generic locus, including a spill of nothing else; how many samples the run
//! has, at one and at six; whether the rate estimate settles, and whether its fallback rate is the
//! same number as the iteration's starting guess; and — in the two cohort-mismatch tests — the
//! contig, which is the one otherwise-constant dimension this pass reads and reports.
//!
//! **And three they hold constant, argued rather than tested here.** GC content (one bin per
//! histogram), the inbreeding coefficient (every sample outbred), and the coverage model (every
//! sample fitted from the same histogram). All three reach the score through
//! `ParalogScoringContext`, which this pass calls once a record and never indexes into, so varying
//! them changes the number and not the path; the fixtures that vary them live in
//! `scoring_context/tests.rs`, where a one-GC-bin fixture was the defect that hid a wrong
//! argument. **Naming only one of them was this file's own review finding** — a list that names a
//! single constant reads as "we checked, there is one", and the contig hole survived under it.

use std::fs;
use std::path::PathBuf;

use super::{
    LrHistogramShape, PassTwoError, TargetFdr, score_the_parked_records_and_resolve_the_cut,
};
use crate::ng::paralog::{
    CalibrationConfig, CoverageFitConfig, DEFAULT_FALLBACK_PARALOG_PRIOR,
    DEFAULT_LR_HISTOGRAM_BINS, DEFAULT_LR_HISTOGRAM_HI, DEFAULT_LR_HISTOGRAM_LO, EmConfig,
    ParalogLrHistogram, ParalogModelParams,
};
use crate::ng::run::paralog_filter::{
    GenericLocusSample, ParalogScoringContext, RepeatTractSample, SpillEntry, SpillFile,
    SpillFileError, SpilledSamples,
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

/// A target false-discovery rate the type admits.
fn target(fdr: f64) -> TargetFdr {
    TargetFdr::try_new(fdr).expect("a fraction between zero and one is a target")
}

/// **The strictest target a run can ask for.** Zero is not a target — it means *do not run the
/// filter* — so the strictest thing that reaches the scoring is a fraction just above it. Records
/// whose tail false-discovery value underflows to exactly zero still meet it, which is the whole
/// reason zero could not be left as a target.
const THE_STRICTEST_TARGET: f64 = 1e-12;

/// **A configuration whose fallback rate is not the iteration's starting guess.**
///
/// `DEFAULT_EM_START` and `DEFAULT_FALLBACK_PARALOG_PRIOR` are both `0.03`, and an empty
/// histogram's estimate *is* the starting guess — so at the shipped defaults a test cannot tell
/// "the fallback was substituted" from "nothing happened", nor "the fallback was read" from "the
/// starting guess was read". Two mutations survived the first version of this suite for exactly
/// that reason. Every test that asserts the fallback uses these two separated numbers; the
/// shipped default stays pinned by `the_shipped_fallback_rate_is_the_documented_one` here, and by
/// the differential against production in `src/ng/paralog/production_parity.rs`.
const A_FALLBACK_THAT_IS_NOT_THE_SEED: f64 = 0.41;
const A_SEED_THAT_IS_NOT_THE_FALLBACK: f64 = 0.17;

fn a_config_whose_fallback_differs_from_its_seed() -> CalibrationConfig {
    CalibrationConfig {
        em: EmConfig {
            start: A_SEED_THAT_IS_NOT_THE_FALLBACK,
            ..EmConfig::default()
        },
        fallback_prior: A_FALLBACK_THAT_IS_NOT_THE_SEED,
    }
}

/// A fit that cannot settle: one iteration at a tolerance nothing reaches.
fn a_config_whose_estimate_cannot_settle() -> CalibrationConfig {
    CalibrationConfig {
        em: EmConfig {
            max_iter: 1,
            tol: 1e-300,
            start: A_SEED_THAT_IS_NOT_THE_FALLBACK,
        },
        fallback_prior: A_FALLBACK_THAT_IS_NOT_THE_SEED,
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

    let verdicts = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
    .expect("the spill is readable");

    assert_eq!(verdicts.ratios.len(), 1, "the record still gets a slot");
    assert!(
        verdicts.ratios[0].is_nan(),
        "an unscored record must carry NaN, not the scorer's neutral 0.0; it carried {}",
        verdicts.ratios[0]
    );
    assert_eq!(
        verdicts.records_in_the_fit, 0,
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

    let verdicts = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
    .expect("the spill is readable");

    assert!(
        verdicts.ratios[0].is_nan(),
        "the scorer weighed no sample, so the record is unscored; it carried {}",
        verdicts.ratios[0]
    );
    assert_eq!(verdicts.records_in_the_fit, 0);
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

    let verdicts = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
    .expect("the spill is readable");

    let finite = verdicts.ratios.iter().filter(|r| r.is_finite()).count();
    assert_eq!(
        verdicts.ratios.len(),
        4,
        "every parked record keeps a slot, scored or not"
    );
    assert_eq!(finite, 3, "one of the four could not be scored");
    assert_eq!(
        verdicts.records_in_the_fit, finite as u64,
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

    let verdicts = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
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

    let verdicts = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
    .expect("the spill is readable");

    assert_eq!(
        verdicts.records_in_the_fit, 2,
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

    let verdicts = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
    .expect("the spill is readable");

    assert_eq!(verdicts.records_in_the_fit, 1);
    assert!(
        verdicts.ratios[0].is_finite(),
        "a one-sample run must produce a number, not NaN; got {}",
        verdicts.ratios[0]
    );
}

// ---------------------------------------------------------------- the rate and the cut

/// **A run with nothing to fit uses the configured fallback rate, and says it did not fit one.**
///
/// **The configuration here deliberately separates two numbers the defaults make equal.** The
/// iteration's starting guess and the documented fallback are both `0.03`, and an empty
/// histogram's estimate *is* the starting guess — so at the defaults this test passes whether the
/// fallback is substituted at all, and whether the substitution reads the right field. Both
/// mistakes survived the first version of this suite.
#[test]
fn an_empty_spill_falls_back_rather_than_fitting_a_rate_from_nothing() {
    let context = a_context_over(6);
    let spill = a_spill_holding(
        "an_empty_spill_falls_back_rather_than_fitting_a_rate_from_nothing",
        vec![],
    );

    let verdicts = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &a_config_whose_fallback_differs_from_its_seed(),
    )
    .expect("an empty spill is a spill");

    assert!(verdicts.ratios.is_empty());
    assert_eq!(verdicts.records_in_the_fit, 0);
    assert!(!verdicts.calibration.prior.converged);
    assert!(
        (verdicts.calibration.prior.prior_probability - A_FALLBACK_THAT_IS_NOT_THE_SEED).abs()
            < 1e-12,
        "an unfitted rate must be the configured fallback of {A_FALLBACK_THAT_IS_NOT_THE_SEED} \
         and not the iteration's starting guess of {A_SEED_THAT_IS_NOT_THE_FALLBACK}; got {}",
        verdicts.calibration.prior.prior_probability
    );
    assert!(verdicts.why_the_paralog_rate_is_not_fitted().is_some());
}

/// **The rate a run falls back to, when nobody configures one, is the documented default.**
///
/// Its sibling above deliberately configures a fallback the starting guess does not share, so
/// that it can tell the two apart; this is what pins the number an ordinary run actually uses.
#[test]
fn the_shipped_fallback_rate_is_the_documented_one() {
    let context = a_context_over(6);
    let spill = a_spill_holding("the_shipped_fallback_rate_is_the_documented_one", vec![]);

    let verdicts = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
    .expect("an empty spill is a spill");

    assert!(
        (verdicts.calibration.prior.prior_probability - DEFAULT_FALLBACK_PARALOG_PRIOR).abs()
            < 1e-12,
        "the shipped fallback is {DEFAULT_FALLBACK_PARALOG_PRIOR}; got {}",
        verdicts.calibration.prior.prior_probability
    );
    assert_eq!(
        verdicts.config.fallback_prior, DEFAULT_FALLBACK_PARALOG_PRIOR,
        "and the verdicts record which rate the run used"
    );
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
        target(0.01),
        &a_config_whose_estimate_cannot_settle(),
    )
    .expect("the spill is readable");

    assert!(!verdicts.calibration.prior.converged);
    assert!(
        (verdicts.calibration.prior.prior_probability - A_FALLBACK_THAT_IS_NOT_THE_SEED).abs()
            < 1e-12,
        "the fallback must be the configured {A_FALLBACK_THAT_IS_NOT_THE_SEED} and not the \
         iteration's last value; got {}",
        verdicts.calibration.prior.prior_probability
    );
    let warning = verdicts
        .why_the_paralog_rate_is_not_fitted()
        .expect("an unfitted rate is worth a warning");
    assert!(
        warning.contains("2 scored record"),
        "the warning must say how much evidence there was: {warning}"
    );
}

/// **The warning counts the records the fit rested on, not the records that were parked.**
///
/// Its sibling above runs on a spill where every record scores, so the two counts are equal there
/// and reporting either one passes — which is the whole reason `records_in_the_fit` exists as a
/// separate number. Here one of the three records cannot be scored, so the two differ.
#[test]
fn the_warning_counts_the_records_in_the_fit_and_not_the_records_parked() {
    let context = a_context_over(6);
    let spill = a_spill_holding(
        "the_warning_counts_the_records_in_the_fit_and_not_the_records_parked",
        vec![
            a_record_shaped_like_a_variant(100, 6),
            a_record_no_sample_covered(200, 6),
            a_record_shaped_like_a_duplication(300, 6),
        ],
    );

    let verdicts = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &a_config_whose_estimate_cannot_settle(),
    )
    .expect("the spill is readable");

    assert_eq!(verdicts.ratios.len(), 3, "three records were parked");
    assert_eq!(
        verdicts.records_in_the_fit, 2,
        "one of them could not be scored"
    );
    let warning = verdicts
        .why_the_paralog_rate_is_not_fitted()
        .expect("an unfitted rate is worth a warning");
    assert!(
        warning.contains("2 scored record"),
        "the warning must count what the fit rested on, not what was parked: {warning}"
    );
    assert!(
        warning.contains("0.410000"),
        "and it must name the rate the run actually used: {warning}"
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

    let verdicts = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
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
            target(target_fdr),
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

    let strict = removed_at(THE_STRICTEST_TARGET);
    let ordinary = removed_at(0.01);
    let loose = removed_at(0.5);
    assert!(
        strict < ordinary && ordinary < loose,
        "a looser false-discovery target must remove more records; got {strict}, {ordinary}, \
         {loose}"
    );
    // **This holds because four of the seven tail false-discovery values underflow to *exactly*
    // zero**, not merely to something tiny — measured `[0.2747, 0.1538, 2.23e-12, 0, 0, 0, 0]` —
    // and the rule is `q <= target`. A fixture whose strongest records merely came close would
    // remove nothing here.
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
///
/// **The fixture is deliberately not on contig 0.** Every other record in this file is, so a
/// failure that reported a hardcoded first contig would be right on all of them — and on a real
/// run, where almost nothing is on the first contig, it would send whoever debugs the wiring
/// error to the wrong chromosome. That mutation survived the first version of this suite.
#[test]
fn a_record_that_is_not_the_runs_cohort_stops_the_pass() {
    let context = a_context_over(3);
    let mut narrower = a_record_shaped_like_a_variant(700, 2);
    narrower.contig = ContigId(7);
    let spill = a_spill_holding(
        "a_record_that_is_not_the_runs_cohort_stops_the_pass",
        vec![narrower],
    );

    let failure = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
    .expect_err("a two-sample record in a three-sample run is a wiring error");

    let PassTwoError::CohortSize(mismatch) = failure else {
        panic!("expected a cohort-size failure, got {failure:?}")
    };
    assert_eq!(mismatch.record, 2);
    assert_eq!(mismatch.cohort, 3);
    assert_eq!(mismatch.position, 700, "the failure names the record");
    assert_eq!(mismatch.contig, 7, "on the record's own contig");
}

/// **A record wider than the cohort stops the pass too.**
///
/// Its sibling above is narrower than the run. The two are different wrong answers — a narrow
/// record would be scored on a prefix of the cohort, a wide one truncated — and the scorer
/// answers both with a neutral score rather than an error, so both have to be refused here.
#[test]
fn a_record_wider_than_the_runs_cohort_stops_the_pass() {
    let context = a_context_over(2);
    let spill = a_spill_holding(
        "a_record_wider_than_the_runs_cohort_stops_the_pass",
        vec![a_record_shaped_like_a_variant(700, 5)],
    );

    let failure = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
    .expect_err("a five-sample record in a two-sample run is a wiring error");

    let PassTwoError::CohortSize(mismatch) = failure else {
        panic!("expected a cohort-size failure, got {failure:?}")
    };
    assert_eq!(mismatch.record, 5);
    assert_eq!(mismatch.cohort, 2);
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

    let failure = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
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

/// **An unscored record is never flagged, however loose the operator's target.**
///
/// "Kept, never flagged, and out of the fit" is the whole contract for a record no sample could
/// speak for, and the other two thirds are covered several times over. This is the third.
///
/// **It is closer to failing than it looks.** The removal rule is `q <= target`, and the
/// false-discovery curve's own answer for a value that is not a number is *not* large — measured
/// on this fixture it is `0.4999999999999852`, which sits **inside** a target of one in two. So
/// the only thing keeping an unscorable record out of the removed set is that the verdict screens
/// on finiteness before it consults the curve. Nothing in this pass would notice if that screen
/// moved: the record would be dropped from the VCF as a hidden duplication, nothing would panic,
/// and `records_in_the_fit` would still say it was never scored.
#[test]
fn an_unscored_record_is_never_flagged_however_loose_the_target() {
    let context = a_context_over(6);
    let spill = a_spill_holding(
        "an_unscored_record_is_never_flagged_however_loose_the_target",
        vec![
            a_record_shaped_like_a_variant(100, 6),
            a_record_no_sample_covered(200, 6),
            a_record_shaped_like_a_duplication(300, 6),
        ],
    );

    for target_fdr in [THE_STRICTEST_TARGET, 0.01, 0.5, 0.999] {
        let verdicts = score_the_parked_records_and_resolve_the_cut(
            &spill,
            &context,
            target(target_fdr),
            &default_config(),
        )
        .expect("the spill is readable");
        assert!(
            verdicts.ratios[1].is_nan(),
            "the middle record is the unscored one"
        );
        assert!(
            !verdicts.calibration.flags(verdicts.ratios[1]),
            "an unscored record must never be removed; at a target of {target_fdr} the curve \
             answers {} for a value that is not a number",
            verdicts.calibration.curve.q_of_lr(f64::NAN)
        );
    }

    // **And at the loosest target the two scored records are removed**, so the assertions above
    // are not passing because nothing is removed at all.
    let loose = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.999),
        &default_config(),
    )
    .expect("the spill is readable");
    assert_eq!(
        loose
            .ratios
            .iter()
            .filter(|ratio| loose.calibration.flags(**ratio))
            .count(),
        2,
        "both scored records are removed at the loosest target, so only the unscored one is \
         held back; the ratios were {:?}",
        loose.ratios
    );
}

/// **The ratios follow the file's order, not the evidence's.**
///
/// Its sibling `the_ratios_are_in_spill_order_and_each_belongs_to_its_own_record` parks its
/// records with the evidence rising, so spill order and ascending order are the same order there
/// and an implementation that returned the ratios *sorted* would satisfy it. Here the file's order
/// is deliberately not sorted — the strongest record is first and the weakest is in the middle —
/// so each ratio must rank exactly as its own record's evidence does. Pass three pairs the *i*th
/// spilled record with the *i*th ratio, and a misalignment there gives every record its
/// neighbour's verdict with nothing about the file looking wrong.
#[test]
fn the_ratios_follow_the_spill_and_not_the_evidence() {
    let context = a_context_over(6);
    let carriers_in_file_order = [6usize, 2, 0, 4, 1];
    let spill = a_spill_holding(
        "the_ratios_follow_the_spill_and_not_the_evidence",
        carriers_in_file_order
            .iter()
            .enumerate()
            .map(|(i, &carriers)| a_record_with_carriers(100 + i as u64 * 10, 6, carriers))
            .collect(),
    );

    let verdicts = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
    .expect("the spill is readable");

    let mut by_ratio: Vec<usize> = (0..carriers_in_file_order.len()).collect();
    by_ratio.sort_by(|&a, &b| {
        verdicts.ratios[a]
            .partial_cmp(&verdicts.ratios[b])
            .expect("every ratio here is a number")
    });
    let mut by_carriers: Vec<usize> = (0..carriers_in_file_order.len()).collect();
    by_carriers.sort_by_key(|&i| carriers_in_file_order[i]);
    assert_eq!(
        by_ratio, by_carriers,
        "the i-th ratio must belong to the i-th record of the file; got {:?} for carrier counts \
         {carriers_in_file_order:?}",
        verdicts.ratios
    );
}

/// **A spill pass one has not finished writing is refused, not read.**
///
/// Whatever has reached the disk is a prefix, and scoring it would calibrate the run on however
/// much happened to be flushed. Calling pass two before pass one has ended is a plausible C4
/// wiring mistake, and this is the only test that takes that branch.
#[test]
fn a_spill_pass_one_has_not_finished_writing_is_refused() {
    let context = a_context_over(6);
    let output = scratch("a_spill_pass_one_has_not_finished_writing_is_refused").join("cohort.vcf");
    let _ = fs::remove_file(SpillFile::beside(&output).path());
    let mut spill = SpillFile::beside(&output);
    spill
        .append(&a_record_shaped_like_a_variant(100, 6))
        .expect("the entry is appended");
    // deliberately no `finish_writing`

    let failure = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
    .expect_err("a spill still being written is not a spill");

    let PassTwoError::Spill(why) = &failure else {
        panic!("an unfinished spill is a read failure, got {failure:?}")
    };
    assert!(
        matches!(why, SpillFileError::PassOneHasNotEnded { .. }),
        "the refusal must say pass one has not ended; it said {why}"
    );
}

/// **A spill holding more records than pass one counted is refused too**, not calibrated on the
/// surplus.
///
/// Its sibling truncates; this is the other direction, and it is a different path in the reader —
/// the count is compared only at end of file, so the extra records are decoded, scored and folded
/// *before* the mismatch is seen. A long file is what two writers on one path, or a run resumed
/// over a leftover, would leave.
#[test]
fn a_spill_holding_more_records_than_pass_one_counted_is_refused() {
    let context = a_context_over(6);
    let spill = a_spill_holding(
        "a_spill_holding_more_records_than_pass_one_counted_is_refused",
        vec![a_record_shaped_like_a_variant(100, 6)],
    );
    let bytes = fs::read(spill.path()).expect("the spill reads");
    let mut doubled = bytes.clone();
    doubled.extend_from_slice(&bytes);
    fs::write(spill.path(), &doubled).expect("the spill is rewritten");

    let failure = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
    .expect_err("a spill longer than its own count is not a spill");

    let PassTwoError::Spill(why) = &failure else {
        panic!("a long spill is a read failure, got {failure:?}")
    };
    assert!(
        why.to_string().contains("paralog-spill"),
        "the failure must name the file it happened on; it said {why}"
    );
}

/// **A target that is not a fraction of one is refused, and both wrong ends are refused.**
///
/// Handed a bare `f64` the calibration would take either without a word: a target of `5` — what
/// someone typing five for "five percent" gets — removes every scored record, and a negative or
/// not-a-number one removes none while reporting no cut, which is also what an unreachable target
/// reports. Neither is distinguishable afterwards, so the type refuses them up front.
#[test]
fn a_target_that_is_not_a_fraction_of_one_is_refused() {
    // **Zero is among the refused, and it is the one worth saying.** It means *do not run the
    // filter*, which is a different thing from a target — and handed to the scoring it would not
    // remove nothing: a strongly duplicated record's tail false-discovery value underflows to
    // exactly zero, and zero is not above zero. `WhatTheOperatorAskedFor::from_the_flags` is
    // where zero becomes "no filter", once, so that this type never has to hold it.
    for refused in [5.0, 1.0, 0.0, -0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            TargetFdr::try_new(refused).is_err(),
            "{refused} is not a target false-discovery rate"
        );
    }
    for admitted in [1e-12, 0.001, 0.01, 0.5, 0.999_999] {
        assert_eq!(
            TargetFdr::try_new(admitted)
                .expect("a fraction below one is a target")
                .get(),
            admitted
        );
    }
}

/// **The shape the verdicts record is the histogram the ratios were actually folded into.**
///
/// The run report reads the range off the verdicts, and `ratios_outside_the_histogram` is counted
/// against it — so a run that folded into one range and reported another would mis-count how many
/// ratios saturated and mis-describe the cut. **The first version of this test could not see
/// that**: it compared the recorded shape against the constants and never against the histogram,
/// so doubling the range where the histogram was built survived it. The pass now builds the
/// histogram *from* the shape, and this asserts the shape's own histogram is the shipped one.
#[test]
fn the_recorded_shape_is_the_histogram_the_ratios_were_folded_into() {
    let context = a_context_over(6);
    let spill = a_spill_holding(
        "the_recorded_shape_is_the_histogram_the_ratios_were_folded_into",
        vec![a_record_shaped_like_a_duplication(100, 6)],
    );

    let verdicts = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
    .expect("the spill is readable");

    assert_eq!(
        verdicts.lr_histogram,
        LrHistogramShape {
            lowest_ratio: DEFAULT_LR_HISTOGRAM_LO,
            highest_ratio: DEFAULT_LR_HISTOGRAM_HI,
            bins: DEFAULT_LR_HISTOGRAM_BINS,
        },
        "the run reports the shipped range"
    );
    // **The same constructor the pass uses**, so a range that drifted inside it is caught here
    // rather than being reported as whatever the shape happens to say.
    assert_eq!(
        verdicts.lr_histogram.histogram(),
        ParalogLrHistogram::with_defaults(),
        "the histogram the ratios were folded into is the shipped one"
    );
}

/// **A ratio past the end of the histogram is counted, because at a large cohort they all will
/// be.**
///
/// The likelihood ratio is a sum over the samples that were weighed, so it grows with the cohort
/// while the histogram's range stays at ±100: measured on one duplication-shaped record, 24.2 at
/// one sample, 156.3 at six. Saturating costs nothing while it happens to one class — its
/// probability is saturated long before — but a cohort large enough to push *both* classes past
/// the same edge leaves them sharing one bin, with nothing left for the target to move. This
/// count is what the plan's D2 and D3 runs report on real data.
#[test]
fn a_ratio_past_the_end_of_the_histogram_is_counted() {
    let context = a_context_over(6);
    let spill = a_spill_holding(
        "a_ratio_past_the_end_of_the_histogram_is_counted",
        vec![
            a_record_shaped_like_a_variant(100, 6),
            a_record_no_sample_covered(200, 6),
            a_record_shaped_like_a_duplication(300, 6),
        ],
    );

    let verdicts = score_the_parked_records_and_resolve_the_cut(
        &spill,
        &context,
        target(0.01),
        &default_config(),
    )
    .expect("the spill is readable");

    assert!(
        verdicts.ratios[2] > verdicts.lr_histogram.highest_ratio,
        "the duplication-shaped record's ratio must be past the top of the range for this test \
         to mean anything; it was {}",
        verdicts.ratios[2]
    );
    assert_eq!(
        verdicts.ratios_outside_the_histogram, 1,
        "one of the three is outside, one is inside, and the unscored one is neither"
    );
    assert!(
        verdicts.ratios[0] > verdicts.lr_histogram.lowest_ratio
            && verdicts.ratios[0] < verdicts.lr_histogram.highest_ratio,
        "the variant-shaped record's ratio must be inside the range; it was {}",
        verdicts.ratios[0]
    );
}
