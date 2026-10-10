//! **What a spilled row hands the scorer.** Every record's rows carry coverage and the allele
//! split, repeat tracts included (spec §3.2), and a sample whose rows show no read at all is absent
//! from the score, as in production — at a tract as at any other record.

use super::{ParalogScoringContext, WhyNoCoverageModel};
use crate::paralog::{CoverageFitConfig, CoverageModelError, ParalogModelParams};
use crate::run::paralog_filter::{SpillEntry, SpilledSample};
use crate::types::{ContigId, InbreedingF, Position};
use crate::window_coverage::{CoverageByGcHistogram, SampleHistogram, WindowCoverage};

/// A histogram a sample's depth model can actually be fitted from: one GC bin, a single-copy
/// peak well inside the regular depth bins.
fn a_fittable_histogram() -> SampleHistogram {
    let gc_bins = 1u32;
    let depth_bins = 40u32;
    let row = depth_bins as usize + 1;
    let mut counts = vec![0u32; gc_bins as usize * row];
    // A peak at bin 10 with shoulders, so the mode is resolvable and the median sits near it.
    for (bin, count) in [(8usize, 200u32), (9, 900), (10, 2000), (11, 900), (12, 200)] {
        counts[bin] = count;
    }
    SampleHistogram::Fitted(CoverageByGcHistogram {
        window_bp: 500,
        gc_bins,
        depth_bin_width: 0.5,
        depth_bins,
        windows_folded: 4200,
        windows_under_the_floor: 0,
        counts,
    })
}

/// The depth this histogram's model calls one copy, near enough — the peak bin's centre.
const ONE_COPY_DEPTH: f32 = 5.25;

fn params() -> ParalogModelParams {
    ParalogModelParams::default()
}

fn fit_config() -> CoverageFitConfig {
    CoverageFitConfig::default()
}

/// `n` outbred coefficients, which is what every fixture here wants.
fn outbred(n: usize) -> Vec<InbreedingF> {
    vec![InbreedingF::try_new(0.0).expect("zero is a coefficient"); n]
}

/// Build a context, expecting the fit configuration to be a usable one.
fn a_context_from(
    histograms: Vec<SampleHistogram>,
    fit: &CoverageFitConfig,
) -> ParalogScoringContext {
    let n = histograms.len();
    ParalogScoringContext::new(histograms, &outbred(n), &params(), fit)
        .expect("the default fit configuration is a usable one")
}

/// A context over `n` samples, every one of which fits.
fn a_context_over(n: usize) -> ParalogScoringContext {
    a_context_from(
        (0..n).map(|_| a_fittable_histogram()).collect(),
        &fit_config(),
    )
}

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

/// A repeat tract whose samples carry windows and read counts — the generic locus's rows, with the
/// tract flag.
fn a_tract(samples: Vec<(WindowCoverage, u32, u32)>) -> SpillEntry {
    SpillEntry {
        position: Position(500),
        is_repeat_tract: true,
        ..a_generic_locus(samples)
    }
}

/// A generic locus whose samples carry windows and read counts.
fn a_generic_locus(samples: Vec<(WindowCoverage, u32, u32)>) -> SpillEntry {
    SpillEntry {
        contig: ContigId(0),
        position: Position(100),
        is_repeat_tract: false,
        line: b"line".to_vec(),
        samples: samples
            .into_iter()
            .map(|(window, ref_reads, alt_reads)| SpilledSample {
                window,
                ref_reads,
                alt_reads,
            })
            .collect(),
    }
}

#[test]
fn a_tract_hands_the_scorer_its_read_split_and_its_zero_read_sample_is_absent() {
    // **Spec §3.2, owner 2026-10-08/09.** A tract is scored as every other record is, so its
    // split reaches the scorer and a sample with no reads says nothing, as at a SNP.
    let context = a_context_over(2);
    let entry = a_tract(vec![
        (a_window(ONE_COPY_DEPTH), 2, 7),
        (a_window(ONE_COPY_DEPTH), 0, 0),
    ]);

    let observation = context
        .observation_of(&entry, 0)
        .expect("a covered sample with reads");
    assert_eq!((observation.alt_reads, observation.total_reads), (7, 9));
    assert!(
        context.observation_of(&entry, 1).is_none(),
        "a tract's sample with no reads entered the score"
    );
}

#[test]
fn a_generic_locus_sample_with_no_reads_is_absent_from_the_score() {
    // The other half of the same rule, and production's: at a locus whose reads report their
    // allele, a sample with none says nothing about allele balance.
    let context = a_context_over(2);
    let entry = a_generic_locus(vec![
        (a_window(ONE_COPY_DEPTH), 5, 5),
        (a_window(ONE_COPY_DEPTH), 0, 0),
    ]);

    assert!(context.observation_of(&entry, 0).is_some());
    assert!(
        context.observation_of(&entry, 1).is_none(),
        "a sample with no reads entered a generic locus's score"
    );
}

#[test]
fn a_generic_locus_hands_the_scorer_its_alternatives_summed() {
    // Pooling is what lets a multiallelic site take this path at all: the collapsed-duplication
    // story asks whether the non-reference share sits at some whole number of copies over the
    // total, and it does not care which alternative each read carried.
    let context = a_context_over(1);
    let entry = a_generic_locus(vec![(a_window(ONE_COPY_DEPTH), 4, 8)]);

    let observation = context
        .observation_of(&entry, 0)
        .expect("a covered sample with reads");

    assert_eq!(observation.alt_reads, 8);
    assert_eq!(observation.total_reads, 12);
}

#[test]
fn a_sample_with_no_usable_window_is_absent_whatever_its_reads_say() {
    // The `NaN` pair is *this sample has no window here*. It has to lose to the reads rather
    // than the other way round: a sample with plenty of reads and no window has no coverage
    // evidence, and coverage is the signal the filter rests on.
    let context = a_context_over(1);
    let entry = a_generic_locus(vec![(no_window(), 20, 20)]);

    assert!(context.observation_of(&entry, 0).is_none());

    let tract = a_tract(vec![(no_window(), 20, 20)]);
    assert!(context.observation_of(&tract, 0).is_none());
}

#[test]
fn a_sample_whose_coverage_model_was_refused_is_absent_from_every_record() {
    // Not scored as neutral — absent. A sample the fit refused has no yardstick, so its depth
    // is a number with no meaning rather than an average one.
    let histograms = vec![
        a_fittable_histogram(),
        SampleHistogram::NoWindowFinalised,
        a_fittable_histogram(),
    ];
    let context = a_context_from(histograms, &fit_config());
    let entry = a_generic_locus(vec![
        (a_window(ONE_COPY_DEPTH), 5, 5),
        (a_window(ONE_COPY_DEPTH), 5, 5),
        (a_window(ONE_COPY_DEPTH), 5, 5),
    ]);

    assert!(context.observation_of(&entry, 0).is_some());
    assert!(
        context.observation_of(&entry, 1).is_none(),
        "a sample with no coverage model entered a score"
    );
    assert!(context.observation_of(&entry, 2).is_some());
}

#[test]
fn a_sample_beyond_the_cohort_is_absent_rather_than_a_panic() {
    // The record and the context both claim the run's sample count, so this cannot happen from
    // inside the run — but `observation_of` takes an index and indexing is where a wiring error
    // turns into a crash on someone's whole-genome run.
    let context = a_context_over(1);
    let entry = a_generic_locus(vec![(a_window(ONE_COPY_DEPTH), 5, 5)]);

    assert!(context.observation_of(&entry, 1).is_none());
    assert!(context.observation_of(&entry, usize::MAX).is_none());
}

#[test]
fn each_way_a_sample_can_lose_its_model_is_kept_apart_from_the_others() {
    // The run report groups these, and "no window was ever finalised" is a different thing to
    // tell an operator than "the fit refused what it had".
    let histograms = vec![
        a_fittable_histogram(),
        SampleHistogram::NoWindowFinalised,
        SampleHistogram::EveryWindowUnderTheFloor,
        SampleHistogram::MedianDepthNotPositive,
        // An empty histogram: it exists, and the fit refuses it.
        SampleHistogram::Fitted(CoverageByGcHistogram {
            window_bp: 500,
            gc_bins: 1,
            depth_bin_width: 0.5,
            depth_bins: 40,
            windows_folded: 0,
            windows_under_the_floor: 0,
            counts: vec![0; 41],
        }),
    ];
    let context = a_context_from(histograms, &fit_config());

    assert_eq!(context.sample_count(), 5);
    assert_eq!(context.how_many_samples_have_a_coverage_model(), 1);

    let why = context.why_no_model();
    assert!(why[0].is_none(), "the fitted sample has no reason");
    assert_eq!(why[1], Some(WhyNoCoverageModel::NoWindowFinalised));
    assert_eq!(why[2], Some(WhyNoCoverageModel::EveryWindowUnderTheFloor));
    assert_eq!(why[3], Some(WhyNoCoverageModel::MedianDepthNotPositive));
    assert!(
        matches!(
            why[4],
            Some(WhyNoCoverageModel::TheFitRefusedTheHistogram(
                CoverageModelError::NoTiles
            ))
        ),
        "the fit's own reason was not kept: {:?}",
        why[4]
    );
}

#[test]
fn the_sigma_slice_is_the_cohorts_length_with_nan_where_a_sample_has_no_model() {
    // The scorer returns a neutral score if this slice's length disagrees with the cohort's, so
    // a short slice would not fail — it would quietly score every record as saying nothing.
    let histograms = vec![
        a_fittable_histogram(),
        SampleHistogram::NoWindowFinalised,
        a_fittable_histogram(),
    ];
    let context = a_context_from(histograms, &fit_config());

    let sigma = &context.single_copy_depth_sd;

    assert_eq!(sigma.len(), 3, "the slice must be the cohort's length");
    assert!(sigma[0].is_finite());
    assert!(
        sigma[1].is_nan(),
        "a sample with no model must be NaN here, not zero: a zero σ₀ is a claim of perfect \
         precision rather than an absence"
    );
    assert!(sigma[2].is_finite());
}

#[test]
fn the_relative_copy_number_doubles_when_the_window_depth_does() {
    // The coverage signal itself, end to end through the fit: a sample covered at twice its
    // one-copy level reads as about two copies, which is the footprint the filter hunts.
    let context = a_context_over(1);

    let one = a_generic_locus(vec![(a_window(ONE_COPY_DEPTH), 5, 5)]);
    let two = a_generic_locus(vec![(a_window(ONE_COPY_DEPTH * 2.0), 5, 5)]);

    let at_one = context
        .observation_of(&one, 0)
        .expect("a covered sample")
        .relative_copy_number;
    let at_two = context
        .observation_of(&two, 0)
        .expect("a covered sample")
        .relative_copy_number;

    assert!(
        (at_two / at_one - 2.0).abs() < 1e-9,
        "doubling the window depth gave {at_one} then {at_two}, a ratio of {}",
        at_two / at_one
    );
    // And the fixture's own claim, checked rather than asserted in a comment: the depth it
    // calls one copy really does read as about one. Without this the test above would pass on
    // a model whose scale was anywhere at all, since only the ratio is compared.
    assert!(
        (at_one - 1.0).abs() < 0.1,
        "ONE_COPY_DEPTH reads as {at_one} copies, so the fixture no longer means what its name          says and the doubled sample is not the two-copy footprint the filter hunts"
    );
}

#[test]
#[should_panic(expected = "one histogram and one inbreeding coefficient per sample")]
fn a_cohort_whose_two_slices_disagree_is_refused_rather_than_truncated() {
    let _ = ParalogScoringContext::new(
        vec![a_fittable_histogram(), a_fittable_histogram()],
        &outbred(1),
        &params(),
        &fit_config(),
    );
}

// ---------------------------------------------------------------------------
// The GC half of the copy number, and the two halves of the finiteness guard
// ---------------------------------------------------------------------------

/// Four GC bins whose one-copy depth rises with GC: bin 0 reads at 9.5, bins 1–2 at 10.5, bin 3
/// at 11.5.
///
/// **The single-bin fixture above cannot test GC at all** — the copied model returns
/// `gc_bias_curve[0]` for every GC value on a one-bin curve, so the GC argument is unobservable
/// there and a wrong one passes every assertion. This one is the fixture for the tests that are
/// about GC; the other stays right for the tests that are not.
fn a_gc_varying_histogram() -> SampleHistogram {
    let gc_bins = 4usize;
    let depth_bins = 30usize;
    let row = depth_bins + 1;
    let mut counts = vec![0u32; gc_bins * row];
    for (gc_bin, depth_bin) in [(0usize, 9usize), (1, 10), (2, 10), (3, 11)] {
        counts[gc_bin * row + depth_bin] = 150;
    }
    SampleHistogram::Fitted(CoverageByGcHistogram {
        window_bp: 500,
        gc_bins: gc_bins as u32,
        depth_bin_width: 1.0,
        depth_bins: depth_bins as u32,
        windows_folded: 600,
        windows_under_the_floor: 0,
        counts,
    })
}

/// A fit whose smoother is the identity, so the fitted curve is the four bins' own medians and
/// the expected values below can be written down rather than solved for.
fn unsmoothed_fit_config() -> CoverageFitConfig {
    CoverageFitConfig {
        smooth_window: 1,
        ..CoverageFitConfig::default()
    }
}

fn a_window_at(gc_fraction: f32, mean_depth: f32) -> WindowCoverage {
    WindowCoverage {
        gc_fraction,
        mean_depth,
    }
}

#[test]
fn the_relative_copy_number_uses_the_windows_own_gc_and_not_a_constant() {
    // Two samples at the *same* depth and different GC must read as different copy numbers,
    // because the depth a single copy produces depends on GC. Without this, passing a constant —
    // or a stale window, or the wrong sample's — is invisible: every other test in this file uses
    // a one-GC-bin histogram, where the curve is flat and the argument cannot be observed.
    let context = a_context_from(
        vec![a_gc_varying_histogram(), a_gc_varying_histogram()],
        &unsmoothed_fit_config(),
    );
    let entry = a_generic_locus(vec![
        (a_window_at(0.125, 9.5), 5, 5),
        (a_window_at(0.875, 9.5), 5, 5),
    ]);

    let low = context
        .observation_of(&entry, 0)
        .expect("a covered sample")
        .relative_copy_number;
    let high = context
        .observation_of(&entry, 1)
        .expect("a covered sample")
        .relative_copy_number;

    assert!(
        (low - 1.0).abs() < 1e-9,
        "the low-GC window sits exactly at its bin's one-copy depth, so it must read as one \
         copy: {low}"
    );
    assert!(
        (high - 9.5 / 11.5).abs() < 1e-9,
        "the high-GC window is below its bin's one-copy depth, so it must read below one: {high}"
    );
}

#[test]
fn a_window_with_only_one_field_absent_is_absent() {
    // Both halves of the guard, separately. Every other fixture sets both fields to `NaN`, so
    // either half alone catches them and one half could be deleted without a red test.
    //
    // The GC half is load-bearing against a *panic*, not a wrong answer: a `NaN` GC fraction
    // reaches `gc_multiplier`, where both range comparisons are false, the floor saturates to
    // zero, and the interpolation indexes one past the end of a one-bin curve — inside a file
    // copied verbatim from production, which ng does not edit.
    let context = a_context_over(2);
    let entry = a_generic_locus(vec![
        (a_window_at(f32::NAN, ONE_COPY_DEPTH), 5, 5),
        (a_window_at(0.4, f32::NAN), 5, 5),
    ]);

    assert!(
        context.observation_of(&entry, 0).is_none(),
        "a window with no GC fraction but a depth was scored"
    );
    assert!(
        context.observation_of(&entry, 1).is_none(),
        "a window with a GC fraction but no depth was scored"
    );
}

#[test]
fn an_infinite_window_depth_is_absent_rather_than_a_winsorised_paralog() {
    // This is where the guard is deliberately stricter than `WindowCoverage::is_absent`, which
    // asks only whether either field is `NaN`. An infinite depth divides to `+∞`, which the
    // scorer winsorises to its maximum relative copy number — so accepting it would turn a
    // measurement that never happened into the filter's strongest possible paralog signal.
    let context = a_context_over(2);
    let entry = a_generic_locus(vec![
        (a_window_at(0.4, f32::INFINITY), 5, 5),
        (a_window_at(f32::NEG_INFINITY, ONE_COPY_DEPTH), 5, 5),
    ]);

    assert!(context.observation_of(&entry, 0).is_none());
    assert!(context.observation_of(&entry, 1).is_none());
}

#[test]
fn a_sample_at_zero_window_depth_reads_as_zero_copies_rather_than_absent() {
    // Zero depth is a measurement — the sample was looked at and nothing was there — and it is
    // the opposite of the paralog signal. Treating it as an absence would drop exactly the
    // samples that argue hardest against a duplication.
    let context = a_context_over(1);
    let entry = a_generic_locus(vec![(a_window_at(0.4, 0.0), 5, 5)]);

    let observation = context
        .observation_of(&entry, 0)
        .expect("zero depth is a measurement, not an absence");

    assert_eq!(observation.relative_copy_number, 0.0);
}

#[test]
fn a_sample_far_above_its_one_copy_depth_is_handed_over_unwinsorised() {
    // The winsor cap belongs to the scorer, which applies it inside its own arithmetic. Capping
    // here as well would be a second cap that nothing states, and it would silently change the
    // input production's copy was validated on.
    let context = a_context_over(1);
    let entry = a_generic_locus(vec![(a_window_at(0.4, ONE_COPY_DEPTH * 20.0), 5, 5)]);

    let observation = context.observation_of(&entry, 0).expect("a covered sample");

    assert!(
        (observation.relative_copy_number - 20.0).abs() < 1e-9,
        "handed over as {} copies, so something capped it on the way",
        observation.relative_copy_number
    );
}

#[test]
fn a_record_mixing_covered_and_absent_samples_hands_over_only_the_covered_ones() {
    // The ordinary cohort shape, which no other test here has: some samples covered, some not.
    let context = a_context_over(4);
    let entry = a_generic_locus(vec![
        (a_window(ONE_COPY_DEPTH), 5, 5),
        (no_window(), 5, 5),
        (a_window(ONE_COPY_DEPTH * 2.0), 6, 6),
        (no_window(), 0, 0),
    ]);

    let scored: Vec<bool> = (0..4)
        .map(|sample| context.observation_of(&entry, sample).is_some())
        .collect();

    assert_eq!(scored, vec![true, false, true, false]);
}

#[test]
fn a_record_whose_width_disagrees_with_the_cohort_is_refused_rather_than_scored_on_a_prefix() {
    // Per sample there is nothing to check against — an index past the record's rows looks
    // exactly like an uncovered sample — so a narrow record would score on a prefix and a wide
    // one would be truncated, and the scorer answers a length mismatch with a *neutral score*
    // rather than an error. Both are wiring errors between the sink and this context.
    let context = a_context_over(3);
    let mut out = Vec::new();

    let two_rows = a_generic_locus(vec![
        (a_window(ONE_COPY_DEPTH), 5, 5),
        (a_window(ONE_COPY_DEPTH), 5, 5),
    ]);
    let refused = context
        .observations_of(&two_rows, &mut out)
        .expect_err("a record narrower than the cohort must be refused");
    assert_eq!(refused.record, 2);
    assert_eq!(refused.cohort, 3);

    let four_rows = a_generic_locus(vec![
        (a_window(ONE_COPY_DEPTH), 5, 5),
        (a_window(ONE_COPY_DEPTH), 5, 5),
        (a_window(ONE_COPY_DEPTH), 5, 5),
        (a_window(ONE_COPY_DEPTH), 5, 5),
    ]);
    assert!(context.observations_of(&four_rows, &mut out).is_err());
}

#[test]
fn a_records_observations_are_the_cohorts_length_and_refill_one_buffer() {
    let context = a_context_over(3);
    let entry = a_generic_locus(vec![
        (a_window(ONE_COPY_DEPTH), 5, 5),
        (no_window(), 5, 5),
        (a_window(ONE_COPY_DEPTH), 0, 0),
    ]);
    let mut out = vec![None; 99];

    context
        .observations_of(&entry, &mut out)
        .expect("the record is the cohort's width");

    assert_eq!(out.len(), 3, "the buffer is refilled, not appended to");
    assert!(out[0].is_some());
    assert!(out[1].is_none(), "no window");
    assert!(out[2].is_none(), "no reads at a generic locus");
}

#[test]
fn a_fit_configuration_that_is_not_one_fails_construction_rather_than_every_sample() {
    // One wrong knob is an operator's mistake, and the copied fit re-validates on every call —
    // so folding it into the per-sample outcome would report that the coverage model rested on
    // 0 of 3 samples, with three identical reasons, which reads as a statement about the cohort.
    let refused = ParalogScoringContext::new(
        vec![a_fittable_histogram(); 3],
        &outbred(3),
        &params(),
        &CoverageFitConfig {
            single_copy_lo: f64::NAN,
            ..CoverageFitConfig::default()
        },
    )
    .expect_err("a NaN band edge is not a configuration");

    assert!(
        !refused.reason.is_empty(),
        "the refusal must say what was wrong with it"
    );
}

#[test]
fn hom_alt_samples_at_a_tract_argue_against_a_duplication() {
    // **The case that decided spec §3.2's rule** (owner, 2026-10-08): a SNP written inside an STR
    // record, in samples whose coverage is raised. Under a collapsed duplication a sample showing
    // no reference read is near-impossible, while a third of the reads non-reference is what one
    // carrier copy in three gives — so the hom-alt tract must score lower. The flag changes
    // nothing: the tract scores as a generic locus with the same rows.
    let context = a_context_over(6);
    let raised = a_window(ONE_COPY_DEPTH * 1.6);
    let mut scratch = Vec::new();
    let mut score = |entry: &SpillEntry| {
        context
            .score(entry, &mut scratch)
            .expect("the cohort's width")
    };

    let hom_alt = score(&a_tract(vec![(raised, 0, 9); 6]));
    let one_copy_in_three = score(&a_tract(vec![(raised, 6, 3); 6]));
    assert!(
        hom_alt < one_copy_in_three,
        "hom-alt samples did not argue against a duplication: {hom_alt} against \
         {one_copy_in_three}"
    );
    assert_eq!(
        hom_alt.to_bits(),
        score(&a_generic_locus(vec![(raised, 0, 9); 6])).to_bits(),
        "the tract flag changed the score"
    );
}
