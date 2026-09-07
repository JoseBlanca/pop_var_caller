//! **The one that carries the step: a tract's samples all score, and a SNP's zero-read sample
//! does not.**
//!
//! Production drops a sample at zero total reads, which is right at a locus where reads report
//! their allele and wrong at a repeat tract, where the filter hands the scorer no read counts at
//! all because slippage has smeared the split. Copied without that distinction, the skip empties
//! every tract's score and the run finishes looking correct (spec §6 trap 1). The type is what
//! forecloses it — a tract's rows have no read counts to be zero — and
//! `every_sample_of_a_tract_scores_even_where_no_read_supports_an_allele` is what says so.

use super::{ParalogScoringContext, WhyNoCoverageModel};
use crate::ng::paralog::{CoverageFitConfig, CoverageModelError, ParalogModelParams};
use crate::ng::run::paralog_filter::{
    GenericLocusSample, RepeatTractSample, SpillEntry, SpilledSamples,
};
use crate::ng::types::{ContigId, Position};
use crate::ng::window_coverage::{CoverageByGcHistogram, SampleHistogram, WindowCoverage};

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

/// A context over `n` samples, every one of which fits.
fn a_context_over(n: usize) -> ParalogScoringContext {
    let histograms: Vec<SampleHistogram> = (0..n).map(|_| a_fittable_histogram()).collect();
    let inbreeding = vec![0.0; n];
    ParalogScoringContext::new(histograms, &inbreeding, &params(), &fit_config())
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

/// A repeat tract whose samples carry windows and, by construction, no read counts.
fn a_tract(windows: Vec<WindowCoverage>) -> SpillEntry {
    SpillEntry {
        contig: ContigId(0),
        position: Position(500),
        is_repeat_tract: true,
        line: b"line".to_vec(),
        samples: SpilledSamples::RepeatTract(
            windows
                .into_iter()
                .map(|window| RepeatTractSample { window })
                .collect(),
        ),
    }
}

/// A generic locus whose samples carry windows and read counts.
fn a_generic_locus(samples: Vec<(WindowCoverage, u32, u32)>) -> SpillEntry {
    SpillEntry {
        contig: ContigId(0),
        position: Position(100),
        is_repeat_tract: false,
        line: b"line".to_vec(),
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

#[test]
fn every_sample_of_a_tract_scores_even_where_no_read_supports_an_allele() {
    // **Spec §6 trap 1.** Production's rule drops a sample at zero total reads. A tract's rows
    // carry no read counts at all, so every sample with a window must score — if this ever
    // fails, every tract in the run scores on nothing and is silently never filtered.
    let context = a_context_over(3);
    let entry = a_tract(vec![
        a_window(ONE_COPY_DEPTH),
        a_window(ONE_COPY_DEPTH * 2.0),
        a_window(ONE_COPY_DEPTH),
    ]);

    let scored: Vec<_> = (0..3)
        .map(|sample| context.observation_of(&entry, sample))
        .collect();

    assert!(
        scored.iter().all(Option::is_some),
        "a sample of a repeat tract was skipped: {scored:?}"
    );
    for observation in scored.into_iter().flatten() {
        assert_eq!(observation.alt_reads, 0);
        assert_eq!(observation.total_reads, 0);
        assert!(observation.relative_copy_number.is_finite());
    }
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

    let tract = a_tract(vec![no_window()]);
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
    let context = ParalogScoringContext::new(histograms, &[0.0; 3], &params(), &fit_config());
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
    let context = ParalogScoringContext::new(histograms, &[0.0; 5], &params(), &fit_config());

    assert_eq!(context.sample_count(), 5);
    assert_eq!(context.samples_with_a_coverage_model(), 1);

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
    let context = ParalogScoringContext::new(histograms, &[0.0; 3], &params(), &fit_config());

    let sigma = context.single_copy_depth_sd();

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
        &[0.0],
        &params(),
        &fit_config(),
    );
}
