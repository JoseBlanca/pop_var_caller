use super::*;
use crate::memory_limit::LimitSource;
use crate::neighbour_evidence::windows::{DEFAULT_MARGIN_DENSITY, MarginThinning};

/// The default window sizes at `budget` bytes, margins thinned to the default density.
fn default_sizes(budget: u64) -> WindowSizes {
    WindowSizes {
        budget_bytes: budget,
        thinning: MarginThinning::parse(DEFAULT_MARGIN_DENSITY).expect("the default parses"),
        ..WindowSizes::default()
    }
}

/// tomato25's chromosome 1, about 90 Mb: longer than any window.
const LONG_CONTIG: u64 = 90_000_000;

fn limit(bytes: u64) -> MemoryGuard {
    MemoryGuard {
        limit: Some(MemoryLimit {
            bytes,
            source: LimitSource::ControlGroupV2,
        }),
        ignore_limit: false,
    }
}

/// **tomato25 at 70 GB, 48 threads**: the budget sets the window. 189,230 walked records (70 GB
/// over 170 bytes a sample for 2,176 samples), 9.54 GB of arrays, 92.3 GB with the slack.
#[test]
fn tomato25_at_70_gb_and_48_threads_expects_about_92_gb() {
    let estimate =
        PeakEstimate::before_any_window(&default_sizes(70_000_000_000), 2_176, 48, LONG_CONTIG);
    assert_eq!(estimate.counted_bytes, 70_000_000_000);
    assert_eq!(estimate.walked_records, 189_230);
    assert_eq!(estimate.walking_threads, 48);
    assert_eq!(estimate.walk_bytes(), 9_537_192_000);
    assert_eq!(estimate.peak_bytes() / 10_000_000, 9_226);
}

/// **The measured replay**, glibc, 48 threads: 103.0 GB at an 80-GB budget and 88.7 GB at 65 GB.
/// The estimates, 105.4 and 85.7 GB, are within 4% of each, either way.
#[test]
fn the_estimates_lie_within_four_percent_of_the_measured_peaks() {
    for (budget, measured) in [(80e9, 103.0e9), (65e9, 88.7e9)] {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let estimate =
            PeakEstimate::before_any_window(&default_sizes(budget as u64), 2_176, 48, LONG_CONTIG);
        #[allow(clippy::cast_precision_loss)]
        let ratio = estimate.peak_bytes() as f64 / measured;
        assert!((0.96..1.04).contains(&ratio), "{budget}: {ratio}");
    }
}

/// **A small cohort's window is set by its length, not the budget**: two samples at the default
/// 32 GB count at most 4 Mb of records at 150 a kb, 0.9 GB, and walk at most 324,500 of them, on two
/// threads however many the run has.
#[test]
fn a_small_cohorts_estimate_is_bounded_by_the_windows_length() {
    let estimate = PeakEstimate::before_any_window(
        &default_sizes(crate::neighbour_evidence::windows::DEFAULT_WINDOW_BYTES),
        2,
        256,
        LONG_CONTIG,
    );
    assert_eq!(estimate.counted_bytes, 600_000 * 1_500);
    // The core and 50 kb either side at 150 a kb, the 1.9 Mb beyond at 5.
    assert_eq!(estimate.walked_records, 315_000 + 9_500);
    assert_eq!(estimate.walking_threads, 2);
    assert!(estimate.peak_bytes() < 2_000_000_000, "{estimate:?}");
}

/// **A short genome bounds the window too**: a 300-kb contig is the whole window.
#[test]
fn a_short_contig_bounds_the_estimate() {
    let estimate = PeakEstimate::before_any_window(&default_sizes(32_000_000_000), 20, 8, 300_000);
    assert_eq!(estimate.walked_records, 45_000);
    assert_eq!(estimate.counted_bytes, 45_000 * (1_000 + 250 * 20));
}

/// **Unthinned margins are walked in full.**
#[test]
fn unthinned_margins_are_walked_in_full() {
    let sizes = WindowSizes {
        thinning: None,
        ..default_sizes(32_000_000_000)
    };
    let estimate = PeakEstimate::before_any_window(&sizes, 2, 4, LONG_CONTIG);
    assert_eq!(estimate.walked_records, 600_000);
}

/// **Where it fits, the run goes on**; above 90% of the limit it stops, unless told to go on.
#[test]
fn the_run_stops_above_ninety_percent_unless_told_to_go_on() {
    let budget_at =
        |bytes: u64| PeakEstimate::before_any_window(&default_sizes(bytes), 2_176, 48, LONG_CONTIG);
    let estimate = budget_at(70_000_000_000);
    let found = check(
        estimate,
        limit(118_000_000_000),
        PeakMoment::BeforeCalling,
        budget_at,
    );
    let Ok(PeakCheck::Fits { share }) = found else {
        panic!("{found:?}");
    };
    assert!((0.78..0.79).contains(&share), "{share}");

    // 100 GB under the same limit: 131.8 GB expected.
    let estimate = budget_at(100_000_000_000);
    let refused = check(
        estimate,
        limit(118_000_000_000),
        PeakMoment::BeforeCalling,
        budget_at,
    )
    .expect_err("131.8 GB is over 90% of 118 GB");
    // 80.57 GB is the largest budget within 106.2 GB, rounded down.
    assert_eq!(refused.budget_that_fits, Some(80_000_000_000));
    assert_eq!(refused.threads_that_fit, None);
    let ignored = check(
        estimate,
        MemoryGuard {
            ignore_limit: true,
            ..limit(118_000_000_000)
        },
        PeakMoment::BeforeCalling,
        budget_at,
    );
    assert_eq!(ignored, Ok(PeakCheck::Ignored(refused)));
}

/// **No limit read, nothing checked**, whatever the estimate.
#[test]
fn without_a_limit_nothing_is_checked() {
    let budget_at = |bytes: u64| PeakEstimate::for_window(bytes, 1_000_000, 100, 64, bytes);
    let found = check(
        budget_at(u64::MAX / 4),
        MemoryGuard::default(),
        PeakMoment::BeforeCalling,
        budget_at,
    );
    assert_eq!(found, Ok(PeakCheck::NoLimit));
}

/// **Where the threads' arrays are most of it, fewer threads are suggested**, and the budget that
/// would fit at this many.
#[test]
fn fewer_threads_are_suggested_where_they_would_fit() {
    // A window of 10 GB counted and 300,000 walked records, 100 samples, 96 threads: 30.2 GB of
    // arrays, 46.7 GB in all, against a 40-GB limit.
    let window =
        |budget: u64| PeakEstimate::for_window(budget / 2, budget * 3 / 200_000, 100, 96, budget);
    let refused = check(
        window(20_000_000_000),
        limit(40_000_000_000),
        PeakMoment::Window {
            contig: 0,
            position: 1,
        },
        window,
    )
    .expect_err("46.7 GB is over 36 GB");
    // (36 / 1.16 − 10) GB over 1.05 kB × 300,000 is 66.8 threads.
    assert_eq!(refused.threads_that_fit, Some(66));
    let message = refused.to_string();
    for expected in [
        "about 46.7 GB in its window from contig 0 position 1",
        "117% of the 40.0 GB",
        "lower --refine-window-bytes to",
        "run on 66 threads or fewer",
        "--refine-ignore-memory-limit",
        "pvc resume-calling",
    ] {
        assert!(message.contains(expected), "{expected}: {message}");
    }
}

/// **The message before calling**, as a person reads it, at tomato25's 100 GB on a 118-GB unit.
#[test]
fn the_message_before_calling_names_the_numbers_and_what_to_change() {
    let budget_at =
        |bytes: u64| PeakEstimate::before_any_window(&default_sizes(bytes), 2_176, 48, LONG_CONTIG);
    let refused = check(
        budget_at(100_000_000_000),
        limit(118_000_000_000),
        PeakMoment::BeforeCalling,
        budget_at,
    )
    .expect_err("over");
    let message = refused.to_string();
    assert_eq!(
        message,
        "genotype refinement is expected to need about 131.8 GB, 112% of the 118.0 GB this run \
         may use (the memory.max of the control group pvc runs in), and pvc stops above 90% so \
         that the system does not kill the run part of the way through. The estimate is 100.0 GB \
         for a window as --refine-window-bytes counts it (100.0 GB allowed), 13.6 GB for 48 \
         threads' working arrays (270,328 walked records at 1.05 kB each a thread), and 16% on \
         top for the memory allocator. To go on, lower --refine-window-bytes to 80000000000 \
         (80.0 GB) or less, or give --refine-ignore-memory-limit to run anyway, if this run's \
         peak is known to be lower."
    );
}
