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

fn limit(bytes: u64) -> Option<MemoryLimit> {
    Some(MemoryLimit {
        bytes,
        source: LimitSource::ControlGroupV2,
    })
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

/// **Within 85% of the limit the estimate fits; above it the run is warned**, and told the largest
/// budget that would fit.
#[test]
fn above_eighty_five_percent_the_run_is_warned() {
    let budget_at =
        |bytes: u64| PeakEstimate::before_any_window(&default_sizes(bytes), 2_176, 48, LONG_CONTIG);
    let found = check(
        budget_at(70_000_000_000),
        limit(118_000_000_000),
        PeakMoment::BeforeCalling,
        budget_at,
    );
    let PeakCheck::Fits { share } = found else {
        panic!("{found:?}");
    };
    assert!((0.78..0.79).contains(&share), "{share}");

    // 80 GB under the same limit: 105.4 GB expected, 89%, above the warning at 100.3 GB.
    let PeakCheck::NearLimit(near) = check(
        budget_at(80_000_000_000),
        limit(118_000_000_000),
        PeakMoment::BeforeCalling,
        budget_at,
    ) else {
        panic!("105.4 GB is over 85% of 118 GB");
    };
    // 76.1 GB is the largest budget within 100.3 GB, rounded down; at 80 GB, (100.3 / 1.16 − 80) GB
    // over 1.05 kB × 216,262 walked records is 28.5 threads.
    assert_eq!(near.budget_that_fits, Some(76_000_000_000));
    assert_eq!(near.threads_that_fit, Some(28));
}

/// **No limit read, nothing checked**, whatever the estimate.
#[test]
fn without_a_limit_nothing_is_checked() {
    let budget_at = |bytes: u64| PeakEstimate::for_window(bytes, 1_000_000, 100, 64, bytes);
    let found = check(
        budget_at(u64::MAX / 4),
        None,
        PeakMoment::BeforeCalling,
        budget_at,
    );
    assert_eq!(found, PeakCheck::NoLimit);
}

/// **Where the threads' arrays are most of it, fewer threads are suggested**, and the budget that
/// would fit at this many; the window's warning names the command that resumes the run.
#[test]
fn fewer_threads_are_suggested_where_they_would_fit() {
    // A window of 10 GB counted and 300,000 walked records, 100 samples, 96 threads: 30.2 GB of
    // arrays, 46.7 GB in all, against a 40-GB limit.
    let window =
        |budget: u64| PeakEstimate::for_window(budget / 2, budget * 3 / 200_000, 100, 96, budget);
    let PeakCheck::NearLimit(near) = check(
        window(20_000_000_000),
        limit(40_000_000_000),
        PeakMoment::Window {
            contig: 0,
            position: 1,
        },
        window,
    ) else {
        panic!("46.7 GB is over 34 GB");
    };
    // (34 / 1.16 − 10) GB over 1.05 kB × 300,000 is 61.3 threads.
    assert_eq!(near.threads_that_fit, Some(61));
    let message = near.warning(Some(Path::new("/scratch/pvc_tmp_ab12")));
    for expected in [
        "refinement memory, window from contig 0 position 1: about 46.7 GB expected",
        "117% of the 40.0 GB",
        "Within 85% would be --refine-window-bytes",
        "or 61 threads or fewer (--threads)",
        "`pvc resume-calling --scratch-dir /scratch/pvc_tmp_ab12 --refine-window-bytes ",
        "Later windows warn again only where their estimate is higher still.",
    ] {
        assert!(message.contains(expected), "{expected}: {message}");
    }
    assert!(!message.contains("upper bound"), "{message}");
    // The whole line, as the implementation report quotes it.
    assert_eq!(
        message,
        "refinement memory, window from contig 0 position 1: about 46.7 GB expected at its peak, \
         117% of the 40.0 GB this run may use (the memory.max of the control group pvc runs in), \
         above the 85% at which pvc warns; the run goes on. The estimate is 10.0 GB for a window \
         as --refine-window-bytes counts it (20.0 GB allowed), 30.2 GB for 96 threads' working \
         arrays (300,000 walked records at 1.05 kB each a thread), and 16% on top for the memory \
         allocator. Within 85% would be --refine-window-bytes 14000000000 (14.0 GB) or less, or \
         61 threads or fewer (--threads). If the system kills the run for memory after calling \
         has ended, its scratch directory keeps the calls, and `pvc resume-calling --scratch-dir \
         /scratch/pvc_tmp_ab12 --refine-window-bytes 14000000000` finishes it. Later windows \
         warn again only where their estimate is higher still."
    );
}

/// **The warning before calling**, as a person reads it, at tomato25's 100 GB on a 118-GB unit:
/// it says the estimate is an upper bound, what would fit, and how a killed run is finished.
#[test]
fn the_warning_before_calling_names_the_numbers_and_what_to_change() {
    let budget_at =
        |bytes: u64| PeakEstimate::before_any_window(&default_sizes(bytes), 2_176, 48, LONG_CONTIG);
    let PeakCheck::NearLimit(near) = check(
        budget_at(100_000_000_000),
        limit(118_000_000_000),
        PeakMoment::BeforeCalling,
        budget_at,
    ) else {
        panic!("over");
    };
    assert_eq!(
        near.warning(Some(Path::new("/scratch/pvc_tmp_ab12"))),
        "refinement memory, before calling: about 131.8 GB expected at its peak, 112% of the \
         118.0 GB this run may use (the memory.max of the control group pvc runs in), above the \
         85% at which pvc warns; the run goes on. This is an upper bound: no window is planned \
         yet, so each is taken to be as large as --refine-window-bytes and the longest contig \
         allow; each window is estimated again once it is planned. The estimate is 100.0 GB for a \
         window as --refine-window-bytes counts it (100.0 GB allowed), 13.6 GB for 48 threads' \
         working arrays (270,328 walked records at 1.05 kB each a thread), and 16% on top for the \
         memory allocator. Within 85% would be --refine-window-bytes 76000000000 (76.0 GB) or \
         less. If the system kills the run for memory after calling has ended, its scratch \
         directory keeps the calls, and `pvc resume-calling --scratch-dir /scratch/pvc_tmp_ab12 \
         --refine-window-bytes 76000000000` finishes it."
    );
    // With the hidden-duplication filter off there is nothing to resume from.
    assert!(
        near.warning(None)
            .ends_with("so a run the system kills for memory has to be started again."),
        "{}",
        near.warning(None)
    );
}
