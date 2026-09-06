# window coverage — A3: a depth axis scaled to the sample it measures

**Date:** 2026-09-06
**Plan step:** [window_coverage.md](../../ng/impl_plan/window_coverage.md) Milestone A, step A3
**Spec:** [window_coverage.md](../../ng/spec/window_coverage.md) §3.4, §3.6
**Review:** [ng_window_coverage_a3_2026-09-06.md](../reviews/ng_window_coverage_a3_2026-09-06.md)
**Branch:** `ng-window-coverage`
**Builds on:** [A1](ng_window_coverage_a1_2026-09-06.md), [A2](ng_window_coverage_a2_2026-09-06.md)

## The answer

**Each sample's histogram now has its own depth axis, cut from its own depth.** Production bins
depth at a fixed 0.5× in 2,000 bins up to 1,000×. That costs 400 kB a sample against a 500 kB
per-open-sample budget, and it serves one end of the depth axis at a time: at three reads a
position a 0.5× bin is a sixth of the single-copy peak's own scatter, and at 300 reads the range
has to reach 1,200× before a four-copy carrier is inside it at all.

So the accumulator holds its first windows back, takes the median of their mean depths, and sets
the bin width so the 400 regular bins span ten times that median. Then it folds those windows and
every later one. **50 GC bins × 401 columns × 4 bytes is 80.2 kB a sample**, and both ends of the
depth range are inside it whatever the sample's coverage.

**A sample the width cannot be fitted from has no histogram at all** — one that finalised no
window, one whose every window was under A2's floor, one whose median depth came out zero. It is
carried absent, which is what the spec asks for and what the filter's scorer already does with an
absent sample.

## The defect the review found, and it is the one this step was warned about

**A fit that failed part-way through a sample silently retried and dropped the windows it failed
on.** Both reviewers found it independently. The fit took the held-back windows out of the
accumulator, found no positive median, and returned — leaving the width unset *and* the buffer
empty, so the next windows started a second scale sample. The sample then got a histogram whose
width was fitted from a later stretch of its genome, with the earlier windows folded nowhere and
the count agreeing with the cells. Measured on six windows at depths `0,0,0,80,30,30`: six
emitted, none absent, **`windows_folded` = 2**. Two doc comments promised the opposite.

The plan names this failure for this step: *"its failure is silent — a width off by a factor is a
histogram whose mode sits in the wrong bin, and every fit anchors on it."*

The fix replaces the width's `Option<f64>` and its separate buffer with one three-state value —
still collecting, fitted, or **unfittable and latched** — so a sample that cannot be scaled stays
unscaled, and the held-back list stays bounded for free. Two more the reviews found: the fitted
width was unguarded against a contig reset (clearing it there left all 42 tests green, while a
probe showed a width of 2.5 becoming 10 on a deeper second contig — one histogram's rows cut on
two axes); and the differential's narrowing had thrown out two properties that have nothing to do
with the depth axis, that a histogram comes back at all and that every emitted window reaches it.

## Assumptions the spec leaves open

1. **A median of zero, or worse, means no histogram** rather than some fallback width. A width of
   zero puts every window in the overflow column, and there is no defensible number to invent
   instead: the sample reported no depth for a yardstick to be built from. Same answer, same
   reason, as for a sample with no windows.
2. **The median is the upper of the two middles on an even count**, not their average. Either is
   a defensible median; this one needs no arithmetic on the values, so the width a sample gets is
   derived from a depth that sample actually reported.
3. **The held-back windows are stored as the exact `f64` pair the fold will use**, so a window
   folded after the fit lands in the cell it would have landed in had the width been known when
   it was finalised. That equality is asserted by a test rather than argued.

## What the differential against production loses, and why that costs nothing

`production_parity.rs` compared both the emitted windows and the histogram, bit for bit. **The
histogram half is gone at A3**: production cuts its depth axis at a fixed width and ng at one
fitted per sample, so the two matrices are not the same object and comparing them would mean
either freezing ng's width or reading production's — both of which would test the wrong thing.

**The windows are still compared, and are still bit-identical over 447,581 of them** — the same
count as before A3, which is what says the fit changed nothing about what a window reports.

The narrowing was paid for in advance. A1's review found four behaviours the differential was
the sole guard of — the overflow column's boundary, the GC clamp, lowercase reference bases, and
the closing frontier — and gave each a unit test at that point, before the histogram comparison
went away.

## Changes made

- **[`mod.rs`](../../../../src/ng/window_coverage/mod.rs)** — `depth_bin_width` leaves
  `WindowCoverageConfig` and appears only on the finished histogram, where its doc says it is
  fitted rather than chosen. `depth_scale_windows` and `depth_range_in_medians` join the
  configuration, with guards. Three new provisional constants — `DEPTH_BINS` 400,
  `DEPTH_SCALE_WINDOWS` 10,000, and `DEPTH_RANGE_IN_MEDIANS` 10.0, a multiple rather than a count
  so that D2 can land between two whole numbers — each marked soft, against the two that are
  production's and settled. `cell_index` takes the fitted width as an argument.
- **[`accumulator.rs`](../../../../src/ng/window_coverage/accumulator.rs)** — the width and the
  windows held back for it as one three-state value, still collecting or fitted or unfittable;
  `WindowMeans`, a named pair rather than a `(f64, f64)` read by index, because transposing the
  two would fit every width from a median GC fraction and leave the differential green;
  `fold_or_hold_back`, which folds when the width is known and buffers when it is not; the fit,
  which takes the median, sets the width and folds everything held back; and a `finish` that fits
  from whatever it has and returns `None` when there is nothing to fit from.
- **[`production_parity.rs`](../../../../src/ng/window_coverage/production_parity.rs)** — the
  histogram's *cells* no longer compared, with the reason in the module's own doc; two properties
  that do not depend on either side's depth axis kept and now asserted explicitly — that a stream
  of real windows produces a histogram at all, and that every window emitted reaches it.

## Tests

**48 in the module now** — six added with the step, six more the reviews required, and one
retired: the `NaN`-width configuration guard, whose hazard moved from a value a caller could type
to a value the fit could produce, and which came back as `new_panics_on_a_nan_depth_range` once
the range multiple became a float.

| test | what it pins |
|---|---|
| `the_depth_width_is_fitted_from_the_median_of_the_sample_own_window_depths` | depths 1 to 5 fed out of order: median 3, range 7, 5 bins → width 4.2, so neither the first nor the last value could have produced it |
| `a_stream_gives_one_histogram_whether_its_windows_are_folded_early_or_late` | one stream run with the width fitted after the first window and with every window held back to the end: identical histograms, over a fixture that fills more than one cell |
| `a_stream_shorter_than_the_scale_sample_fits_its_width_at_the_end` | three windows against a scale sample of a million: the width is fitted at `finish` from what there is |
| `sliding_empty_stream_yields_no_window_and_no_histogram` | nothing measured, so no invented yardstick |
| `a_sample_whose_windows_all_report_zero_depth_gets_no_histogram` | a median of zero is not a scale |
| `new_panics_on_a_zero_depth_scale_sample`, `new_panics_on_a_zero_depth_range`, `new_panics_on_a_nan_depth_range` | the new configuration guards |
| `a_failed_fit_does_not_start_a_second_scale_sample` | the defect above: a sample that cannot be scaled stays unscaled |
| `the_fitted_depth_width_survives_a_contig_change` | the axis is the sample's, not the contig's |
| `the_scale_sample_is_exactly_as_many_windows_as_it_says` | `>=` slipping to `>` moves a fitted width from 2.5 to 1.0 |
| `median_depth_returns_the_upper_middle_on_an_even_count` | empty, one, even, odd, all-equal, and a transposed pair |
| `a_held_back_window_lands_in_the_same_cell_as_one_folded_immediately` | held-back and live folds agree **at depths that differ**, which the uniform-depth test above cannot see |
| `sliding_uniform_all_one_cell` (strengthened) | the cell is named, not just counted |

Four fixtures whose windows are all under A2's floor changed their claim rather than their
arithmetic: they used to assert `windows_folded == 0` and now assert that there is no histogram
at all, which is the same fact one step further on.

**Nine mutations, each caught, run in the container on this tree** — four against the step as
first written and five against it after the reviews' fixes:

| mutation | the test that fails |
|---|---|
| the median replaced by the smallest of the held-back depths | `the_depth_width_is_fitted_from_the_median_of_the_sample_own_window_depths` |
| the held-back windows never folded after the fit | `a_stream_gives_one_histogram_whether_its_windows_are_folded_early_or_late`, and 15 others |
| the non-positive-width guard removed | `a_sample_whose_windows_all_report_zero_depth_gets_no_histogram` |
| `finish` not fitting from what it has | `a_stream_shorter_than_the_scale_sample_fits_its_width_at_the_end` |
| the unfittable latch removed, so a failed fit retries | `a_failed_fit_does_not_start_a_second_scale_sample`, alone |
| the fitted width cleared at a contig change | `the_fitted_depth_width_survives_a_contig_change`, and the differential's restored assertion |
| the scale sample's `>=` slipped to `>` | `the_scale_sample_is_exactly_as_many_windows_as_it_says`, alone |
| the median taken on the GC fraction | `median_depth_returns_the_upper_middle_on_an_even_count`, and 12 others |
| held-back windows folded against the median rather than the width | `a_held_back_window_lands_in_the_same_cell_as_one_folded_immediately`, and 4 others |

## Validation results

In the container (`./scripts/dev.sh`):

- `cargo test --lib --all-features` — **6,323 passed, 0 failed, 15 ignored** (6,312 at A2; 6,275
  without this module). 48 in this module.
- `cargo clippy --lib --all-features` — no diagnostic in `src/ng/window_coverage/`.
  `rustfmt --check` clean on all three module files.
- **The standing calling oracle, re-run on this tree because the plan lands A3 alone:** six
  tomato accessions over the first two 100 kb intervals of `benchmarks/tomato1/regions.bed`,
  called from the CRAMs and from stored files — **2,311 records, sha256
  `84ad19c22dd14de583cd85805dcd2e5169e799d7a63691c979b7fa43d400590d` on both sides, identical to
  the baseline taken before A1**. Nothing in Milestone A moved a byte of a VCF.

## Tradeoffs and follow-ups

- **Three of the four numbers this step introduces are provisional**, and the code says so in
  each one's first sentence: how many bins the depth axis has, how many windows the scale is
  fitted from, and how many medians the range spans. Plan step D2 measures the overflow fraction
  on both benchmarks against the fit's own rejection guard and confirms or moves them.
- **The held-back list costs more than the spec says.** Sixteen bytes a window is right, but the
  list grows by doubling, so at the configured 10,000 windows its capacity reaches 16,384 and the
  transient peaks at **262 kB**, not the 120 kB spec §3.4 quotes. That paragraph's budget still
  closes — 123 + 80 + 262 = 465 kB against 500 — but with about 35 kB of headroom rather than the
  ~180 the quoted figure implies. Plan step D3 measures it.
- **The absent-sample answer is one bit, not a reason.** A caller at plan step C5 is told there
  is no histogram, not whether the sample had no windows, had them all silenced by the floor, or
  reported no depth. The review implemented a four-variant answer and I deferred it: spec §3.6
  fixes `finish`'s signature and §3.5 its prose, so it is the owner's call at Checkpoint A. At
  one low-coverage sample — the project's declared hardest corner — "every window under the
  floor" is the likely case, and it is a reading on `MIN_WINDOW_POSITIONS` that plan step D1
  exists to settle.
