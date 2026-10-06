# Fit precision, after checkpoint E — the standard errors pinned across platforms

**Date:** 2026-10-06. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md), "Decided at checkpoint
E", item 2. **Branch:** `fit-precision`. Review:
[fit_precision_errors_pinned_cross_platform_2026-10-06.md](../reviews/fit_precision_errors_pinned_cross_platform_2026-10-06.md),
fixes: [fixes_applied_fit_precision_errors_pinned_cross_platform_2026-10-06.md](../reviews/fixes_applied_fit_precision_errors_pinned_cross_platform_2026-10-06.md).

## 1. What was built

The cross-platform check ([cross_platform_digests.rs](../../../../src/cli/cross_platform_digests.rs)) records the
checksums of a small cohort's fitted parameters file and of the calls made with it, so that a change of platform or of
thread count that moves a written number fails a test. Its two-sample fixture has no sequencing error in any read, so
every library's error rate sits on its bound and the file carries no fit standard error: the errors this milestone
writes were not pinned. **Two cohorts are added, each with its own pair of checksums:**

| test | cohort | what its file carries | which computation of the errors |
|---|---|---|---|
| `a_fit_with_standard_errors_writes_the_same_bytes_on_every_platform` | the same two samples, a sequencing error in about one read in four | a standard error on each of 3 multipliers (values 2.97–6.43, errors 1.42–17.71) and on the tract's substitution rate | the whole information matrix (cohorts of at most 20 samples) |
| `errors_computed_a_block_at_a_time_write_the_same_bytes_on_every_platform` | grown to 21 samples, each with a substitution of its own | errors on all 22 multipliers (values 0.35–7.97, errors 1.66–3.01) | a sample's block at a time (above 20 samples, every large cohort's path) |

- **The errors are placed where reads are hashed**: from the sample's name and the read, at an offset and to a base
  that vary, never on the read's own sample's site. Few positions collect one wrong base from two reads, and the calls
  are only the records designed into each cohort (3, and 22).
- **The thread-count test runs all three cohorts** at one, four and seven threads.
- **Recorded in the Linux container (arm64, glibc) and the same on macOS (arm64)**, run natively; no x86_64 machine was
  compared. The first fixture's two checksums, re-recorded in the Linux container at the previous commit, are the same on
  macOS too.
- The fixture's `AVaryingCohort` now carries its samples' names, and the digest run reads its psps by them.

## 2. What it is shown to catch

| change planted | caught by |
|---|---|
| the per-sample block sums rounded differently (`(s + ½·x) + ½·x` for `s + x`) | only the 21-sample test — the other three pass |
| (review) the multiplier's error, or a square root of a variance, multiplied by 1 + ε; the whole matrix's sums rounded differently | the two-sample errors test, the old test passing |
| (review) a chunk size taken from the pool's width | every cohort of the thread-count test |

**Not pinned, measured in the review:** the cohort's score for the share of mismapped positions — computed with the
platform's `exp`, or doubled — moves no byte in the two-sample cohort; a join order chosen by the pool cannot show at
three chunks (held by `fit`'s many-chunk unit test); the inbreeding coefficients carry no error in either cohort.

## 3. Validation

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- The four digest tests: pass in the Linux container and natively on macOS.
- Full suite (`cargo test --all-targets --all-features --no-fail-fast`): 5,043 passed, 3 failed (the pre-existing
  `examples/ng_generic_loci_dump.rs` × 2 and `examples/ng_ssr_loci_dump.rs` × 1), 9 ignored.
