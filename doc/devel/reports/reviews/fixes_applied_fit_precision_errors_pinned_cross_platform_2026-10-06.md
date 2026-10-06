# Fix Application Report: fit_precision_errors_pinned_cross_platform_2026-10-06.md

**Date:** 2026-10-06
**Source review:** `doc/devel/reports/reviews/fit_precision_errors_pinned_cross_platform_2026-10-06.md`
**Source state reviewed against:** review object `755d29b6`
**Execution mode:** non-interactive
**Overall status:** Completed

---

| ID | Final status | What was done |
|---|---|---|
| M1 | Applied | `errors_computed_a_block_at_a_time_write_the_same_bytes_on_every_platform`: the fixture grown to 21 samples (`a_larger_varying_cohort_with_sequencing_errors_on_disk`), its own two checksums, and in the thread-count test. The reviewer's surviving mutation — the block sums rounded differently — now fails this test and only this test (planted and restored). |
| Mi1 | Applied | The constant's doc says the cohort's own scores, the mismapped share's among them, are not pinned by that pair. |
| Mi2 | Applied | "never on the base that read's own sample varies at (another sample's site, and the tract, can take one)". |
| Mi3 | Applied | "few positions collect the same wrong base from more than one read". |
| Mi4 | Applied | The ranges re-measured after the hash change and written as measured. |
| Mi5 | Applied | The comment says what the thread-count test guards and points to `fit`'s many-chunk test for the join order. |
| Mi6 | Applied | "recorded on arm64 Linux (glibc) and arm64 macOS; x86_64 not yet compared". |
| Nits | Applied | The whole name is hashed; the guard's coordinates explained; the constants named (splitmix64's); an array for the three bases. The boolean parameter kept: its two callers are named wrappers, and a third argument now says how many further samples. |

Changing the hash moved the two-sample error fixture's file, so its checksums were recorded again, with the new
21-sample pair, in the Linux container and confirmed natively on macOS. Validation in the implementation report.
