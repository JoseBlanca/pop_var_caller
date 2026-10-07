# Fix Application Report: fit_precision_a6a_2026-09-28.md

**Date:** 2026-09-28
**Source review:** `doc/devel/reports/reviews/fit_precision_a6a_2026-09-28.md`
**Source state reviewed against:** `e529839d` (review-only commit on `fit-precision`)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 4
- Minors: 10
- Nits: 3 (grouped)

### Outcome totals
- Applied: 8 (B1, M1, M2, M3, M4, Mi4, Mi5, Mi6)
- Applied with adaptation: 3 (Mi1, Mi2, Mi7)
- Already fixed: 0
- Deferred: 3 (Mi3, Mi9, Mi10)
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 1 (Mi8, by the report's §3 item 8)
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`,
  `cargo test --all-targets --all-features --no-fail-fast`: see §11 (`tmp/fit_precision/suite_a6a_fixes.log`).
- `cargo test --release --lib parameter_estimation::joint::fit` → 63 passed, 1 ignored.
- Mutation re-check: the six survivors below re-run against the fixed tests, all six fail
  (`tmp/fit_precision/a6f_mutations.log`), tree restored and compared byte for byte.
- `cargo doc`, `cargo audit` → not run (no public API or dependency change).
- Performance check → not triggered by the fixes (tests and the interim error marking only; the
  likelihood's cost was measured in review, recorded in the implementation report).

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| B1 | Blocker | three branches' per-library crediting untested | Apply | Applied | information.rs (test) | mutation `seg_first`, `dup_first` fail |
| M1 | Major | other errors computed as if a multi-library sample's rates were known | Apply | Applied | standard_errors.rs, information.rs, fit.rs (docs) | tests |
| M2 | Major | information test does not rule out one library's rates for all | Apply | Applied | information.rs (test) | mutation `last_rates` fails |
| M3 | Major | per-library depth handling untested | Apply | Applied | information.rs (test) | mutation `sample_cov` fails |
| M4 | Major | a sample's positions with reads untested for several libraries | Apply | Applied | fit.rs (generator, test) | mutation `first_depth` fails |
| Mi1 | Minor | "up to 482" figure | Apply | Applied with adaptation | impl report | — |
| Mi2 | Minor | "above the census cap" imprecise | Apply | Applied with adaptation | impl report | — |
| Mi3 | Minor | empty-library branch reachable, unreported | Defer | Deferred | impl report (recorded) | — |
| Mi4 | Minor | the marking call untested | Apply | Applied | fit.rs (test) | mutation `no_marker` fails |
| Mi5 | Minor | interim rule decided in four places | Apply | Applied | standard_errors.rs | tests |
| Mi6 | Minor | stale doc on `sample::CLEAN_ERROR_RATE` | Apply | Applied | information.rs | — |
| Mi7 | Minor | log wording | Apply | Applied with adaptation | standard_errors.rs | test |
| Mi8 | Minor | cost not stated | Apply | Superseded | impl report | — |
| Mi9 | Minor | a library's `observations` is the position count | Defer | Deferred | — | — |
| Mi10 | Minor | `Scratch::new` takes three counts by position | Defer | Deferred | — | — |

## 3. Questions asked and answers

None.

## 4. Per-finding log

### B1 — three branches' per-library crediting untested
- **Final status:** Applied. `each_librarys_tallies_hold_its_own_reads_in_every_branch`: over each
  of 600 positions, with sample 0 split over two libraries at different rates and the duplicated
  class fitted, each library's tallies summed over classes, branches and genotypes hold its own
  non-reference reads and its own reference reads (worst relative gap 8.5 × 10⁻¹², asserted below
  10⁻⁹); every branch carries mass (invariant 531.8, fixed 5.4, segregating 59.3, duplicated 3.5).
  The review's sketch checked non-reference reads only; the reference reads are added. The
  whole-fit test's fixture was not changed (the reviewer suggested a lower invariant share); the
  conservation test is exact and does not need it.
- **Verification:** crediting the segregating branch (`seg_first`) or the duplicated branch
  (`dup_first`) to a sample's first library fails it.
- **Residual risk:** under a depth range, the rate each branch's expected reference reads is taken at
  can be the wrong library's without this test seeing it (reference reads are exact at a single
  depth, and the new range test checks weights, not tallies). No range occurs at the shipped cap.

### M1 — errors computed as if a multi-library sample's rates were known
- **Final status:** Applied. `with_errors_awaiting_each_librarys_scores` now marks **every** error of
  a cohort that holds a sample of several libraries (the variant is renamed
  `AwaitingEachLibrarysScores`, its text "not computed yet where a sample has several libraries"),
  and leaves a cohort of one library a sample untouched (asserted). This follows the plan's split
  clause ("a multi-library sample's errors reported absent") with the review's measurement: the
  sample's own rates are not the only errors that go wrong. Docs updated in all three files.

### M2 — the information test does not rule out one library's rates for all
- **Final status:** Applied. `a_librarys_reads_are_scored_under_its_own_rates_and_no_other`: an
  empty second library leaves the likelihood bit-for-bit the first library's alone, at two rate
  pairs; exchanging the libraries' reads with their rates leaves it bit-for-bit unchanged.
- **Verification:** scoring every library under its sample's last library's rates (`last_rates`)
  fails it.

### M3 — per-library depth handling untested
- **Final status:** Applied. `each_library_is_weighted_around_its_own_depth` (two samples, a 40-read
  and a 132-read library each, cap 140; 597 ranged library-positions).
- **Verification:** centring each range on the sample's first library (`sample_cov`) fails it.
- **Residual risk:** as B1's; and no fixture has a library unwalked where its sibling is walked.

### M4 — positions with reads untested
- **Final status:** Applied. The generator counts each sample's drawn positions with one and two
  reads (its libraries added, no random number drawn); the two-library whole-fit test asserts them.
- **Verification:** counting only the first library (`first_depth`) fails it.

### Mi1, Mi2 — the report's figures
- **Final status:** Applied with adaptation: the report now says 483 to 964 read groups on kimura,
  and that per-library coverage moves nothing at the shipped cap of 124. The plan's own "up to 482"
  is the owner's text and is raised at checkpoint A′ rather than edited.

### Mi3 — empty-library branch
- **Final status:** Deferred, recorded in the implementation report's §3 item 7. A sample holding no
  ordinary-position section panicked before this step; now it contributes no reads. A test needs a
  census of repeat tracts only, which the drawing generator does not build.

### Mi4 — the marking call untested
- **Final status:** Applied (the whole-fit test asserts every error is marked); `no_marker` fails it.

### Mi5 — the interim rule in four places
- **Final status:** Applied. `named` reads each read group's sample row and no longer re-derives
  the rule; the rule is decided in `with_errors_awaiting_each_librarys_scores` alone.

### Mi6, Mi7
- **Final status:** Mi6 applied. Mi7 applied with adaptation: the reason reads "not computed yet where
  a sample has several libraries", and the closing clause says how many samples and why, without
  "scored".

### Mi8 — cost
- **Final status:** Superseded by the implementation report's §3 item 8, which quotes the review's
  measurement (2.3 to 3.6% slower for one library a sample).

### Mi9, Mi10
- **Final status:** Deferred. `observations` of a library's noise estimate is the position count for
  every library, as before this step; what it should count is step E1's, which rewrites `Estimate`.
  `Scratch::new`'s three counts are private and called from one production site.

## 5. Deferred findings to carry forward
- Mi3 — a census of repeat tracts only, untested.
- Mi9 — a library's `observations`, to step E1.
- Mi10 — `Scratch::new`'s argument order.

## 6. Disputed findings to return to reviewer

None.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Performance check

- **Triggered:** no — the fixes add tests and change which errors are reported absent; the
  likelihood's cost was measured in review.

## 10. Commands run
- `scripts/dev.sh bash -c 'cargo fmt && cargo test --release --lib parameter_estimation::joint::fit …'`
- six mutations, each `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit`
  with the file restored and compared after (`tmp/fit_precision/a6f_mutations.log`)
- `scripts/dev.sh bash -c 'cargo fmt --check; cargo clippy --all-targets --all-features -- -D warnings; cargo test --all-targets --all-features --no-fail-fast'`

## 11. Command results
- module tests → 63 passed, 1 ignored.
- mutations → 6 run, 6 fail a named test.
- full gate → recorded in the commit message (`tmp/fit_precision/suite_a6a_fixes.log`).

## 12. Notes

The one-library arithmetic is untouched by the fixes, so the oracle's checksums (measured on the
reviewed code) still apply: the fixes change tests, the generator's counts (no random draw), and
which errors a multi-library cohort reports.
