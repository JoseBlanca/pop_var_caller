# Fix Application Report: fit_precision_a4_2026-09-27.md

**Date:** 2026-09-27
**Source review:** `doc/devel/reports/reviews/fit_precision_a4_2026-09-27.md`
**Source state reviewed against:** `42400835` (review-only commit), branch `fit-precision`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 6
- Minors: 11
- Nits: 7

### Outcome totals
- Applied: 18 (five of them nits)
- Applied with adaptation: 3
- Deferred: 3 (two of them nits)
- Disputed: 0
- Failed validation: 0
- Awaiting user answer: 0 (two items carried to checkpoint A)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --all-targets --all-features --no-fail-fast` → counts in the commit message
  (`tmp/fit_precision/suite_a4b.log`); the three pre-existing failures only; both
  `cli::cross_platform_digests` tests pass unchanged
- `cargo test --release --lib the_errors_mean_what_they_say -- --ignored --nocapture` → 0, 1 passed,
  every check in `check_the_coverage` held (`tmp/fit_precision/a4_coverage4.log`)
- Performance check → not applicable (tests only)

### Unresolved high-priority findings
- None in this step's code. Two readings go to checkpoint A: spec question 3 (M3, M5) and the
  Beta-shape update, which the review measured should be fixed before plan step B1 (M4).

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files changed | Follow-up |
|---|---|---|---|---|---|---|
| M1 | Major | the coverage test cannot fail | Apply | Applied | `information.rs` | No |
| M2 | Major | missing kinds and an empty run pass | Apply | Applied | `information.rs` | No |
| M3 | Major | "the full matrix holds" overstates | Apply | Applied with adaptation | report, `information.rs` doc, PROJECT_STATUS | Checkpoint A |
| M4 | Major | the gap grows with the positions | Apply | Applied | report | Checkpoint A |
| M5 | Major | the shapes' errors hold at the maximum; the invariant share's at 4 samples do not; no spread | Apply | Applied | `information.rs`, report | Checkpoint A |
| M6 | Major | `full_matrix` untested | Apply | Applied | `information.rs` | No |
| Mi1 | Minor | a bad cohort count passes silently | Apply | Applied | `information.rs` | No |
| Mi2 | Minor | the variable never reaches the container | Apply | Applied | `information.rs`, report | No |
| Mi3 | Minor | full-matrix count hidden | Apply | Applied | `information.rs` | No |
| Mi4 | Minor | every sample drawn at the same values | Defer | Deferred | None | Later |
| Mi5 | Minor | truths and estimates read by position | Apply | Applied | `information.rs` | No |
| Mi6 | Minor | the 48-node probe not like for like | Apply | Applied | report | No |
| Mi7 | Minor | fewer than 4 samples unmeasured | Apply | Applied with adaptation | `information.rs`, report | No |
| Mi8 | Minor | the stop rule not ruled out | Apply | Applied | report | No |
| Mi9 | Minor | the coupling stated as fact | Apply | Applied | report | No |
| Mi10 | Minor | "1 to 2 errors" | Apply | Applied | report | No |
| Mi11 | Minor | the digamma cause stated as fact in the doc | Apply | Applied with adaptation | `information.rs` doc | No |
| Nits | Nit | seven | Apply (5), Defer (2) | Applied / Deferred | `information.rs`, report | No |

## 3. Questions asked and answers

None asked. Carried to checkpoint A: spec question 3 per kind (M3, M5), and the Beta-shape update
before B1 (M4).

## 4. Per-finding log

### M1 — the coverage test cannot fail
- **Final status:** Applied. `check_the_coverage` holds each regime, at the default 200 cohorts
  only (the draws are fixed by their seeds, so the numbers move only with the code): at most five
  fits stop at the pass limit (4 and 20 samples); where the errors cover — 20 samples at both depths,
  4 samples at 30 reads — the error rates, the homozygote excess and the mismapped share, and at 20
  samples the fixed share, land between 0.58 and 0.78 within one of the blocks' errors and between
  0.90 and 0.99 within two; at 4 samples and 3 reads the full matrix puts at least 0.08 more of the
  error rates within one error than the blocks; at 20 samples the Newton-corrected estimate of each
  density shape sits within half an error of the truth on average.
- **Which reviewed mutations each check kills** (by construction, from the reviewers' printed
  numbers; not re-run as mutations): a wrong truth or error slot, a within-one threshold of 1.5, and
  doubled block errors move the within-one share out of its band; blocks and full swapped, or the
  full matrix reduced to the blocks, break the ordering at 4 samples and 3 reads; a flipped Newton
  sign puts the corrected shapes about 2.6 errors off at 20 samples; shapes `a` and `b` swapped move
  the corrected shapes several errors off.

### M2 — missing kinds
- **Final status:** Applied. The tally must hold exactly the eight kinds, each estimated by both
  matrices (`full.count == blocks.count`), at any cohort count; the count must be a positive whole
  number (Mi1).

### M3 — the full matrix is too wide where it differs
- **Final status:** Applied with adaptation: the reading is corrected, but not to the review's. The report's §2.2 and §6, the A3 comparison's doc and PROJECT_STATUS
  say which side is off and by how much, now from the spread the test prints (M5): at 4 samples and
  3 reads the estimates scatter 1.23 to 1.67 times the blocks' errors and 0.89 to 1.02 times the
  full matrix's — close, and on the wide side, most for the clean rates (about 12%); at 20 samples
  the two are about 10% apart, either side of right. The reviewers' "three to four spreads too
  wide" read the within-one share, which for the mismapped share is high because its distances are
  peaked, not because its errors are wide (spread 1.02).

### M4 — the gap grows with the positions
- **Final status:** Applied. The report quotes the correctness reviewer's measurement (in `a`'s own
  units the gap shrinks with the samples and not with the positions, so in errors it grows as the
  square root of the positions) and recommends fixing the shape update before B1. Not re-measured
  here: the reviewer's probe was a test in its own worktree, and the figures are quoted as the
  review's.

### M5 — the spread, and the coverage at the likelihood's maximum
- **Final status:** Applied. `Tally` now keeps the sum of squares, so each kind prints its spread;
  `Coverage::at_flat_point` tallies the Newton-corrected estimate against the truth. The report
  splits the reading per kind and answers spec question 3 per kind.

### M6 — `full_matrix` untested
- **Final status:** Applied. `the_full_matrix_sums_the_products_between_two_samples`: hand-built
  scores for two samples at two positions; a cross-sample entry (45), a cohort–sample entry (31) and a
  diagonal entry (25), each the sum of its products. A full matrix without the products between two
  samples gives 0 for the first.

### Mi1–Mi3, Mi5
- Applied. Mi1: a value that is not a positive whole number panics with its text. Mi2: renamed
  `NG_FIT_PRECISION_COVERAGE_COHORTS`, which `dev.sh` forwards; the doc gives the `dev.sh env …`
  form too. Mi3: each line prints both counts. Mi5: `TrueValues` (one named struct, which also feeds
  the draw) and `value_at`, which reads a slot by its named constant for both the truth and the
  estimate.

### Mi4 — every sample at the same values
- **Final status:** Deferred. Drawing each sample at its own rates changes every draw, so every
  number of the measurement, and the reviewer confirmed the present indexing is correct (read groups
  are numbered in sample order). Recorded at `value_at`'s doc.

### Mi6, Mi8–Mi10
- Applied in the report: the 48-node probe against the same 40 seeds at 16 nodes
  (`tmp/fit_precision/a4_probe16_40.log`), with the estimates' own shift and the node-count law the
  correctness reviewer measured; the stop thresholds 10⁴ times tighter, and 300 further passes, leave
  the shapes where they were; the coupling's direction, marked as not measured; "1 to 2.3 errors".

### Mi7 — fewer than 4 samples
- **Final status:** Applied with adaptation. A regime of 2 samples over 30,000 positions at 3 reads —
  the corner the A3 review named. Not at 30 reads: at 4 samples both matrices already hold at 30
  reads, and the 2-sample regime costs about 11 minutes. Its fits stop at the pass limit 34 times in
  200, so it is printed and held only to the checks that apply at any count.

### Mi11 — the cause stated as fact
- **Final status:** Applied with adaptation. The correctness reviewer's probes established the cause
  (the fit is an exact fixed point of the update, where the digamma form is zero to 0.012 and the
  likelihood's slope in `a` is −8 to −20), so the doc now states it with that evidence rather than
  hedging it.

### Nits
- Applied: `Coverage` holds three named `Tally`s instead of `[_; 2]` arrays; `label_of(slot)`; "fitted
  by the fit's own rule"; "floors … at 0.2"; "0.938 to 0.965".
- Deferred: the true values as literals in the A3 tests (touches step A3's tests, no defect); the
  seeds' overlap between depths (harmless, and changing it would redraw every 30-read cohort).

## 5. Deferred findings to carry forward
- Mi4 — per-sample truths in the coverage draw.
- Nits — shared true values across the A3 tests; seed stride.

## 6.–8.
None disputed, failed or blocked.

## 9. Performance check
- Not triggered: tests only.

## 10. Commands run
- `scripts/dev.sh bash -c 'cargo fmt && NG_FIT_PRECISION_COVERAGE_COHORTS=3 cargo test --release --lib parameter_estimation::joint::fit::information -- --include-ignored --nocapture'`
- `scripts/dev.sh cargo test --release --lib the_errors_mean_what_they_say -- --ignored --nocapture`
  (twice: the first, `a4_coverage3.log`, failed its own convergence check on the new 2-sample regime,
  which was then scoped to 4 and 20 samples; the second is `a4_coverage4.log`)
- `scripts/dev.sh bash -c 'FIT_PRECISION_COVERAGE_COHORTS=40 cargo test --release --lib the_errors_mean_what_they_say -- --ignored --nocapture'` at 16 nodes, before the rename, for the like-for-like probe
- `scripts/dev.sh bash -c 'cargo fmt --check && cargo clippy --all-targets --all-features -- -D warnings && cargo test --all-targets --all-features --no-fail-fast'`

## 11. Command results
- module tests with the short coverage run → 0, 17 passed
- the full coverage run → see §1
- fmt, clippy, full suite → see §1 and the commit message

## 12. Notes
- The 48-node and 16-node probes were run by editing the test from a backup copy and restoring it,
  checked byte for byte.
