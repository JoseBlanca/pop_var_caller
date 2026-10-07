# Fix Application Report: fit_precision_a5_2026-09-27.md

**Date:** 2026-09-27
**Source review:** `doc/devel/reports/reviews/fit_precision_a5_2026-09-27.md`
**Source state reviewed against:** `c032eab0` (review-only commit), branch `fit-precision`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 2
- Minors: 7
- Nits: 8

### Outcome totals
- Applied: 16 (eight of them nits)
- Applied with adaptation: 1
- Disputed: 0
- Deferred: 0 (one of the review's missing tests is deferred, §4)
- Failed validation: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --all-targets --all-features --no-fail-fast` → counts in the commit message
  (`tmp/fit_precision/suite_a5b.log`); the three pre-existing failures only; both
  `cli::cross_platform_digests` tests pass unchanged
- `cargo test --release --lib parameter_estimation::joint::fit` → 0, 59 passed, 1 ignored
- Performance check → not applicable: the fixes change what is printed and traced, not a pass

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files changed | Follow-up |
|---|---|---|---|---|---|---|
| M1 | Major | "no error" reads as "error-free" | Apply | Applied | `standard_errors.rs` | No |
| M2 | Major | the trace's error rows untested | Apply | Applied | `fit_trace.rs`, `fit.rs`, `information.rs` | No |
| Mi1 | Minor | the log line's arguments untested | Apply | Applied | `standard_errors.rs`, `information.rs` | No |
| Mi2 | Minor | `named` repeats the trace's names | Apply | Applied | `standard_errors.rs` | No |
| Mi3 | Minor | the `expect` a function away from its guarantee | Apply | Applied | `fit.rs` | No |
| Mi4 | Minor | the ratio untested; "inf", "(0ms)" | Apply | Applied | `fit.rs` | No |
| Mi5 | Minor | "the ratio measures the information's cost" | Apply | Applied with adaptation | `fit.rs` doc, report | No |
| Mi6 | Minor | wrong and loose numbers in the report | Apply | Applied | report | No |
| Mi7 | Minor | `fit_jointly`'s eight-value tuple | Apply | Applied | `fit.rs` | No |
| Nits | Nit | eight | Apply | Applied | `standard_errors.rs`, `fit.rs`, `fit_trace.rs` | No |

## 3. Questions asked and answers

None.

## 4. Per-finding log

### M1 — the log's wording
- **Final status:** Applied. Every string says "standard error": "median standard error …, largest …,
  none missing", "2 missing (1 no information, 1 not identified)", "none has a standard error (4
  held fixed)", "no standard error (it came out at …, wider than the parameter's whole range)"; the
  count reads "(4 samples)". Each cohort-level parameter now prints its value beside its error
  (`1.9811e-2 ± 4.89e-3`), the design reviewer's open question: an error alone does not say whether
  it is tight. The names are the reader's: "fixed non-reference share", "allele-frequency shape a",
  "carrier-frequency shape a", "error rates at ordinary positions", "error rates at mismapped
  positions"; the start line says "collects what the standard errors are computed from" rather than
  "sums the information". The line's head says "at the returned values", which holds also when no
  pass ran.
- **Tests:** `the_logged_line_summarises_each_kind` pins the new line exactly.

### M2 — where the trace's rows land
- **Final status:** Applied. `fit_trace::captured` (test-only) runs a fit with the calling thread's
  rows captured instead of written; the fit writes every row from that thread.
  `the_trace_files_the_returned_fit_after_the_winning_starts_last_pass` checks the rows after the
  winner's last pass: start 2 (the fixture's winner, neither first nor last), pass `passes + 1`, the
  returned values and their errors under the value rows' names, the final pass's log-likelihood.
  Writing them at `passes` fails it (run, then restored and checked byte for byte).
- **Adopted from the correctness reviewer's out-of-scope note:** the trace now also writes the
  returned values beside their errors, since a start whose last accelerated step is refused at the
  pass limit returns the step before it, and the last pass's rows are then not what was returned.

### Mi1 — the line's arguments
- **Final status:** Applied. `described` takes the parameters and the read-group map and works out
  the duplicated class and the later read groups itself, so no caller can pass them wrong; tested at
  one later read group, and through a whole fit with the class off
  (`a_fit_without_the_duplicated_class_says_so`).

### Mi2 — `named`
- **Final status:** Applied. It takes `&Parameters` and zips its errors with `Parameters::named`,
  asserting one error a name; the literal list and the flag are gone.

### Mi3 — the `expect`
- **Final status:** Applied. `maximise` takes the information out of the pass that asked for it, with a
  `PANIC-FREE` note, and `StartOutcome` carries it; `fit_jointly` has no `expect` for it.

### Mi4 — the ratio
- **Final status:** Applied. `cost_of_the_last_pass(whole, last, passes)` and `finely` (tenths of a
  millisecond below a second, hundredths of a second below a minute); no ratio when there is no pass,
  or none measurable, before the last; `the_last_passes_cost_reads_as_a_ratio`.

### Mi5 — what the ratio measures
- **Final status:** Applied with adaptation. The design reviewer held that the ratio does not isolate
  the information, since the final pass has no maximisation and keeps the per-position posteriors;
  the correctness reviewer measured both at nothing measurable (0.436 s against 0.435 s), and the
  ratio with against without the information at 1.43, 1.72 and 1.75. So the claim stays, with the
  measurement in `cost_of_the_last_pass`'s doc and the report, and the one unmeasured caveat (a run
  keeping genotype posteriors).

### Mi6 — the report's numbers
- **Final status:** Applied. 2.11 errors at pass 90 (not 2.0); shape `b` first at 50 on pass 101,
  there from pass 105; shape `a` 20.13, 19.87, 19.74 over passes 198–200, not "still climbing"; the
  "lower bounds" claim withdrawn; the cross-reference is to §6.

### Mi7 — the tuple
- **Final status:** Applied. `FittedCohort`; the score is the winning statistics' log-likelihood, so
  `best` no longer duplicates it.

### Nits
- Applied: `reason()` per variant and named reason constants in place of the parallel index; the name
  constants private; `winning_start`; `JointFit::standard_errors`' doc says the duplicated class's
  three read no information when not fitted; `fit_trace.rs` says the returned fit's rows carry the
  final pass's log-likelihood, and that a process fitting twice interleaves both; the ratio's doc
  counts SQUAREM's jump passes; the log vocabulary (M1).
- The 'standard errors at the fitted values' heading, printed also when no pass ran, now reads 'at the returned values' (M1).

### The review's missing tests
- Added: every one but `a_sample_without_reads_has_no_error_in_the_fit`, which is deferred. The
  pass-level test from step A2 (`a_sample_without_reads_carries_no_information`) pins the blocks'
  zeros, and step A3's unit tests pin that zero information gives `NoInformation`; the reliability
  reviewer ran the whole-fit case as a probe and it held. Building a census with a silent sample in
  the tests' fixtures is the cost.

## 5. Deferred findings to carry forward
- The whole-fit test of a sample without reads (above).

## 6.–8.
None disputed, failed or blocked.

## 9. Performance check
- Not triggered.

## 10. Commands run
- `scripts/dev.sh bash -c 'cargo fmt && cargo clippy --all-targets --all-features -- -D warnings; cargo test --release --lib parameter_estimation::joint::fit'`
- `scripts/dev.sh bash -c 'cargo test --release --lib the_trace_files_the_returned_fit'` (with the
  rows filed at `passes`; failed as expected; restored)
- `scripts/dev.sh bash -c 'cargo test --release --lib a_fit_carries_the_errors_at_the_parameters_it_returns -- --nocapture'` (the report's example lines)
- `scripts/dev.sh bash -c 'cargo fmt --check && cargo clippy --all-targets --all-features -- -D warnings; cargo test --all-targets --all-features --no-fail-fast'`

## 11. Command results
- fit module → 0, 59 passed, 1 ignored
- the planted slip → 1 failed; after the restore the module passes again
- fmt, clippy, full suite → see §1 and the commit message

## 12. Notes
- The oracle was not re-run after the fixes: they change the log's wording and the trace's last
  rows, not a computed number, and the full suite's two cross-platform checksum tests are the check
  on that.
