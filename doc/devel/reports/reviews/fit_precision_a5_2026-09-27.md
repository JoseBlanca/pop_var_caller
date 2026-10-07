# Code Review: fit_precision_a5
**Date:** 2026-09-27
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step A5 of `fit_precision.md` — the fit computes, prints and carries its standard errors
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the diff `8bcb9cbb..c032eab0` (`c032eab0`, a review-only commit on no branch).
- **In-scope files:** [src/parameter_estimation/joint/fit.rs](../../../../src/parameter_estimation/joint/fit.rs)
  (`maximise` → `StartOutcome`, the final pass, `fit_jointly`, `JointFit::standard_errors`),
  [fit/standard_errors.rs](../../../../src/parameter_estimation/joint/fit/standard_errors.rs)
  (`described`, `named`, `summary_of`), [fit_trace.rs](../../../../src/parameter_estimation/joint/fit_trace.rs),
  the tests in [fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs),
  the A5 implementation report, `PROJECT_STATUS.md`.
- **Categories dispatched**, three reviewers, each in its own worktree:
  - *design and claims* — naming, idiomatic, module structure, refactor safety, smells, defaults;
    conformance to plan A5 and spec §3.3, §3.5; the prose and the printed log line (18 mutations);
  - *reliability and errors* — tests, mutation, edge cases (17 mutations; one sample, a sample
    without reads, `max_passes` 0 to 4, the duplicated class off, no starting point);
  - *correctness, float portability, concurrency, cost and memory* — whether anything moves, the
    final pass's memory at kimura's size, its cost measured directly, determinism, the winning start
    (7 mutations).

## 2. Verdict

**Approve-with-changes.** Nothing moves: the correctness reviewer ran the final pass with and without
the information at 1 and 4 threads on 1, 4 and 64 samples, and every other field came out
bit-identical; the oracle's 18 checksums equal main's; both cross-platform checksum tests pass. The
errors come structurally from the start whose parameters are returned. Every mutation to *which*
errors the fit carries was caught. What no test reads is the wiring from `fit_jointly` to the log
line and the trace, and the log line's wording can be misread.

## 3. Execution status

- Each reviewer: `cargo test --release --lib parameter_estimation::joint::fit` — 55 passed, 1 ignored,
  at `c032eab0`.
- Mutations: reliability 17 run, 5 survived with a behaviour change; design 18 run, 5 survived, 1 no
  behaviour change; correctness 7 run, 3 survived, 1 no behaviour change (`>` to `>=` on the winning
  score: no test produces a tie). The survivors are the same five places in all three, all in
  `fit_jointly`'s calls to the log and the trace.
- The correctness reviewer measured the final pass's cost (one thread, median of five interleaved
  runs): with the information against a plain pass, 1.42 at 4 samples (200,000 positions), 1.72 at 64
  (20,000), 1.75 at 64 with the duplicated class off; the maximisation and the per-position lists add
  nothing measurable (0.436 s against 0.435 s).
- Findings labelled *Needs verification*: 0.

## 4. Open questions and assumptions

1. **Printing the fitted values beside their errors** (design, Nit): the line prints errors alone,
   so a reader cannot tell whether 3.62e-4 on the mismapped share is tight or loose. Adopted for the
   cohort-level parameters, where one value each is short enough to print.

## 5. Top 3 priorities

1. **M1** — the log line says "no error" where it means "no standard error".
2. **M2** — test where the trace's error rows land, since the checkpoint's `SETTLED_FRACTION` reading
   comes from them.
3. **Mi1** — take the log line's and the trace's arguments out of `fit_jointly` into functions a test
   can call.

## 6. Findings

### Major

#### M1: standard_errors.rs `described`, `summary_of` — "no error" reads as "error-free"
- **Categories:** design
- **Confidence:** High
- **Problem:** at one sample the fit printed `homozygote excesses (1): no error, 1 without one (1 held
  fixed)`; beside kinds called "error rates", "every one has an error" is ambiguous; the `(4)` is an
  unlabelled count of samples. Spec §3.2: an absent error is reported as absent, never as zero.
- **Suggested fix:** "standard error" in every such string; "(4 samples)".

#### M2: fit.rs `fit_jointly` — the trace's error rows are not tested to belong to the winning start or to follow its last pass
- **Categories:** reliability, design, correctness (convergent)
- **Confidence:** High
- **Problem:** writing them under another start, or under `passes` instead of `passes + 1` (where they
  collide with the last value rows), passes every test. The report gives the trace as the source of
  the `SETTLED_FRACTION` reading.
- **Suggested fix:** a test-only capture of the trace, and a test on where the rows land.

### Minor

#### Mi1: fit.rs `fit_jointly` — the log line's arguments are untested
- **Categories:** reliability, design, correctness (convergent) · **Confidence:** High
- **Problem:** the later-read-group count without its − 1, and the duplicated-class flag forced
  either way, survive; the unit test uses 0 and 3 later read groups, never exactly 1; no fixture has
  a sample with two read groups.
- **Suggested fix:** functions for the count and the line; tests at one later read group and with
  the class off.

#### Mi2: standard_errors.rs `named` — the trace's names repeated, and a flag each caller derives
- **Categories:** design · **Confidence:** High
- **Suggested fix:** take `&Parameters` and zip with `Parameters::named`.

#### Mi3: fit.rs `fit_jointly` — the `expect` on the final pass's information sits a function away from the flag that guarantees it
- **Categories:** reliability · **Confidence:** High
- **Suggested fix:** take the sums out in `maximise`, where they are asked for, into `StartOutcome`.

#### Mi4: fit.rs — the start line's ratio untested, and a zero average prints "inf" or "(0ms)"
- **Categories:** reliability, design, correctness (convergent) · **Confidence:** High

#### Mi5: fit.rs `StartOutcome` doc, report §1, §6 — "the ratio measures what the information costs"
- **Categories:** design · **Confidence:** High
- **Problem (design):** the final pass has no maximisation and keeps the per-position lists.
- **Resolution:** the correctness reviewer measured both at under 1% of a pass, so the ratio is, in
  practice, the information's cost (1.42, 1.72, 1.75 against A2's 1.41, 1.72, 1.75). The doc says so
  with the measurement; the one caveat is a run keeping genotype posteriors, unmeasured.

#### Mi6: report §6 — wrong and loose numbers
- **Categories:** design, correctness (convergent) · **Confidence:** High
- **Problem:** the invariant share was 2.11 errors (not 2.0) from its pass-200 value at pass 90;
  shape `a` is not "still climbing at pass 200" (20.13, 19.87, 19.74 at passes 198–200) and the
  invariant share falls over those passes too, so the "lower bounds" claim does not hold; shape `b`
  first reaches 50 at pass 101 and stays there from pass 105; the §2 item 1 cross-reference is to §6.

#### Mi7: fit.rs `fit_jointly` — about 340 lines, and an eight-value tuple
- **Categories:** design · **Confidence:** High
- **Suggested fix:** a named struct for the closure's result; the two reports as functions (Mi1).

### Nits

- `absent_reason` indexes a parallel array, and `printed` hard-codes `ABSENT_REASONS[3]`.
- The two name constants are `pub(super)` and used only in their module.
- Log vocabulary: "fixed share" (for which allele?), "density/carrier shape a" (code objects),
  "wider than its range, 5.22e1" (does not say the number is the error), "sums the information"
  (undefined).
- `winner` holds a start's number: `winning_start`.
- `JointFit::standard_errors`' doc does not say the duplicated class's three read `NoInformation`
  when the run does not fit it.
- `fit_trace.rs`: the error rows' log-likelihood is the final pass's; one process that fits twice
  interleaves both fits' rows.
- The ratio's denominator includes the SQUAREM jump passes; say so.
- "standard errors at the fitted values" is printed also when no pass ran (`max_passes` below 3).

## 7. Out of scope observations

- `fit.rs` `maximise` (pre-existing): when a start's last jump is rejected at the pass limit it returns
  the second plain step's values, so the trace's last value row is not the returned parameters. The
  correctness reviewer suggests writing the returned parameters beside the error rows — adopted.
- Memory (correctness): at 2,169 samples and 128 chunks the information adds 573,128 bytes a chunk,
  73.4 MB; the final pass's transient hold goes from about 174 MB to about 247 MB (16 nodes, one read
  group a sample assumed), one pass, three times a fit.

## 8. Missing tests to add now

- `the_trace_files_the_returned_fit_after_the_winning_starts_last_pass` (M2).
- `later_read_groups_counts_all_but_each_samples_first`, and the line at one later read group (Mi1).
- `a_fit_without_the_duplicated_class_has_no_information_on_it` (Mi1).
- `a_fit_with_no_pass_reports_errors_at_its_start` (`max_passes` 0).
- `the_last_passes_cost_reads_as_a_ratio` (Mi4).
- `a_sample_without_reads_has_no_error_in_the_fit` (the pass-level test exists since step A2).

## 9. What's good

- The errors and the returned parameters come out of one stored start, so they cannot disagree, on a
  tie included (correctness).
- `a_fit_carries_the_errors_at_the_parameters_it_returns` compares bit for bit against an independent
  pass, and its fixture's winner is start 2 of 3, so first-start and last-start slips both fail it.
- The start line prints a cost every run measures, which the kimura run will report at scale.

## 10. Commands to re-verify

- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit -- --nocapture`
- `scripts/dev.sh bash -c 'cargo fmt --check && cargo clippy --all-targets --all-features -- -D warnings && cargo test --all-targets --all-features --no-fail-fast'`

Per-category files: `tmp/review_2026-09-27_fit_precision_a5/` (`design.md`, `reliability.md`,
`correctness.md`).
