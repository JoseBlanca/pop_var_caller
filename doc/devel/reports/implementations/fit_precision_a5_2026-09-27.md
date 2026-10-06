# Fit precision, step A5 — the fit computes and prints its errors

**Date:** 2026-09-27. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step A5.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §3.3 (the final pass), §3.5.
**Branch:** `fit-precision`. **Review:** [fit_precision_a5](../reviews/fit_precision_a5_2026-09-27.md);
this report is the version after its fixes.

## 1. What was built

Every SNP/indel fit now ends with its parameters' standard errors. The pass the fit already makes at
the parameters it returns — to keep each position's probability of being mismapped — also sums the
information (step A2), on every start; the winning start's sums become the errors (step A3). No
fitted number moves: the pass computes the same statistics as before, and the errors are only
printed and carried.

**What the run's log gains**, on the `estimating:` lines of `estimate-parameters` (examples from a
test's four-sample fixture):

- **each start's line** says how long its final pass took and how many times the average pass before
  it that is. That ratio is what collecting the information costs (§6):

  ```
  SNP/indel fit, start 2 of 3: stopped at the pass limit after 200 pass(es), log-likelihood
  -6.461300e3; its last pass, which also collects what the standard errors are computed from, took
  1.6ms, 1.43 times the average pass before it (1.1ms); …
  ```

- **one line of errors** once the best start is chosen: each cohort-level parameter's value and
  standard error, or why it has none; for each kind of a sample's own parameters, the median and the
  largest standard error and how many are missing, by reason; and how many read groups have none
  because they are a sample's second or later (the likelihood reads only a sample's first, A1 report
  §2 item 1):

  ```
  SNP/indel fit: standard errors at the returned values: mismapped share 1.9811e-2 ± 4.89e-3,
  invariant share 8.8968e-1 ± 4.81e-2, fixed non-reference share 1.2487e-2 ± 2.17e-3,
  allele-frequency shape a 1.1259e0 ± 1.22e0, allele-frequency shape b 3.6846e0 ± 2.13e0,
  duplicated share 3.0474e-9 ± 2.45e-2, carrier-frequency shape a 3.1712e1, no standard error (it
  came out at 5.45e9, wider than the parameter's whole range), carrier-frequency shape b 4.9206e1,
  no standard error (it came out at 8.16e9, wider than the parameter's whole range); error rates at
  ordinary positions (4 samples): median standard error 3.78e-4, largest 4.08e-4, none missing;
  error rates at mismapped positions (4 samples): median standard error 1.51e-2, largest 2.08e-2,
  none missing; homozygote excesses (4 samples): median standard error 9.75e-2, largest 1.00e-1,
  none missing
  ```

**What `JointFit` gains:** a private field `standard_errors`, in the fit's own layout. Each
`Estimate` gains its own error only at plan step E1, which will read this field; until then only the
tests read it.

**What the per-pass trace gains** (`PVC_JOINT_FIT_TRACE`): once the best start is chosen, under that
start and the pass after its last, one row a parameter with the value the fit returns and one with
its standard error (`standard_error:` and the parameter's name, NaN where there is none). Each pass's
move can then be read in units of the parameter's error, against the values the fit returned — which
is what checkpoint A's `SETTLED_FRACTION` decision needs, and what spec §1.1 says is unknown on
kimura ("which parameter was moving on kimura is not known").

## 2. Assumptions and deviations

1. **Every start's final pass sums the information, not only the winner's.** The winner is known
   only once every start has run, and the pass that sums the information is the one each start
   already makes at its returned parameters; summing only on the winner would cost one more pass.
   The price is the collecting pass's extra time on the two losing starts (§6).
2. **`JointFit` carries the errors in a private field.** The plan says `JointFit` gains the errors
   and `Estimate` gains its own only at step E1. A public field would need public types for the
   errors' layout, which E1 replaces; a private one is read by the tests and by E1, both inside the
   module.
3. **A sample's second and later read groups are counted, not listed.** Their errors do not exist
   (the likelihood never reads their rates), so the log gives their number and the trace gives each
   one NaN. How the parameters file writes them is step E2's, with the owner's ruling at checkpoint A.
4. **The median of an even count is the lower of the two middle values**, as the tests' comparisons
   take it.
5. **The trace also gains the returned values and their errors**, which the plan does not name. It is
   the only way a run shows each parameter's move in units of its own error, which checkpoint A's
   `SETTLED_FRACTION` decision is to be taken from; it costs nothing unless `PVC_JOINT_FIT_TRACE` is
   set. The returned values are written because they are not always the last pass's: a start whose
   last accelerated step is refused at the pass limit returns the step before it (found in review).
6. **The line prints each cohort-level parameter's value beside its error** (a review suggestion): an
   error alone does not say whether it is tight or loose.

## 3. Changes

- [fit.rs](../../../../src/parameter_estimation/joint/fit.rs): the final pass of each start keeps the
  information and is timed; `maximise` returns a `StartOutcome` — the parameters, the final pass's
  statistics, the information taken out of them, the passes, whether it converged, its trace and the
  final pass's duration; `fit_jointly` keeps the winning start's number, computes
  `StandardErrors::of`, prints the line of errors, writes the trace rows
  (`trace_the_returned_fit`), puts the errors on `JointFit`, and hands its result out of the
  evidence's loan as a `FittedCohort` in place of an eight-value tuple; `cost_of_the_last_pass` and
  `finely` write the start line's cost.
- [fit/standard_errors.rs](../../../../src/parameter_estimation/joint/fit/standard_errors.rs): the
  module's `expect(dead_code)` is gone; `StandardErrors::described` (the line), `StandardErrors::named`
  (the trace's error rows, zipped with `Parameters::named`), `summary_of`, the reason and name
  constants; `inverse_of_positive_definite` is now test-only, which is all that uses it; a misplaced
  test doc paragraph removed.
- [fit_trace.rs](../../../../src/parameter_estimation/joint/fit_trace.rs): its doc describes the rows
  of the returned fit; a test-only capture of a thread's rows (`captured`).
- [fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs), tests only:
  `fitted` returns the whole fit beside its parameters.

## 4. Tests

| test | what it shows |
|---|---|
| `a_fit_carries_the_errors_at_the_parameters_it_returns` (information.rs) | The errors on the `JointFit` are, bit for bit, those of one more pass at its returned parameters that sums the information — on four samples with the duplicated class on (three starts; the winner is start 2) and on one sample. At four samples every error rate and homozygote excess has an error; at one, the excess says it is held fixed and the two rates keep theirs. |
| `the_trace_files_the_returned_fit_after_the_winning_starts_last_pass` (information.rs) | With the trace captured: the rows after the winner's last pass are under start 2 and pass `passes + 1`, and are exactly the returned values and their errors, under the value rows' names, with the final pass's log-likelihood. Writing them at `passes` fails it (checked). |
| `a_fit_without_the_duplicated_class_says_so` (information.rs) | With the class off, its three errors say no information, and the log line names neither the duplicated share nor a carrier shape. |
| `a_fit_with_no_pass_reports_errors_at_its_start` (information.rs) | With a pass limit below one cycle, the fit runs no pass and still gives every sample's two error rates an error. |
| `the_logged_line_summarises_each_kind` (standard_errors.rs) | The line, exactly: each cohort parameter's value and error or its reason; each sample kind's median (the lower middle of an even count) and largest, and the count missing by reason; values chosen so a mean, the upper middle, or the last value listed in place of the largest print a different line; the duplicated class only when fitted, and one later read group counted as one. |
| `the_trace_names_each_error_after_its_value` (standard_errors.rs) | The trace's error rows carry the value rows' names in their order, with and without the duplicated class; on two samples over three read groups, a sample's second read group gets NaN and the others their sample's errors. |
| `the_last_passes_cost_reads_as_a_ratio` (fit.rs) | The start line's cost: the last pass left out of the average (an average with it in gives a different ratio); tenths of a millisecond below a second; no ratio with no pass, or none measurable, before it. |

The existing `a_census_of_many_chunks_fits_to_the_same_bits_at_any_pool_width` compares the whole
`JointFit` by its debug print at one and four threads, so it now covers the errors too (a review
mutation that changed them only above one thread failed it).

Measured with `cargo test --release --lib parameter_estimation::joint::fit` in the container: 59
passed, 1 ignored (the coverage measurement).

## 5. Validation

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- `cargo test --all-targets --all-features --no-fail-fast` (`tmp/fit_precision/suite_a5b.log`):
  counts in the commit message; the three failures are the pre-existing ones
  (`examples/ng_generic_loci_dump.rs` × 2, `examples/ng_ssr_loci_dump.rs` × 1). **Both
  cross-platform checksum tests pass unchanged.**
- **The identity oracle** (`scripts/promote_ng_oracle.sh`, four tomato accessions at about 3 reads
  over tomato2's first 20 regions; `tmp/fit_precision/oracle_a5/`), run on the reviewed code: every
  checksum equals main's, the two fit lines included (`fitted.parameters.toml`
  `c590689d40eb66359b2973c2bd1cbfcd`, `from_fit.comparable` `5cdeafd0aaa805cfd57b79e3bbb4e3d2`). The
  script reports those two as differing from `scripts/promote_ng_oracle.baseline`, whose fit lines
  predate the fit's acceleration (A1 report §6); main gives the same two values. The review's fixes
  change only what is printed and traced.
- The review ran the final pass with and without the information at 1 and 4 threads, on 1, 4 and 64
  samples: every other field of the pass came out bit-identical.

## 6. On the oracle cohort

Measured on the run above (`tmp/fit_precision/oracle_a5/fit.log`), before the review's fixes changed
the lines' wording; it ran beside the full test suite, so its times are rough:

```
SNP/indel fit, start 1 of 3: stopped at the pass limit after 200 pass(es), log-likelihood -7.991599e5; its last pass, which also sums the information, took 1s against 776ms a pass before it
SNP/indel fit, start 2 of 3: stopped at the pass limit after 200 pass(es), log-likelihood -7.991602e5; its last pass, which also sums the information, took 1s against 710ms a pass before it
SNP/indel fit, start 3 of 3: stopped at the pass limit after 200 pass(es), log-likelihood -7.991633e5; its last pass, which also sums the information, took 988ms against 674ms a pass before it
SNP/indel fit: standard errors at the fitted values: mismapped share 3.62e-4, invariant share 9.90e-5, fixed share 8.15e-6, density shape a none (wider than its range, 5.22e1), density shape b none (wider than its range, 1.31e2), duplicated share 1.98e-4, carrier shape a 7.04e-1, carrier shape b 1.27e0; clean error rates (4): median 1.59e-5, largest 2.17e-5, every one has an error; mismapped error rates (4): median 1.16e-3, largest 1.57e-3, every one has an error; homozygote excesses (4): median 1.00e-2, largest 1.37e-2, every one has an error
```

- **What collecting the information costs.** On the oracle, 1.47 times an ordinary pass on the one
  start whose time the log then gave in milliseconds (988 ms against 674 ms); the other two printed
  "1s", which that form floors from anything between one and two seconds, so their ratios lie
  between 1.29 and 2.8 (the start line now prints the ratio itself). Measured directly in review
  (one thread, median of five interleaved runs), a final pass with the information against one
  without: 1.43 at 4 samples (200,000 positions), 1.72 at 64 (20,000), 1.75 at 64 with the duplicated
  class off; the maximisation and keeping the per-position posteriors each add nothing measurable,
  so the printed ratio is that cost. Three such passes a fit, against 600 ordinary ones on the oracle.
- **What it holds.** At kimura's 2,169 samples and 128 chunks the information adds 573,128 bytes a
  chunk, 73.4 MB, to the final pass's transient hold, which goes from about 174 MB to about 247 MB
  (review; 16 quadrature nodes and one read group a sample assumed). It lasts one pass, three times a
  fit.
- **The density's two shapes have no error**: they come out wider than their whole range. The trace
  shows why: the second shape first reaches its upper bound, 50, at pass 101 and stays there from
  pass 105, and the first ends near 20 (14.75 at pass 100, 18.37 at 150, then 20.13, 19.87 and 19.74
  over passes 198 to 200). The data of four samples do not place them.
- **How far the fit's last passes move each parameter, in its own errors**, from the trace
  (`tmp/fit_precision/oracle_a5_trace.tsv`, winning start 1; `settled_realised.py`): the invariant
  share is 2.01 errors from its pass-200 value at pass 100, 0.51 at pass 150 and 0.13 at pass 190;
  every parameter with an error stays within 1 error of its pass-200 value from pass 123 on, within
  0.5 from pass 152, within 0.2 from pass 184 and within 0.1 only from pass 199. The invariant share
  is the furthest at 198 of the 199 passes before pass 200. Two caveats: the pass-200 values are not
  the likelihood's maximum (the invariant share rose over the run, 0.99740 at pass 25, 0.99833 at
  100, 0.99853 at 200, with small reversals within cycles); and this trace predates the rows of the
  returned values, so whether pass 200's values are the ones returned — they are not when the last
  accelerated step is refused at the pass limit — cannot be told from it.
- **The spec's projected distance cannot be replayed faithfully from the trace.** The trace records
  every plain pass, but not where the fit's accelerated cycles begin and end, so moves taken over
  three passes straddle cycles. Replayed that way (contractions capped at 0.95), it called every
  parameter settled at 0.1 of its error by pass 90, when the invariant share was still 2.11 errors
  from its pass-200 value (`settled_replay.py`). Step B1's comparison against a 1,000-pass fit is
  where the projection is measured.
