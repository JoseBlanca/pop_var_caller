# Fit precision, step B3 — each start of the SNP/indel fit says how it ended, and why when it did not converge

**Date:** 2026-09-29. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step B3.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §3.5.
**Branch:** `fit-precision`. **Review:** [fit_precision_b3_2026-09-29.md](../reviews/fit_precision_b3_2026-09-29.md);
**fixes:** [fixes_applied_fit_precision_b3_2026-09-29.md](../reviews/fixes_applied_fit_precision_b3_2026-09-29.md).
This report describes the step as committed, after the review's fixes.

## 1. What was built

Spec §3.5 asks that a fit which did not converge say which parameter held it, "or a correct number invites
being discarded". **Each start now records how it ended, and when it did not converge, how far from settled
it was at the values it returned.**

- **`JointFit::starts`**, one `StartRecord` a starting point, in order: its number (from one), how it ended
  (`StartEnding`: `Converged`, `AtTheLimit`, or `Agreed { with_start }`), its passes and log-likelihood, and
  `furthest_from_settled` — for a start that did not converge, how many parameters were not yet within a
  tenth of their standard error of where the likelihood peaks, of how many, and the furthest by name with
  its distance in standard errors (`FurthestFromSettled`). Both types are public; the returned fit is the
  winner's record.
- **Names a reader knows** (`Parameters::plain_names`). The cohort's parameters in the words the
  standard-error line already uses ("allele-frequency shape a", "duplicated share"); a sample's excess by the
  sample's name ("homozygote excess of LA1589"); a library's rates by its sample, with its place among the
  sample's libraries where it has several ("error rate at ordinary positions of TS-1's library 2"). The
  per-cycle progress line uses the same names. At kimura's 2,169 samples the trace's keys (`hom_excess_1922`)
  could not be mapped back to a sample from the log.
- **Where the reason comes from.** The information and slope the start's last pass summed, at the values it
  returned — its final pass, or for a start that agreed, the pass that judged it — solved as step B1's
  settled test solves them (`how_far_from_settled`). A start at the pass limit is judged at the point it
  returns, not at its last judged cycle.
- **The start's log line** gains, when it did not converge (`describe_short_of_settled`): "at the values it
  returned 4 of 17 parameter(s) still more than 0.1 standard errors from where the likelihood peaks;
  furthest, the allele-frequency shape a, at 0.23 standard errors". A start that ran out of passes before
  any cycle judged it, with every parameter already within the fraction, says so; one with no parameter to
  judge says that nothing can be said.
- **The per-cycle progress line** already gave the count not settled and the furthest (step B1).
- **The per-pass trace** (`PassSummary::not_settled`, when the run asks for it) records, on each pass that
  judged the parameters, how many were not settled — found by the judging pass's number. This is the hook
  B1's review asked for, to test at the level of the loop that the fit stops only on a pass that found every
  parameter settled.
- `StartEnd` is renamed `StartEnding`, and its `with` field `with_start` (B2 review's naming nits).

Nothing fitted moves: the record and the log are read off passes the fit already runs. A start that did not
converge gains one inversion for the errors and one Newton solve at the values it returned: 5.45 ms a start
on a drawn cohort of kimura's size (2,169 samples, 2,653 read groups; the correctness reviewer's measurement),
where one pass over 2,000 positions took 524 ms.

## 2. What was measured

- **The four-accession oracle cohort** (`tmp/fit_precision/oracle_b3f/`): every checksum
  matches `scripts/promote_ng_oracle.baseline`. Every start still runs to the limit, and now says why: start 1
  (198 passes) "at the values it returned 4 of 20 parameter(s) still more than 0.1 standard errors from where
  the likelihood peaks; furthest, the duplicated share, at 0.39 standard errors"; start 2 (198) 6 of 20, the
  duplicated share at 0.40; start 3 (200) 4 of 20, at 0.34.
- **The cross-platform checksums** pass unchanged.
- On 4 drawn samples at 3 reads where every start runs to the limit (`each_start_is_recorded_as_it_ended`),
  each start reports 4 of 17 parameters not settled, the density's first shape furthest at 0.23, 0.24 and
  0.20 standard errors.

## 3. Tests

| test | what it shows |
|---|---|
| `information::tests::each_start_is_recorded_as_it_ended` | 20 samples: start 1 converged with no reason; starts 2 and 3 agreed with 1, each with one. 4 samples: every start at the limit with its reason. The count not settled is non-zero exactly when the furthest is at least 0.1 errors; the returned fit is the winner's record |
| `information::tests::a_fit_stops_only_on_a_pass_that_found_every_parameter_settled` | the winner's trace on 20 samples whose judged passes find 3, then 1, then 0 unsettled (passes 10, 13, 16), so stopping with one unsettled would fail: the last judged pass found none, each earlier one some; one verdict a cycle, on the pass that judged; none on the first pass. A fit at the limit never found none |
| `information::tests::a_start_at_the_limit_reports_the_settled_test_at_the_values_it_returned` | at a settled fraction of 0.15, the recorded count (3), total (17), furthest parameter and distance equal those recomputed from a fresh pass at the returned parameters, the distance to 10⁻⁹ relative |
| `information::tests::the_log_names_each_parameter_by_its_sample` | the plain names, in the fit's order, for a one-library and a two-library sample |
| `settled::tests::a_start_that_did_not_converge_says_how_far_it_stopped_from_settled` | the log clause: empty when converged; count, name and distance; every parameter within; nothing to judge |

- Fit module and cross-platform (`cargo test --release --lib -- parameter_estimation::joint::fit
  cli::cross_platform_digests`): 113 passed, 1 ignored.
- **Six mutations the review found surviving now fail a test** (`tmp/fit_precision/b3/mut/`): stopping with
  one parameter unsettled; each verdict written on pass 1; the first parameter's name; the distance doubled; a
  fraction of zero; a total five too high.
- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- Full suite (`cargo test --all-targets --all-features --no-fail-fast`): 4,984 passed, 3 failed (the
  pre-existing `examples/ng_generic_loci_dump.rs` × 2 and `examples/ng_ssr_loci_dump.rs` × 1), 5 ignored.

## 4. Deviations and what is left

1. **The trace's per-pass verdict** is not in the plan's step text; B1's review deferred its loop-level
   test to this step, and the verdict is what makes it testable.
2. **Not done:** the parameters file recording each start's outcome (spec §5.2) is step E2's.
