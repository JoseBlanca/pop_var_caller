# Fix Application Report: fit_precision_b2_2026-09-29.md

**Date:** 2026-09-29
**Source review:** `doc/devel/reports/reviews/fit_precision_b2_2026-09-29.md`
**Source state reviewed against:** `d9f2adab` (a review object of the working tree, parent `ca11a2d8`)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 5
- Minors: 12
- Nits: grouped

### Outcome totals
- Applied: 11 (B1, M1, M4, M5, Mi1, Mi2, Mi6, Mi9, Mi10, Mi11, Mi12)
- Applied with adaptation: 3 (M2, M3, Mi3)
- Superseded: 2 (Mi4 and Mi5, by the new measurement and by M1)
- Deferred: 2 (Mi7 in part, Mi8), and most Nits
- Disputed: 1 Nit (the oracle baseline's commit)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --release --lib parameter_estimation::joint::fit` → 0, 106 passed, 1 ignored
- `cargo test --release --lib cli::cross_platform_digests` → 0, 2 passed, both checksums as recorded at B1
- `cargo test --all-targets --all-features --no-fail-fast` → 101, 4,979 passed, 3 failed (the pre-existing
  `examples/ng_generic_loci_dump.rs` × 2, `examples/ng_ssr_loci_dump.rs` × 1), 5 ignored
- `cargo doc --no-deps`, `cargo audit` → not run (no public API or dependency change)
- Performance check → not applicable (no bench covers the fit's stopping; the step saves passes)

### Unresolved high-priority findings
- None. M3's cost is measured and reported; whether to accept it is the owner's call at checkpoint B.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| B1 | Blocker | spec §3.4's two traps untested | Apply | Applied | `fit/information.rs` | Pass; mutations below |
| M1 | Major | agreement before the settled test | Apply | Applied | `fit.rs` | Pass; fixture back to B1's bytes |
| M2 | Major | wrong mechanism for starts stopping above the answer | Apply | Applied with adaptation | impl report | N/A |
| M3 | Major | the cost of not competing understated | Apply | Applied with adaptation | impl report, `fit.rs` docs and log | Pass |
| M4 | Major | an agreed start's statistics lack posteriors | Apply | Applied | `fit.rs` | Pass |
| M5 | Major | nothing tests `fit_jointly`'s choices | Apply | Applied | `fit/settled.rs`, `fit.rs`, `fit/information.rs` | Pass |
| Mi1 | Minor | one of three lengths asserted | Apply | Applied | `fit/settled.rs` | Pass |
| Mi2 | Minor | vacuous agreement with no errors | Apply | Applied | `fit/settled.rs`, `fit.rs` | Pass |
| Mi3 | Minor | choice 4 against the trap's wording | Apply | Applied with adaptation | impl report | N/A |
| Mi4 | Minor | passes column omits cohorts | Apply | Superseded | impl report | N/A |
| Mi5 | Minor | "fifth significant figure" | Apply | Superseded (by M1) | — | N/A |
| Mi6 | Minor | earlier answer built twice | Apply | Applied | `fit.rs`, `fit/settled.rs` | Pass |
| Mi7 | Minor | three locals kept in step by hand | Apply in part | Applied in part / Deferred | `fit.rs` | Pass |
| Mi8 | Minor | start number as a bare `usize` | Defer | Deferred | None | N/A |
| Mi9 | Minor | module doc | Apply | Applied | `fit/settled.rs` | N/A |
| Mi10 | Minor | nothing at one sample | Apply | Applied | impl report | N/A |
| Mi11 | Minor | `best.expect`'s invariant | Apply | Applied | `fit.rs` | Pass |
| Mi12 | Minor | "1, 4 and 8 threads" | Apply | Applied | impl report | N/A |
| Nits | Nit | names, magic index, NaN | Apply in part | Applied in part / Deferred / Disputed | `fit/information.rs`, `fit/settled.rs` | Pass |

## 3. Questions asked and answers

None during the run. M3's cost goes to the owner at checkpoint B with a recommendation.

## 4. Per-finding log

### B1 — the traps untested
- **Implementation:** `a_start_at_the_pass_limit_is_no_yardstick` runs every start alone on 4 samples × 3 reads,
  where each reaches the limit, and requires `fit_jointly` to return the best of them pass for pass.
  `a_later_start_heading_where_the_first_converged_stops_there` moved to 4 samples × 3 reads, where the second
  start agrees after 20 passes (58 alone) with its value still 1.09 errors from the answer; it now asserts that
  distance exceeds the agreement fraction, so only the endpoint could have agreed.
- **Verification:** the mutation that makes a start at the limit the yardstick fails
  `a_start_at_the_pass_limit_is_no_yardstick` (and `the_yardstick_is_the_best_converged_start`); the one that
  compares the value fails `a_later_start_heading_where_the_first_converged_stops_there`
  (`tmp/fit_precision/b2/mut/`).
- **Review suggestion used verbatim?** The first test adapted from the tests reviewer's code.

### M1 — agreement before the settled test
- **Implementation:** `maximise` asks `agrees` only when the judging pass finds some parameter not settled; a
  settled start finishes its cycle, converges and competes. `a_start_both_settled_and_agreeing_converges`.
- **Verification:** the cross-platform fixture writes B1's bytes again (`22dda760…`, `2c0a8946…`): every start
  converges after 18 passes. The checksum re-record is withdrawn.

### M2 — the mechanism
- **Adaptation:** with M1 the fixture has no agreeing start, so the wrong explanation is removed rather than
  corrected; the non-monotone last cycle is recorded as out of scope in the review.

### M3 — the cost of not competing
- **Adaptation:** the probe re-run on 54 cohorts, nine regimes including one sample and the correctness
  reviewer's 2 × 20 regime; the report's §4.2 tabulates every cohort losing more than 0.01 units (six, the
  largest 0.165, all at one to four samples) with the cause of each, and says why competing would not recover
  them. The doc of `agreement_fraction` and the log line now say "on every parameter the earlier answer gives an
  error". The reviewer's guard (refuse agreement when the stop already scores above the answer) was measured —
  it would apply to 1 of 78 agreeing starts — and not built.

### M4 — posteriors of an agreed start
- **Implementation:** `StartOutcome::statistics` and `last_pass_took` documented for an agreed start; `fit_jointly`
  asserts the winner's statistics come from a final pass (`collect_noisy_posterior`).

### M5 — `fit_jointly`'s choices
- **Implementation:** `settled::is_the_new_yardstick` and `is_the_new_best`, used by the loop, with
  `the_yardstick_is_the_best_converged_start` and `the_answer_is_the_best_start_that_did_not_agree`;
  `later_starts_that_agree_leave_the_first_converged_answer` through `fit_jointly`.

### Mi1, Mi2 — `agrees`
- All four lengths asserted (`a_distance_list_of_another_length_is_refused`); no agreement when no parameter is
  compared (`nothing_agrees_with_an_answer_that_gives_no_error_or_at_a_nan_endpoint`); docs say zero, negative
  or NaN never agrees.

### Mi3 — choice 4
- The report now states the case (start 1 at the limit, start 2 converged, start 3 judged against start 2) and
  quotes both sentences of spec §3.4.

### Mi4, Mi5 — superseded
- The new table gives passes over all six cohorts of each regime, and the totals over the cohorts with an
  agreement separately. With M1 no fixture number moves.

### Mi6, Mi7 — shape
- `Parameters::values`, `EarlierAnswer::of` (its log-likelihood a field). The loop's verdict is one tuple with
  the start being headed to; `ended` and the kept statistics are set together in one place. **Deferred:** a single
  enum for the verdict.

### Mi8 — start numbers
- **Deferred:** B3 publishes the per-start record and is the place to type a start's number.

### Mi9, Mi10, Mi11, Mi12
- Module doc section "Agreement between starts"; one-sample rows in the report; `best.expect`'s invariant named
  in a PANIC-FREE comment; the threads corrected.

### Nits
- `cohort::P_INVARIANT` for the magic index and `pub(super)` fields: applied. A NaN endpoint: tested. Names
  (`StartEnd`, `with`): deferred to B3, which makes the type public. The oracle baseline's commit: disputed —
  it is already its own commit, `ca11a2d8`, which the reviewed range spanned.

## 5. Deferred findings to carry forward
- Mi7 (one verdict enum), Mi8 (start-number type), the names — to B3.

## 6. Disputed findings to return to reviewer
- The oracle baseline Nit — above.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — no bench covers the fit's stopping rule; the step removes passes (47% fewer in the cohorts where a
start agrees) and adds an O(parameters) comparison a judging pass.

## 10. Commands run
- `scripts/dev.sh cargo fmt`, `cargo clippy --all-targets --all-features -- -D warnings`
- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit`, `… cli::cross_platform_digests`
- `scripts/dev.sh cargo test --all-targets --all-features --no-fail-fast`
- the probe (`tmp/fit_precision/b2/probe2.rs`), the oracle (`tmp/fit_precision/b2/run_oracle_b2f.sh`)

## 11. Command results
In the implementation report, [fit_precision_b2_2026-09-29.md](../implementations/fit_precision_b2_2026-09-29.md) §5.

## 12. Notes
- The reviewed commit re-recorded the fit checksum; after M1 it is B1's again, and the re-record is withdrawn.
