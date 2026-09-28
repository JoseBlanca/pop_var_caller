# Fix Application Report: fit_precision_a6b_2026-09-28.md

**Date:** 2026-09-28
**Source review:** `doc/devel/reports/reviews/fit_precision_a6b_2026-09-28.md`
**Source state reviewed against:** `09f58fd6` (review-only commit on `fit-precision`)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 6
- Minors: 5
- Nits: 1

### Outcome totals
- Applied: 10 (M1–M6, Mi1–Mi4)
- Applied with adaptation: 1 (Mi5)
- Already fixed: 0
- Deferred: 0
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0
- Nit: applied (the report names the log the figure is in).

### Validation summary
- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`,
  `cargo test --all-targets --all-features --no-fail-fast`: recorded in the commit message
  (`tmp/fit_precision/suite_a6b_fixes.log`).
- `cargo test --release --lib parameter_estimation::joint::fit` → 68 passed, 1 ignored; the four
  tests with new `debug_assert`s on their path also pass in a debug build.
- Mutation re-check: six survivors from the review re-run against the fixed tests, all six fail
  (`tmp/fit_precision/a6b_mutations.log`), files restored and compared byte for byte.
- Performance check → not triggered (tests, `debug_assert`s, a doc and one constant named).

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Final status | Validation |
|---|---|---|---|---|
| M1 | Major | `information.rs` double-encoded | Applied | valid UTF-8; no garbled sequence left; diff against the parent 736 → 510 changed lines |
| M2 | Major | whole-fit test cannot tell a sibling's error | Applied | own-error ordering assertion, both rates |
| M3 | Major | no pass runs a three-library sample | Applied | `third_dropped` fails |
| M4 | Major | slots' intervals untested | Applied | `class_swapped`, `excess_as_clean`, `call_site_capped` fail |
| M5 | Major | a sample holding no library untested | Applied | `phantom_group` fails |
| M6 | Major | row-skip rule untested for several libraries | Applied | `any_zero_skip` fails |
| Mi1 | Minor | `TOLERANCE` doc | Applied | — |
| Mi2 | Minor | rows zipped unchecked | Applied | debug-build run |
| Mi3 | Minor | `JointFit::standard_errors` doc | Applied | — |
| Mi4 | Minor | ranged test's coordinate count | Applied | test passes |
| Mi5 | Minor | report's cross-reference and cost | Applied with adaptation | — |

## 3. Questions asked and answers

None.

## 4. Per-finding log

### M1 — the double-encoded file
- **Final status:** Applied. The file was decoded from UTF-8 once and every character at or below
  U+00FF written back as the byte it stood for, the one character above (the minus sign at line
  1644, written after the damage) re-encoded as UTF-8. Checked: valid UTF-8 (`iconv`), no garbled
  sequence left, line 2 reads "a time — the raw material…", line 1644 keeps its "−0.038". The first
  attempt emptied the file (`Encode::decode` with a croak flag clears the string it is handed, and the
  verification step handed it the result); it was restored from the backup taken before it. **Cause
  and guard:** a perl edit whose replacement held a wide character made perl write the whole buffer
  as UTF-8; later edits with non-ASCII text went through the byte-safe block-replacement script, and
  every edited file was checked with `iconv` after.

### M2 — the whole-fit test cannot tell a sibling's error
- **Final status:** Applied with the reviewer's suggestion made specific: within each sample the
  library drawn at the higher clean rate must have the larger clean-rate error, and the library drawn
  at mismapped rate 0.10 the larger mismapped-rate error than the one at 0.04 — at equal depth an
  error grows with its rate, and an error read from the sibling's slot gives the two one error.

### M3 — three libraries
- **Final status:** Applied. `a_sample_of_three_libraries_has_every_slope_right`: one sample of one
  library beside one of three, the duplicated class fitted; every library's two slopes and each
  excess against a central difference, largest disagreement 4.4 × 10⁻⁸.

### M4 — the slots' intervals
- **Final status:** Applied. `each_slot_of_a_samples_row_is_its_own_parameter_with_its_own_bounds`
  (samples of one to four libraries: distinct slots filling the row, each naming its class and
  carrying its class's bounds), and in `the_blocks_give_the_errors_of_an_arrow_of_mixed_libraries` a
  second library's rate rescaled to an error of 0.3: an error as a mismapped rate, wider than its range
  as a clean one — through the inversion, which is where the call-site mutation lives. `class_of` now
  names `ONE_LIBRARY_SAMPLE_PARAMETERS` where it wrote 3.

### M5 — a sample holding no library
- **Final status:** Applied. `a_sample_holding_no_library_carries_no_information`: the likelihood is
  bit-for-bit the one where the sample's library holds no read, its row is zeros, its three errors say
  no information, the log counts three read groups.

### M6 — the row-skip rule
- **Final status:** Applied. `a_pass_sums_the_scores_multiplied_pairwise` also runs on samples of one,
  two and three libraries at four reads a library (a library with no read at a position is common
  there); worst entry within 1.2 × 10⁻¹⁵ of its bound.

### Mi1–Mi5
- **Mi1:** the doc gives 6.8 × 10⁻⁷ as the largest disagreement, and where.
- **Mi2:** `debug_assert`s: a sample's scores against its block's size in `add_position`; two chunks'
  block sizes in `absorb`; `named` and `described` check each row's length against the library map
  (`laid_out_for`).
- **Mi3:** the doc now lists a sample's own parameters once, in their order.
- **Mi4:** the ranged test reads the cohort's coordinate count from `cohort_coordinates`.
- **Mi5:** the report's cross-reference is §6; the cost measured in review is in its §3 item 4, and
  the encoding repair in item 5.

## 5. Deferred findings to carry forward

None.

## 6. Disputed findings to return to reviewer

None.

## 7. Failed-validation findings

None.

## 8. Blocked-by-context-mismatch findings

None.

## 9. Performance check

- **Triggered:** no.

## 10. Commands run
- `scripts/dev.sh bash -c 'cargo fmt && cargo clippy … && cargo test --release --lib parameter_estimation::joint::fit && cargo test --lib -- <the four new tests>'`
- six mutations with restore and `cmp` (`tmp/fit_precision/a6b_mutations.log`)
- `scripts/dev.sh bash -c 'cargo fmt --check; cargo clippy --all-targets --all-features -- -D warnings; cargo test --all-targets --all-features --no-fail-fast'`

## 11. Command results
- module tests → 68 passed, 1 ignored; mutations → 6 of 6 fail a named test; full gate → in the
  commit message.

## 12. Notes

No fitted number and no one-library error changes with these fixes: they add tests, `debug_assert`s,
docs, and one constant named where a literal stood. The oracle and bit-identity measurements of the
review commit stand.
