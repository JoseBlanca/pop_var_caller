# Fix Application Report: fit_precision_a2_2026-09-27.md

**Date:** 2026-09-27
**Source review:** `doc/devel/reports/reviews/fit_precision_a2_2026-09-27.md`
**Source state reviewed against:** `582c8b4c` (review-only commit), branch `fit-precision`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 1
- Minors: 6
- Nits: 14

### Outcome totals
- Applied: 6
- Applied with adaptation: 0
- Already fixed: 0
- Deferred: 1
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0 (one owner item carried to checkpoint A)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean (two needless borrows in the
  new test helper fixed on the way)
- `cargo test --all-targets --all-features --no-fail-fast` → counts in the commit message
  (`tmp/fit_precision/suite_a2b.log`); the three pre-existing failures only; both
  `cli::cross_platform_digests` tests pass unchanged
- `cargo test --release --lib parameter_estimation::joint::fit` → 0, 36 passed
  (`tmp/fit_precision/module_a2c.log`)
- `cargo doc`, `cargo audit` → not run (no public API or dependency change)
- Performance check → not applicable (no pass collects the information yet); the review measured
  the cost of a collecting pass, recorded in the A2 report §6

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | final pass's join untested | Apply | Applied | No | `information.rs` | Pass | No |
| Mi1 | Minor | `absorb` drops a mismatch silently | Apply | Applied | No | `fit.rs` | Pass | No |
| Mi2 | Minor | two positional bools | Apply | Applied | No | `fit.rs`, `information.rs` | Pass | No |
| Mi3 | Minor | "only on a run's last pass" | Apply | Applied | No | `fit.rs` | N/A | No |
| Mi4 | Minor | read-slope index written twice | Apply | Applied | No | `information.rs` | Pass | No |
| Mi5 | Minor | spec and plan describe another block shape | Defer | Deferred | No | None | N/A | Owner, checkpoint A |
| Mi6 | Minor | debug checks never exercised | Apply | Applied | No | `information.rs` | Pass | No |
| Nits | Nit | fourteen items | Apply | Applied (13), not applied (1) | No | `information.rs`, `fit.rs`, A2 report | Pass | No |

## 3. Questions asked and answers

None asked. Mi5 is carried to checkpoint A: the spec's and the plan's description of a sample's block
is the owner's to amend, together with the ruling on later read groups.

## 4. Per-finding log

### M1 — the final pass's join
- **Final status:** Applied. `a_pass_sums_the_scores_multiplied_pairwise` checks the values on both
  joins (worst 1.43 × 10⁻¹⁵ on the tree, 1.37 × 10⁻¹⁵ in position order, of each entry's bound);
  `the_information_is_the_same_bits_at_any_pool_width` `expect`s the information before comparing.
  The reviewers' surviving mutants — the ordered join setting the information to `None`, and keeping
  only its first chunk — now fail the first test (by the `expect`, and by a disagreement of 0.73 of the
  bound, as the reviewers measured on their probes).
- **Files changed:** `information.rs` · **Residual risk:** None.

### Mi1 — `absorb`'s silent drop
- **Final status:** Applied. `Statistics::absorb` matches both sides and is `unreachable!` on a
  mismatch; in the pass, one `Option` (the scoring tables and the row) decides whether a position is
  scored, and the sums are `expect`ed when it is. `InformationSums::absorb` states its equal-length
  precondition in a `debug_assert_eq!`.

### Mi2 — the two bools
- **Final status:** Applied. `PassKeeps { per_position_posteriors, information }` with
  `PassKeeps::SUMS_ONLY` for the iterating passes; the `too_many_arguments` allowance is gone.

### Mi3 — the doc comments
- **Final status:** Applied. `Statistics.information` and `PassKeeps` say the information is kept only
  when a pass asks for it; the report says which steps will ask (A5, B1).

### Mi4 — the table's index
- **Final status:** Applied. `read_slope_index` and `read_slope_table_len`, each written once.

### Mi5 — the design text
- **Final status:** Deferred to the owner at checkpoint A (review open question 1).

### Mi6 — the debug checks
- **Final status:** Applied. `tables_built_for_another_rule_are_refused`, debug builds only, hands the
  scorer tables built at a density shape `a` 0.01 higher than the pass's and expects the panic. It
  passes, and would fail if the check always returned true.

### Nits
- Applied: `describes`; the branch shares built once a pass inside `ScoringTables` (which replaces
  `PassRuleSlopes`); fields `sample_blocks` and `sample_cohort_blocks`; `NoiseClassTerms` bound as
  `class_terms`; `Parameters` imported; `PartialEq` dropped; `absorb`'s length precondition;
  "every finite `f64`"; the transposition claim limited to the sample-with-cohort block; 1, 4 and 8
  threads; the report's review codes removed and its 450 MB tied to kimura's 7,478 parameters; the
  genotype-0 comment says "in exact arithmetic"; the zero-row skip described as a saving in time.
- Also applied, from a cross-category note: the value test runs on a one-sample cohort too (6.40 ×
  10⁻¹⁵ and 6.54 × 10⁻¹⁵ of each entry's bound).
- Not applied: comparing every count through `Statistics`' `Debug` — `Statistics` does not derive
  `Debug`, and deriving it for a test is more than the nit asks; the four fields compared are the ones
  every maximisation reads first.

## 5. Deferred findings to carry forward
- Mi5 — the spec's and the plan's block description (owner).
- Out of scope, recorded: the underflow early return has no fixture.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no — nothing collects the information yet. The review's measurement of a collecting
  pass (1.41× a plain pass at 4 samples, 1.72–1.75× at 64) is recorded in the A2 report for A5.

## 10. Commands run
- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit -- --nocapture`
- `scripts/dev.sh cargo test --lib tables_built_for_another_rule_are_refused`
- `scripts/dev.sh bash -c 'cargo fmt --check && cargo clippy --all-targets --all-features -- -D warnings && cargo test --all-targets --all-features --no-fail-fast'`

## 11. Command results
- fit module → 0, 36 passed
- debug-check test → 0, 1 passed
- fmt, clippy, full suite → see §1 and the commit message

## 12. Notes
- A perl substitution while applying the one-sample change pasted the test into the `RuleSlopes`
  struct (a `|` delimiter turned the pattern's escaped pipes into alternation); it was reverted by
  hand and the non-test code compared line by line against the copy taken before the edits — every
  difference is formatting or the two intended comment changes.
