# Fix Application Report: snp_fit_convergence_CD_2026-10-17.md

**Date:** 2026-10-17
**Source review:** `doc/devel/reports/reviews/snp_fit_convergence_CD_2026-10-17.md` (committed with these fixes)
**Source state reviewed against:** branch `snp-fit-convergence`, HEAD `0682d6f5`
**Execution mode:** non-interactive (a subagent run; the lead set each finding's treatment)
**Overall status:** Completed

---

## 1. Executive summary

The review found no Blocker and no Major: four Minor findings, two Nits and one spec-level
observation. All six findings are applied; the observation is left as the review recorded it.

- **The errors of `S`, `Hexp` and `E[f]` now use exact slopes along the shares.** The three numbers
  are linear in the invariant share, the fixed non-reference share and the duplicated share, so
  their slopes along those shares are known in closed form; only the shapes are still differenced.
  On a fixed non-reference share at its lower bound of 10⁻¹², the difference had given −0.0867 for
  an exact −0.0533.
- **A trade the cancellation leaves at rounding level is no trade** (`+∞`), with one wording in the
  run log, where one unit in the last place used to give 4.66 × 10¹⁶ and "over a million".
- **A cost below two never prints as two**, and the decimals of a printed cost are chosen after
  rounding.
- **The condition for writing the density's table is tested**, and the review's mutation of it now
  fails a test.
- **Either format-4 table in a file stamped before version 4 is refused**; both specs carry a dated
  amendment.

**What moved.** On the 1-sample oracle the parameters file is byte-identical. On the 4-sample oracle
four standard errors move, by at most 1.1 × 10⁻¹⁰ of themselves (§10), and the run log's verdict on
the density changes wording; no value and no verdict moves. **No digest checksum moves.**

### Review totals
- Blockers: 0
- Majors: 0
- Minors: 4
- Nits: 2

### Outcome totals
- Applied: 5
- Applied with adaptation: 1 (Minor 2)
- Already fixed: 0
- Deferred: 0
- Disputed: 0
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --release --features bench-fixtures --no-fail-fast` → 101; the library 5,207 passed, 3 failed — the three debug-only should-panic tests
  (`left_align_repeated` and `left_align_structured`'s `an_offset_past_the_reference_panics`,
  `lgamma_panics_on_non_positive_in_debug`), as before; every other target passed
- `cargo doc --no-deps` → not run
- `cargo audit` → not run (no dependency changed)
- Performance check → not applicable (§9)

### Unresolved high-priority findings
None.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| Minor 1 | Minor | Central difference wrong near a share's lower bound | Apply | Applied | No (lead's instruction) | `fit/determined.rs`, `parameter_prepass_joint_fit.md` | Pass | No |
| Minor 2 | Minor | No test of the condition for writing the density's table | Apply | Applied with adaptation | No | `run/census_fit.rs` | Pass | No |
| Minor 3 | Minor | "No trade" gets two wordings | Apply | Applied | No | `fit/determined.rs`, `fit/standard_errors.rs`, `parameter_prepass_joint_fit.md` | Pass | No |
| Minor 4 | Minor | A cost just under 2 prints "2.00" | Apply | Applied | No | `fit/standard_errors.rs` | Pass | No |
| Nit 1 | Nit | Format-4 tables accepted in a version-3 file | Apply | Applied | No (lead's instruction) | `parameters_file/validate.rs`, `parameters_file.md`, `parameter_prepass_joint_fit.md` | Pass | No |
| Nit 2 | Nit | The unknown-key test edits both tables at once | Apply | Applied | No | `parameters_file/mod.rs` | Pass | No |
| Spec-level observation | — | One sample's `fixed_non_reference_share` keeps its error | Record | Deferred (left as the review recorded it, by the lead's instruction) | No | None | N/A | Owner's call |

## 3. Questions asked and answers

None. The lead's brief fixed each finding's treatment.

## 4. Per-finding log

### Minor 1 — central difference wrong near a share's lower bound
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** `Group::gradient` took a central difference with a step of 10⁻⁵ of each number's
  value. For `p_fixed_alt` on its lower bound of 10⁻¹² that step, 10⁻¹⁷, is below the last place of
  the segregating share `σ = 1 − p_invariant − p_fixed_alt`, so the difference measured rounding.
- **Implementation summary:** `gradient`, `reported` and `determined_share` take the slopes known
  exactly, `[Option<f64>; N]`; a `None` is differenced as before. The callers supply, from the shapes
  alone: `∂S/∂p_invariant = ∂S/∂p_fixed_alt = −P(f > 0.01)`; `∂Hexp/∂p_invariant = ∂Hexp/∂p_fixed_alt
  = −2ab/((a+b)(a+b+1))`; `∂E[f]/∂p_invariant = −a/(a+b)`, `∂E[f]/∂p_fixed_alt = 1 − a/(a+b)`; and for
  the duplicated class `∂S/∂δ = P(c > 0.01)`. The per-segregating-position heterozygosity and mean
  frequency come from `FrequencyDensity`'s own functions at `p_invariant = p_fixed_alt = 0`, so the
  formulas are not written twice. Module doc and the spec's gradient bullet amended (dated).
- **Review suggestion used verbatim?:** Yes, in substance.
- **Adaptation:** None.
- **Verification performed:** the new test failed before the change with the review's figure
  (`S … 8.673617379884035e-5 against 5.334440746162261e-5`, i.e. −0.0867 against −0.0533) and passes
  after. The three 59-, 4- and 1-sample reference tests still pass at their tolerances (10⁻⁶ on
  errors).
- **Files changed:** `src/parameter_estimation/joint/fit/determined.rs`,
  `doc/devel/ng/spec/parameter_prepass_joint_fit.md`
- **Tests added or modified:** `a_fixed_non_reference_share_on_its_bound_gets_exact_slopes` — the
  review's three measured cases, `(p_invariant, a, b) = (0.84, 0.02, 4.28)`, `(0.9, 0.6, 2.2)`,
  `(0.99, 0.3, 14.5)` at `p_fixed_alt = 10⁻¹²`, with only that share uncertain, so each error is the
  slope times 10⁻³. Exact `P(f > 0.01)` from SciPy (`tmp/review_CD/py/pfa_bound.py`): 0.0533444,
  0.893036, 0.398873; `S` to 10⁻⁹, `Hexp` and `E[f]` to 10⁻¹².
- **Validation:** `cargo test --release --lib determined` → 0, 14 passed.
- **Residual risk:** None.

### Minor 2 — no test of the condition for writing the density's table
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied with adaptation
- **Reasoning:** `parameters_file_of` wrote the density's table only on the `fitted_curve` rung, and
  nothing but the macOS digests covered that or the field-by-field copy.
- **Implementation summary:** the decision and the two copies moved from `parameters_file_of` into
  `the_informational_tables_as_written(rung, &DensityReport, Option<&DuplicatedReport>)`, which
  `parameters_file_of` calls. Two tests on hand-built reports whose every value and error differs.
- **Adaptation:** the review asked for a test "on a `JointFit` with known reports". A `JointFit`
  cannot be built by hand outside its module (a private field), and getting the rung other than
  `fitted_curve` would need a whole `RunParameters`; the test therefore reaches the same decision
  through the function it now lives in, taking the rung and the two reports directly.
- **Verification performed (mutations, each applied with the Edit tool and reverted, diff re-read):**
  the review's `rung == SeedRung::FittedCurve` → `true` fails
  `on_any_other_rung_only_the_duplicated_class_is_written` (`tmp/cd_fix/mut_rung.log`); swapping the
  written `Hexp` and `E[f]` fails `on_the_fitted_curve_every_field_is_its_report_field`
  (`tmp/cd_fix/mut_swap.log`).
- **Files changed:** `src/run/census_fit.rs`
- **Tests added or modified:** `the_informational_tables_in_the_file::on_the_fitted_curve_every_field_is_its_report_field`,
  `…::on_any_other_rung_only_the_duplicated_class_is_written` (also: no duplicated table when the
  class was not fitted).
- **Validation:** `cargo test --release --lib the_informational_tables_in_the_file` → 0, 2 passed.
- **Residual risk:** None.

### Minor 3 — "no trade" gets two wordings
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Reasoning:** where `σ` and `S` move together, `V_trade` is a cancellation that lands on exactly
  zero (`+∞`, "cannot change without changing it") or one unit in the last place above (4.66 × 10¹⁶,
  "over a million"): one state, two wordings.
- **Implementation summary:** `trade_cost_of(share, variance, trade_variance)` returns `+∞` when
  `trade_variance ≤ NO_TRADE_SHARE · variance`, `NO_TRADE_SHARE = 10⁻¹²`. Such a trade, were it real,
  would cost about 2.4 × 10¹¹ · σ²/Var(σ), past two unless σ's error is over 300,000 times σ.
  `a_cost`'s doc reworded: a cost over a million is rounding or a genuine trade far past the line
  (drawn cohorts gave 1.4 and 4.3 million). The spec's trap amended (dated).
- **Review suggestion used verbatim?:** Yes.
- **Verification performed:** on the 4-sample oracle the density's verdict now reads "the
  segregating share cannot change without changing it" where it read "costs over a million
  log-likelihood units"; the 4-sample reference test now asserts the cost is `+∞`.
- **Files changed:** `src/parameter_estimation/joint/fit/determined.rs`,
  `src/parameter_estimation/joint/fit/standard_errors.rs`,
  `doc/devel/ng/spec/parameter_prepass_joint_fit.md`
- **Tests added or modified:** `a_trade_variance_at_or_below_zero_is_no_trade` renamed
  `a_trade_variance_at_rounding_level_is_no_trade` (zero, negative, one ulp, exactly 10⁻¹² of the
  variance → `+∞`; the next number up and an ordinary value → the formula);
  `a_share_that_moves_only_with_what_is_seen_has_no_trade` and the 4-sample test tightened from
  "cost > 10¹²" / "> 10⁶" to `+∞`.
- **Validation:** `cargo test --release --lib parameter_estimation::joint::fit::` → 0, 156 passed.
- **Residual risk:** None.

### Minor 4 — a cost just under 2 prints "2.00"
- **Severity:** Minor
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** `in_three_figures` rounds to three significant figures first and
  chooses the decimals (and scientific notation from a thousand) from the rounded value: 0.9996 →
  `1.00`, 9.996 → `10.0`, 999.6 → `1.00e3`. `a_cost` prints a cost below
  `SEPARATION_LOG_LIKELIHOOD` with more decimals until it shows below two: 1.996 → `1.996`, 1.99996 →
  `1.99996`; the largest number below two prints below two.
- **Review suggestion used verbatim?:** Yes.
- **Files changed:** `src/parameter_estimation/joint/fit/standard_errors.rs`
- **Tests added or modified:** `a_cost_that_rounds_up_a_power_of_ten_keeps_three_figures`,
  `a_cost_below_the_line_never_prints_as_reaching_it`. Both fail on the old code (0.9996 printed
  `1.000`; 1.996 printed `2.00`).
- **Validation:** as Minor 3.
- **Residual risk:** None.

### Nit 1 — format-4 tables accepted in a version-3 file
- **Severity:** Nit
- **Initial decision:** Apply
- **Final status:** Applied
- **Implementation summary:** `the_informational_tables_hold_numbers` refuses either table in a file
  whose `format_version` is below `FORMAT_THAT_ADDED_THE_INFORMATIONAL_TABLES` (4), naming the
  table, worded as `from_toml.rs`'s refusal of a format-3 key in a file stamped 1 or 2. A version-1
  or -2 file holding a table is migrated to version 3 and refused the same way. The other two cases
  the review listed (`fitted_density` on a rung other than `fitted_curve`; a pair whose halves
  disagree on `separately_determined`) are left accepted, as the lead's brief named only this one.
  Amended with dated lines: `parameters_file.md` §3.6 and `parameter_prepass_joint_fit.md`'s
  validation bullet, which said "nothing else is checked".
- **Files changed:** `src/calling/parameters_file/validate.rs`,
  `doc/devel/ng/spec/parameters_file.md`, `doc/devel/ng/spec/parameter_prepass_joint_fit.md`
- **Tests added or modified:** `a_format_four_table_in_a_version_three_file_is_refused` (each table
  alone refused, naming it; version 3 with neither and version 4 with both accepted).
  `a_version_this_build_does_not_read_is_refused` and `a_version_one_file_is_accepted` now take both
  tables out of the fixture before stamping it 3 or 1 — both failed until they did, which is the
  check working on files that never had the tables.
- **Validation:** `cargo test --release --lib calling::parameters_file` → 0, 290 passed.
- **Residual risk:** None.

### Nit 2 — the unknown-key test edits both tables at once
- **Severity:** Nit
- **Final status:** Applied
- **Implementation summary:** the test inserts `later_key = 1.0` under one table header at a time
  (asserting each header occurs once in the written text) and asserts the refusal names `later_key`.
- **Files changed:** `src/calling/parameters_file/mod.rs`
- **Validation:** as Nit 1.

### Spec-level observation — one sample's `fixed_non_reference_share` keeps its error
- **Final status:** Deferred — left as the review recorded it, by the lead's instruction. At one
  sample the error (2.08 × 10⁻⁵) is computed with the unidentified shapes held fixed; the code does
  what the spec says. A spec decision for the owner if it is to change.

## 5. Deferred findings to carry forward
- Spec-level observation — whether one sample's `fixed_non_reference_share` should keep its error.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no. The changed code runs once per fit, after the winning start's final pass (the
  reports) or once per file (writing and validating); no bench in `benches/` reaches it.

## 10. Identity

- **Oracle:** the release `pvc` of the final tree (`tmp/cd_fix/pvc`, run by `tmp/cd_fix/run.sh` at 8
  threads) on the 1 and 4 oracle psps, reference `S_lycopersicum_chromosomes.4.00.fa`, the tomato1
  catalog, against `tmp/checkpoint_D/{s1,o4}.toml`:
  - **1 sample: byte-identical.** At one sample the derived numbers have no error to move.
  - **4 samples: four standard errors move, nothing else.**

    | field | before | after | moved by |
    |---|---|---|---|
    | `share_segregating_above_stated_frequency.standard_error` | 3.2707211572639323e-4 | 3.270721157265594e-4 | 5.1 × 10⁻¹³ of itself |
    | `expected_heterozygosity.standard_error` | 1.0858529460580352e-4 | 1.085852946074162e-4 | 1.5 × 10⁻¹¹ |
    | `expected_alternative_frequency.standard_error` | 8.361351901696355e-5 | 8.361351901778871e-5 | 9.9 × 10⁻¹² |
    | `share_with_carriers_above_stated_frequency.standard_error` | 1.5948924139317435e-4 | 1.5948924141061228e-4 | 1.1 × 10⁻¹⁰ |

    All from Minor 1 (exact slopes in place of differences of a share at 1.7 × 10⁻³ and 0.014).
    No value moves and neither verdict: the density's pair separately determined, the duplicated
    pair not.
  - **Run log:** the 4-sample density verdict's wording (Minor 3): "halving or doubling the
    segregating share with it held costs over a million log-likelihood units" → "the segregating
    share cannot change without changing it". The duplicated pair's printed cost, 0.00550, does not
    move; nothing else in either log moves but timings and the output path.
- **Digests:** none moves, so none is re-recorded. All five `cli::cross_platform_digests` tests pass
  with their recorded checksums — the fitted files' and the calls' alike — in the full suite. The
  fixture cohorts' written errors evidently do not change in any printed bit; the VCF checksums
  could not have moved in any case, since calling reads neither table.

## 11. Commands run
- `cargo test --release --lib a_fixed_non_reference_share_on_its_bound` (before the fix: failed as
  expected)
- `cargo test --release --lib determined`
- `cargo test --release --lib parameter_estimation::joint::fit::`
- `cargo test --release --lib the_informational_tables_in_the_file` (and under two mutations)
- `cargo test --release --lib calling::parameters_file`
- `cargo fmt`, `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo build --release --bin pvc`; `tmp/cd_fix/run.sh` for `s1` and `o4`
- `cargo test --release --features bench-fixtures --no-fail-fast`

## 12. Notes
- `PROJECT_STATUS.md` has no block for this feature, and the branch's two earlier fix runs did not
  add one; left unchanged, as they did.
- Native `cargo` was used (the light dev loop), so builds are under `target/`, not
  `target-container/`.
