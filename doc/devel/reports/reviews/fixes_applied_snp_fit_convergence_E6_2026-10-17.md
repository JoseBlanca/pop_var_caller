# Fix Application Report: snp_fit_convergence_E6_2026-10-17.md

**Date:** 2026-10-17
**Source review:** `doc/devel/reports/reviews/snp_fit_convergence_E6_2026-10-17.md` (committed as 6d43ab97)
**Source state reviewed against:** branch `snp-fit-convergence`, HEAD `b73be771`
**Execution mode:** non-interactive (a subagent run; the lead set each finding's treatment)
**Overall status:** Completed for M1, Mi1–Mi5 and the nits. M2, the merge with main, is done in its
own commit after this one and reported there and in the merge's commit message.

---

## 1. Executive summary

- **M1 — a pair is marked whenever a number of its group has no standard error, and a derived
  number's error wider than 1 is withheld.** On the reviewer's two-sample cohort the report had
  said the share segregating above 1% was 0.073 ± 11.0, and called the density's pair separately
  determined at a trade cost of `+∞`, while every number of the density's group had an error wider
  than its range. Now both pairs there are marked, with no trade cost, and `S`, `Hexp` and `E[f]`
  have no error. Spec §2.1.3 amended (dated).
- **Mi1 — not fixed; recorded as a known limit** in `fit_precision.md` §2 and the plan's
  follow-ups. The few-lines option (fix every parameter within a tolerance of the shortest room)
  does not help the probe's case, whose crossings are spread over the whole step.
- **Mi2 — a whole-fit test** that fails with the bound fix switched off (8 of 10 excesses land on
  exactly 1) and passes with it.
- **Mi3–Mi5 and the nits** — docs, two spec passages, one constant moved, the whole Newton step
  taken as the target itself, two test docs' counts, two link targets.

**What moved, against checkpoint D (`tmp/checkpoint_D/`) and the E3 re-run (`tmp/f/runs/`):**

| input | parameters file | verdicts |
|---|---|---|
| 1 oracle sample | byte-identical to checkpoint D | both pairs marked, as before |
| 4 oracle samples | the density's pair now marked: `invariant_share` and `segregating_shape_a` lose their errors (were 3.36 × 10⁻⁴ and 36.8) and say `separately_determined = false`. Nothing else moves against the C–D fixes' run (`tmp/cd_fix/runs/o4.toml`) | **density: separately determined → marked**, because its shape `b` sits on its bound of 50 with an error of 115 against a range of 49.98. Duplicated pair: still marked; its cost of 0.0055 is no longer computed, since its share has no error |
| 59 tomato2 samples | byte-identical to the E3 re-run with the bound fix (`tmp/f/runs/t59.toml`) | both pairs marked, as at checkpoint D |

The 4-sample change is the lead's conservative rule doing what it says; the run log now reads
"what halving or doubling the segregating share costs cannot be computed, as the data do not
determine every number it is computed from".

### Review totals
- Blockers: 0
- Majors: 2 (M1, M2)
- Minors: 5
- Nits: 3

### Outcome totals
- Applied: 7 (M1, Mi2, Mi3, Mi4, Mi5, nit 2, nit 3)
- Applied with adaptation: 1 (nit 1)
- Deferred: 1 (Mi1, recorded as a known limit, on the lead's decision)
- Handled in the following merge commit: 1 (M2)
- Already fixed, Disputed, Failed validation, Blocked, Superseded, Awaiting answer: 0

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean (after replacing an `f64::powi` the crate's lint refuses in a new test)
- `cargo test --release --features bench-fixtures --no-fail-fast` → 101: only the three debug-only should-panic tests fail (§11)
- `cargo doc --no-deps` → not run; the two redundant link targets the review named are removed
- `cargo audit` → not run (no dependency changed)
- Performance check → not applicable (§9)

### Unresolved high-priority findings
- None. M2 is resolved in the merge commit that follows.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | Impossible errors and a misleading verdict at two samples | Apply (lead's rule) | Applied | Lead's decision | `fit/determined.rs`, `fit.rs` (test), spec §2.1.3 | Pass | No |
| M2 | Major | The branch does not merge cleanly into main | Apply | Handled in the merge commit | Lead's instructions | — | — | Merge commit |
| Mi1 | Minor | Box solve cost grows with the square of the samples | Defer (lead) | Deferred | Lead's decision | `fit_precision.md` §2, plan | N/A | Yes |
| Mi2 | Minor | No whole-fit test for the E3 defect | Apply | Applied | No | `fit.rs` (test) | Pass; fails under the mutation | No |
| Mi3 | Minor | Docs describing removed behaviour | Apply | Applied | No | `standard_errors.rs`, `settled.rs`, `information.rs` | Pass | No |
| Mi4 | Minor | Specs lag the code | Apply | Applied | No | spec §3.3.1, §2.1.3 | N/A (docs) | No |
| Mi5 | Minor | A constant between `refuse`'s doc comment and the function | Apply | Applied | No | `validate.rs` | Pass | No |
| N1 | Nit | A whole Newton try can land one ulp outside a bound | Apply | Applied with adaptation | No | `fit.rs` | Pass | No |
| N2 | Nit | Two test docs misstate counts | Apply | Applied | No | `fit/information.rs` | Pass | No |
| N3 | Nit | Two redundant explicit link targets | Apply | Applied | No | `determined.rs`, `standard_errors.rs` | Pass | No |

## 3. Questions asked and answers

None asked; the lead decided M1's rule (mark whenever any number of the group has no standard
error), Mi1's treatment, and the rest before the run.

## 4. Per-finding log

### M1 — impossible errors and a misleading verdict at two samples
- **Implementation:** `determined::error_from_variance` withholds an error wider than 1, the
  range of every number derived there (`S`, `Hexp`, `E[f]`, the duplicated class's `S`). `Group`
  records whether every number of the group has a standard error; when one has none,
  `determined_share` computes no trade cost, so the pair is marked. *Choice made here:* the cost is
  not computed in that case (`None`), rather than computed and overridden, because a covariance
  whose diagonal exceeds the range is not a usable local picture, and the run log's existing
  wording for `None` ("cannot be computed, as the data do not determine every number it is computed
  from") says exactly that. `DeterminedShare::trade_cost`'s doc and the module doc follow.
- **Range for `Hexp`:** taken as [0, 1] as the lead decided, although `Hexp` cannot exceed 0.5.
- **Spec:** §2.1.3's formula block ("no standard error" in place of "no covariance"), a new
  dated bullet "A number with no standard error", the derived errors' bullet, the `DeterminedShare`
  sketch, and a dated note under the measured table saying which cells change (the 4-sample row)
  and which rows were not re-run (the two read in part).
- **Tests:** `a_pair_whose_numbers_have_no_error_is_marked_and_no_derived_error_exceeds_its_range`
  (fit.rs, the reviewer's fixture: 2 samples, 3 reads, 10,000 positions, duplicated class on, seed
  77) — failed before the fix with `share_above` 0.0731 ± 11.03 and `separately_determined = true`,
  passes after; `a_number_with_a_covariance_but_no_error_marks_the_pair` and
  `a_derived_error_wider_than_its_range_is_withheld` (determined.rs). Changed:
  `at_four_oracle_samples_…` renamed to `…both_pairs_are_marked_as_a_number_of_each_has_no_error`,
  which also checks that the trades themselves are as before when the missing errors are supplied
  (`+∞` for the density, 0.0055 for the duplicated class); `a_share_that_moves_only_with_what_is_seen_has_no_trade`
  now gives the carrier shapes a variance of 10⁻³⁰ instead of 0, so that they have errors.
- **Verification:** the 1-, 4- and 59-sample runs above (`tmp/e6_fix/run.sh`, release `pvc` of
  this tree).

### M2 — the merge with main
Done in the next commit; see that commit's message.

### Mi1 — box solve cost
- **Deferred on the lead's decision.** The few-lines option was considered and not built: the
  probe's excesses cross 0 at fractions spread over the whole step (`0.5 + draw` slopes), so
  grouping crossings within a tolerance of the shortest would still take about one round each.
- Recorded in `fit_precision.md` §2 ("A known limit") with the measurements (10 ms a solve at 59
  samples, 5.4 to 5.9 s at 2,000 in the contrived case; 1.3 to 1.9 ms at 2,000 with no parameter
  near an end) and the three options; added to the plan's follow-ups.

### Mi2 — whole-fit test of the E3 defect
- `fit_jointly_returns_a_coefficient_for_every_sample_when_excesses_sit_at_one`: drawn cohorts of
  10 (seed 5) and 20 (seed 7) samples, 3 reads, 20,000 positions, excess drawn at 1.0; asserts no
  excess at 1 and a coefficient for every sample.
- **Mutation:** `Coordinate::newton_range` made to return the whole interval — the test fails with
  8 of 10 samples on 1 (`tmp/e6_fix/mi2_mutation.log`); reverted.

### Mi3 — docs describing removed behaviour
- `standard_errors.rs` module doc: stops on the decrement, the per-parameter distances reported only.
- `settled.rs`: the bound `√(2λ)` covers parameters that have an error; those without are profiled
  out (plan step B4).
- `information.rs`: the ignored coverage measurement's floors are marked as measured under the
  B1 rule and not re-measured.

### Mi4 — specs lag the code
- `parameter_prepass_joint_fit.md` §3.3.1 step 1: dated note on the golden-section numbers' target
  stopping half the search's width inside each end, and the whole step being `θ*` itself.
- §2.1.3: `V_trade ≤ 10⁻¹² · Var(σ)` in the formula block.

### Mi5 — `validate.rs`
- `FORMAT_THAT_ADDED_THE_INFORMATIONAL_TABLES` moved above `refuse`'s doc comment.

### N1 — whole Newton try one ulp outside a bound
- **Adaptation:** of the review's two options, the one the lead chose: at fraction 1 the try is the
  target itself; other fractions are as before. Test `a_whole_newton_step_lands_exactly_on_its_target`
  (200 starts of `a` stepped to its bound 0.02; the fixture asserts at least one start whose
  `a + 1·(0.02 − a)` rounds below 0.02). No output moved: the 1- and 59-sample files are
  byte-identical to before, and the 4-sample file moves only by M1.

### N2 — test docs' counts
- The accepting fixture: eight steps accepted (a quarter five times, a half, the whole twice),
  converged at its ninth judging pass, pass 25. The refusing fixture: stops at the 200-pass limit
  with its outcome counting 199. Both re-measured on this tree with a temporary probe (removed).

### N3 — redundant link targets
- `determined.rs` and `standard_errors.rs`: `[`cohort`](super::information::cohort)` → `[`cohort`]`.

## 5. Deferred findings to carry forward
- Mi1 — the box solve's cost when many parameters cross an end at different fractions; measure a
  judging pass's share of wall time on a real cohort of thousands first.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no. M1 runs once per fit, after the final pass; N1 removes one multiply-add from
  the whole step's try. No bench in `benches/` reaches either.

## 10. Commands run
- `cargo test --release --lib a_pair_whose_numbers_have_no_error` (before the fix: failed as expected)
- `cargo test --release --lib parameter_estimation::joint` → 478 passed, 1 ignored
- `cargo test --release --lib fit_jointly_returns_a_coefficient_for_every_sample` under the mutation → failed as expected
- `cargo fmt`, `cargo fmt --check`
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo build --release --bin pvc`; `tmp/e6_fix/run.sh` for `s1`, `o4` (8 threads) and `t59` (18 threads, alone)
- `cargo test --release --features bench-fixtures --no-fail-fast`

## 11. Command results
- Full suite → 101; the library 5,213 passed, 3 failed — the three debug-only should-panic tests
  (`left_align_repeated` and `left_align_structured`'s `an_offset_past_the_reference_panics`,
  `lgamma_panics_on_non_positive_in_debug`), 7 ignored; every other target passed. No
  `cross_platform_digests` checksum moved, so none is re-recorded here.
- `tmp/e6_fix/run.sh`: s1 107 s, o4 258 s, t59 1,508 s, all exit 0.

## 12. Notes
- `PROJECT_STATUS.md` has no block for this feature, and the earlier fix runs of this branch did
  not add one; left unchanged.
