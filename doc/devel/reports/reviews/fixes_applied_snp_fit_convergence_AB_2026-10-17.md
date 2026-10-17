# Fix Application Report: snp_fit_convergence_AB_2026-10-16.md

**Date:** 2026-10-17
**Source review:** `doc/devel/reports/reviews/snp_fit_convergence_AB_2026-10-16.md` (with its addendum of 2026-10-17)
**Source state reviewed against:** `snp-fit-convergence` at 93f22f6c
**Execution mode:** non-interactive (the owner's decisions on the three open questions came with the task)
**Overall status:** Completed

---

## 1. Executive summary

The review found why every start of the cross-platform checks' small cohorts — two samples without
sequencing errors, two with, and 21 — ran to the 200-pass limit after the decrement stopping rule was
built. Two defects stacked (the box solve that failed to finish, and an error-rate swap that made the
alternation lose tens of log-likelihood units every other cycle), and a third cause was a design
question the owner settled: the decrement now counts only the parameters that have a standard error.
**With the three fixes every start of all three fixtures converges** (16, 16, 16; 22, 13, 13; 84, 38,
44 passes), and on the 1, 4 and 59 tomato samples nothing about the fit's outcome changed: the 1- and
4-sample parameters files match the review's prototype run byte for byte, and the 59 samples take the
same passes to the same log-likelihood.

### Review totals
- Blockers: 2 (B1, B2)
- Majors: 2 (M1, M2)
- Minors: 2 (Mi1, Mi2)
- Nits: 2 (the threshold's doc; the self-swap)

### Outcome totals
- Applied: 6 (B1, B2, M2, Mi1, nit 1, nit 2)
- Applied with adaptation: 2 (M1, Mi2)
- Already fixed, Deferred, Disputed, Failed validation, Blocked, Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → 0, clean, at every commit
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean, at every commit
- `cargo test --release --lib parameter_estimation::joint` → 0 at every code commit; 441, 442, 444,
  450 passed, 1 ignored
- `cargo test --release --lib cli::cross_platform_digests` → 0 at every code commit, re-recorded
- `cargo test --release --features bench-fixtures --no-fail-fast` (full suite, at the last commit) →
  lib 5,173 passed, 3 failed (the three debug-only panic tests), 7 ignored; all other targets pass
- `cargo doc --no-deps`, `cargo audit` → not run: no public doc or dependency change
- Performance check → not run with criterion (§9); timings measured directly instead

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| B1 | Blocker | box solve fails on small cohorts, returns a step outside the box | Ask (replace or fall back) | Applied | Yes: replace, no fallback | `fit/standard_errors.rs`, `fit/settled.rs` (docs) | Pass | No |
| B2 | Blocker | error-rate swap makes the alternation lose up to 909 units | Ask (keep the order?) | Applied | Yes: delete; note if inverted | `fit.rs` | Pass | No |
| M1 | Major | decrement counts directions the data do not determine | Ask | Applied with adaptation | Yes: profile them out, via the arrow | `fit/standard_errors.rs`, `fit.rs`, `fit/information.rs` | Pass | No |
| M2 | Major | four of five control-flow mutations survive | Apply | Applied | No | `fit.rs`, `fit/information.rs` | Pass | No |
| Mi1 | Minor | progress line prints the size, not the sign | Apply | Applied | No | `fit.rs` | Pass | No |
| Mi2 | Minor | an accepted try's pass left unjudged at the limit | Apply | Applied with adaptation | No | `fit.rs` | Pass | No |
| Nit 1 | Nit | `decrement_threshold` doesn't say what zero does | Apply | Applied | No | `fit.rs` | Pass | No |
| Nit 2 | Nit | `swap(group, group)` | Apply | Applied (deleted with B2) | No | `fit.rs` | Pass | No |

## 3. Questions asked and answers

1. **M1** — Should the decrement count only directions the data determine?
   - **Answer (owner, 2026-10-17):** yes; only parameters with a standard error enter it, the others
     adjust freely (profiled); it must equal today's decrement when every parameter has an error,
     and scale to thousands of samples through the arrow, not a dense inverse.
2. **B2** — Must a library's mismapped rate stay above its ordinary rate?
   - **Answer:** no; the swap is deleted, the 0.5 cap on the mismapped share defines the class;
     verify nothing outside the fit reads the mismapped rate; a note only if a library ends inverted.
3. **B1** — Replace the box solve, or fall back to the new one only on failure?
   - **Answer:** replace, no fallback.

All three are recorded in `fit_precision.md` §2 (amendment of 2026-10-17) and
`parameter_prepass_joint_fit.md` §3.3 and §11 question 13 parts 5–7 (commit 06f57ef4).

## 4. Per-finding log

### B1 — the box solve
- **Final status:** Applied (commit ebfebf7b, plan step B2).
- **Implementation:** `bounded_step` is the primal active-set method: from `d = 0`, solve the free
  parameters' Newton step from the current point, move to the first end met and fix every parameter
  that meets an end at that point, free after a whole step the one fixed parameter whose slope points
  furthest back in (above 10⁻⁹ of its terms, `SLOPE_POINTING_IN`), stop when none does. It returns
  `BoundedStep { moves, finished }`; after `10 × side + 10` rounds it gives up and `newton_step`
  reports the decrement as NaN. Docs on `NewtonStep`, `has_converged` and `settled.rs` rewritten.
- **Adaptation from the prototype:** every parameter meeting an end at the same step length is fixed
  in one round (the prototype fixed one), so a cohort whose homozygote excesses sit on 0 does not
  spend a round per sample — measured: 1,061 parameters on an end at 2,000 samples, 3.0 ms a solve.
  The 4-sample parameters file is nonetheless byte-identical to the prototype's run.
- **Test:** `the_box_step_stays_inside_the_box_and_is_worth_at_least_nothing` (400 random cases, whole
  matrix and arrow, condition numbers to about 10¹², positions on, beside and between the ends; inside
  the box, worth ≥ 0, finished, and the KKT conditions on well-conditioned cases). It fails on the old
  solve at its second case (checked by putting the old solve back in the test).
- **What moved:** 1 oracle sample, file byte-identical; 4 samples, byte-identical to the review's
  `o4_primal` run — against B1, same passes, multipliers within 4 × 10⁻¹², inbreeding coefficients
  within 5 × 10⁻¹⁰, concentrations within 3 × 10⁻¹¹, repeat-tract rows within 1.9 × 10⁻⁵ of
  themselves. Digests: two samples without errors converge (37, 41, 41); the other two still at the
  limit.

### B2 — the error-rate swap
- **Final status:** Applied (commit b043ee0c, plan step B3).
- **Verification that nothing reads the mismapped rate:** calling reads `JointFit::sequencing_error_rates`
  (ordinary rate only); contamination reads the per-position mismapped posterior; `parameters_file/`
  contains no mismapped rate; the run log summarises its standard errors; the fit trace prints it
  (diagnostic); `census_file.rs` compares it only in a test.
- **The note:** a **run-log** line after the standard errors (`note_on_inverted_rates`), printed only
  when a library's mismapped rate ends below its ordinary one; **no parameters-file change**. It never
  printed on the 1, 4 or 59 tomato samples; it prints on the digest fixtures with errors (2 of 3
  libraries; 12 of 22).
- **Tests:** `plain_passes_never_lower_the_log_likelihood` gains three small cohorts where the swap
  fires (239 falling plain passes with it, up to 0.45 units; none without);
  `the_log_names_the_libraries_whose_mismapped_rate_ended_below_the_ordinary_one`.
- **Also:** the comment that `θ₂` "never needs a check" now states the guarantee and its limits.
- **What moved:** oracle nothing; digests: with-errors starts end within 0.001 of each other, 21
  samples −5,053.47 → −5,045.64.

### M1 — the decrement over determined parameters
- **Final status:** Applied with adaptation (commit c0929246, plan step B4).
- **Implementation:** `determined_decrement` computes
  `λ_K = min over e_U ½ (d_K, e_U)ᵀ I (d_K, e_U) = ½ d_Kᵀ (I_KK − I_KU I_UU⁻¹ I_UK) d_K`
  with one `arrow_solve`/`whole_matrix_solve` over `U` and two products with the information;
  `newton_step` (now taking the pass's `StandardErrors`) uses it when some parameter the step moves
  lacks an error, and the old `g·d − ½ dᵀ I d` otherwise, bit for bit.
- **Adaptation:** the prototype inverted the dense matrix; this goes through the arrow's block
  elimination, as the owner required. Cost on drawn arrows: **7 µs a judging pass at 59 samples, 160
  µs at 2,000** (box solve 0.12 and 3.0 ms; standard errors 0.08 and 1.8 ms). On the oracle and the
  59 tomato2 samples every parameter the step moves has an error, so the new path is not taken there.
- **Tests:** against a dense `½ d_Kᵀ (Σ_KK)⁻¹ d_K` on three matrices and three choices of profiled
  parameters (holding them still instead fails it — checked); bit-identity with the old decrement
  when every parameter has an error; exactly zero when none has.
- **What moved:** oracle nothing; digests: every start converges (16 ×3; 22, 13, 13; 84, 38, 44).

### M2 — mutation coverage
- **Final status:** Applied (commit of step B5).
- **`PassSummary::newton`** (`NewtonTrial`: first fraction, fraction kept, last pass the tries ran).
- **Each mutation, run against `parameter_estimation::joint` and reverted** (logs `tmp/ab_fix/mut_*.log`):

  | mutation | caught by |
  |---|---|
  | negative decrement counted converged, in `has_converged` | `a_start_has_converged_below_the_threshold_and_never_at_a_negative_or_nan_decrement` |
  | the same at the call site in `maximise` | **nothing — an equivalent mutant now**: the box solve cannot give a negative decrement (NaN fails `< 0` too); `every_judged_decrement_is_at_least_zero` pins that on 28 judging passes |
  | the "end at the judged point" guard removed | `refused_tries_that_leave_fewer_than_two_passes_end_the_start_at_the_judged_point` |
  | an accepted try's pass re-run | `an_accepted_newton_step_costs_one_pass`, `a_start_at_the_pass_limit_is_no_yardstick` |
  | first try always the whole step | `the_first_newton_try_after_an_accepted_quarter_is_half_the_step`, and two others |
  | a start at the limit used as yardstick | `a_start_at_the_pass_limit_is_no_yardstick` (now at 23 passes, with a check that start 2 would agree with an answer at the limit), `the_yardstick_is_the_best_converged_start` |
  | Mi2 reverted (old loop condition) | `an_accepted_trys_pass_is_judged_at_the_pass_limit` |

### Mi1 — the progress line's sign
- **Final status:** Applied. `how_the_log_likelihood_moved` prints `rose (+2.28e0)`, `fell (-5.00e1)`
  or `did not move`; unit-tested.

### Mi2 — an accepted try's pass unjudged at the limit
- **Final status:** Applied with adaptation. Rather than `passes + 2 <= max_passes` for such a cycle,
  the loop always enters with an accepted try's pass in hand (judging it costs no pass), and the
  guard before the SQUAREM cycle (two passes) now applies to every cycle, not only after refused
  tries — equivalent on every other path, since the loop's three-pass condition already left two.
  Test: `an_accepted_trys_pass_is_judged_at_the_pass_limit`.

### Nits
- `JointFitConfig::decrement_threshold` now says that zero or below never converges.
- The no-op `swap(group, group)` went with B2's loop.

## 5–8. Deferred, disputed, failed, blocked

None.

## 9. Performance check

- **Triggered:** the box solve and the decrement run every judging pass; they are not reachable from
  `benches/` as a separate bench, so criterion was not run.
- **Measured instead:** the micro-timings in M1 above, and whole runs: 1 oracle sample 75 s, 4 oracle
  samples 253–261 s at 4 threads (B1's final runs: 61 s and 239 s, recorded by the implementer on
  another day); **59 tomato2 samples, full, alone, 18 threads: 942 s** against B1's 993 s (the
  review's run), SNP/indel stage 15 min 20 s against 16 min 10 s, peak memory 414 MB against 454 MB.
  Same passes as B1, so the difference is run-to-run variation, not the fix.

## 10. Commands run

- `cargo fmt --check`; `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --release --lib parameter_estimation::joint`
- `cargo test --release --lib cli::cross_platform_digests` (and with `--test-threads 1` to split the
  failing checksums' files out: `tmp/ab_fix/digests.sh`, `split.py`)
- `tmp/ab_fix/run.sh <binary> <name> {s1|o4|t59} [threads]` (estimate-parameters on the oracle psps
  and tomato2's 59)
- `cargo test --release --features bench-fixtures --no-fail-fast`

## 11. Command results

- Full suite at the last commit, `cargo test --release --features bench-fixtures --no-fail-fast` →
  exit 101: lib 5,173 passed, **3 failed — the three debug-only panic tests**
  (`left_align_repeated`/`left_align_structured` `an_offset_past_the_reference_panics`,
  `genetics::lgamma_panics_on_non_positive_in_debug`), 7 ignored; every other target passed
  (`tmp/ab_fix/full_suite.log`).

## 12. Notes

- `PROJECT_STATUS.md` was not updated: it has no block for this feature, and the lead owns it.
- The 59 tomato2 run's B1 reference is the review's `t59_committed` log (log-likelihood to seven
  figures); no B1 trace exists, so the comparison is to seven significant figures.
