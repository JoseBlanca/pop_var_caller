# Fix Application Report: snp_fit_convergence_B_2026-10-17.md

**Date:** 2026-10-17
**Source review:** `doc/devel/reports/reviews/snp_fit_convergence_B_2026-10-17.md`
**Source state reviewed against:** `snp-fit-convergence` at dda9e929 (the review read f2ef6574; the
code it names is unchanged between the two except for the C and D steps' additions)
**Execution mode:** non-interactive (the lead's task named how to fix each finding)
**Overall status:** Completed

---

## 1. Executive summary

The review found three ways the box solve of the Newton step and the decrement read from it could
keep a start from stopping where it should, and three parts of the solve no test guarded. All eight
findings and both nits are applied. **The box solve now frees no parameter twice, so it always
finishes** — the review's suggestion (block only a freed parameter that will not move back in) does
not end every cycle, as its own case 573 shows; **a parameter the golden section left beside an end
is on it**, which takes the solve at 2,000 outbred samples at their floor from 1.83 s to 1.9 ms; **a
rounding-negative determined decrement is raised to zero**. On the 1- and 4-sample oracle both
parameters files are byte for byte those of the build before. Of the digest checksums only the
21-sample fit's moves, by the beside-an-end rule: same starts and passes, numbers within 6 × 10⁻⁶
of themselves (§10).

### Review totals
- Majors: 3 (findings 1, 2, 3)
- Minors: 5 (findings 4 to 8)
- Nits: 2

### Outcome totals
- Applied: 7 (2, 3, 5, 6, 7, 8, nit 1)
- Applied with adaptation: 2 (1, 4)
- Already fixed, Deferred, Disputed, Failed validation, Blocked, Superseded: 0
- Kept in view, no change asked: 1 (nit 2)
- Awaiting user answer: 0

### Validation summary
- `cargo fmt --check` → 0
- `cargo clippy --release --all-targets --all-features -- -D warnings` → 0
- `cargo test --release --lib parameter_estimation::joint` → 0, 473 tests (469 run, 4 ignored)
- `cargo test --release --features bench-fixtures --no-fail-fast` → lib 5,201 passed, 3 failed (the
  three debug-only panic tests: two `an_offset_past_the_reference_panics`,
  `lgamma_panics_on_non_positive_in_debug`), 7 ignored; every other target passes
- `cargo doc --no-deps`, `cargo audit` → not run: no public doc or dependency change
- Oracle identity (`estimate-parameters` on the 1- and 4-sample oracle psps) → both files
  byte-identical to `tmp/checkpoint_D/{s1,o4}.toml`

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| 1 | Major | the box solve can cycle | Apply | Applied with adaptation | `fit/standard_errors.rs` | Pass |
| 2 | Major | a rounding-negative determined decrement never converges | Apply | Applied | `fit/standard_errors.rs` | Pass |
| 3 | Major | the NaN for an unfinished solve is untested | Apply | Applied | `fit/standard_errors.rs` | Pass |
| 4 | Minor | quadratic cost beside an end | Apply | Applied with adaptation | `fit.rs`, `fit/standard_errors.rs` | Pass |
| 5 | Minor | S1, S3, S13 survive | Apply | Applied | `fit/standard_errors.rs` | Pass |
| 6 | Minor | B4 pinned only by checksums | Apply | Applied | `cli/cross_platform_digests.rs` | Pass |
| 7 | Minor | the note does not say which rate calling uses | Apply | Applied | `fit.rs`, `parameter_prepass_joint_fit.md` | Pass |
| 8 | Minor | errors tied to their sums by a length only | Apply | Applied | `fit/standard_errors.rs` | Pass |
| Nit 1 | Nit | "the smallest that length can be" overstates | Apply | Applied | `fit/standard_errors.rs`, `fit_precision.md` | — |
| Nit 2 | Nit | with no parameter determined a start stops at its first judging pass | — | No change (accepted by the owner) | None | — |

## 3. Questions asked and answers

None. The lead's task said how to fix each finding.

## 4. Per-finding log

### 1 — the box solve can cycle
- **Final status:** Applied with adaptation.
- **Reproduced first:** the review's stress test ported to this tree reproduces its numbers exactly:
  28 of 6,000 cases unfinished (18, 28, 108, 203, 573, …), 125 of 3,546 determined decrements
  negative.
- **Why the suggested fix was not enough:** blocking a freed parameter whose next step points out
  (the suggestion) left case 5,398 unfinished. Traced: two parameters, 9 and 11, take turns — 9
  freed, a whole step; 11 freed, moving until 9 meets its end again (reach 7 × 10⁻⁷), then 11 meets
  its own (8 × 10⁻⁷); repeat. Each cycle gains about 10⁻⁹ of a step worth 10⁻⁶, so the model's value
  rises steadily and the limit comes first. Every freed parameter moves inwards, so blocking on an
  outward step never fires. Case 573, in the review, progresses the same way.
- **Implementation:** **no parameter is freed twice in one solve** (`freed_before`). A freed parameter
  whose step then points out has no room (it sits on its end) and is fixed again in the next round.
  Each parameter is then fixed at most twice and freed at most once, and every round fixes one or
  frees one, so the solve takes at most `3 × side + 1` rounds: it always finishes, and the limit of
  `10 × side + 10` is a guard.
- **What it costs:** on information where a second freeing is needed the point the solve finishes at
  need not be the box's maximum. Measured: on information well conditioned no parameter needed a
  second freeing — 0 of the 2,400 harsh cases summed from twice or three times as many positions as
  parameters, 0 of the property test's 134 well-conditioned cases — while 705 of the 3,600 nearly
  singular harsh cases (and 31 of the property test's nearly singular ones) did. On the 1,200 best
  conditioned harsh cases the step passes a check of the box's maximum (the slope left at each
  parameter is zero inside its interval and points out on an end); freeing every inward-pointing
  parameter at once fails it.
- **Test:** `the_box_solve_finishes_on_nearly_singular_information` — all 6,000 harsh cases finish
  within `3 × side + 1` rounds, every move inside its interval, worth ≥ 0, and the box's maximum on
  the 1,200 best conditioned. Allowing a second freeing fails it (and two other tests).
- **Fixture ported:** the review's stress sequence, as `harsh_cases()` (the same draws, so case 203
  is case 203).

### 2 — a rounding-negative determined decrement
- **Final status:** Applied.
- **Implementation:** `determined_decrement` returns 0 for a negative length and leaves a NaN a NaN
  (written as `if length < 0.0 { 0.0 } else { length }`, since `f64::max` would turn a NaN into 0).
  `has_converged` keeps refusing negatives and NaN.
- **Test:** `on_exactly_singular_information_the_determined_decrement_is_not_below_zero` — the 100
  floor-0 arrows among the first 1,000 harsh cases; without the raise it fails on case 311
  (−2.1 × 10⁻¹¹, one of the review's listed cases). After the fix 0 of 3,574 profiled harsh cases
  are negative.

### 3 — the NaN for an unfinished solve is untested
- **Final status:** Applied.
- **Implementation:** `bounded_step` takes its round limit; `newton_step` passes
  `box_solve_round_limit(side)` through a private `newton_step_within`, which a test calls with a
  short limit. With finding 1's fix the limit cannot be reached otherwise.
- **Test:** `a_box_solve_that_runs_out_of_rounds_gives_a_decrement_that_is_not_a_number` — a bound
  binding on both matrices; given one round the decrement is NaN and `has_converged` refuses it;
  given two, it equals `newton_step`'s bit for bit. Mutation S4 fails it.

### 4 — quadratic cost beside an end
- **Final status:** Applied with adaptation (the second of the two fixes offered).
- **Why not the first:** fixing together every parameter whose room is within `reach × (1 + 10⁻¹²)`
  of the shortest does nothing for the timing case: the excesses sit 1 to 2 × 10⁻¹⁰ above 0 at
  different distances, so their rooms differ by factors, not by 10⁻¹².
- **Implementation:** a parameter the solve has not yet moved that sits within `BESIDE_AN_END` of
  the end its step heads to has no room (`room_to_an_end`). `BESIDE_AN_END` is the golden section's
  stopping width, now one constant in `fit.rs` (`GOLDEN_SECTION_WIDTH`, 10⁻⁹ of one plus the end's
  size) used by both. **Only before the first move**: applying it every round undid the small moves
  of freed parameters and made 205 of the 6,000 harsh cases cycle (measured before finding 1's rule
  was in; with the rule in, the same variant fails the box-maximum check).
- **Timing** (the review's case, `bounded_step` alone, release build, this machine):

  | samples (parameters) | excesses in the middle | beside the floor, before | beside the floor, after |
  |---|---|---|---|
  | 59 (185) | 40 µs | 1.66 ms | 58 µs |
  | 500 (1,508) | 0.32 ms | 114 ms | 0.45 ms |
  | 2,000 (6,008) | 1.29 ms | 1.83 s | 1.9 ms |

- **Test:** `parameters_beside_an_end_are_fixed_there_in_one_round` — the 2,000-sample case finishes
  within two rounds, every excess moved exactly its gap to 0. Without the rule it fails.

### 5 — S1, S3, S13 survive
- **Final status:** Applied.
- **S1** (free every inward-pointing fixed parameter at once): caught by
  `the_box_solve_finishes_on_nearly_singular_information` (case 5,475: a parameter on its ceiling
  with its slope pointing in).
- **S3** (`SLOPE_POINTING_IN` = 0): caught by
  `a_parameter_on_an_end_whose_slope_is_rounding_is_not_freed` — eight parameters on their ends
  where the box's maximum has them exactly, so the slope there is rounding; the solve finishes in two
  rounds on the whole matrix and three on the arrow, and needs five on the whole matrix under S3.
- **S13** (room without `.max(0.0)`): **no test of the solve's answer can catch it**, because after a
  negative room the next round's solve returns to the same maximum — so the room is now its own
  function, `room_to_an_end`, and `a_parameters_room_to_an_end_is_never_below_zero` pins its
  contract (no room for a parameter past its end; the gap's share inside; none beside an end before
  the first move). S13 fails it. An end-to-end test with a value handed over outside its interval
  was tried first and passed under S13.

### 6 — B4 pinned only by checksums
- **Final status:** Applied.
- **Test:** `every_start_on_the_cohort_with_sequencing_errors_converges` (in
  `cli::cross_platform_digests`, which has the walk-and-fit plumbing; `walked_and_fitted` split out
  of `walked_fitted_and_called`) — reads the parameters file's `snp_indel_fit_starts`: three starts,
  all converged, each under 50 passes (22, 13, 13). Under S7 (every parameter given an error) all
  three run to 200 passes and it fails.

### 7 — calling's rate on an inverted library
- **Final status:** Applied.
- **Implementation:** the run-log note ends "calling uses the ordinary-position rate as each
  library's sequencing error rate, which for these is the higher of the two"; its doc says the same;
  its test checks the clause. `parameter_prepass_joint_fit.md` §3.3 gains a bullet recording it as
  accepted.

### 8 — errors from other sums
- **Final status:** Applied.
- **Implementation:** `newton_step_within` begins with
  `debug_assert_eq!(errors, &StandardErrors::of(sums))`. Every caller passes errors computed from the
  same sums one line earlier; the check costs one more error computation a judging pass in debug
  builds only.

### Nit 1 — "the smallest that length can be" overstates
- **Applied:** `determined_decrement`'s doc and `fit_precision.md` §2 say that when the solve over the
  profiled parameters leaves one still, the length is a minimum over fewer directions and can exceed
  `λ_K` (12 of 3,546 cases, up to 3.3 times, never above 2.5 × 10⁻⁵).

### Nit 2 — no parameter determined
- No change; the review lists it as accepted by the owner.

## 5. Deferred findings to carry forward
None.

## 6. Disputed findings to return to reviewer
- Finding 1's suggested rule, not the finding: blocking only a freed parameter whose step points out
  does not make the solve terminate (case 5,398 above). Applied with the stronger rule instead.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no criterion run. No bench in `benches/` covers the box solve; the change is timed
  directly instead (finding 4's table) and on the oracle runs.
- The two oracle runs at 8 threads: 1 sample 38 s, 4 samples 150 s (not comparable with
  `tmp/checkpoint_D`'s, which were run at another width).

## 10. Identity
- **Oracle:** the release `pvc` of the final tree, run on the 1 and 4 oracle psps (reference
  `S_lycopersicum_chromosomes.4.00.fa`, the tomato1 catalog), writes parameters files byte-identical
  to `tmp/checkpoint_D/s1.toml` and `o4.toml` (`cmp`). Checked twice: after the first version of the
  solve change and after the final one.
- **Digests:** the two-sample checksums (with and without sequencing errors) do not move, nor any
  calls checksum. **The 21-sample fit's file moves** (`728256994c0c90109ae7f4a15eeb00a9` →
  `c548c669d12d20097217ff03e1299e3a`), re-recorded with its reason in the constant's doc. Cause, by
  switching each change off: the beside-an-end rule alone (with it off the old checksum returns).
  The solve reaches the same maximum by another path, so its steps differ in their last bits and
  so do the points the fit tries them at. The starts' rows are unchanged (same endings and passes);
  the fitted numbers move by at most 6 × 10⁻⁶ of themselves and 4 × 10⁻⁶ of their standard errors,
  the errors by at most 4.4 × 10⁻⁵ of themselves.

## 11. Mutation results
Each applied with the Edit tool, `cargo test --release --lib parameter_estimation::joint` run, then
reverted (logs `tmp/b_fix/mut_*.log`); S7 run against the new fit-level test.

| mutation | caught by |
|---|---|
| S1 free every inward-pointing parameter at once | `the_box_solve_finishes_on_nearly_singular_information` |
| S3 `SLOPE_POINTING_IN` = 0 | `a_parameter_on_an_end_whose_slope_is_rounding_is_not_freed` |
| S4 a decrement for an unfinished solve | `a_box_solve_that_runs_out_of_rounds_gives_a_decrement_that_is_not_a_number` |
| S7 every parameter given an error | `every_start_on_the_cohort_with_sequencing_errors_converges` |
| S13 room without `.max(0.0)` | `a_parameters_room_to_an_end_is_never_below_zero` |
| no raise of a negative determined decrement | `on_exactly_singular_information_the_determined_decrement_is_not_below_zero` |
| a parameter may be freed again | that test, `the_box_solve_finishes_on_nearly_singular_information`, `the_box_step_stays_inside_the_box_and_is_worth_at_least_nothing` |
| no beside-an-end rule | `parameters_beside_an_end_are_fixed_there_in_one_round`, `a_parameters_room_to_an_end_is_never_below_zero` |
| the rule applied every round | `the_box_solve_finishes_on_nearly_singular_information`, `a_parameters_room_to_an_end_is_never_below_zero` |

## 12. Notes
- One slip: the S3 mutation was first applied with `sed -i` (against the task's rules); it was
  reverted with the Edit tool and its backup moved to `tmp/b_fix/`, and every other mutation used
  Edit.
- Scratch: `tmp/b_fix/` (probe output, traces, mutation logs, oracle runs).
