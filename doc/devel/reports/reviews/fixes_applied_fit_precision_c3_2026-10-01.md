# Fix Application Report: fit_precision_c3_2026-10-01.md

**Date:** 2026-10-01 (applied); 2026-10-02 (round limit, prose, oracle)
**Source review:** `doc/devel/reports/reviews/fit_precision_c3_2026-10-01.md`
**Source state reviewed against:** `f7aa3174` (a review object of the working tree, parent `91995ed5`)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 7
- Minors: 21
- Nits: grouped

### Outcome totals
- Applied: 17 (M1, M2, M4, M6, M7; Mi2–Mi7, Mi12, Mi15, Mi17–Mi20), and six of the Nits
- Applied with adaptation: 2 (M3, M5)
- Applied in part: 1 (Mi8)
- Deferred: 8 (Mi1, Mi9, Mi10, Mi11, Mi13, Mi14, Mi16, Mi21), and four of the Nits
- Disputed: 1 (part of M3: who chose the put-off and the stop on a losing round)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --release --lib parameter_estimation::joint::ssr_fit` → 0, 60 passed, 3 ignored
- `cargo test --all-targets --all-features --no-fail-fast` → 101, 5,010 passed, 3 failed (the pre-existing
  `examples/ng_generic_loci_dump.rs` × 2 and `examples/ng_ssr_loci_dump.rs` × 1), 8 ignored; the cross-platform
  checksum tests pass unchanged
- Performance check → not applicable: no file under `benches/` reaches the repeat-tract climb. The step's cost is
  measured on the oracle cohort instead (implementation report §2.2).

### Unresolved high-priority findings
None. Mi8's stale plan and spec text needs the owner's approval to change (checkpoint C).

## 2. Findings table

| ID | Severity | Title | Final status | Files changed |
|---|---|---|---|---|
| M1 | Major | Newton carriage and indexing untested | Applied | `ssr_fit.rs` tests |
| M2 | Major | the round limit's reason contradicted | Applied | `ssr_fit.rs`, spec §4.3, impl report |
| M3 | Major | the §4.3 amendment disagrees with the code | Applied with adaptation; attribution Disputed | `ssr_fit.rs`, spec §4.3 |
| M4 | Major | report §2.2 miscounted | Applied | impl report |
| M5 | Major | the trigger's units untested at the call | Applied with adaptation | `ssr_fit.rs` |
| M6 | Major | the held point's score unchecked | Applied | `ssr_fit.rs` test |
| M7 | Major | the comparison asserts nothing | Applied | `ssr_fit.rs` test |
| Mi1 | Minor | no error at all counts as settled | Deferred | — |
| Mi2 | Minor | a `NaN` start never stops | Applied | `ssr_fit.rs` |
| Mi3 | Minor | GQ moved by up to 2 | Applied | baseline, impl report |
| Mi4 | Minor | 307 evaluations a round, 1.7 rounds a judgement | Applied | spec §4.3, impl report |
| Mi5 | Minor | "with the same answers" | Applied | impl report |
| Mi6 | Minor | "no walk settled" | Applied | impl report |
| Mi7 | Minor | the 17.6 blamed on the old limit | Applied | `ssr_fit.rs` |
| Mi8 | Minor | stale plan, spec and doc text | Applied in part | `ssr_fit.rs`, spec §4.3 |
| Mi9 | Minor | `ending` and `walks` can disagree | Deferred | — |
| Mi10 | Minor | one triple in two tuple orders | Deferred | — |
| Mi11 | Minor | `Walked` repeats `Climb` | Deferred | — |
| Mi12 | Minor | distances paired with errors by position | Applied | `ssr_fit.rs` |
| Mi13 | Minor | "climb" and "walk" for one thing | Deferred | — |
| Mi14 | Minor | `PeriodLengthSpectrum::converged` a `bool` | Deferred | — |
| Mi15 | Minor | `settled_fraction`'s default unstated | Applied | `ssr_fit.rs` |
| Mi16 | Minor | empty `starting_points` panics | Deferred | — |
| Mi17 | Minor | the climb summary untested | Applied | `ssr_fit.rs` test |
| Mi18 | Minor | no test sets `settled_fraction` | Applied | `ssr_fit.rs` test |
| Mi19 | Minor | schedule tests blind to endings and reuse | Applied | `ssr_fit.rs` tests |
| Mi20 | Minor | no two-group fixture | Applied | `ssr_fit.rs` test |
| Mi21 | Minor | the gain after a losing round unspecified | Deferred | — |
| Nits | Nit | | Applied in part | `ssr_fit.rs`, spec §4.3 |

## 3. Questions asked and answers

None during the run. Mi1 and Mi8 go to the owner at checkpoint C.

## 4. Per-finding log

### M1 — the Newton distance's carriage and indexing untested
- **Applied.** `newton_distances_are_carried_to_each_numbers_scale_like_the_errors` builds a coupled inverse and a
  slope by hand and checks every distance against `I⁻¹ g` carried as the errors are: `p(1 − p)` for each slippage
  number, the concentration itself, and each share through its slopes in the log-ratios, the largest class's
  included. It runs once with every coordinate kept and once with one dropped, which pins the indexing.
  `the_furthest_number_counts_every_kind_and_a_distance_not_a_number_is_unsettled` checks that a share and the
  concentration can each be the furthest number.
- **Verification:** the mutations `share_cross_terms` (the cross terms dropped) and `concentration_unscaled` (the
  concentration left on the log scale) are killed, each by the first test (`tmp/fit_precision/c3/mut/`).

### M2 — the round limit's reason contradicted
- **Applied.** The comparison was re-run with room to 60 rounds
  (`tmp/fit_precision/c3/fixes_{3,13}class_60.log`). Walks that settled took 2 to 6 rounds at three classes and
  6 to 36 at thirteen; one thirteen-class walk had not settled by 60, and none settled between 37 and 60.
  `DEFAULT_MAX_ROUNDS` is 40, its doc states that range and that no measurement says whether a walk settles after
  40. The spec and the report say the same, and the report gives the put-off's lengthening (17 rounds to 29 on one
  walk).

### M3 — the amendment disagrees with the code
- **Applied with adaptation.** The code was changed, not only the text: a losing round with nothing better than
  its start — or whose gain is not a number — is now judged whatever the put-off says, so a walk never ends as
  having lost a round at a point it did not judge. A losing round with a better point part-way through still goes
  back to that point unjudged while a judgement is put off; the spec now says so. The spec, `Climb`'s and the
  judgement's docs now say the errors are reused whenever the winning walk's last round was judged at the point it
  returns, settled or not.
- **Disputed in part:** the put-off and the stop at a stuck losing round are the owner's. The owner approved the
  stopping rule on 2026-10-01 with both in it: an unsettled verdict puts the next judgement off by one round more
  each time, and a losing round is undone and, if stuck, judged and the walk stopped.
- **Tests:** `an_unsettled_verdict_puts_the_next_judgement_off` now has a losing round inside the put-off judged
  (judgements at 1, 2 and 2 again).

### M4 — report §2.2 miscounted
- **Applied.** §2.2 rewritten from the final run and from the 40-round comparison's per-stratum lines
  (`tmp/fit_precision/c3/oracle_newton_lines.txt`), every count re-derived: five strata, 0.000 to 0.054 units
  short; the near-zero shares named with their values; 2:9's share of 0.031 against 0.028 stated as the rule's
  verdict being wrong there; 1:15's walks as recorded; 1:8 as equal by construction.

### M5 — the trigger's units untested at the call
- **Applied with adaptation.** Rather than a judgement-count test, the curvature and the walk now read the tract
  count from one function, `tracts_a_mean_is_over`, so the walk's units cannot drift from the errors' without the
  error tests seeing it: the mutation making that function return 1 is killed by three tests
  (`a_number_the_climb_left_at_the_end_of_its_reach_is_not_placed`, `a_stratums_errors_shrink_as_its_tracts_grow`,
  `tracts_without_reads_change_no_error`).
- **Residual risk:** a constant written at the call in `climb_from` itself, in place of the function, is not
  caught by a test.

### M6 — the held point's score unchecked
- **Applied.** `the_held_point_scores_what_it_holds`: over four rounds from each starting point of a drawn
  stratum, the held point scored afresh equals its held score to the bit, and is at least the round's start and
  end.

### M7 — the comparison asserts nothing
- **Applied.** At the default draws, every fit whose winning walk settled must be within a tenth of an error of
  the reference on every number and within ½ · p · 0.1² of its log-likelihood. Both 60-round runs passed them.

### Mi1 — no error at all counts as settled
- **Deferred.** It is the spec's rule ("settled by definition"); changing it is the owner's. Raised at
  checkpoint C.

### Mi2 — a `NaN` start never stops
- **Applied.** A gain that is not a number counts as stuck, so the walk is judged and stops as having lost a
  round. `a_walk_from_a_start_that_is_not_a_number_lost_a_round`; `a_point_that_scores_not_a_number_is_not_judged_settled`.

### Mi3 — GQ moved by up to 2
- **Applied**, re-measured on the final oracle run (report §2.4, baseline comment).

### Mi4 — 307 evaluations a round
- **Applied.** Counted in `climb_one_round`: 17 golden sections of 18 evaluations and the score where the round
  ends; a judgement, 513, is about 1.7 rounds.

### Mi5, Mi6, Mi7
- **Applied.** Report §4 item 3 gives the 22 rounds and draw 0's moved answer; §2.1 says two of the fifth draw's
  walks settled; `DEFAULT_MAX_ROUNDS`'s doc drops the 17.6.

### Mi8 — stale text
- **Applied in part.** `DEFAULT_REFUSAL_FLOOR`'s 83% is dated to the rule before; the §4.3 amendment says §4.1's
  description is the rule before; a test doc's "the default five" corrected. **For the owner:** the plan's C3
  entry and checkpoint C, spec §4.1, §4.2, §4.5 item 2 and §6 — design text this run may not edit.

### Mi9, Mi10, Mi11, Mi13 — the walk's types
- **Deferred.** Recording the winner instead of duplicating its ending, named structs for the two tuples, one
  type for `Walked` and `Climb`, and one word for a walk are one refactor of the climb's types with no behaviour
  change; done together, after step D, which reshapes the climb again (the warm start).

### Mi12 — distances paired by position
- **Applied.** `SlippageDistances { level, shorter_share, fall_off }`, paired with `SlippageErrors` by name.

### Mi14 — `PeriodLengthSpectrum::converged`
- **Deferred.** A public field's type change, which belongs with Mi9.

### Mi15 — the default unstated
- **Applied.** The field's doc says "[`SETTLED_FRACTION`], 0.1, by default".

### Mi16 — empty `starting_points`
- **Deferred.** It predates this step.

### Mi17, Mi18, Mi19, Mi20 — missing tests
- **Applied.** `the_climb_summary_counts_each_ending`; `the_settled_fraction_is_read_from_the_config`;
  `the_winning_walks_errors_are_the_errors_at_its_answer_in_both_schedules`, and the schedule tests' flattened
  numbers now carry each walk's ending, rounds and judgements; `a_walk_at_two_slippage_groups_is_judged_on_both`.

### Mi21 — the gain after a losing round
- **Deferred.** The spec does not say which gain the next projection uses; the projection is only a trigger, so
  either choice moves when a walk is judged, not whether it has settled.

### Nits
- **Applied:** `WalkRecord`'s fields documented; "one per starting point"; `LostARound`'s doc leads with what it
  is; `remaining_gain_target` names the stopped report's file; "not placed" defined in the spec; "the default
  limit" in a test doc.
- **Deferred:** `is_settled()` (with Mi9); the one re-score a round, 1 evaluation in 307; the comparison test's
  knobs (a measurement tool); `HeldPoint::offer` on exact ties (equivalent).

## 5. Deferred findings to carry forward
- Mi1 — a judgement with no error at all counts as settled (the owner's rule).
- Mi9, Mi10, Mi11, Mi13, Mi14 — one refactor of the walk's types, after step D.
- Mi16 — empty `starting_points` (code that predates the step).
- Mi21 — the gain after a losing round.

## 6. Disputed findings to return to reviewer
- M3, in part — the put-off and the stop at a stuck losing round were approved by the owner on 2026-10-01.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
- **Triggered:** no — nothing under `benches/` reaches the repeat-tract climb. Its cost on the oracle cohort is in
  the implementation report §2.2.

## 10. Commands run
- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::ssr_fit`
- `scripts/dev.sh cargo test --release --lib the_climbs_stop_loses_nothing_against_a_longer_climb -- --ignored
  --nocapture`, at three and at thirteen classes, limit 60
- `tmp/fit_precision/c3/mut/run.sh` (three mutations)
- `scripts/dev.sh sh tmp/fit_precision/c3/run_oracle_c3g.sh tmp/fit_precision/oracle_c3g`
- `scripts/dev.sh cargo fmt --check`; `cargo clippy --all-targets --all-features -- -D warnings`;
  `cargo test --all-targets --all-features --no-fail-fast`

## 11. Command results
- Module: 60 passed, 3 ignored
- Comparison, three classes: passed, 2,422 s; thirteen classes: passed, 1,839 s
- Mutations: all three killed
- Oracle: the two fit lines differ from the old baseline, as expected, and are re-recorded; the other five match
- fmt, clippy, suite: as in §1

## 12. Notes
- **Two defects found after the review.** The baseline comment and report §2.4 said every one of the 50 moved
  lines was a slippage row, and the review's number check agreed; counted on the final run, they are 37 slippage
  rows and 13 strata's length spectra. Corrected in both. And clippy rejected Mi2's `!(gain > 0.0)`
  (`neg_cmp_op_on_partial_ord`); it is now `!matches!(gain.partial_cmp(&0.0), Some(Ordering::Greater))`, true for
  the same values, a `NaN` included, so the oracle run made before it stands.
- The two 60-round comparison runs used `DEFAULT_MAX_ROUNDS` = 60, a placeholder, before it was set to 40. At 40 the
  drawn strata's answers are the same: the one walk that ran past 40 (thirteen classes, draw 3) did not win, and a
  walk's score never falls from one round to the next, so cutting it at 40 cannot make it win.
