# Fix Application Report: str_bundles_hard_simulator_2026-10-13.md

**Date:** 2026-10-13
**Source review:** `doc/devel/reports/reviews/str_bundles_hard_simulator_2026-10-13.md`
**Source state reviewed against:** b2714614 plus the uncommitted diff, branch `str-bundles-hard` (worktree `pop_var_caller-str-bundles-hard`); nothing committed
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 5
- Minors: 22
- Nits: 7 (grouped)

### Outcome totals
- Applied: 17 (M1–M5, Mi1, Mi2, Mi4–Mi11, Mi21, Mi22; M5 in a second pass, on top of 90659177)
- Applied with adaptation: 2 (Mi3, Mi18)
- Already fixed: 0
- Deferred: 7 (Mi12–Mi17, Mi20)
- Disputed: 1 (Mi19)
- Failed validation: 0
- Blocked by context mismatch: 0
- Superseded: 0
- Awaiting user answer: 0

### Validation summary
- `cargo fmt -- --check examples/ng_str_stutter_benchmark.rs` → exit 0.
- `cargo clippy --release -j 2 --example ng_str_stutter_benchmark` → exit 0. No warning in the example. With `-D warnings` it exits 101 on a pre-existing library lint, `src/run/edge_deletion_gathering.rs:201`, outside the diff.
- `cargo test --release -j 2 --example ng_str_stutter_benchmark` → exit 0, `22 passed; 0 failed`.
- `cargo test --all-targets`, `cargo doc`, `cargo audit` → not run. Only the example and two scripts changed, and the request limited builds on the loaded host.
- Performance check → not applicable: no `benches/` path touched.

### Unresolved high-priority findings
None.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | User input | Files changed | Validation | Follow-up |
|---|---|---|---|---|---|---|---|---|
| M1 | Major | Run-rule test misses allele windows | Apply | Applied | No | `examples/ng_str_stutter_benchmark.rs` | Pass | No |
| M2 | Major | Hard-site test checks properties, not the rule | Apply | Applied | No | same | Pass | No |
| M3 | Major | `ng_str_bundle_runs` failure is silent | Apply | Applied | No | `scripts/ng_str_stutter_benchmark.sh` | Pass (success path in 4 runs) | No |
| M4 | Major | Bundles with no reference run almost never laid | Apply | Applied | No | example, both specs | Pass | No |
| M5 | Major | Scanner copy blind to the flank | Defer, then Apply (coordinator) | Applied | No | example, `str_stutter_benchmark.md` | Pass | No |
| Mi1 | Minor | Junction after a tail abuts at 0.4 | Apply | Applied | No | example, spec | Pass | No |
| Mi2 | Minor | Zero-share test doc comment | Apply | Applied | No | example | Pass | No |
| Mi3 | Minor | Kinds test cannot tell both-kinds shapes | Apply | Applied with adaptation | No | example | Pass | No |
| Mi4 | Minor | Unchecked flank subtraction | Apply | Applied | No | example | Pass | No |
| Mi5 | Minor | Missing caller row shown as `other_runs` | Apply | Applied | No | check script | Pass | No |
| Mi6 | Minor | Summary wording reads as a pass | Apply | Applied | No | check script | Pass | No |
| Mi7 | Minor | Bare `KeyError` | Apply | Applied | No | check script | Pass | No |
| Mi8 | Minor | `small_cohort` silently on hard shapes | Apply | Applied | No | example | Pass | No |
| Mi9 | Minor | Default cites a gitignored script | Apply | Applied | No | example | Pass | No |
| Mi10 | Minor | Unnamed `30` and `3` | Apply | Applied | No | example | Pass | No |
| Mi11 | Minor | Wildcard arm absorbs `(false, false)` | Apply | Applied | No | example | Pass | No |
| Mi12 | Minor | Tuple triples | Defer | Deferred | No | None | N/A | Yes |
| Mi13 | Minor | Run-rule copy's structure | Defer | Deferred | No | None | N/A | Yes |
| Mi14 | Minor | Probe built three times | Defer | Deferred | No | None | N/A | Yes |
| Mi15 | Minor | Interruption kinds as `below(4)` | Defer | Deferred | No | None | N/A | Yes |
| Mi16 | Minor | `repeats` shadowing, long function | Defer | Deferred | No | None | N/A | Yes |
| Mi17 | Minor | Anonymous pair return | Defer | Deferred | No | None | N/A | Yes |
| Mi18 | Minor | `read_runs` rows by column number | Apply | Applied with adaptation | No | check script | Pass | No |
| Mi19 | Minor | Share not in the run's log | Dispute | Disputed | No | None | N/A | No |
| Mi20 | Minor | Helpers without direct tests | Defer | Deferred | No | None | N/A | Yes |
| Mi21 | Minor | "Level" overstated in step 5 | Apply | Applied | No | `str_bundle_loci.md` | N/A (prose, numbers re-derived) | No |
| Mi22 | Minor | Flank claim too broad; rule-placement unstated | Apply | Applied | No | example, `str_stutter_benchmark.md` | Pass | No |

## 3. Questions asked and answers

None. Open question 1 of the review (should stutter ever be placed by anything but the caller's rule) is recorded in the spec as a limit of the benchmark, not acted on.

## 4. Per-finding log

### M1: run-rule test misses allele windows
- **Severity:** Major. **Initial decision:** Apply. **Final status:** Applied.
- **Implementation summary:** `the_run_rule_copy_agrees_with_the_library` now also compares the two on every allele window of `small_cohort(4, 150)`: flanks, length changes, substitutions, SNPs beside the hull and edge deletions. It asserts more than 300 windows were seen. No disagreement was found.
- **Tests added or modified:** `the_run_rule_copy_agrees_with_the_library`.
- **Validation:** `cargo test --release -j 2 --example ng_str_stutter_benchmark` → 0, 22 passed.

### M2: hard-site test checks properties, not the rule
- **Final status:** Applied.
- **Implementation summary:** `a_hard_bundle_stutters_in_pure_stretches_at_the_floor` now asserts that the reference allele's sites equal `find_reference_runs` on its hull, moved by the 30-base flank. It also asserts that there are more sites than hard bundles, and at least one hard bundle with none. The earlier property checks on every allele stay.
- **Validation:** as M1.

### M3: silent failure of `ng_str_bundle_runs`
- **Final status:** Applied.
- **Implementation summary:** on failure the driver prints `ng_str_bundle_runs failed; its log, …:` and the log's last 20 lines to stderr, then exits 1.
- **Validation:** the success path ran in all four 1-sample re-runs. The failure path was not provoked.

### M4: bundles with no reference run
- **Final status:** Applied.
- **Implementation summary:** new constant `NO_RUN_HARD_BUNDLE_SHARE = 123/1447`. For each hard shape, `draw_hard_bundle_shape` draws once whether this bundle has no run. If so, every repeat is interrupted and the shape is redrawn until `reference_runs(hull)` is empty. The draw happens inside the hard draw, so a zero share still spends none.
- **Tests:** `hard_shapes_keep_their_rules` requires 1 in 20 to 1 in 8 of 2,000 shapes to have no run. `a_hard_bundle_stutters_in_pure_stretches_at_the_floor` requires at least one.
- **Verification performed:**
  - **Byte identity**, `tmp/identity.sh` against the build before hard shapes: 23, 15 and 17 files equal on the three configurations, 0 differing; the only new file is `bundle_runs.tsv`.
  - **Calibration** on 3,000 hard-only shapes: 243 have no run (8.3 in 100, HG002 8.5); impure 79, nested or overlapping 35, abutting runs 16, touching 19 in 100; hull median and 90th percentile 48 and 80.
  - **Default-share runs**: 297 of 500 hard by HG002's classification at the default seed and 300 at seed 7. 31 no-run bundles in a default run.
  - **Truth scored as a caller** on the new `pcr_1s_5x`: all 1,202 variant and 1,798 reference segment cells right, and all 203 variant and 297 reference bundle cells.
  - **1-sample comparisons re-run** (PCR and PCR-free, 5× and 30×, `tmp/bench/run_norun.sh`, results in `tmp/bench/runs_norun/`, summary `tmp/bench/summary_norun.tsv`). Step 5 of `str_bundle_loci.md` gains the new table. The 2026-10-12 table stays, labelled as measured before the change.
- **Follow-up:** none; the 50-sample rows were re-run after the M5 fix (see Notes).

### M5: scanner copy blind to the flank
- **Final status:** Applied (second pass, at the coordinator's request, on top of commit 90659177).
- **Cause:** shapes were accepted by scanning the hull with a spacer either side, never the 150-base flanks laid after. A repeat can read on through a one-base spacer: `(AAT)6 A8`, spacer `C`, flank `AACAG…` was typed 5 bases longer than laid.
- **Implementation summary:** `build_reference` now lays each bundle through `lay_bundle`, which re-runs the scanner copy over flank, hull and flank and redraws the flanks (truncating back) until the typing sees the bundle as accepted: the drawn members for a separated shape, one bundle over exactly the hull for a hard one. `scanner_types_one_bundle_over` and `scanner_finds_exactly_the_members` were split out of the two existing checks so both take any probe.
- **Separated shapes were affected too**: 2 of 20,000 failed with their first flanks (a `TC` repeat read back into the left flank; a `GGA`-like repeat read on into the right). Under the old build those runs would have stopped at the typing check.
- **Tests added:** `a_bundle_is_never_typed_past_its_flanks` (the literal M5 sequence fails the flank-aware check though it passed the old one; `lay_bundle` gives that shape passing flanks over 50 seeds; every bundle of a 400-bundle hard-only cohort is typed over exactly its hull with its flanks) and `separated_bundles_rarely_need_new_flanks` (pins 2 of 20,000 and that `lay_bundle` fixes each).
- **Validation:** tests 24 passed; fmt clean; clippy no example warnings. Byte identity at share 0 and at no bundles: `tmp/identity.sh`, 0 of 23, 15 and 17 files differ, and checksums over each configuration's files equal (760193a0bec5092c, ef2c3c6bb4c53991, e453af587b0d1937). Typing check: 3,000 of 3,000 on hard-only runs at the default seed and seed 11; 500 of 500 at the default seed and seed 7.

### Mi1: junction after a tail
- **Final status:** Applied. Junctions are now drawn in order, each after any tail rewrite at the junction before it. The spec states the rule.

### Mi2: zero-share doc comment
- **Final status:** Applied. The comment now says the truth VCF and reference are hashed, and that the value was taken from this build after the full-file comparison. The hash is unchanged and the test passes.

### Mi3: kinds test
- **Final status:** Applied with adaptation. `AbuttingAndInterrupted` must now both abut and be interrupted; `Separated` is `unreachable!`. The suggested `interrupted == 0` for `Abutting` failed on correct shapes: two repeats of one motif a gap apart look interrupted to the test's count. That arm asks only for an abutting junction, and a comment says why.

### Mi4: unchecked subtraction
- **Final status:** Applied. `checked_sub(WINDOW_FLANK_BP).expect("a reference stutter site lies inside the hull")`.

### Mi5, Mi6, Mi7, Mi18: check script
- **Final status:**
  - Mi5, Mi6 and Mi7: Applied.
  - Mi18: Applied with adaptation — named locals and a corrected docstring, not a `DictReader` rewrite.
- **What changed:**
  - A missing caller row is now its own outcome, `no_caller_runs_row`.
  - A bundle missing from `bundle_runs.tsv` exits with its name and position.
  - The summary says "matched means the caller's reference runs equal the simulator's stutter sites too".
  - The two row lookups are bound once, as `laid` and `caller`.
- **Validation:** the check ran on the 3 typing runs and the 4 bench runs. It passed 500 of 500 bundles on each default-share run, and flagged the one M5 bundle on the hard-only run.

### Mi8, Mi9, Mi10, Mi11
- **Final status:** Applied.
  - Mi8: `small_cohort` sets `hard_bundle_share: DEFAULT_HARD_BUNDLE_SHARE` explicitly.
  - Mi9: the default's comment points at §11.4 and says the measuring script was scratch.
  - Mi10: `INTERRUPTED_MEMBER_MOST_EXTRA_COPIES` and `TAIL_LEAD_LEAST_COPIES` replace the literals.
  - Mi11: an exhaustive match with `(false, false) => unreachable!`.

### Mi12–Mi17, Mi20
- **Final status:** Deferred.
- **Reasoning:** these are readability and structure changes. Several touch pre-existing scanner code (`scanner_repeats`, `scanner_reports_exactly_the_members`). Mi15 would change the draw stream again for no change in meaning. None changes a result, and the run-rule and scanner agreement tests guard the behaviour they concern. Worth a follow-up cleanup pass before merge if the owner wants it.

### Mi19: share not recorded
- **Final status:** Disputed. The simulator already prints its whole `Settings`, `hard_bundle_share` included, into `simulate.log`. Seen in `tmp/calibrate/hard6/simulate.log`: `settings : Settings { …, hard_bundle_share: 1.0, … }`.

### Mi21, Mi22: spec and doc claims
- **Final status:** Applied.
  - Step 5 now gives the hull-length numbers: 13 against 13 and 15 against 18 of 29 cells in hulls of at most 30 bases, against 18 against 25 and 24 against 30 of 52 in hulls over 50.
  - `stutter_sites_of`'s doc restricts the no-run-in-the-flank claim to the reference allele.
  - §11.4 now says stutter is placed by the caller's own rule, so this benchmark cannot show that rule wrong, and adds the overlap-cut stretch to what is left out.

## 5. Deferred findings to carry forward
- Mi12–Mi17, Mi20: structure and naming cleanups, and direct tests of three helpers.

## 6. Disputed findings to return to reviewer
- Mi19: the share is already in `simulate.log`.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped: no `Apply` touched code reachable from `benches/`.

## 10. Commands run
- `scripts/dev.sh bash tmp/review_2026-10-13_str_bundles_hard/fmt_test_build.sh` (rustfmt, fmt check, the example's tests, clippy, release build)
- `scripts/dev.sh bash tmp/review_2026-10-13_str_bundles_hard/clippy_build.sh`
- `scripts/dev.sh bash tmp/identity.sh`
- `scripts/dev.sh bash tmp/review_2026-10-13_str_bundles_hard/verify.sh` (typing check: default seed, seed 7, 3,000 hard-only)
- `scripts/dev.sh bash tmp/bench/run_norun.sh pcr 1 5 pcr 1 30 pcr_free 1 5 pcr_free 1 30`
- `bash tmp/review_2026-10-13_str_bundles_hard/truth_as_caller.sh` on `runs/pcr_1s_30x`, `runs/pcr_50s_30x` and `runs_norun/pcr_1s_5x`
- `uv run --no-project python tmp/hg002/compare.py …`, `tmp/bench/summary_norun.py`, `nocall_by_hull.py`, `lag.py`, `norun_split.py`

## 11. Command results
- Tests → exit 0, 22 passed. The first run after the edits failed 3 tests: `small_cohort(_, 300)` did not fit its 0.3 Mb genome, and Mi3's suggested `Abutting` assertion was wrong. Both were fixed as above.
- Identity → 0 files differing in 3 configurations.
- Typing check:
  - default seed: 500 of 500 (94 abutting, 53 both, 158 interrupted, 195 separated);
  - seed 7: 500 of 500;
  - hard-only: 2,999 of 3,000 (M5); after the M5 fix 3,000 of 3,000 at two seeds.
- The four 1-sample runs → typing 500 of 500 each, scored. STR segment cells: 0 of 3,000 differ between the bundle caller and the uncalled build in each.

## 12. Notes
- `PROJECT_STATUS.md` has no block for STR bundles, so it was not edited.
- After the M5 fix all eight runs were re-run (`tmp/bench/runs_m5/`, `tmp/bench/summary_m5.py`); `str_bundle_loci.md` §10 step 5 now holds one current table dated 2026-10-13, the earlier two kept as superseded history.
