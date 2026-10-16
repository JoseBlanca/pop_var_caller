# Code Review: str_penalised_fit_D8
**Date:** 2026-10-16
**Reviewer:** rust-code-review skill (single reviewer, every category covered directly, no sub-agents)
**Scope:** the uncommitted diff on top of `c50d2ce6` (branch `str-penalised-fit`): plan step D8, where a thin class's errors start
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** `git diff` against `HEAD c50d2ce6` in the worktree `pop_var_caller-str-penalised-fit`, 11 files.
- **Reviewed against:** HEAD `c50d2ce6` plus the working-tree changes, exported to `tmp/review_d8/src_copy/` for the probes.
- **In-scope files:**
  - `src/parameter_estimation/joint/ssr_fit/penalised_fit.rs` (constants, `errors_of_a_class`, tests)
  - `src/parameter_estimation/joint/ssr_fit/penalised_fit/run_log.rs` (`withheld_errors_line`, tests)
  - `src/parameter_estimation/joint/ssr_fit.rs` (`StratumError::TooFewLoci` doc)
  - `src/parameter_estimation/joint/stratum_fits.rs` (`median_concentration` doc)
  - `src/calling/parameters_file/to_toml.rs`, `src/calling/parameters_file/testdata/every_shape_as_written.toml`
  - `src/cli/cross_platform_digests.rs`
  - `doc/devel/ng/spec/str_penalised_fit.md` (header, §2.6, §2.7, §2.10.4, §7), `doc/devel/ng/spec/str_fit_to_calling.md` §3.5, `doc/devel/ng/arch/str_penalised_fit.md` (header, §1.1, §1.4), `doc/devel/ng/impl_plan/str_penalised_fit.md` (D8)
- **Deliberately out of scope:** the fit's arithmetic (D8 changes no fitted number); D7's changes (already in HEAD).
- **Categories covered:** reliability (rule-by-rule comparison with HEAD, mutation runs, test challenge); errors (the `UndeterminedErrors` and absent-error paths); naming; defaults (the three constants and where the file and log state them); idiomatic; refactor safety (closure restructure); smells; float_portability (no new arithmetic reaches output: the change is comparisons of `usize` counts and one existing `f64` product, so nothing new is platform-dependent); extras/"diff matches stated intent" (the owner's rule in spec §7). `unsafe_concurrency` and `tooling` do not apply.

## 2. Verdict

**Approve-with-changes.** `errors_of_a_class` does exactly what the owner set on 2026-10-16, and nothing else moved: over 37,800 hand-built cases it agrees with HEAD's function in every outcome except the two the owner changed. The run log counts correctly on four drawn cohorts. The digests moved by note text alone. What needs changing is one wrong explanation in a doc comment (Mi1) and two boundaries that no test pins (Mi2, Mi3).

## 3. Execution status

All runs native cargo 1.98.0 on macOS arm64, `CARGO_TARGET_DIR=tmp/review_d8/target`, `-j 8`.

| command | exit | result |
|---|---|---|
| `cargo test --release --lib review_d8` (the review's probes, in `src_copy`) | 0 | 3 passed; `review_d8 grid: 37800 cases agree` |
| `cargo test --release --lib -- penalised_fit stratum_fits parameters_file cross_platform_digests ssr_fit` (D8 tree) | 0 | `451 passed; 0 failed; 2 ignored` |
| `cargo test --release --lib cross_platform_digests` in `tmp/review_d8/src_undo` (D8's code, HEAD's `to_toml.rs`, HEAD's three fitted checksums) | 0 | `4 passed; 0 failed` |
| 8 mutation runs, `tmp/review_d8/mut/run.sh` | — | 6 killed, 2 survived (§6, Mi2, Mi3) |
| `cargo fmt --check` on the five changed `.rs` files outside `penalised_fit.rs` | 1 | every reported diff is in the review's own `probe.rs`; the reviewed files are clean |

Not run: `cargo clippy` and the full suite (the implementer's logs `tmp/d8/clippy.log`, `tmp/d8/test_full.log` report them clean apart from the three re-recorded checksums and the three debug-only panic tests); the reference fixture; anything at 2,000 samples.

Findings labelled "Needs verification": 0.

## 4. Open questions and assumptions

None that block. One assumption, used in Mi1: the varying cohort's only repeat tract is one locus. `src/cli/test_fixtures.rs:126-136` builds a 600-base contig with one tract (`TRACT = (200, 220)`), and the 21-sample cohort uses the same contig.

## 5. Top 3 priorities

1. **Mi1** — the digest comment gives the wrong reason the fixture's row carries no error. Fix the sentence; it will mislead whoever re-records next.
2. **Mi2** — the share rule's boundary (share × loci exactly 1) is not tested; a `<` → `<=` mutation survives.
3. **Mi3** — no run-log test has a class at exactly 3 loci; a `<` → `<=` mutation in the "fewer than 3" count survives.

## 6. Findings

### What was checked and holds

**The rule matches spec §7 (both 2026-10-16 entries), branch by branch.** I pasted HEAD's `errors_of_a_class` into a probe (`tmp/review_d8/probe.rs`, `head_errors_of_a_class`). I built the expected answer as HEAD with two edits: the share line set to 1, and the shorter share and fall-off at 3–7 loci judged as HEAD judges them at 8 or more. The grid covered:
- loci 0 to 20;
- slipped reads and slipped reads with a direction, each from {absent, 0, 0.5, 1, 5};
- errors absent, 0, NaN, ∞, 0.01 and 10⁶, and the fit giving no errors at all;
- prior weight 0, 1 and 3;
- concentration present and absent;
- shares including 0, NaN and exactly one locus's worth (0.125 × 8).

All 37,800 cases agree. So, against HEAD:
- every `NotPlaced` path (absent concentration, slipped reads under the prior, an absent, zero, non-finite or too-wide error) is unchanged;
- F's error, which comes from `errors.concentration` whether F was read off pairs or in the counting, still has its line at 8;
- the slip rate keeps 8, and the slipped-reads rule still comes after the loci rule for every stutter number.

The `UndeterminedErrors` paths reach `errors_of_a_class` as `None` errors, so they fall under the "absent error → `NotPlaced`" row. That is the one visible change from D8 on those paths: a class of 3–7 loci whose `H` failed now shows its shorter share and fall-off as `NotPlaced` rather than `TooFewLoci`. This follows the table's order and is correct.

**The order is sensible.** The loci rules come before the slipped-reads rule, so a number below its line is always `TooFewLoci`, as before. The share line applies only from 8 loci, because the `too_few_loci ||` short-circuit at `penalised_fit.rs:1148` catches thinner classes first. That keeps the share rule's meaning when its constant drops from 3 to 1.

**The C9 median still uses 8.** `stratum_fits.rs:1290` filters on `LOCI_FOR_AN_ERROR`, which is still 8, and the line is unchanged. Its test (`stratum_fits.rs:2068`) passes.

**The run-log counts are right.** `review_d8_withheld_line_counts_match_the_rule` fits four drawn cohorts with classes of 1, 2, 3, 4, 6, 7, 8, 9, 12, 20, 25, 30 and 40 loci. For each, it recomputes the withheld counts from the fitted values by the owner's rule, not by counting `TooFewLoci`:
- below 8 loci: 1 for the slip rate, plus the number of shares, plus 1 if F is present, plus 2 if the class has fewer than 3 loci;
- from 8 loci: the shares with share × loci < 1.

All four lines match in total, in the thin-class part, in the "in K of them with fewer than 3" clause, and in the rare-share part. One example line (seed 100): `32 in all; 32 in the 4 classes with fewer than 8 loci with reads (each one's slip rate, F and shares of allele lengths, and, in 2 of them with fewer than 3, the shorter share and fall-off too)`.

**The parameters file's notes say what the code does.** The own-fit-error paragraph says no error "below 8 loci" for the first number and "below 3" for the other two. Both numbers come from the constants through `format!` (`to_toml.rs:340-343`). The fallback-F paragraph says 8 is where the between-locus share and the slip share lose their error, which is true. The file does not carry the loci count, as both notes say: only test fixtures set `tracts_of_its_own` in `parameters_file/`.

**The digests moved by note text alone — verified independently.** I did not rerun `undo_notes.py`, which edits the file text pulled from a failing test's log. Instead I built a tree with D8's code, HEAD's `to_toml.rs`, and the three HEAD checksums restored (`tmp/review_d8/src_undo`). All 4 digest tests pass there. So D8's code, with HEAD's notes, writes byte for byte the files HEAD wrote.

**The oracle-psp comparison.** `tmp/d8/compare_head_vs_d8.txt` lists exactly 7 classes that gained `shorter_share_own_fit_standard_error` and `fall_off_own_fit_standard_error`, and "lines differing beyond those keys: 0". `diff` of the two files shows 11 changed lines outside the comments. Four of them are period-2 rows that `diff` printed as one block, so the count is consistent.

**The docs agree with the code and with each other.** The following all state the same two lines, the share line of 1, and the same order as `errors_of_a_class`:
- spec §2.6 (its table and the "Below 3 loci no number gets an error" sentence);
- spec §2.7's row and §2.10.4's bullet;
- §7's annotations;
- arch §1.1 (constants with values 1.0 / 8 / 3) and §1.4 (enum doc and table);
- companion §3.5.

Spec §2.10.4's bullet names only the stutter numbers' lines, which is right for a bullet about the stutter prior.

**Mutation runs (8; 6 killed, 2 survived, 0 changed-no-behaviour):**

| mutation | killed by |
|---|---|
| M1 direction numbers back on the line of 8 | 6 tests |
| M2 share rule ignores the class's 8-loci line | 5 tests |
| M3 slipped-reads rule before the loci rule | `a_direction_number_resting_on_the_prior_is_not_placed`, `the_lines_of_loci_fall_where_the_owner_set_them` |
| M4 share `<` → `<=` | **survived** (Mi2) |
| M5 run log "fewer than 3" `<` → `<=` | **survived** (Mi3) |
| M6 F's error given at any loci count | 5 tests |
| M7 direction line `<` → `<=` | 2 tests |
| M8 run log drops the "fewer than 3" clause | `the_withheld_errors_line_counts_the_classes_below_three_loci` |

Both survivors change behaviour. M4 changes the answer at share × loci exactly 1; the review's grid includes 0.125 at 8 loci. M5 changes the clause's count for any class of exactly 3 loci; probe cohorts seed 100 and seed 200 each have one.

### Blocker

None.

### Major

None.

### Minor

- `src/cli/cross_platform_digests.rs:160-165` — **[Minor]** Mi1: The digest note gives the wrong reason the row carries no error
- **Confidence:** High
- **Assumptions:** the varying cohort's one repeat tract is one locus (§4).
- **Problem:** The note says no number or key moved because "the file's one repeat-tract row has 0.03 expected slipped reads, under the stutter prior's one read, so it carries no error before or after". The slipped-reads rule is never reached for that row. The fixture's 600-base contig holds a single tract, so the class has 1 locus with reads. That is below both lines: below 8 at HEAD, and below 3 after D8. Every error is `TooFewLoci` before the slipped-reads rule is consulted (`penalised_fit.rs:1113-1114`). The same wrong reason is repeated for the 21-sample file at `:324` ("this tract has 0.0023 expected slipped reads"). The conclusion is right; the mechanism is not.
- **Why it matters:** The next person to re-record will look for a slipped-reads effect that does not occur. They may also believe a better-slipping fixture would exercise the error keys, when the fixture's single locus stops that whatever the slips.
- **Suggested fix:**
  ```rust
  /// moved: the file's one repeat-tract row is one tract, a class of 1 locus with reads, which is
  /// below both lines of loci (8 before D8, 3 for the direction numbers after), so it carries no
  /// error before or after; putting the two paragraphs' earlier wording back into the new file
  /// gives `b673dc9611c7126b5c35e9d7629493b3` exactly. ...
  ```
  Drop the "(this tract has 0.0023 expected slipped reads)" parenthesis at `:324`, or replace it with "one locus, as above".

- `src/parameter_estimation/joint/ssr_fit/penalised_fit.rs:3950-3962` — **[Minor]** Mi2: The share line's boundary is not pinned
- **Confidence:** High
- **Problem:** `the_lines_of_loci_fall_where_the_owner_set_them` tests shares worth 0.9 and 1.1 loci, not exactly 1. Mutation M4 (`<` → `<=` at `:1148`) passes every test in the module. The spec says "less than one locus" (§2.6) and the run log says "expected at less than one locus", so a share worth exactly 1 locus must keep its error.
- **Why it matters:** Each class is a small integer count of loci, and shares such as 1/8 or 0.5 are exactly representable, so the boundary can be reached in practice. Nothing stops a later edit from moving it silently.
- **Suggested fix:** in the same test, add a share worth exactly one locus (0.125 × 8 = 1.0 exactly in `f64`):
  ```rust
  // 8 loci: a share of 0.125 is expected at exactly one locus and keeps its error.
  let mut at_one = a_fitted_class(Some(2.0));
  at_one.allele_length_variation = vec![0.125, 0.875];
  at_one.errors = Some(PenalisedClassErrors {
      allele_length_variation: vec![0.01, 0.01],
      ..a_fitted_class(Some(2.0)).errors.expect("every error given")
  });
  let errors = errors_of_a_class(&at_one, Some(2.0), 8, STUTTER_PRIOR_READS);
  assert_eq!(errors.allele_length_variation, vec![Estimated(0.01), Estimated(0.01)]);
  ```

- `src/parameter_estimation/joint/ssr_fit/penalised_fit/run_log.rs:845-875` — **[Minor]** Mi3: No run-log test has a class at exactly 3 loci
- **Confidence:** High
- **Problem:** `the_withheld_errors_line_counts_the_classes_below_three_loci` uses classes of 40, 2, 5 and 40 loci. Mutation M5 (`fit.tracts_of_its_own < LOCI_FOR_A_DIRECTION_ERROR` → `<=` at `run_log.rs:590`) passes every test. That mutant would report a class of 3 loci as one "with fewer than 3", while its shorter share and fall-off were not withheld.
- **Why it matters:** The clause's count and the counted errors would then disagree in the log, and no test would notice.
- **Suggested fix:** add a class of exactly 3 loci to the fixture, e.g. `a_drawn_class(2, 12, 3, 12, 8, 45)`. Add `withheld_in(&outcomes[4], 1)` to `in_thin_classes`, and assert "in the 3 classes with fewer than 8 … in 1 of them with fewer than 3". The review's probe shows drawn classes of 3 loci at 12 samples × 8 reads get `Estimated` direction errors: seed 100 got 0.031 and 0.071, seed 200 got 0.045 and 0.100.

### Nits

- `penalised_fit.rs:3624-3626` — the class of 5 loci asserts `assert_ne!(error, StratumError::TooFewLoci)` for the shorter share and fall-off. The probe shows both come out `Estimated` (0.043 and 0.098). `NotPlaced` would also pass, so a mutant that sent direction numbers at 3–7 loci to `NotPlaced` would survive this test, though `the_error_mapping_follows_its_table` would catch it. Tighten to `assert!(matches!(error, StratumError::Estimated(_)), "{error:?}")`.
- `penalised_fit.rs:3605` (fixture; and `run_log.rs:699`) — the class of 8 loci (seed 35) passes on one share of 0.1113, worth 0.891 loci: a margin of 0.11 loci. The run-log fixture's class of 8 (seed 38) has two shares worth 0.646 and 0.624 loci. Both seeds were chosen by scanning. A change to the fit's arithmetic can tip the first. A one-line comment naming the margin would save the next person a scan.
- `ssr_fit.rs:597-598` — "held the truth about as often as honest errors do". Give the numbers, as the constant's doc does: 89 to 93 times in 100 against 95 (CLAUDE.md, "never assert a property without its size").
- `stratum_fits.rs:1278-1281` — the edit left a short line ("…allele-length shares. A thin stratum's") followed by the old wrap. Re-wrap the paragraph.
- `run_log.rs:606-614` — "(each one's slip rate, F and shares of allele lengths…)": F is withheld only where the class has one, since an absent F's error is `NotPlaced`, not counted. The counts are right; the wording slightly over-claims. HEAD's "every number fitted" had the same looseness.

## 7. Out of scope observations

- `src/cli/cross_platform_digests.rs` — because the varying cohort's only stratum is one locus, `a_fit_with_standard_errors_writes_the_same_bytes_on_every_platform` and `errors_computed_a_block_at_a_time_write_the_same_bytes_on_every_platform` cannot pin any repeat-tract own-fit error in the written file. Every one is `TooFewLoci` at 1 locus. The repeat-tract half of those checksums covers values only, not errors. Follow-up: a fixture with at least 3 tracts of one class (8 for the slip rate) would put direction errors into the file.
- `run_log.rs:211-233` (pre-existing) — `classes_whose_directions_rest_on_the_prior` counts classes under the prior's weight whatever their loci count. A class under 3 loci (under 8 before D8) has its direction numbers withheld for too few loci, not for the prior, yet it is counted here. The log's sentence ("which carry no error") stays true. D8 narrowed the overlap.

## 8. Missing tests to add now

- `errors_of_a_class_keeps_a_share_expected_at_exactly_one_locus` — boundary of `LOCI_FOR_A_SHARE_ERROR`; catches `<` → `<=`. Code in Mi2.
- `withheld_errors_line_does_not_count_a_class_of_three_loci_as_below_three` — boundary of the run log's direction-line clause; catches `<` → `<=` at `run_log.rs:590`. Specification in Mi3.
- Optionally promote the review's grid (`tmp/review_d8/probe.rs`, `review_d8_grid_matches_head_but_for_the_owners_changes`) to a permanent table-driven test of `errors_of_a_class`. Its cost is under a millisecond, and it pins every branch's outcome at once.

## 9. What's good

- `to_toml.rs:340-343` builds the own-fit-error note with `format!` from `LOCI_FOR_AN_ERROR` and `LOCI_FOR_A_DIRECTION_ERROR`, so the file cannot state a different line from the one the fit applies.
- `penalised_fit.rs:1112-1121` gives `of_a_stutter_number` the line as an argument (`loci_needed`), which keeps the three stutter numbers in one closure with each line visible at the call site.
- `the_lines_of_loci_fall_where_the_owner_set_them` (`penalised_fit.rs:3912`) asserts the three constants' values first, so a change to the owner's numbers fails one named test rather than several indirect ones.
- `run_log.rs:585` counts `TooFewLoci` outcomes instead of recomputing the rule, so the log cannot disagree with what the file and calling receive.

## 10. Commands to re-verify

```
# from tmp/review_d8/src_copy (D8 tree plus the probe module appended to penalised_fit.rs)
CARGO_TARGET_DIR=$PWD/../target cargo test --release -j 8 --lib review_d8 -- --nocapture
CARGO_TARGET_DIR=$PWD/../target cargo test --release -j 8 --lib -- penalised_fit stratum_fits parameters_file cross_platform_digests ssr_fit
# from tmp/review_d8/src_undo (D8 code, HEAD's to_toml.rs and HEAD's three fitted checksums)
CARGO_TARGET_DIR=$PWD/../target cargo test --release -j 8 --lib cross_platform_digests
# mutations
tmp/review_d8/mut/run.sh
```

`PROJECT_STATUS.md` was not updated: the review was dispatched read-only on the worktree apart from this report.
