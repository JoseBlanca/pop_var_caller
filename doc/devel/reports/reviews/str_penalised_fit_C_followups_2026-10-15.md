# Code Review: str_penalised_fit_C_followups
**Date:** 2026-10-15
**Reviewer:** rust-code-review skill (orchestrator, single reviewer: no sub-agents were spawned, at the lead's instruction)
**Scope:** the four follow-up commits of milestone C on `str-penalised-fit`, read together: `git diff 0b12941d^..53a3e4f9`
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the combined diff of `0b12941d` (milestone C review fixes: errors that
  could not be determined are left out, and a fit that could not be factorised is never chosen),
  `fdfbd0ad` (the owner's answers to C8's review, docs only), `adbc8bd0` (C8, the stutter prior),
  `2cc85d75` (C9, the stated F as a median over strata with 8 or more loci) and `53a3e4f9` (C10,
  the cohort fitted in one fixed order). Each was reviewed on its own before
  (`tmp/C_fix_errors_review.md`, `tmp/C8_review.md`, `tmp/C9_review.md`, `tmp/C10_review.md`).
  This review looks for what shows only when they are read together.
- **Reviewed against:** commit `53a3e4f9`, exported with `git archive` to
  `tmp/review_C_followups/src_copy/` (the working tree is being edited by other agents).
- **In-scope files:** `src/parameter_estimation/joint/ssr_fit/penalised_fit.rs`,
  `penalised_fit/run_log.rs`, `src/parameter_estimation/joint/stratum_fits.rs`,
  `src/parameter_estimation/joint/census.rs`, `src/cli/estimate_parameters.rs` and its `tests.rs`,
  `src/cli/test_fixtures.rs`, `src/run/psp_caller.rs`, `src/calling/parameters_file/to_toml.rs` and
  its golden file, `src/calling/likelihood/stutter_rates.rs`, `examples/ng_str_fit_bench.rs`,
  `scripts/research/str_fit_fixture.py`; spec `doc/devel/ng/spec/str_penalised_fit.md`, arch
  `doc/devel/ng/arch/str_penalised_fit.md`, plan `doc/devel/ng/impl_plan/str_penalised_fit.md`.
- **Out of scope:** the known defect that the error matrix leaves out the allele-length prior's
  curvature (`penalised_fit.rs` ~1428; diagnosed in
  `doc/devel/reports/implementations/str_penalised_fit_thin_errors_diagnosis_2026-10-15.md`), and
  the stale comment at `ssr_fit.rs` ~1700. Finding M1 says how the combined diff interacts with the
  first.
- **Categories covered** (by one reviewer, not dispatched): reliability (tests that still test what
  they say), errors, refactor safety, smells, float portability (order of summation), and
  "diff matches stated intent" for the docs. Naming, idiomatic, tooling and module structure were
  read for but gave nothing new beyond the single reviews.

## 2. Verdict

**Approve-with-changes.** The four steps fit together: C10's reordering touches nothing the fit
reads by position except through the evidence it reorders, C9's median does not depend on order,
and C8's prior and `0b12941d`'s error handling compose as the spec says (761 tests pass). Two
things need attention: three tests and two pieces of user-facing text pin the known
allele-length defect as expected behaviour and name the wrong cause (M1), and the rule for the
fall-off's error does not match what the parameters file and the owner's decision say (Mi1).

## 3. Execution status

- `cargo test --release -j 6 --lib -- ssr_fit stratum_fits estimate_parameters census parameters_file stutter_rates psp_caller --test-threads 6`
  on the exported commit plus one probe test: exit 0, `test result: ok. 761 passed; 0 failed; 7 ignored`
  (`tmp/review_C_followups/tests_scoped.log`).
- Probe 1, the fall-off rule (`review_probe_fall_off_error_below_one_slipped_read`, appended to the
  copy's test module): `tmp/review_C_followups/probe1.log`, `probe2.log`, `probe3.log`.
- Probe 2, the stiff three-base cohort's error matrix: the diagnosis report's dump block inserted in
  the copy, the test `errors_that_need_a_ridge_are_left_undetermined` run with `PVC_DUMP_H`, and the
  diagnosis script `tmp/diag_errors/remedy.py diag` run on the dump
  (`tmp/review_C_followups/dump/stiff.p3.json`).
- **Not run:** `cargo fmt --check`, `cargo clippy`, `cargo doc`, `cargo audit`, the full test suite,
  the bench grid and the oracle. Each commit's report records clippy and its oracle; nothing in this
  review depends on them.
- Findings labelled "Needs verification": none.

## 4. Open questions and assumptions

1. **Which rule withholds the fall-off's error** (Mi1): fewer *repeats slipped* than the prior's
   weight, as the code and the first sentence of spec §2.10.4 say, or fewer *slipped reads*, as the
   owner's decision (§7), the defaults table (§2.7), §2.10.5 and the parameters file's note say?
   *Recommended:* fewer slipped reads, for all three stutter numbers. It is one rule a reader can
   check against the row's `expected_slipped_reads`, which is the one count the file writes.
2. **When the allele-length curvature fix lands** (M1), should the tests of the "errors left
   undetermined" path keep a real cohort that still fails, or use the existing test seam? I assume
   the seam, since the diagnosis found no real failure that the fix does not remove.

## 5. Top 3 priorities

1. **M1** — rewrite the doc of `a_stiff_three_base_cohort` and the cause given in `UndeterminedErrors`
   and in the run-log line now, and plan the two tests' replacement with the curvature fix: they
   assert the defect.
2. **Mi1** — settle question 1 and make the code, the file's note and the spec say the same.
3. **Nit, the arch's stale pointers** — cheap, and the arch is what D-step agents read first.

## 6. Findings

### Major

- `src/parameter_estimation/joint/ssr_fit/penalised_fit.rs:3774` (`a_stiff_three_base_cohort`), with the tests at `:3827` (`errors_that_need_a_ridge_are_left_undetermined`) and `penalised_fit/run_log.rs:1136` (`undetermined_errors_are_logged_and_written_as_not_placed`); the text at `run_log.rs:429` (`undetermined_errors_line`) — **[Major]** M1: Three tests and two texts pin the known allele-length defect as the expected behaviour, and name the wrong cause
- **Confidence:** High
- **Categories:** reliability, refactor safety, "diff matches stated intent"
- **Problem.** `0b12941d` added a test cohort, "a stiff three-base cohort whose curvature cannot be
  inverted". Its doc says the cause is the stutter curves: they come out near a strength of 10¹⁰,
  "so each is close to a straight line, and with so few slipped reads some such line has next to
  no reads behind it". Two tests assert that this cohort's errors are left undetermined. C8 then
  built on it (the run-log test reuses the cohort). I dumped this cohort's error matrix H, at the
  committed code, and ran the diagnosis report's script on it:
  - the factorisation fails at the last row, the share two repeats longer at 11 repeats, pivot −0.90;
  - the one negative direction of H (eigenvalue −0.143) is the trade between the shares one and two
    repeats longer, at 11, 6 and 10 repeats;
  - **the 18 stutter numbers alone have positive curvature**, smallest eigenvalue 0.044;
  - **with the allele-length prior's curvature added, H inverts.** Its smallest eigenvalue is then
    0.044, on the fall-off.

  So the cohort fails for the known reason (the allele-length prior left out of H), not the one its
  doc gives. The diagnosis reached the same answer on 54 of 56 failing bench runs.
- **Why it matters.** When the owner approves the curvature fix,
  `errors_that_need_a_ridge_are_left_undetermined` and
  `undetermined_errors_are_logged_and_written_as_not_placed` will fail. Their docs will then send
  whoever fixes them looking for a stiff stutter line that is not there, or invite changing the
  cohort until it fails again. Meanwhile the doc of `UndeterminedErrors` and the run-log line
  (`undetermined_errors_line`, printed in every user's log when it fires) offer "the slope of a stiff
  curve through a single repeat count's slipped reads" as the example cause; on the cases measured
  the cause is a trade between two rare allele lengths. The line's second clause ("the curvature is
  not positive in some direction") is the true one.
- **Suggested fix.**
  1. Now: rewrite the cohort's doc to say what was measured: H fails along a trade between two
     allele-length shares at the edge of the window, because H leaves out the allele-length prior;
     the stutter numbers alone are well determined.
  2. Now: in `UndeterminedErrors`' doc and in `undetermined_errors_line`, replace the stiff-curve
     example with "a trade between two rare allele lengths at a thin class, for one". (The diagnosis
     report asks for the doc change too.)
  3. With the curvature fix: invert the two tests' expectation on this cohort (errors determined;
     slip-rate logit errors finite and positive), and test the "undetermined" path through the
     existing seam (`failing_factorisation_above` is for Newton; add a matching seam for
     `penalised_sandwich`, or call it with a hand-built non-positive matrix as
     `the_sandwich_tries_no_ridge_larger_than_the_fits_own` already does).

### Minor

- `src/parameter_estimation/joint/ssr_fit/penalised_fit.rs:1016` (the `fall_off` arm of `errors_of_a_class`, `:982`), with `src/calling/parameters_file/to_toml.rs:348` and `penalised_fit/run_log.rs:219` (`classes_whose_directions_rest_on_the_prior`) — **[Minor]** Mi1: A row whose expected slipped reads are below one can still carry a fall-off error, which the file's note and the owner's decision say it cannot
- **Confidence:** High
- **Problem.** The code withholds the fall-off's error when the class's **repeats slipped** are fewer
  than the prior's weight. A slipped read slips one repeat or more, so the repeats slipped can reach
  one while the slipped reads are still below one. Probe, at the committed code: the one-class motif
  length of the C8 tests (10 loci, one sample, 3 reads each), with one read two repeats short, at a
  homozygote excess of 0.95 (0.6 to 0.9 give the same):

  | | value |
  |---|---|
  | expected slipped reads (written as `expected_slipped_reads`) | 0.99999986 |
  | slips with a direction | 0.99999986 |
  | repeats slipped | 3.11 |
  | slip rate, shorter share | 0.035, 0.75: errors `NotPlaced` |
  | fall-off | **0.525, error `Estimated(0.345)`** |

  Three texts say this row has no error on any of its three numbers:
  - the parameters file's note (C8): "a row with `expected_slipped_reads` below one has … no error on
    any of its three numbers";
  - spec §2.10.5 ("by §2.10.4, no error on any of its three stutter numbers"), the defaults table of
    §2.7 ("a slip rate, shorter share or fall-off under one expected slipped read") and the owner's
    answer in §7 ("no error on a shorter share or fall-off whose class has fewer than one expected
    slipped read");
  - the run log counts this class and says its "shorter share and fall-off come from their
    neighbours on the curve and the stutter prior, not from their own reads". Its fall-off of 0.525
    came from the read: the prior's centre is 0.05, and the class has no neighbours.

  Only the first sentence of spec §2.10.4 and the arch's §1.4 row say "repeats slipped", which is
  what the code does. The spec's parenthesis there ("its repeats slipped, never fewer, the
  fall-off's") reads the inequality the wrong way round for the conclusion it draws.
- **Why it matters.** A reader who trusts the note and treats every row with
  `expected_slipped_reads` below one as prior-led will read an `Estimated` fall-off error on such a
  row as a contradiction, or worse, will not look. Where it happens: two- and three-base classes
  whose fall-off is high (a mean slip of 1.4 repeats at a fall-off of 0.3), with 0.7 to 1 slipped
  read. That is a thin class at any cohort size, the corner where the prior matters most.
- **Suggested fix.** Depends on question 1.
  - If the owner's wording stands (recommended): judge all three numbers on `expected_slipped_reads`.
    ```rust
    let few_slipped = fitted.expected_slipped_reads.unwrap_or(0.0) < stutter_prior_reads;
    // level, shorter_share and fall_off: NotPlaced when !too_few_loci && few_slipped
    ```
    `expected_slips_with_a_direction` and `expected_repeats_slipped` then have no reader and can go,
    and `classes_whose_directions_rest_on_the_prior` counts on the slipped reads too. Add the probe
    above as a test: a lone class with one read two repeats short, homozygote excess 0.9, all three
    errors `NotPlaced`.
  - If the code's rule stands: correct the file note, §2.7, §2.10.5 and §7 to "a fall-off whose
    class's repeats slipped are below one", and the run-log clause to name the fall-off only where
    it is withheld.

## Nits

- `doc/devel/ng/arch/str_penalised_fit.md` §1.1a cites `opening_line (:67)` and
  `how_the_fit_ran (:156)`; they are at `run_log.rs:79` and `:235`. It also names a closure
  `of_a_direction`; the code's is `of_a_stutter_number`. The design citations to `penalised_fit.rs`
  line numbers in spec §2.10.4 (`:1384-1408`, `:2301-2334`, `:2180`) predate the code by about 240
  lines. Refresh them once, before milestone D's agents read them.
- `src/calling/parameters_file/to_toml.rs:348`: the stutter-prior sentence writes "one extra read's
  worth" and the centre's three numbers out by hand, while C9's paragraph beside it builds its
  threshold from `LOCI_FOR_AN_ERROR`. Build the stutter sentence from `STUTTER_PRIOR_READS` and
  `STUTTER_PRIOR_CENTRE` the same way, so a change of either cannot leave the file's note wrong.
- `doc/devel/ng/impl_plan/str_penalised_fit.md:36-39` and `:319`: the reflowed "Isolate the silent steps"
  bullet has a line well past the file's width; two blank lines follow "Every D step runs after
  C10".

## 7. Out of scope observations

- **C10 makes the output one file, but the amplification behind it is still there.** The C10 commit
  records that a difference in the last digits of the SNP/indel fit became a two-base fall-off
  strength of 5.7×10⁹ against 1.2 in the repeat-tract fit, because the marginal likelihood is almost
  flat near the top of the strength search. Sorting the cohort removes one source of such
  differences; adding or removing one sample, or a change in how a sum is taken, will move the
  strength as far again. Worth a line in D1 or D3: how much the *fitted values* (not the strengths)
  move between two strengths that close in score.
- **The thresholds of C8 sit on a knife's edge for one read**, as any threshold must: in the probe,
  the same one slipped read gives 1.00002 expected slipped reads at a homozygote excess of 0.98 (all
  three errors `Estimated`) and 0.99999986 at 0.95 (two or three `NotPlaced`). Not a defect; the
  owner may want to know the rule moves on that little.

## 8. Missing tests to add now

Grouped by the function under test.

- `errors_of_a_class`:
  `errors_of_a_class_withholds_fall_off_below_one_slipped_read_even_when_repeats_slipped_reach_one`
  — a lone class of 10 loci at 3 reads with one read two repeats short, homozygote excess 0.9; after
  question 1 is settled, all three stutter errors `NotPlaced` (or, if the code's rule stands, the
  fall-off `Estimated` and the note's wording changed). Catches Mi1.
- `standard_errors` (with the curvature fix): `a_stiff_three_base_cohort_has_its_errors_determined`
  — the cohort of M1, `undetermined_errors` `None`, every class with a finite positive slip-rate
  error. Replaces the current assertion that it fails.
- `penalised_sandwich` / `undetermined_errors_line`: keep the "undetermined" path covered through a
  hand-built non-positive matrix, independent of any cohort.

## 9. What's good

- `SampleCensusEvidence::with_read_groups_in_declared_id_order` (`census.rs:1680`) moves sections,
  declarations and claimed-error totals with one map, so the evidence cannot be relabelled without
  moving with its name, and its test checks exactly that.
- `fit_robustly_at` scores an unfactorised fit as minus infinity rather than as a number
  (`penalised_fit.rs`, `marginal = ln_determinant.map_or(f64::NEG_INFINITY, …)`): an absent answer
  is typed as absent, and each count of what was rejected reaches the log.
- C9's note builds its threshold from the constant (`to_toml.rs:333-336`), so the file cannot state a
  different count from the one the fit used.
- Weight 0 reproduces the fit without the prior by construction: the prior enters only as
  `+ 0.0` to counts and penalty and a `< 0.0` test that never fires.

## 10. Commands to re-verify

```
mkdir -p tmp/review_C_followups/src_copy && git archive 53a3e4f9 | tar -x -C tmp/review_C_followups/src_copy
cd tmp/review_C_followups/src_copy
CARGO_TARGET_DIR=$PWD/../target RAYON_NUM_THREADS=6 cargo test --release -j 6 --lib -- \
  ssr_fit stratum_fits estimate_parameters census parameters_file stutter_rates psp_caller --test-threads 6
# Probe 1 (the copy already carries it): review_probe_fall_off_error_below_one_slipped_read
CARGO_TARGET_DIR=$PWD/../target cargo test --release -j 6 --lib -- review_probe --nocapture --test-threads 1
# Probe 2 (the copy carries the dump block from tmp/diag_errors):
PVC_DUMP_H=$PWD/../dump/stiff CARGO_TARGET_DIR=$PWD/../target cargo test --release -j 6 --lib -- \
  errors_that_need_a_ridge_are_left_undetermined --test-threads 1
cd ../../..
uv run --no-project --with numpy python -I tmp/diag_errors/remedy.py diag tmp/review_C_followups/dump/stiff.p3.json
```

`PROJECT_STATUS.md` has no block for this feature and was not edited; the lead can add one.

### Author response convention

Address each finding by its identifier (M1, Mi1) with `fixed in <commit>` / `disputed because …` /
`deferred to …` / `won't fix because …`. Answer the open questions in section 4 first.
