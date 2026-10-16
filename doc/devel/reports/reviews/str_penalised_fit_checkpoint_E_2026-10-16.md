# Code Review: str_penalised_fit checkpoint E
**Date:** 2026-10-16
**Reviewer:** rust-code-review skill (orchestrator; every category covered by the orchestrator itself, no sub-agents, at the caller's request)
**Scope:** the commits since checkpoint D read together — `7cbc5f78..2276c93c` (D7, E1, E2 and its doc notes, D8, the curve-code deletion) — and the merge of main, `909266ee`
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** a branch diff, read as one change. Each step already has its own review
  (`str_penalised_fit_{D7,E1,D8}_2026-10-16.md`); this one looks for what shows only when the steps
  are read together: interactions between D7's errors and D8's thresholds, code the combination left
  dead, deleted code still named, tests that no longer test what they say, and whether the spec,
  companion, architecture and plan match the final code.
- **Reviewed against:** `909266ee` (branch `str-penalised-fit`), exported with `git archive` to
  `tmp/review_E/src_copy/`; the worktree was not touched except for this file.
- **In-scope files:** the 69 files of `git diff --stat 7cbc5f78..2276c93c`, chiefly
  `src/parameter_estimation/joint/{ssr_fit.rs, ssr_fit/penalised_fit.rs, ssr_fit/penalised_fit/run_log.rs, stratum_fits.rs, slippage_curve.rs, share_curve.rs}`,
  `src/calling/parameters_file/{to_toml.rs, from_run_parameters.rs, to_run_parameters.rs}`,
  `src/calling/inference/repeat_tract_parameters.rs`, `src/cli/cross_platform_digests.rs`, the
  examples and scripts the deletion touched; `doc/devel/ng/spec/str_penalised_fit.md`,
  `doc/devel/ng/spec/str_fit_to_calling.md`, `doc/devel/ng/arch/str_penalised_fit.md`,
  `doc/devel/ng/impl_plan/str_penalised_fit.md`. For the merge: the 16 files main changed in
  `e4bd1e94..98230b01`.
- **Deliberately out of scope:** the arithmetic of D7's block and D8's rule table, each verified by
  its own review with mutations; the research reports' numbers (quoted, not re-derived, except where
  §8a says).
- **Categories covered:** reliability, errors, naming, defaults, idiomatic, refactor_safety,
  module_structure, float_portability (D7 adds arithmetic that reaches the file's errors), smells,
  tooling, extras (diff matches stated intent). `unsafe_concurrency`: nothing in scope uses `unsafe`,
  locks or channels; the parallel sums were made ordered at C3 and D7 adds none (checked).

## 2. Verdict

**Approve-with-changes.** No blocker. The full suite passes except the three debug-only panic tests
known to fail in release; clippy and fmt are clean. D7's errors and D8's thresholds fit together: the
coverage D8 rests on was measured with D7's errors, and D7 narrows the errors D8 newly releases by a
median of under 3 in 100. The merge of main has no semantic conflict. What is left: the oracle
baseline still predates D7 and D8 (M1), and the deletion left one error variant and the curves'
evaluation code with no production caller, documented as if reading needed them (Mi1, Mi2), plus
stale descriptions in docs and type comments (Mi3, Mi4).

## 3. Execution status

All on the exported copy, native `cargo` (macOS arm64), `CARGO_TARGET_DIR` under `tmp/review_E/`, 6
jobs and 6 test threads.

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | no output |
| `cargo clippy --all-targets --all-features -j 6 -- -D warnings` | 0 | `Finished`; no warning in `pop_var_caller` |
| `cargo test --release --all-features --no-fail-fast` | 101 | **5,244 passed, 3 failed, 7 ignored**. The three: `alignment::left_align_repeated::tests::an_offset_past_the_reference_panics`, `alignment::left_align_structured::tests::an_offset_past_the_reference_panics`, `genetics::tests::lgamma_panics_on_non_positive_in_debug` — the debug-only panic tests, as in every recent report |
| `cargo doc --no-deps` | 101 | `error: unresolved link to super::cohort_merge::LocusSummary::over_depth_ceiling` at `src/run/depth_ceiling.rs:14:47`. **Pre-existing on main**: the file and `src/run/cohort_merge/` are identical at `98230b01` and `909266ee`, and `38487520`, which last touched it, is on main. Out of scope (§7) |
| `cargo check --lib` and `--all-targets`, default features | 0 | no warning from the library; the only project warnings are in `examples/dhat_ng_joint_fit.rs` and `examples/dhat_ng_merge.rs` without the `dhat-heap` feature (§7) |
| mutation: `SlippageCurve::level_at` and `ShareCurve::share_at` panic unless an env var is set; `cargo test --release --all-features --lib --tests` | 101 | 5,233 passed, 13 failed: the three above, 9 tests of `slippage_curve.rs`/`share_curve.rs` about the evaluation itself, and `stratum_fits::tests::the_level_is_the_one_the_fit_emitted_and_not_the_curves_own`, which uses `level_at` as its yardstick. **No calling or reading test reaches either function** (Mi2) |

Not run: `cargo audit` (not installed in this setup); the oracle (`scripts/promote_ng_oracle.sh`),
which needs the tomato data and was being run by another agent; nothing at 2,000 samples.

Findings labelled "Needs verification": 0.

## 4. Open questions and assumptions

1. Is the oracle run in progress meant to re-record the baseline after D7 and D8 (plan step D8's
   last sentence)? If so, M1 closes with it. (M1)
2. Should the curves' evaluation code go, or stay with an honest comment? The owner chose deletion
   for the curve fitting on 2026-10-16; the reviewer recommends the same here. (Mi2)

## 5. Top 3 priorities

1. **M1** — re-record `scripts/promote_ng_oracle.baseline` after D7 and D8, as plan step D8
   requires; until then checkpoint E's oracle reports `DIFFERS`.
2. **Mi2** — delete `level_at`/`share_at` and what only they use, or stop calling them "what reading
   an older file needs"; nothing reads a curve back.
3. **Mi1** — delete `StratumError::NotIdentified`: since E1 nothing builds it, and its doc describes a
   case the penalised fit maps to `NotPlaced`.

## 6. Findings

### Blocker

None.

### Major

#### M1: [scripts/promote_ng_oracle.baseline](../../../../scripts/promote_ng_oracle.baseline) — the oracle baseline predates D7 and D8

- **Confidence:** High
- **Problem:** Plan step D8 ends "D4's re-recording is repeated after D7 and D8"
  ([impl_plan/str_penalised_fit.md:401-402](../../ng/impl_plan/str_penalised_fit.md#L401-L402)). The
  baseline was last changed by `c4a425a6` (D4), before D7: `git log -- scripts/promote_ng_oracle.baseline`
  lists nothing after it. Both later steps move the fitted file on the oracle's four psps: D8's own
  commit message says "on the oracle's four tomato psps the file gains 14 errors", and D7 moved
  59,443 of the reference fixture's 62,538 numeric error fields (recomputed here, §8a), which the
  file writes for every class of 8 or more loci. So `fitted.parameters.toml`'s line cannot match. The
  plan ticks D8 without mentioning the re-recording.
- **Why it matters:** checkpoint E is the merge gate; a branch merged with this baseline makes the
  next oracle run on main report a difference nobody explains, which is exactly what the baseline's
  header asks never to happen ("say in the commit which lines moved and why").
- **Suggested fix:**
  1. Run `scripts/promote_ng_oracle.sh` at `909266ee`.
  2. Expect `fitted.parameters.toml` to move (errors from D7; 14 new errors and two note paragraphs
     from D8) and nothing else — errors are written to the file only, and calling reads none of them
     (`own_fit_standard_errors` has no reader outside `parameters_file/` and its tests). Any other
     moved line is investigated first.
  3. Re-record with a header line naming both steps, as D4 did, and add one sentence to plan step D8.

### Minor

#### Mi1: [src/parameter_estimation/joint/ssr_fit.rs:574-576](../../../../src/parameter_estimation/joint/ssr_fit.rs#L574-L576) — `StratumError::NotIdentified` has had no producer since E1

- **Confidence:** High
- **Problem:** Before E1, only the climb's errors built it (`7cbc5f78:ssr_fit.rs:2237`, `:2240`, `:2252`,
  `:2302`). At `909266ee`, `git grep` finds no constructor in `src/`; the only other mentions are
  `.value()`'s match (`:624`) and an example's match arm
  ([examples/ng_str_counting_fit_errors.rs:384](../../../../examples/ng_str_counting_fit_errors.rs#L384)).
  Its doc says the variant means "flat, or a saddle", but the penalised fit sends that case to
  `NotPlaced` (`errors_of_a_class`, table row "the fit gave no error … not finite"). Arch §1.4 still
  lists it ([arch/str_penalised_fit.md:229](../../ng/arch/str_penalised_fit.md#L229)).
- **Why it matters:** a reader matching on `StratumError` handles a state that cannot occur, and the
  doc points them at the wrong variant for the case it describes. The type is not read from a file
  (the file holds `Option<f64>`), so nothing older needs it.
- **Suggested fix:** delete the variant, its arm in `value()` and the example's arm, and the line in
  arch §1.4. If kept for the research example, say in its doc that only `counting_fit_errors` could
  produce it and that it does not.

#### Mi2: [src/parameter_estimation/joint/slippage_curve.rs:127-174](../../../../src/parameter_estimation/joint/slippage_curve.rs#L127-L174), [share_curve.rs:140-171](../../../../src/parameter_estimation/joint/share_curve.rs#L140-L171) — the curves' evaluation code is dead in production, and documented as needed

- **Confidence:** High (verified by mutation, §3)
- **Problem:** E1 removed the climb's calls to `level_at`/`share_at` (`7cbc5f78:ssr_fit.rs:3673-3675`),
  and `2276c93c` removed the blend's. What is left — `SlippageCurve::level_at`, `level_on_the_line`,
  `RiseShape::is_multiplying`, `LEVEL_FLOOR`, `ShareCurve::share_at`, `share_on_the_curve`,
  `SHARE_FLOOR`, `SHARE_CEILING` — has no caller outside its own module's tests (plus `LEVEL_CEILING`
  in one test of `stutter_rates.rs`, and one `stratum_fits` test that uses `level_at` as a yardstick).
  Reading an older file takes the numbers the file wrote and calls only `reach`
  (`to_run_parameters.rs:1704`, `:1711`, `:1729`). Yet both module docs say they hold "what reading
  and calling such a file needs" and list the evaluation and the bounds
  ([slippage_curve.rs:10-18](../../../../src/parameter_estimation/joint/slippage_curve.rs#L10-L18),
  [share_curve.rs:16-18](../../../../src/parameter_estimation/joint/share_curve.rs#L16-L18)), and arch
  §4.1's last row keeps "the level and share floors and ceilings" for the same reason
  ([arch/str_penalised_fit.md:496](../../ng/arch/str_penalised_fit.md#L496)). With both functions made
  to panic, 5,233 tests pass and none of the 13 failures is a calling or file-reading test.
- **Why it matters:** about 150 lines of code and nine tests are kept on a justification that is not
  true, and `cargo` cannot flag them because they are `pub`. The owner's ruling on 2026-10-16 was to
  delete the curve code nothing in a run calls.
- **Suggested fix:** delete the items listed, their tests, and
  `the_level_is_the_one_the_fit_emitted_and_not_the_curves_own` (it guards against recomputing a level
  from the curve, which needs the function it would delete); correct the two module docs and arch
  §4.1 to "`SlippageCurve`/`ShareCurve` as written, `reach`, …". If the owner prefers to keep them,
  change the docs to "kept to read a written curve back; nothing in a run does".

#### Mi3: [src/parameter_estimation/joint/ssr_fit.rs:526-528](../../../../src/parameter_estimation/joint/ssr_fit.rs#L526-L528) and four more — type docs still describe the climb's blend as current

- **Confidence:** High
- **Problem:** the comments on types every fit and every file read pass through still describe the
  climb as what produces them:
  - `StratumFit::level_provenance` (`ssr_fit.rs:526-528`): "The level is the only one of the four
    numbers a curve supplies" — the penalised fit's level and both shares are all points on curves.
  - `SharesProvenance` (`:673-676`): "The two shares are smoothed exactly as the level is — each gets
    its period's curve, and each stratum departs from that curve by how much evidence it has".
  - `LevelProvenance::slipped_reads` (`:658-668`) and `SharesProvenance::slipped_reads` (`:684-691`):
    "None where the stratum has no level of its own because it borrowed", "the weight the blend gave
    that answer"; the penalised fit fills it from the E-step's expected slipped reads.
  - `Library::stutter_multiplier` (`:203-206`): "Held while a stratum's numbers are fitted; one until
    the fit sets it" — the fit that set it (`library_stutter.rs`) went with E1; nothing sets it now.
  - `warrant_of` ([repeat_tract_parameters.rs:1134](../../../../src/calling/inference/repeat_tract_parameters.rs#L1134)):
    "The climb reports three sources per number".
- **Why it matters:** these are the definitions a reader goes to when a row of the file looks odd;
  each one now explains the row by a mechanism that no longer runs.
- **Suggested fix:** one sentence each in the present tense for the penalised fit, and "an older file
  may carry …" for the climb's shapes. E.g. for `SharesProvenance::slipped_reads`: "The class's
  expected slipped reads from the fit's last E-step (`PenalisedClass::expected_slipped_reads`); in a
  file the climb wrote, `None` where the stratum borrowed."

#### Mi4: [doc/devel/ng/spec/str_penalised_fit.md:5](../../ng/spec/str_penalised_fit.md#L5) and [arch/str_penalised_fit.md:31, 167, 254, 474](../../ng/arch/str_penalised_fit.md#L167) — the spec's status line and the architecture predate E1

- **Confidence:** High
- **Problem:**
  - spec line 5: "A research prototype exists and is measured; nothing is wired into production" —
    `estimate-parameters` has run the penalised fit since C7;
  - arch line 31: "The climb's code in it is deleted at plan milestone E" (future tense);
  - arch §1.2's type block (line 167) still shows `pub climb: Option<ClimbRecord>` as a field, and
    the reuse table (line 474) lists it;
  - arch §1.4 (line 254) judges `NOT_PLACED` "as the climb judges it (`ssr_fit.rs:2172`)", a line
    E1 deleted;
  - plan step D7 (line 391) is still headed "Implemented (2026-10-16), for review", though reviewed
    and ticked.
  The arch header says its line numbers are at `e0b41ad4`, which covers the citations but not the
  claims about what exists.
- **Why it matters:** the next plan (`str_library_stutter.md`, rebuilt on this fit) starts from these
  documents; a reader of arch §1.2 looks for a `climb` field and a reader of the spec's status line
  thinks the fit is not in production.
- **Suggested fix:** a dated line in each header ("2026-10-16, plan E1/E2 done: the climb, its record,
  `NotIdentified`'s producer and the curve fitting are deleted (`3198d55d`, `2276c93c`); line numbers
  below are at `e0b41ad4`"), the status sentence replaced, and the `climb` field struck in §1.2 and
  the table.

### Nits

- [to_toml.rs:342](../../../../src/calling/parameters_file/to_toml.rs#L342) and `:367` — every new
  file's notes explain `samples_fitted_on` subsets and numbers "taken from a parametric curve, or
  blended", which only a file the climb wrote can contain. Say so: "(written only by the fit before
  2026-10-16)".
- [penalised_fit.rs:399-403](../../../../src/parameter_estimation/joint/ssr_fit/penalised_fit.rs#L399-L403),
  [ssr_fit.rs:591-596](../../../../src/parameter_estimation/joint/ssr_fit.rs#L591-L596) and spec
  §2.6 quote the direction numbers' coverage at 3 to 7 loci at 2,000 samples only (89 to 93 in 100).
  The same file the figures come from (`tmp/d7/coverage_by_number.md`, table (d)) gives 89 to 90 in
  100 over every cohort size, and D1's table shows the thinnest corners lower (fall-off at 1 sample,
  3–4 loci, 0.86; shorter share at 60 samples, 5–7 loci, 0.82). CLAUDE.md asks for the corner a
  figure covers; quote table (d) and the lowest cell beside it.
- [stratum_fits.rs:1789-1791](../../../../src/parameter_estimation/joint/stratum_fits.rs#L1789-L1791),
  [repeat_tract_parameters.rs:2667](../../../../src/calling/inference/repeat_tract_parameters.rs#L2667)
  and `:2720` — test docs mention `blend_level` and "no climb record" as current.
- Commit `680435b2` says "Fitted numbers: unchanged to the bit (152,442 values in the fixture)". Zero
  moved is confirmed; the count is not reproducible: comparing `tmp/d7/fixture_head.jsonl.gz` with
  the committed fixture gives 146,689 numeric fields outside `errors` and `seconds`. Name what was
  counted, or drop the figure.

## 7. Out of scope observations

- `src/run/depth_ceiling.rs:14` — a broken intra-doc link to `cohort_merge::LocusSummary`, which
  lives in `cohort_merge::observation_cache`. Fails `cargo doc --no-deps` on main as well. One-line
  fix on main.
- `examples/dhat_ng_joint_fit.rs` (touched by E1) and `examples/dhat_ng_merge.rs` (not touched) emit
  9 and 12 unused-item warnings in a default-feature build, because their bodies run only under
  `dhat-heap`. Clippy's `--all-features` gate hides it. A `#![cfg_attr(not(feature = "dhat-heap"), allow(unused))]`
  per file would make it explicit.
- Carried over from the D8 review, still true: no cross-platform digest pins a repeat-tract error
  (the fixture's one class has 1 locus), so D7's new arithmetic in the errors is pinned on no second
  platform. It is plain `+ − × ÷` in a fixed order (checked), so the risk is low; a 3-tract fixture
  would close it.
- 39 occurrences of "length spectrum" remain in `src/` at `909266ee`, all from before the branch's
  rename; main's merge added none.

## 8. Missing tests to add now

None that the findings need: Mi1 and Mi2 delete code rather than add behaviour. If Mi2's code is kept
instead, the guard it would need is
`reading_an_older_file_never_evaluates_a_curve` — read `every_shape_as_written.toml` with
`level_at`/`share_at` instrumented to panic, as in §3's mutation.

## 8a. The diff's own numbers

| claim | where | result |
|---|---|---|
| `ssr_fit.rs` 9,999 → 2,099 lines | plan E1, commit `3198d55d` | CHECKED-CORRECT (`7cbc5f78`: 9,999; `3198d55d`: 2,099) |
| three modules, 688 lines | plan E1 | CHECKED-CORRECT (521 + 63 + 104) |
| seven examples, 1,600 lines | plan E1 | CHECKED-CORRECT (193 + 324 + 222 + 202 + 164 + 279 + 216) |
| "about 12,000 lines removed" | `3198d55d` | CHECKED-CORRECT (14,146 deleted, 2,281 added) |
| "about 2,700 lines" | `2276c93c` | CHECKED-CORRECT (2,683 deleted, 116 added) |
| 59,443 of ~63,000 error fields move; no fitted number moves | plan D7, `680435b2` | CHECKED-CORRECT (59,443 of 62,538 numeric; 63 null errors became numbers; 0 values moved) |
| 152,442 fitted values | `680435b2` | NOT REPRODUCED (146,689 by this count; Nits) |
| full suite 5,230 passed | `2276c93c` | consistent: 5,244 at `909266ee`, which adds main's tests |

**Interactions checked and found sound:**

- **D8's direction line rests on D7's errors.** The 89-to-93 coverage was measured after D7 on all
  96 runs at 2,000 samples, the 33 D7 recovered included (`tmp/d7/coverage_by_number.md`, table (c)).
- **D8's share line (3 → 1 expected locus) on D7's errors.** D1 measured 96 in 100 coverage for
  shares at 1 to 3 expected loci before D7. D7 narrows exactly those errors by a median of 2.8 in 100
  (5th percentile 7.2) — the most of any bin, against 0.7 to 1.9 elsewhere
  (`tmp/review_E/scripts/share_err_ratio.py`, on the fixtures before and after D7, classes of 8 or
  more loci). That moves 96 by well under one point. And share errors reach only the run log, not the
  parameters file.
- **No error reaches calling.** `own_fit_standard_errors` is read only by the file's writer,
  validator and tests, so D7 and D8 change the file and the log, never a call.
- **The merge of main.** Only `Cargo.toml` (different sections) and
  `src/calling/allele_candidates/ssr.rs` (the branch changed two doc words, main a lookup) were
  changed on both sides; main's 16 files read nothing the branch renamed or removed, and main's
  additions use none of the renamed vocabulary.

## 9. What's good

- `to_toml.rs:333-344` builds both threshold paragraphs from `LOCI_FOR_AN_ERROR` and
  `LOCI_FOR_A_DIRECTION_ERROR`, so the file can never state a different line from the one applied.
- `StratumFits::with_rows_of_a_file` (`stratum_fits.rs:826-874`) sends the tests' climb-shaped rows
  through `of_gathered_rows`, the real door an older file uses, instead of keeping a
  production-dead `Derived` variant alive for them.
- `a_cohort_whose_reads_alone_curve_upward_gets_its_errors` (`penalised_fit.rs:4334`) asserts the
  same matrices fail without the prior's block before asserting they succeed with it, so the fixture
  provably exercises D7.
- `refuse_removed_options` (`cli/estimate_parameters.rs:386-396`) refuses the removed flag with the
  reason, and its test checks the help no longer lists it.

## 10. Commands to re-verify

```
# copy
git archive 909266ee | tar -x -C tmp/review_E/src_copy
cd tmp/review_E/src_copy
export CARGO_TARGET_DIR=$PWD/../target
cargo fmt --check
cargo clippy --all-targets --all-features -j 6 -- -D warnings
RUST_TEST_THREADS=6 cargo test --release --all-features --no-fail-fast
cargo doc --no-deps                          # fails on src/run/depth_ceiling.rs:14, as on main

# Mi2: make level_at / share_at panic, then
cargo test --release --all-features --lib --tests   # 13 failures, none in calling or file reading

# §8a and the interaction checks
uv run --no-project python tmp/review_E/scripts/d7_claims.py tmp/d7/fixture_head.jsonl.gz \
  src/parameter_estimation/joint/ssr_fit/bench_fixtures/d_tract_prior_reference.jsonl.gz
uv run --no-project python tmp/review_E/scripts/share_err_ratio.py tmp/d7/fixture_head.jsonl.gz \
  src/parameter_estimation/joint/ssr_fit/bench_fixtures/d_tract_prior_reference.jsonl.gz

# M1
scripts/promote_ng_oracle.sh …               # expect fitted.parameters.toml to differ until re-recorded
```

Scratch: `tmp/review_E/` — `src_copy/` (restored after the mutation), `clippy.log`, `fmt.log`,
`test.log`, `doc.log`, `mutation_curve_eval.log`, `scripts/`.
