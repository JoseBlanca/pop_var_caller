# Code Review: ng_paralog_filter_c1

**Date:** 2026-09-07
**Reviewer:** rust-code-review skill (orchestrator), three sub-agents in isolated worktrees
**Scope:** step C1 of the hidden-duplication filter plan — the scoring context
**Status:** Request-changes

---

### 1. Scope

- **What was reviewed:** commit `0ceb463b` on branch `ng-paralog-filter` — one new module and its
  tests, plus a declaration and two re-exports. Agents detached to `ca3ebb1c`, which is that
  commit plus a `PROJECT_STATUS.md` update.
- **In-scope files:**
  - [scoring_context.rs](../../../../src/ng/run/paralog_filter/scoring_context.rs)
  - [tests.rs](../../../../src/ng/run/paralog_filter/scoring_context/tests.rs)
  - [mod.rs](../../../../src/ng/run/paralog_filter/mod.rs) — the declaration and re-exports
  - [ng_paralog_filter_c1_2026-09-07.md](../implementations/ng_paralog_filter_c1_2026-09-07.md) —
    its quantitative and mechanism claims
- **Deliberately out of scope:** `spill.rs`, `spill_file.rs`, `patch.rs`, the VCF writer (B1–B3,
  reviewed); `src/ng/window_coverage/` (reviewed on the branch it merged from); `src/ng/paralog/`
  (Milestone A — and its five files are byte-for-byte copies of production behind a textual guard,
  so findings there are raised, never edited); `src/paralog/`, `src/var_calling/`,
  `src/sample_summary/` (production, frozen).
- **Categories dispatched:** reliability + extras (a per-sample, per-record path on the filter's
  inner loop, and a numerical boundary); errors + naming + idiomatic + defaults + the diff's own
  numbers; refactor_safety + module_structure + smells. `tooling` skipped — `Cargo.toml` untouched.
  `unsafe_concurrency` skipped — no `unsafe`, no threads, no shared state.

### 2. Verdict

**Request-changes.** One Blocker, four Major, thirteen Minor.

**The Blocker is a fixture, not a defect in the code**, and it is the sharpest kind: the one
histogram every coverage test is built on has a **single GC bin**, and the copied model returns
`gc_bias_curve[0]` for every GC value on a one-bin curve. So the GC half of the relative copy
number — half of what turns a window depth into a copy number — is not exercised by any of the ten
tests, and replacing the GC argument with a literal `0.9` leaves all ten green. A wrong GC there
does not panic and does not give `NaN`; it gives a plausible copy number, wrong by the bias
multiplier, on every sample of every record.

The Majors are of one kind too: **three separate defences whose two halves the suite cannot tell
apart**, and one length invariant that is documented and unenforced.

### 3. Execution status

Run by the orchestrator in the dev container at `ca3ebb1c`, and re-run independently by two agents:

| command | result |
|---|---|
| `cargo test --all-features --lib "ng::run::paralog_filter::scoring_context"` | `ok. 10 passed; 0 failed; 6559 filtered out` |
| `cargo test --all-features --lib --bins --tests` | lib `ok. 6554 passed; 0 failed; 15 ignored`; one integration test failed |
| `cargo clippy --lib --bins --tests --all-features -- -D warnings` | 9 errors, **none in `src/ng/run/paralog_filter/`** |
| `cargo fmt --check` | dirty on 9 unique files, none of them this step's |

The integration failure (`a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`) and
three `needless_lifetimes` in `cohort_merge` are `main`'s at `a33ada0f`; six further lints and the
nine `fmt` files arrive with the merged `ng-window-coverage` branch and fire there identically.
`--all-targets`, `cargo doc --no-deps` and `cargo audit` were not run, for the reasons B1 recorded.

**Findings labelled "Needs verification": 0.** Every finding was produced by a mutation or a probe
run in a reviewer's own worktree, with output quoted in the per-category files under
[tmp/review_2026-09-07_ng_paralog_filter_c1/](../../../../tmp/review_2026-09-07_ng_paralog_filter_c1/).

**Mutation totals: 5 run, 3 survived, 0 changed no behaviour.** Every survivor was proved to change
behaviour on a fixture the file does not contain — which is the finding in each case.

### 4. Open questions and assumptions

1. **Is a half-absent window pair representable in practice?** M1 and M2 below rest on it. The
   owning type says yes in as many words — `WindowCoverage::is_absent`'s doc notes both fields are
   `pub`, "so a depth without the GC fraction that says what depth to expect there is
   representable, and it is not a usable measurement" — and the spill codec round-trips both floats
   by bit pattern. Affects **M1**. No decision needed; recorded because the fix is justified by it.
2. **Should `observations_of` introduce an error type before C3 exists?** The record-level entry
   point that fixes **M2** wants one. The alternative is a `debug_assert_eq!`, which is loud in
   tests and silent in release. Affects **M2**; resolved in favour of the error type, see §6.

### 5. Top 3 priorities

1. **B1** — the GC curve is never exercised, so a wrong GC passes the whole suite.
2. **M1** — neither half of the finiteness guard is tested, and deleting the GC half lets a
   `NaN` GC reach an out-of-bounds index inside a file the copy guard forbids editing.
3. **M2** — a record narrower than the cohort is scored on a prefix, silently.

### 6. Findings

#### Blocker

- **B1: [tests.rs:21-42](../../../../src/ng/run/paralog_filter/scoring_context/tests.rs#L21-L42) — the fixture every coverage test rests on has one GC bin, so a wrong GC argument passes all ten tests**
- **Categories:** reliability
- **Confidence:** High
- **Problem:** `a_fittable_histogram()` sets `gc_bins = 1`, and `gc_multiplier`
  ([coverage_model.rs:349-366](../../../../src/ng/paralog/coverage_model.rs#L349-L366)) returns
  `gc_bias_curve[0]` for every GC value on a one-bin curve. So
  `f64::from(window.gc_fraction)` at
  [scoring_context.rs:172](../../../../src/ng/run/paralog_filter/scoring_context.rs#L172) cannot
  affect any assertion in the file. Measured: replacing it with the literal `0.9` — standing in for
  every way the wrong GC could arrive, a stale window, a constant, the wrong sample's window, a
  percentage where a fraction was wanted — leaves `test result: ok. 10 passed`. The same mutant
  fails at once on a four-GC-bin fixture whose one-copy depth rises with GC.
  **`ONE_COPY_DEPTH` does not rescue it**: on this fixture the fitted scale is exactly 5.25 and the
  curve exactly `[1.0]`, so the doubling test's `at_one` is 1.0 whatever GC is passed.
- **Why it matters:** GC normalisation is half of what turns a window depth into a copy number, and
  a wrong GC neither panics nor produces `NaN` — it produces a plausible relative copy number,
  wrong by the bias multiplier, for every sample of every record, and the run finishes having
  dropped and kept the wrong loci. The rubric names exactly this: no test for a path that, if
  broken, is silently wrong.
- **Suggested fix:** add a four-GC-bin fixture and
  `the_relative_copy_number_uses_the_windows_own_gc_and_not_a_constant` (*Missing tests* 1),
  verified to pass at `ca3ebb1c` and to fail under the mutant. Keep the one-bin fixture for the
  tests that are not about GC — it is the right fixture for those.

#### Major

- **M1: [scoring_context.rs:168-170](../../../../src/ng/run/paralog_filter/scoring_context.rs#L168-L170) — neither half of the finiteness guard has a test that fails when that half is deleted, and one half prevents a panic inside a frozen copy**
- **Categories:** reliability
- **Confidence:** High
- **Problem:** both fixtures that exercise the guard set **both** fields to `NaN`, so either half
  alone still catches them. Deleting the GC half leaves all ten green and lets a
  `(NaN GC, finite depth)` window reach `relative_copy_number`, where `gc_multiplier`'s `x <= 0.0`
  and `x >= (gc_bins - 1) as f64` are both false for `NaN`, `x.floor() as usize` saturates to `0`,
  and the interpolation indexes `gc_bias_curve[1]` on a one-bin curve:
  `panicked at src/ng/paralog/coverage_model.rs:365:67: index out of bounds: the len is 1 but the
  index is 1`. Deleting the depth half also leaves all ten green and returns `Some` with a `NaN`
  relative copy number; a probe measured what that costs — *one `NaN` sample of three: lr = NaN,
  samples_used = 3*, so one sample unscores the whole record while still counting itself.
- **Why it matters:** the guard is the only thing between a half-absent pair and either a
  process-killing panic in code that cannot be edited, or a silently unscored record. Its two
  halves are indistinguishable from one half to the suite, so a refactor that simplifies it loses
  one without a red test.
- **Suggested fix:** *Missing tests* 2 and 3, and replace the comment — which today reads as taste
  ("any non-finite value is unusable, however it arrived") — with what it is actually preventing.

- **M2: [scoring_context.rs:161-196](../../../../src/ng/run/paralog_filter/scoring_context.rs#L161-L196) — a record narrower than the cohort is scored on a prefix, and nothing notices**
- **Categories:** reliability, and convergent with the errors agent's reading of the same lines
- **Confidence:** High
- **Problem:** the only length check is `entry.samples.window(sample)?`, which answers *absent* for
  an index past the record's own rows — so a short record scores as though the missing samples were
  merely uncovered. Verified with a three-sample context and a two-sample record. The doc at
  [:83-85](../../../../src/ng/run/paralog_filter/scoring_context.rs#L83-L85) states the invariant —
  "that length is the cohort size every scored record must also have" — and nothing enforces or
  tests it. The complementary case is the same shape: `score_locus_for_paralogy` answers a length
  mismatch with a **neutral score** rather than an error
  ([locus_score.rs:276-281](../../../../src/ng/paralog/locus_score.rs#L276-L281)), so a long
  record's extra rows would be dropped just as quietly.
- **Why it matters:** the disagreement would be a wiring error between C2's sink and this context,
  and its symptom is a weaker score on some records — the defect shape that survives a whole run
  because nothing about the output looks wrong. It is spec §6 trap 3 in a place the σ₀ slice
  already defends by construction and this does not.
- **Suggested fix:** a record-level entry point that checks once and fills a caller-owned scratch
  buffer, which removes the per-sample second lookup at the same time and matches the project's
  scratch-buffer preference. Body in the reliability findings file; adopted in §8's *Missing tests*
  11.

- **M3: [scoring_context.rs:109-120](../../../../src/ng/run/paralog_filter/scoring_context.rs#L109-L120) — a bad fit configuration is reported once per sample as a verdict about the data**
- **Categories:** errors
- **Confidence:** High
- **Problem:** `SingleCopyCoverageModel::fit` calls `validate_config` on every call, so one wrong
  knob makes a 63-sample run report that the coverage model rested on 0 of 63 samples, with 63
  identical `InvalidConfig` reasons — a configuration mistake dressed as a statement about every
  sample's coverage. Latent today, because §3.6 keeps these knobs as compiled defaults, and live
  the moment they reach the command line.
- **Why it matters:** the run report's line is "samples whose coverage model was rejected, each
  with its reason". A reader of that line would conclude the cohort was uncovered.
- **Suggested fix:** validate once before the loop and fail construction, rather than folding a
  config error into a per-sample verdict.

- **M4: [scoring_context.rs:161-196](../../../../src/ng/run/paralog_filter/scoring_context.rs#L161-L196) — `observation_of` reads the sample's row twice and reads its counts by field, so a field added to either row type compiles here unchanged**
- **Categories:** refactor_safety, idiomatic — convergent
- **Confidence:** High
- **Problem:** `samples.window(sample)` and then the match arm's `samples.get(sample)`, and the
  counts read as `sample.ref_reads` / `sample.alt_reads` rather than destructured. The sibling
  codec in the same module destructures exhaustively with a comment saying why — and the spec's own
  deferred work (§8's tract-aware allele term) is precisely a field arriving on
  `RepeatTractSample`. The second `?` is also unreachable once one `match` yields both halves.
- **Why it matters:** this is the file's only reader of those rows, so a field added to them and
  forgotten here is silently dropped from every score. The enum's *variants* are compiler-checked
  (verified: adding a third stops the build here); its *fields* are not.
- **Suggested fix:** one `match` producing the window and the counts together, destructuring both
  row types exhaustively.

#### Minor

- **Mi1: four different reasons a sample is absent reach the caller as the same `None`, and one of
  them is corruption** — no model, no usable window, no reads, an index past the cohort, and a read
  pair that overflows `u32` via `checked_add`. C3 folds all of them into `samples_used` and cannot
  separate a corrupt spill from a sparse cohort. The record-level entry point removes the index
  case; the overflow wants a counted total or a doc sentence.
- **Mi2: the guard silently diverges from `WindowCoverage::is_absent` and nothing says so.** The
  owning type has this predicate (`is_absent`: either field `NaN`); this hand-rolls `!is_finite()`
  on both. The divergence is **right** — `is_absent` accepts `±∞`, and an infinite depth would
  divide to `+∞`, be winsorised to `max_relative_copy_number` (4.0) and read as a confident
  four-copy paralog — but no comment says it is deliberate and no test pins it, so the next reader
  is as likely to "simplify" it as to keep it.
- **Mi3: the tolerance guarding `ONE_COPY_DEPTH`'s meaning is 0.1 where the fixture is exact.** The
  peak is symmetric and the curve is `[1.0]`, so the value is 1.0 to within `1e-9`; `< 0.1` would
  pass with the fitted scale 9% wrong. Tighten it and say why the fixture is exact.
- **Mi4: `coverage_models` and `why_no_model` are two parallel `Option` vectors spelling one
  per-sample sum type**, so "both" and "neither" are representable. The loop already *builds*
  `Result<SingleCopyCoverageModel, WhyNoCoverageModel>` and then takes it apart. Length drift is
  the weaker half — no `pub` field and no `&mut self` method, so all four are filled in one
  function. `single_copy_depth_sd` should stay separate: the scorer takes it as a contiguous
  `&[f64]`, which is also how production splits it.
- **Mi5: `new` takes `&[f64]` where the accessor the spec points at returns `&[InbreedingF]`** — a
  newtype validated to `[0, 1)`. C2 would strip it at the call site; taking the newtype moves the
  strip inside, to the one boundary with the frozen copied code that wants `f64`.
- **Mi6: `WhyNoCoverageModel` has no `Display`**, so spec §3.5's "each with its reason" has nothing
  to print but `{:?}`.
- **Mi7: `samples_with_a_coverage_model` returns a count but reads as a collection**, five lines
  below `sample_count`, which gets it right.
- **Mi8: `precompute()` is a verb naming an accessor** that computes nothing.
- **Mi9: the by-value `Vec<SampleHistogram>` against three borrows is correct and unexplained.** It
  frees each sample's bins as the fit proceeds, which is the memory lever at three thousand
  samples — but `fit` only borrows, so a reader will find the ownership unnecessary and may
  simplify it away. It needs one line in the doc.
- **Mi10: the fit configuration is dropped after construction**, so nothing can answer "which
  single-copy band produced these models?" — and spec §3.1 flags those constants as inherited from
  a prototype and never re-measured, which is exactly the kind a run report should state.
- **Mi11: five `pub` accessors, three of which exist so the unwritten C3 can reassemble
  `score_locus_for_paralogy`'s arguments from outside.** Production's context of the same name keeps
  those three private behind one `score` method — which is what makes spec §6 trap 3 unspellable by
  a caller. Filed at Medium confidence as a recommendation for C3's shape rather than a change now.
- **Mi12: [mod.rs](../../../../src/ng/run/paralog_filter/mod.rs)'s "what exists so far" paragraph
  still lists the scoring context among "later steps"**, in the file that now declares and
  re-exports it.
- **Mi13: the implementation report defers a cost to step D2 that is now measured**, and the
  measurement says it is not worth deferring — see §6a.

#### Nits

- `new` is `#[must_use]` on a constructor that panics rather than returning a `Result`; harmless,
  but the attribute is carrying no weight.
- This is production's `build_observation` under a different name in the same crate; a line naming
  the correspondence would help the next reader compare them.

### 6a. The diff's own numbers

Fourteen claims re-derived by running something. **Twelve CHECKED-CORRECT, two WRONG**, and one
correct but weaker than it reads.

| claim | verdict | the real value |
|---|---|---|
| 231 / 325 lines, 10 tests | CHECKED-CORRECT | as stated |
| the lib count moves **6,544 → 6,554** | CHECKED-CORRECT | full suite run at both ends |
| both mutations kill exactly the named test, 9 passed 1 failed | CHECKED-CORRECT | both re-run independently |
| the three dropped histogram fields occur **zero times** inside the guarded fit | CHECKED-CORRECT | whole-file counts are 2/1/2, all in header prose or the test fixture |
| a short σ₀ slice gives a **neutral score**, not an error | CHECKED-CORRECT | `locus_score.rs:276-281` |
| the allele term is zero at zero reads in **both** the H1 and H2 arms, so it cancels | CHECKED-CORRECT | both arms read |
| "not `Clone` because `CoverageModelError` is not" | CHECKED-CORRECT | |
| **"694 lines of the fit stay compared byte for byte"** | **WRONG** | **644.** The header is 54 `//!` lines, `#[cfg(test)]` is at line 700, and the guard drops the trailing blank. Production's side computes 644 too, which is why the guard is green while the number is wrong. In two places: `coverage_model.rs:52-54` and `copy_fidelity.rs:241`, both introduced in `dfc6a40c` |
| the module header's list of the fit's refusal reasons | **WRONG** | it names **three of six**, and omits `NoTiles` — the one the file's own test asserts on |
| `copy_fidelity.rs:7`'s table: "1,158 past its module header" | **WRONG** (incidental, pre-existing) | 1,157 |
| `ONE_COPY_DEPTH` is one copy "within 0.1" | true but weak | it is **exactly** 1.0 — the tolerance was set to `0.0` and it still passed. See **Mi3** |

### 7. Out of scope observations

- **[spill.rs:401](../../../../src/ng/run/paralog_filter/spill.rs#L401) and
  [:607](../../../../src/ng/run/paralog_filter/spill.rs#L607) — a third `SpilledSamples` variant
  would compile there and be silently treated as "not a repeat tract".** The writer's check is
  `matches!(entry.samples, SpilledSamples::RepeatTract(_))` and the decoder branches on
  `if is_repeat_tract`; neither appeared among the four errors when a third variant was added.
  C1 was asked the same question and answers it the opposite way — its `match` is exhaustive by
  construction. `spill.rs` was reviewed at B1 and reshaped since, so this is the reshape's, and
  worth closing before a third variant exists rather than after.
- `src/ng/window_coverage/mod.rs` — `WindowCoverage::is_absent` accepts `±∞`. Correct for its own
  purpose; noted because **Mi2** turns on the difference.

### 8. Missing tests to add now

Eleven plus two probes, all compiled and run by the reliability agent against `ca3ebb1c` — final
run `23 passed; 0 failed`. Bodies in
[reliability_extras.md](../../../../tmp/review_2026-09-07_ng_paralog_filter_c1/reliability_extras.md)
and in the salvaged diff beside it.

**`observation_of`**

1. `the_relative_copy_number_uses_the_windows_own_gc_and_not_a_constant` — **B1's fix.** Needs a
   four-GC-bin fixture whose one-copy depth rises with GC, and a fit config whose smoother is the
   identity so the curve is the bins' own medians.
2. `a_window_with_only_one_field_absent_is_absent` — **M1's first half**, the one standing between
   a `NaN` GC and an index panic inside the frozen copy.
3. `an_infinite_window_depth_is_absent_rather_than_a_winsorised_paralog` — **M1's second half and
   Mi2's pin**: `+∞` must not become a confident four-copy paralog.
4. `a_sample_at_zero_window_depth_reads_as_zero_copies_rather_than_absent` — zero depth is a
   measurement, not an absence.
5. `a_sample_far_above_its_one_copy_depth_is_handed_over_unwinsorised` — the winsor cap is the
   scorer's, not this module's, and the raw value must reach it.
6. `the_relative_copy_number_is_the_same_at_three_reads_and_at_three_hundred` — `CLAUDE.md`'s depth
   range, on the quantity that is supposed to be depth-invariant.
7. `a_cohort_of_one_scores_through_the_scorer` — `N = 1`, the hardest case and the one the spec
   calls out.
8. `a_cohort_of_three_thousand_hands_the_scorer_every_sample` — the other end of the declared range.
9. `a_cohort_with_no_samples_at_all_is_built_and_scores_nothing` — the degenerate end.
10. `a_record_mixing_covered_and_absent_samples_hands_over_only_the_covered_ones` — the ordinary
    cohort shape, which no existing test has.
11. `a_record_shorter_than_the_cohort_is_refused_rather_than_scored_on_a_prefix` — **M2's fix**,
    which requires the record-level entry point.

### 9. What's good

- **The two halves of the tract rule are both tested and both mutation-killed**
  ([tests.rs](../../../../src/ng/run/paralog_filter/scoring_context/tests.rs)) — the step's whole
  reason for existing is pinned, and the mutation that reproduces spec §6 trap 1 exactly fails the
  test named for it.
- **The `SpilledSamples` match is exhaustive by construction**: adding a third variant stops the
  build here, and nothing in the file reads `entry.is_repeat_tract`, so the variant is the single
  source of truth for which signals a record carries. Verified by compile probe.
- **`WhyNoCoverageModel` keeps four causes apart where a count would have merged two opposite
  problems** — a sample never covered and one whose depth ran off the top of the histogram.
- **σ₀ is `NaN` and not zero for a sample with no model**, with a test whose message says why: a
  zero σ₀ is a claim of perfect precision rather than an absence.
- **The by-value `Vec<SampleHistogram>`** frees each sample's bins as the fit proceeds — the right
  call at three thousand samples, and the sort of thing usually got wrong in the other direction.

### 10. Commands to re-verify

    ./scripts/dev.sh cargo test --all-features --lib "ng::run::paralog_filter::scoring_context"
    ./scripts/dev.sh cargo test --all-features --lib --bins --tests
    ./scripts/dev.sh cargo clippy --lib --bins --tests --all-features -- -D warnings
    ./scripts/dev.sh cargo fmt --check

### Author response convention

Address each finding by its identifier (`B1`, `M3`, `Mi7`) with one of: `fixed in <commit>` /
`disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer §4's two open
questions first.
