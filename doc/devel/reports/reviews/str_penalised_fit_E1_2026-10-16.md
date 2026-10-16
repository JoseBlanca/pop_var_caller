# Code Review: str_penalised_fit E1 — the climb deleted
**Date:** 2026-10-16
**Reviewer:** rust-code-review skill, run by one agent (no sub-agents, at the dispatcher's request)
**Scope:** the uncommitted change in worktree `pop_var_caller-str-penalised-fit` on HEAD 680435b2 —
plan step E1 of `doc/devel/ng/impl_plan/str_penalised_fit.md`, which deletes the old repeat-tract
fit (*the climb*: a per-stratum coordinate climb over a 256-point average of each tract's
likelihood, with parametric curves drawn across repeat count afterwards) as
`doc/devel/ng/arch/str_penalised_fit.md` §4.1 lists.
**Status:** Approve-with-changes

---

## 1. Scope

- **Reviewed:** `git diff HEAD --diff-algorithm=histogram` — 38 files, 1,550 lines added and
  13,289 removed. The working tree was copied with `rsync` to `tmp/review_e1/src_copy/` and HEAD
  exported with `git archive` to `tmp/review_e1/head_copy/`; every build and run below used those
  copies, with `CARGO_TARGET_DIR` under `tmp/review_e1/`. Nothing in the worktree was changed
  except this report.
- **In scope:** the four files the deletion rewrote (`ssr_fit.rs`, `bench_fixtures/fit_bench.rs`,
  `expected_counts.rs`, `counting_fit_errors.rs`), the fixture edits in `calling/` and `tests/`,
  `stratum_fits.rs`, the eleven touched or deleted examples and the bench, `Cargo.toml`, and the
  scripts and documents that name deleted items.
- **Out of scope:** the penalised fit itself (reviewed at D2a/D6/D7); `slippage_curve.rs` and
  `share_curve.rs`, which the change does not touch, except where something in them lost its last
  caller.
- **Categories covered** (all by this agent): reliability, errors, naming, defaults, idiomatic,
  refactor safety, module structure, float portability (only to confirm the deletion moved no
  number), smells, tooling, extras (diff matches stated intent).

## 2. Verdict

**Approve-with-changes.** Nothing that calls or fits behaves differently, and that was checked
end to end, not read off the diff: the same psps and climb-written parameters files call to
byte-identical VCFs on the HEAD and E1 binaries, the penalised fit's reference fixture is exact,
and the bench's counting-fit rows are unchanged apart from the method's name. The four findings
are all Minor: code that no production path calls anymore but which still compiles into the
production build, two research scripts left speaking of methods that are gone, and two kept types
whose documented reason for being kept is not true.

## 3. Execution status

All on macOS arm64, native cargo (no container runtime is needed for these; CLAUDE.md allows native
builds in the light dev loop).

| command | exit | result |
|---|---|---|
| `cargo build --release --all-targets` (no features), E1 copy | 0 | **8 new `dead_code` warnings in the library** (Mi1); the example warnings are the same as at HEAD |
| `cargo fmt --check`, E1 copy | 0 | clean |
| `cargo clippy --all-targets --all-features -- -D warnings`, E1 and HEAD copies | 101 / 101 | the same 5 errors in both, all in code E1 did not write (§7) |
| `cargo clippy --lib -- -D warnings` (no features), E1 and HEAD | 101 / 101 | HEAD: 1 error (§7). E1: that one plus 8 `dead_code` errors (Mi1) |
| `scripts/check_float_ban.sh`, E1 and HEAD copies | 1 / 1 | fails at HEAD already (§7); E1 adds 8 more compile errors to the same failure |
| `cargo doc --no-deps --document-private-items --features bench-fixtures`, E1 and HEAD | 101 / 101 | the same 56 warning and error lines in both: no broken link added |
| `cargo test --release --features bench-fixtures --lib` | 101 | 5,201 passed, 3 failed, 7 ignored. The 3 are `should_panic` tests on `debug_assert`s, which cannot panic in release; they fail identically in the implementer's own release log at HEAD (`tmp/e1/test_release.log`) |
| `cargo test --release --features bench-fixtures --test '*' --no-fail-fast` | 0 | all six integration test files pass, both `*_on_real_cells` tests included |
| targeted: `cross_platform_digests`, `parameters_file`, `stratum_fits`, `ssr_fit` | 0 | 446 passed, 2 ignored. The digest constants are not in the diff, so the four digest tests pass on unchanged checksums |
| `str_fit_fixture.py run` + `compare --require exact --partial`, the 144 runs below 2,000 samples, E1 bench binary | 0 | **144 of 144 exact**; the 9 runs at 2,000 samples were not run, as instructed |
| calls with climb-written parameters files, HEAD binary against E1 binary | 0 | **9 of 9 VCFs identical** outside `##commandline`/`##parametersFile`, and the parameters file each run writes beside its VCF identical too (table below) |
| `ng_str_fit_bench`, 6 cells, HEAD (`B_raw`) against E1 (`B`) | 0 | **every row, call line and locus line identical** once the HEAD-only `_curves` rows and timing rows are dropped and `_raw` stripped |

**Calls with older parameters files.** Each file was written by the climb. *Format* is the file's
`format_version`; *blend rows* counts lines naming a blended level or share, which only the climb
wrote.

| file | format | blend rows | input | records | HEAD vs E1 |
|---|---|---|---|---|---|
| `tmp/oracle_C6/fitted.parameters.toml` | 3 | 17 | 4 tomato psps, `call-from-psps` | 6,683 | identical; also identical to the recorded `tmp/oracle_C6/from_fit.vcf` |
| `tmp/oracle_A4/fitted.parameters.toml` | 2 | 17 | the same 4 psps | 6,683 | identical |
| `tmp/oracle_B1/fitted.parameters.toml` | 2 | 17 | the same 4 psps | 6,683 | identical |
| `tmp/d5/out/climb/params/hg002_30x`, `_50x` | 3 | 79 each | HG002 tandem-repeat tier BAM, 50,000 regions | 19,330 / 19,583 | identical |
| `tmp/d5/out/climb/params/HG003_30x`, `_50x` | 3 | 10 each | HG003 per-sample BAM, its 100-region BED | 721 / 722 | identical |
| `tmp/d5/out/climb/params/HG004_30x`, `_50x` | 3 | 10 each | HG004 per-sample BAM, its 100-region BED | 711 / 715 | identical |

**The bench before and after.** *Cell* is one draw of the research bench
(`examples/ng_str_fit_bench.rs`). Three of the six were chosen to make the counting fit's floor of 8
loci with reads decide something, since that floor moved from `ssr_fit::DEFAULT_REFUSAL_FLOOR` to
`expected_counts::COUNTING_FIT_FLOOR`: at 7 loci every class is left unfitted, and at 9 loci each
held-out fold trains on 7 or 8.

| cell | methods | rows compared | what it exercises |
|---|---|---|---|
| period 2, 60 samples, 3 reads, PCR, `--loci-at-most 100` | B, B_within, both fallbacks | 1,451 rows, 176 call lines, 3,512 locus lines | the ordinary case |
| period 1, 1 sample, 3 reads, PCR | B, B_within, both fallbacks | 975 / 120 / 3,700 | one sample |
| period 2, 10 samples, 30 reads, PCR-free | B, B_bootstrap, both fallbacks | 1,451 / 176 / 3,512 | the bootstrap errors |
| period 2, 1 sample, 3 reads, 10 loci a class | B, B_within | 1,431 / 176 / 620 | thin classes |
| period 2, 3 samples, 3 reads, 9 loci a class | B, B_within | 1,451 / 176 / 84 | the floor deciding held-out folds |
| period 1, 3 samples, 3 reads, 7 loci a class | B | 339 / 18 / 0 | every class below the floor: 14 of 14 classes have no level, in both |

**Needs verification:** 0.

## 4. Open questions and assumptions

1. **Should production code compile what only the bench and tests call?** Mi1's fix assumes no:
   gate it, as the module already gates `bench_fixtures`. If the owner prefers to keep the
   library-scoring code visible for the library-multiplier work that will rebuild on the new fit
   (arch §4.1, `library_stutter.rs` row), an `#[allow(dead_code, reason = …)]` naming that plan is
   the alternative. Affects Mi1.
2. **Is a type kept only so tests can manufacture older-file shapes worth keeping in the
   production enum?** Mi4 asks the author to choose between documenting that reason and building
   those fixtures through the file-reading path instead. Affects Mi4.

## 5. Top 3 priorities

1. **Mi1** — gate the eight items that the default build reports as dead, so the production build
   is warning-free again.
2. **Mi2** — the research summariser now prints `—` where it printed the penalised fit's
   improvement over the counting fit, because it still looks for `B_raw`; one line fixes it.
3. **Mi4** — the doc comments on `DerivedStratum` and `PeriodAlleleLengthVariation` give a reason
   for keeping them that is not true; a reader trusting them would believe older files are read
   through types they never touch.

## 6. Findings

### Blocker

None.

### Major

None.

### Minor

#### Mi1: `src/parameter_estimation/joint/ssr_fit.rs:850-1160` — eight items are now dead in the default build

- **Confidence:** High
- **Problem:** With the climb gone, `TractLikelihoods` (`:921`, and its `of` at `:1068`),
  `add_ln_reads` (`:1050`), `ReadDistributionsByLevel`'s `distribution_of_library` field (`:952`)
  and its `of_library` / `rows_by_distribution` methods (`:1007`, `:1016`), `ln_sum_exp` (`:1128`),
  `ln_gamma` (`:1136`) and `normalise` (`:850`) are called only from `bench_fixtures` (behind the
  `bench-fixtures` feature) and from `#[cfg(test)]` code. The production
  build therefore reports them, verbatim from `tmp/review_e1/build_new.log`:

  ```
  warning: function `normalise` is never used
  warning: struct `TractLikelihoods` is never constructed
  warning: field `distribution_of_library` is never read
  warning: methods `of_library` and `rows_by_distribution` are never used
  warning: function `add_ln_reads` is never used
  warning: associated function `of` is never used
  warning: function `ln_sum_exp` is never used
  warning: function `ln_gamma` is never used
  warning: `pop_var_caller` (lib) generated 8 warnings
  ```

  At HEAD the library built without a warning (`tmp/e1/baseline_build.log`). The implementer's
  clippy ran with `--all-features`, which turns `bench-fixtures` on and hides all eight.
- **Why it matters:** this is §4.1's "dead without a reason" in its plainest form: the production
  binary compiles a likelihood table nothing in it builds. Every featureless `cargo clippy -- -D
  warnings` now fails on them; `scripts/check_float_ban.sh` runs exactly that and was already
  failing at HEAD for another reason (§7), so E1 adds a second, independent cause to that gate's
  failure.
- **Suggested fix:** gate each item on the same condition as the module that uses it. The
  production fit uses `ReadDistributionsByLevel::of` and `read_distribution_by_allele`
  (`expected_counts.rs:652`, `:1103`, `:1724`), so the struct stays and only the per-library parts
  are gated.

  ```rust
  #[cfg(any(test, feature = "bench-fixtures"))]
  fn normalise(weights: &mut [f64]) { … }

  #[cfg(any(test, feature = "bench-fixtures"))]
  struct TractLikelihoods { … }
  #[cfg(any(test, feature = "bench-fixtures"))]
  impl TractLikelihoods { … }
  #[cfg(any(test, feature = "bench-fixtures"))]
  #[inline(always)]
  fn add_ln_reads(…) { … }
  #[cfg(any(test, feature = "bench-fixtures"))]
  fn ln_sum_exp(…) -> f64 { … }
  #[cfg(any(test, feature = "bench-fixtures"))]
  fn ln_gamma(…) -> f64 { … }

  struct ReadDistributionsByLevel {
      read_distribution_by_allele: Vec<Vec<Vec<f64>>>,
      // read only by the library-aware scoring the bench and tests use
      #[cfg_attr(not(any(test, feature = "bench-fixtures")), allow(dead_code))]
      distribution_of_library: Vec<usize>,
  }
  // and the two methods `of_library`, `rows_by_distribution` under the same cfg
  ```

  Simpler still: move `TractLikelihoods`, `add_ln_reads`, `ln_sum_exp`, `ln_gamma` and `normalise`
  into `bench_fixtures` (which already does `use super::*`), next to their only callers. Re-run
  `cargo build --release` without features and `cargo clippy --lib -- -D warnings`: the eight lines
  should be gone and only the pre-existing `single_match` in §7 should remain.

#### Mi2: `scripts/research/str_fit_summarise.py:288` and `:321-336` — the side-by-side tables still look for `B_raw` and `B_curves`

- **Confidence:** High (run, not read)
- **Problem:** the bench now reports the counting fit as `B` (`fit_bench.rs` module doc), but
  `comparison()` pairs D with a hard-coded `"B_raw"` and the other methods with `"B_curves"`, and
  the one-sample table reads `B_raw`, `B_curves` and `B_curves_own_shares` by name. Run on the E1
  bench's own rows (`tmp/review_e1/summarise_new/`, `--compare B,D`), the paired-difference columns
  come out empty, with no error:

  ```
  | cell | seeds | loci a seed | B | D | D − B_raw | D − B_curves |
  | period 2, 60 samples, 3 reads, pcr | 1 | 878 | -0.342 | -0.587 | — | — |
  ```

  The `D − B` difference is the column the D2 conclusion was read from. The same file's usage
  text (`:19`) and `--compare` help (`:424`) still give `B_raw,B_curves,D` as the example.
  Relatedly, `scripts/research/str_fit_grid.sh:13` shows `METHODS="A B"` as its usage example and
  `:18` lists `A` among the bench's methods; run as written, the bench now stops with
  `no method is called A`.
- **Why it matters:** a dash reads as "not measured", so a reader re-running the comparison on new
  rows is told nothing rather than told the script is out of date.
- **Suggested fix:**

  ```python
  # str_fit_summarise.py, in comparison()
  base = "B" if any((cell, "B") in scores for cell in cells) else "B_raw"   # rows before 2026-10-16
  pairs = [("D", base)] + [(m, "B_curves") for m in methods
                           if not m.startswith("B") and any((c, "B_curves") in scores for c in cells)]
  ```

  and in the one-sample table read `pair(base)`, dropping the `B_curves` columns when no row has
  them. Update `:19` and `:424` to `--compare B,D,D_working_F`, and `str_fit_grid.sh:13, 18` to
  `METHODS="B D_tract_prior"` and "(B, B_within, D, D_tract_prior, …)".

#### Mi3: `src/parameter_estimation/joint/share_curve.rs:533` and `:643` — the climb's share-curve ladder and its blend have no caller left

- **Confidence:** High
- **Problem:** `share_curve_for_a_period` (the four-rung fallback that gave every period a share
  curve) and `blend_share` were called only by the climb's curve step. After E1 the only callers
  are their own unit tests in `share_curve.rs`; `ssr_fit.rs:144-147` still re-exports both. The
  implementer kept the curve-fitting code because "`tests/*_on_real_cells.rs` and
  `str_fit_penalised_em` do" use it — true of `choose_rise_shape`, `choose_share_shape`, `blend_level`
  and the curve types, not of these two (searched `src`, `tests`, `examples`, `benches`:
  `share_curve_for_a_period` has 8 mentions, all in `share_curve.rs`; `blend_share` 11 there, plus
  one comment in `to_run_parameters.rs:1715`). Arch §4.1 lists `draw_and_apply_curves`'s blend
  among what goes; `blend_share` is half of that blend.
- **Why it matters:** about 180 lines (`share_curve.rs:518-700`, with `BlendedShare`, used nowhere
  else) plus their tests are maintained for a caller that no longer exists; they are `pub`, so the
  compiler cannot say so.
- **Suggested fix:** delete `share_curve_for_a_period`, `blend_share`, `BlendedShare`, their tests, and the two names from the re-export at `ssr_fit.rs:144-147`; reword
  the comment at `to_run_parameters.rs:1715` to "what the climb's blend wrote in that arm". If the
  owner wants them kept for reference, say so in their doc comments instead.

#### Mi4: `src/parameter_estimation/joint/ssr_fit.rs:719-724` and `:871-877` — `DerivedStratum` and `PeriodAlleleLengthVariation` are documented as kept for reading older files, but reading never touches them

- **Confidence:** High
- **Problem:** a parameters file is read into calling's tables by `StratumFits::of_gathered_rows`
  (`stratum_fits.rs:722`), which takes `FittedAlleleLengthVariation` per period and slippage rows
  per stratum; it never builds a `StratumOutcome`, a `DerivedStratum` or a
  `PeriodAlleleLengthVariation`. In non-test code, `StratumOutcome::Derived` appears only as match
  arms (`ssr_fit.rs:773-815`, `stratum_fits.rs:591`), and `PeriodAlleleLengthVariation` only as the
  parameter of `with_period_allele_length_variations`, which only tests call. Their three
  bookkeeping fields `tracts_fitted`, `strata_pooled` and `converged` (`ssr_fit.rs:893-898`) are
  read nowhere. Yet the doc comments say the variant "stays … because the climb's parameters files
  hold rows of this shape and the file's writer and calling's tables still take them", and that
  `PeriodAlleleLengthVariation` "is how a caller hands such a table" to calling. The plan's done
  note (`str_penalised_fit.md:425`) has it right — "(tests only now)" — the code comments do not.
- **Why it matters:** the next person deciding what older files need will look here first, and be
  told that two test-only types are on the reading path. The true reason they survive — tests
  push older-file shapes through the writer to round-trip them (`to_run_parameters.rs:1802`,
  `from_run_parameters.rs:1353`, `:1379`) — is a reasonable one, but it is not written down.
- **Suggested fix:** pick one.
  1. Keep both, and say why: replace the sentence at `:723-724` with "No fit produces one and no
     file is read into one; the variant stays because the parameters-file tests write the climb's
     derived rows through it, and calling reads those rows with `StratumFits::of_gathered_rows`."
     Do the same at `:874-877`, and mark `tracts_fitted`, `strata_pooled`, `converged` as written
     by the climb and read by nothing.
  2. Or move those tests onto `of_gathered_rows` and the file writer's row types, then delete
     `StratumOutcome::Derived`, `DerivedStratum`, `PeriodAlleleLengthVariation` and
     `with_period_allele_length_variations`. Larger, and it removes a variant every match in
     calling has to carry.

### Nits

- `ssr_fit.rs:2074` — the test module `the_largest_table` lost the test its name described
  (`the_table_is_rows_times_genotypes_times_eight_bytes_for_the_largest_stratum`); what is left
  tests `rows_with_reads` and `heap_bytes`. Rename to `evidence_size`.
- `ssr_fit.rs:379` — `tracts_with_reads` is "the count the per-stratum floor is measured in"; no
  production fit has a floor now. Say "the count the bench's `COUNTING_FIT_FLOOR` is measured in".
- `examples/ng_str_guard_census.rs:49` — `CLASS_FLOOR: u64 = 8` restates
  `expected_counts::COUNTING_FIT_FLOOR`; use that constant (`as u64`) so the two cannot drift.
- `ssr_fit.rs:748-751` — `StratumRefusal` has one variant, and `StratumOutcome::Refused::tracts`
  is now always 0 (`penalised_fit.rs` asserts it). Fine to keep for the file format's sake; worth
  a sentence saying so.
- Documents outside E2's list that name deleted items: `doc/devel/ng/arch/parameter_prepass_joint_fit.md:325`
  (`fit_strata`), `doc/devel/ng/spec/population_diversity.md:390` (`fit_strata`),
  `doc/devel/ng/spec/str_fit_to_calling.md:184` (`fit_period_allele_length_variations`), and
  `doc/devel/ng/impl_plan/str_slippage_level_curve.md:289`, whose Markdown link to
  `examples/ng_ssr_thin_stratum_gate.rs` is now broken. E2 names only `parameter_prepass_ssr.md` and
  `str_slippage_level_curve.md`; widen it.
- `scripts/portable_float_baseline.sh:139-144` still runs `ng_joint_fit_perf`, which no longer has
  its three repeat-tract groups (`ng_joint_fit/stratum_by_tracts`, `…/stratum_by_samples`,
  `…/strata`); a before/after comparison against a log from before E1 will find them missing
  rather than slower.
- The plan's done note says `ssr_fit.rs` is 2,120 lines; `wc -l` gives 2,119. The other numbers in
  it were checked and are right: the three deleted modules are 521 + 63 + 104 = 688 lines, the
  seven deleted examples 1,600 lines, and the `climb:` field was removed from 11 struct literals
  outside the files whose climb code was deleted (15 in all).

## 7. Out of scope observations

- **The float-ban gate fails at HEAD already.** `scripts/check_float_ban.sh` runs a featureless
  `cargo clippy -- -D warnings`, which stops on `clippy::single_match` at
  `src/run/edge_deletion_gathering.rs:201`; the script then reports every banned method as "listed
  but not refused" (`tmp/review_e1/float_ban_head.log`). Since the gate cannot compile the
  library, it currently proves nothing about the ban. Fix the `match` (an `if let`), then Mi1, and
  re-run the script.
- `cargo clippy --all-targets --all-features -- -D warnings` fails identically at HEAD and E1 on 5
  errors: `expected_counts.rs:2261` (`useless_conversion`), `penalised_fit.rs:3116-3117`
  (`needless_range_loop`, twice), `fit_bench.rs:392` (disallowed `f64::powi`), and the
  `single_match` above. Separate follow-up.
- `cargo doc --document-private-items` fails at HEAD and E1 on the same 56 lines of unresolved links
  (e.g. `src/calling/likelihood/ssr_emission.rs:1061`). E1 adds none.
- Three `should_panic` tests (`left_align_repeated` and `left_align_structured`
  `an_offset_past_the_reference_panics`, `genetics::lgamma_panics_on_non_positive_in_debug`) fail
  in a release-profile test run because the panics are `debug_assert`s. Pre-existing; mark them
  `#[cfg(debug_assertions)]`.

## 8. Missing tests to add now

None that E1 owes. Its deletions removed tests together with the code they tested. The one property
of kept code that lost its test is the claim at `ssr_fit.rs:941-945`, that libraries kept apart at a
multiplier of one score exactly as their reads summed. That test
(`libraries_kept_apart_at_a_multiplier_of_one_fit_as_their_sum_to_the_bit`) ran through the climb.
The production fit collapses every library to one distribution (`one_stutter_for_every_library`,
`expected_counts.rs:1248`), so the claim matters only to the bench's held-out score. If Mi1 keeps
that code, add:

- `tract_likelihoods_of_two_libraries_at_multiplier_one_equal_their_summed_reads` — one sample, two
  libraries of slippage group 0 at multiplier 1, counts `[3, 7, 2]` and `[5, 1, 9]`; `TractLikelihoods::of`
  must give the same bits as one library holding `[8, 8, 11]`. It catches a `rows_by_distribution`
  that fails to sum rows sharing a distribution.

## 9. What's good

- The bench rename was done at the source (`held_out_run`'s variant naming, `fit_bench.rs:862-871`)
  rather than by string-patching outputs, so held-out rows, call lines and locus lines all moved
  together; the comparison above found no row where they disagree.
- `COUNTING_FIT_FLOOR` (`expected_counts.rs:89-93`) says what it is, where its value came from, and
  that the fit's own drivers do not read it — exactly the information the next reader needs.
- The examples that lost a mode say so in place, with the date and the arch section
  (`ng_joint_str_harness.rs:1365-1369`, `ng_joint_records_walk.rs:310-314`,
  `tests/slippage_curve_on_real_cells.rs:19-22`), instead of disappearing silently.
- The deleted `BelowTheFloor` fixtures were replaced with the refusal the new fit actually produces
  (`tracts: 0, NoSpanningReads`), not with a different made-up shape.

## 10. Commands to re-verify

Run from `/Users/jose/devel/pop_var_caller-str-penalised-fit`, scratch under `tmp/review_e1/`:

```
# copies
rsync -a --exclude target --exclude tmp ./ tmp/review_e1/src_copy/
git archive HEAD | tar -x -C tmp/review_e1/head_copy

# Mi1: should show no warning from pop_var_caller (lib) once fixed
(cd tmp/review_e1/src_copy && CARGO_TARGET_DIR=$PWD/../target_new cargo build --release)
(cd tmp/review_e1/src_copy && CARGO_TARGET_DIR=$PWD/../target_clippy_new cargo clippy --lib -- -D warnings)

# older parameters files, HEAD binary against E1 binary
bash tmp/review_e1/call_older_files.sh            # c6, hg002 x2, HG003 x2, HG004 x2

# the bench's counting fit before and after
bash tmp/review_e1/bench_b.sh
uv run --no-project python tmp/review_e1/compare_bench.py tmp/review_e1/bench_b

# reference fixture, the 144 runs below 2,000 samples
uv run --no-project python scripts/research/str_fit_fixture.py run --out tmp/review_e1/fixture_runs \
  --bench tmp/review_e1/target_new/release/examples/ng_str_fit_bench --jobs 6 \
  --only $(grep -v n2000 tmp/review_e1/fixture_runs_all.txt)
uv run --no-project python scripts/research/str_fit_fixture.py compare --require exact --partial \
  src/parameter_estimation/joint/ssr_fit/bench_fixtures/d_tract_prior_reference.jsonl.gz \
  tmp/review_e1/fixture_runs

# Mi2: the D − B column should hold a number once fixed
uv run --no-project --with numpy python scripts/research/str_fit_summarise.py \
  tmp/review_e1/summarise_new --compare B,D
```

`PROJECT_STATUS.md` was not updated: the dispatcher asked for the worktree to stay read-only apart
from this report.
