# Review: D2a, the fitting bench's psp mode (2026-10-16)

**Scope.** The uncommitted change in `/Users/jose/devel/pop_var_caller-str-penalised-fit` on branch
`str-penalised-fit`, HEAD `8c6ce6f5`: `examples/ng_str_fit_bench.rs`,
`src/parameter_estimation/joint/ssr_fit/bench_fixtures/fit_bench.rs`, `src/run/census_fit.rs`,
`src/run/mod.rs`. It was judged against plan step D2a (`doc/devel/ng/impl_plan/str_penalised_fit.md`)
and spec `str_penalised_fit.md` §4 item 4. The allocator lines at the top of the bench were excluded,
as asked. One reviewer, no sub-agents. Every run used 4 threads at most. Scratch files are in
`tmp/review_d2a/`.

## Verdict

**The psp mode fits on exactly the evidence `pvc estimate-parameters` fits on, and the drawn mode is
unchanged.** No Blocker or Major findings. There are four Minor findings and two Nits. None of them
moves the 59-sample tomato2 numbers already produced in `tmp/d2a/`.

## What was verified, and how

### 1. Same evidence as `estimate-parameters`: verified bit for bit

- **Code path.** The psp mode calls the same functions as `fit_and_assemble` → `fit_a_cohort`:
  `open_the_censuses`, which puts samples in name order through `in_sample_name_order`;
  `refuse_another_selection`; `every_read_group_pooled`; and `the_tract_strata_of_a_cohort`, which
  applies the repeat-tract guard. The SNP/indel config also matches: `JointFitConfig::default()` is
  production's `ordinary_position_config` at ploidy 2 with contamination estimated, field by field.
  The penalised fit matches as well: `ThePenalisedFit::production()` calls
  `fit_period_by_penalised_counting`, which calls `fit_one_period`, the same function
  `fit_repeat_tracts_logging` calls for each motif length. Its `with_allele_span(ALLELE_SPAN)` changes
  nothing, because production already uses `ALLELE_SPAN`.
- **Probe.** I added an environment-gated block to a copy of the bench
  (`tmp/review_d2a/probe_copy/examples/ng_str_fit_bench.rs`). On the bench's own evidence and excess,
  it runs `ThePenalisedFit::production()` over every locus. I ran it on the 4 oracle psps with
  `--parameters tmp/oracle_C10/fitted.parameters.toml`. **All 49 fitted classes' `level`,
  `shorter_share` and `fall_off` are bit-identical to that file's `slippage_by_stratum_and_group`
  rows**: 49 of 49, compared as f64 (`tmp/review_d2a/compare_full.py`, `probe_full.tsv`).
- **Sample order.** I ran the psps once by directory and once named in reverse order. The two locus
  files are byte-identical (`smoke.loci.tsv`, `smoke_rev.loci.tsv`). They are also byte-identical to
  the implementer's `tmp/d2a/smoke4_file.loci.tsv`.
- **`--parameters` against the default (excess from the SNP/indel fit).** The implementer's two
  locus files, `smoke4_fit` and `smoke4_file`, are byte-identical (`cmp`). The excesses are joined by
  sample name, a sample missing from the file is refused, and the file writes shortest-round-trip
  f64s. **They are equivalent for a file whose coefficients are `fitted_here`.** All 59 rows of
  `tmp/perf_ep/t2_59.parameters.toml` are `fitted_here`, and the largest value is 0.9999999992, below
  1. Finding 1 covers the files for which the two are not equivalent.

### 2. Held-out correctness

- **Whole loci, each held out exactly once.** `FoldPlan::of` gives every locus exactly one fold, and
  `without_fold` drops whole `TractReads`. The new test checks both. **Mutation:** making
  `without_fold` keep the held-out loci fails the test at `fit_bench.rs:1888`.
- **Leaks.** I found none that favours one method:
  - **F, the allele-length shares, the stutter curves and the fallback concentration** are fitted
    from `training` only. That covers `fit_period(&training, ...)`, B's `with_the_curves(&raw,
    &training, ...)`, and the median that `numbers_for_scoring` takes over the fold's fit.
  - **The homozygote excess** is one number per sample, from roughly 1.9 million ordinary positions
    over the whole cohort, held-out loci included. It is fixed before folding and is the same for
    both methods, which is also what production does. This is not a bias between the methods.
  - **The repeat-tract guard and the selection** are applied per tract, before folding, to all loci.
    They choose which loci exist; they fit nothing.
  - **The class-wide totals that `without_fold` copies from the full class**
    (`fit_bench.rs:851–860`: guard counts, `bases_compared`, `mismatching_bases`) include held-out
    loci. Neither fit nor the scorer reads them; grep finds no reader in `penalised_fit.rs`,
    `fit_bench.rs`'s scorer or the counting fit. They are not a leak today, but see Nit 2.
- **Monte Carlo draws** are shared across methods. The seed is `locus_seed(settings,
  class_seed_across_motif_lengths, locus)`, with the class seed shifted left by 40 bits, so motif
  length and repeat count cannot collide below 2¹⁶ repeats.

### 3. The pairing script, `tmp/d2a/heldout_real.py`

- **Correct.** It keys loci by (motif length, repeat count, locus index in its class). Both methods
  share one `classes` vector, so the indices agree. No key repeats within a variant: I checked the
  59-sample file. The ± is the sample standard deviation (n − 1) divided by √n.
  `tmp/review_d2a/pair_check.py` recomputes it independently and gets the same values (+0.1364 ±
  0.0410 against `B_curves`, +0.2130 ± 0.0488 against `B_raw`).
- Finding 3 covers what the script does not show.

### 4. The drawn mode is unchanged

- **Cell two-base, 60 samples, depth 3, seed 1, methods `B,D_tract_prior`, `--both-fallbacks`.** I
  built the working tree and a `git archive HEAD` copy and ran both. Every non-timing row is identical
  (2,185 lines). The truth files are identical. The locus files are identical in their first eight
  columns; the new ninth column is `period`. The `D_tract_prior` rows are also identical to
  `tmp/d1/grid/p2_n60_d3_pcr_s1.tsv`. That D1 file holds no `B` rows, so `B` was compared against
  HEAD only.
- `scripts/research/str_fit_fixture.py selftest`: **15 checks passed**. The fixture line and its
  flags (`--fixture-out`, `--period`) were not touched by the refactor.

### Commands run (in the exported copy, 4 jobs)

- `cargo fmt --check`: clean.
- `cargo clippy --release --all-targets --all-features -- -D warnings`: the five known errors and
  nothing new. They are `fit_bench.rs:460` (useless conversion), `edge_deletion_gathering.rs:201`
  (single-pattern match), `expected_counts.rs:2138` (`f64::powi` disallowed), and
  `penalised_fit.rs:2915` and `2916` (needless range loop). **Because the library fails clippy,
  clippy never reaches the example**, so the bench itself has not been linted.
- `cargo test --release --features bench-fixtures --lib -- fit_bench:: run::census_fit`: 36 passed.

## Findings

### Minor

- `examples/ng_str_fit_bench.rs:1068–1090` — **[Minor] `--parameters` accepts a file whose
  coefficients are not the excess production's tract fit used**
  - **Confidence:** High (read in code). No file triggering it has been made.
  - **Problem:**
    - `fit_a_cohort` gives the tract fit the raw `hom_excess`, whatever `--inbreeding` said.
    - The file's coefficient is `DeclaredInbreeding::of_each_sample_over(...)`, which differs in
      three cases:
      - **(a)** The run was given `--inbreeding X`. The file then holds X, `supplied`.
      - **(b)** A sample's excess came out at exactly 1.0. `HomozygoteExcess` allows [0, 1] and
        `InbreedingF` only [0, 1), so `fitted_inbreeding_of` drops the sample and the file writes
        0.0, `defaulted`.
      - **(c)** The file was fitted on another cohort, for example a superset of the psps. Extra
        names are ignored and the values are that cohort's.
    - `homozygote_excess_from_file` reads `value` and ignores `warrant` and `[fitted_from]`.
  - **Why it matters:** In each case the held-out fits weigh genotypes differently from
    `estimate-parameters`, with no message. That breaks "the production fit" in D2's comparison.
    tomato2's file is unaffected: all 59 rows are `fitted_here` and the largest is 0.9999999992.
  - **Suggested fix:** Refuse any row whose `inbreeding_coefficient.warrant` is not `FittedHere` or
    `Defaulted`. Refuse a file whose sample list differs from the cohort's. A `defaulted` 0.0 from
    case (b) cannot be told apart from a real default, so the bench's documentation should name that
    case; the robust fix is to keep the default path for D2.

- `examples/ng_str_fit_bench.rs:885–995` — **[Minor] The bench's held-out wiring has no test**
  - **Confidence:** High.
  - **Problem:** The new test covers `fold_plans_across_motif_lengths` and `without_fold`. It does
    not check that `held_out_run` scores fold *f*'s loci with the fit made without fold *f*. A
    mutation such as `*locus_fold == (fold + 1) % settings.folds` at line 919 compiles, and no test
    fails, since the example has no tests. Only a by-hand comparison against recorded outputs would
    catch it.
  - **Why it matters:** That one line decides whether the held-out score is held out at all.
  - **Suggested fix:** Move `HeldOutClasses` and `held_out_run` into `fit_bench.rs`, next to
    `FoldPlan`. Add a test with a stub `PeriodFit` that records which loci it was given; assert that
    no locus scored in fold *f* was among them.

- `tmp/d2a/heldout_real.py:33–35` and `tmp/d2a/impl_report.md` — **[Minor] The pairing hides
  unscored loci in one direction, and the ± assumes independent loci**
  - **Confidence:** High (measured).
  - **Problem 1: missing counts.** The script counts loci that only `D_tract_prior` scored, but not
    loci that only the other method scored. On tomato2 the second count is 0 for both `B` variants
    (`tmp/review_d2a/coverage.py`), but nothing would show it otherwise. Separately, 19 of the 2,392
    loci folded for scoring were scored by no method: `D_tract_prior` scored 2,373. Their classes
    have no reads in the training folds, so `slippage` is `None` and `numbers_for_scoring` skips
    them. The implementation report says "2,392 held-out loci scored"; the right figure is 2,373 for
    the production fit.
  - **Problem 2: the standard error.** The ± treats loci as independent. Loci in one class share a
    fitted set of numbers, so they are not. A standard error clustered by class is about 30% wider:
    0.0527 against 0.0410 for `B_curves`, and 0.0615 against 0.0488 for `B_raw`
    (`tmp/review_d2a/pair_check.py`). The sign holds either way: the mean difference is 2.6 and 3.5
    clustered standard errors from zero.
  - **Suggested fix:** Print the "only old scored" count and the count scored by none. Give D2 the
    class-clustered standard error, or state that the ± is over loci treated as independent.

- `examples/ng_str_fit_bench.rs:1255–1258` — **[Minor] The psp mode silently ignores the drawn
  mode's flags**
  - **Confidence:** High.
  - **Problem:** `--psp` conflicts only with `--fixture-out` and `--truth-out`. `--excess`,
    `--samples`, `--depth`, `--truth`, `--loci`, `--repeats`, `--between-locus-share`,
    `--library-multipliers` and `--multipliers-told` are accepted and do nothing. `--excess` is the
    one that misleads: a user who sets it believes the fit used it.
  - **Suggested fix:** Add these flags to `conflicts_with_all` on `psps`.

## Nits

- **Nit 1.** `examples/ng_str_fit_bench.rs:1249`: `peak_rss_mib` is written under the last motif
  length's columns. This is documented, but a reader grouping rows by `period` will attribute the
  whole process's peak to motif length 6. A `period` of `all` would be clearer.
- **Nit 2.** `fit_bench.rs:851–860`: `without_fold` copies class-wide totals that include the
  held-out loci. They are harmless today (see §2), but a comment there would stop a future fit that
  reads `bases_compared` from leaking unnoticed.

## Scratch

All under `/Users/jose/devel/pop_var_caller-str-penalised-fit/tmp/review_d2a/`:

| What | Files |
|---|---|
| Exported trees | `src_copy/` (working tree), `head_copy/` (HEAD) |
| Probe tree | `probe_copy/`: full-fit probe, plus the leak mutation at `fit_bench.rs:850` gated on `MUT_LEAK` |
| Build directories | `target/`, `target_head/` |
| Drawn-mode check | `new_p2.*`, `head_p2.*` |
| Psp-mode runs | `smoke.*`, `smoke_rev.*`, `probe_full.*` |
| Scripts | `compare_full.py`, `coverage.py`, `pair_check.py` |
