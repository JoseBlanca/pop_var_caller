# Code Review: estimation_memory_issue2
**Date:** 2026-09-26
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** issue 2 of `doc/devel/implementation_plans/estimation_memory.md` — `--str-param-estimates-at-once N`, the memory guess removed, `Scorer::refresh` drops old tables first
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** commit `4178f4b0` against `1d523425` (main), branch `strata-at-once`.
- **In-scope files:** [ssr_fit.rs](../../../../src/parameter_estimation/joint/ssr_fit.rs),
  [progress.rs](../../../../src/parameter_estimation/progress.rs),
  [estimate_parameters.rs](../../../../src/cli/estimate_parameters.rs) and its tests,
  `cross_platform_digests.rs` (the args literal), `examples/ng_ssr_strata_threads.rs`,
  `benches/ng_joint_fit_perf.rs`, the implementation report, the status block.
- **Categories dispatched:** three agents in their own worktrees — reliability + errors +
  refactor_safety; naming + idiomatic + smells + defaults; unsafe_concurrency + float_portability +
  module_structure + tooling + extras. The same grouping as issue 1's review, for the same reason
  (a few hundred lines in one module).

## 2. Verdict

Approve-with-changes. **No concurrency or float-portability defect**: the several-at-once schedule
gives the same bits as one at a time at 2, 4 and 7 at once on a six-stratum fixture; nothing in
`parameter_estimation` reads `rayon::current_num_threads`; and a probe counting live likelihood
tables measured a peak of exactly N for N = 1 to 4. The two Majors are about what the tests can
see.

## 3. Execution status

- clippy `-D warnings` clean; `cargo test --all-targets --all-features --no-fail-fast` 4,867
  passed, 3 failed (the two pre-existing probe examples).
- Reliability agent: 8 mutations, 7 survived, 1 of those changed no behaviour — **6 real
  survivors**, listed under M1 and Mi1.
- `cargo doc`, `cargo audit`: not run.

## 4. Open questions and assumptions

1. The library's progress line names the command-line flag. The plan's own example line does, so
   this is kept (Mi6).

## 5. Top 3 priorities

1. **M1** — nothing observes how many tables are held at once.
2. **M2** — the checksum test cannot reach the several-at-once schedule.
3. **Mi1** — nothing checks that the flag reaches the fit.

## 6. Findings

### Major

**M1: ssr_fit.rs, `fit_strata` — no test observes how many tables the fit holds at once, which
is the only thing `strata_at_once` changes.** **Categories:** reliability. Confidence: High. Every
test of the option compares fitted bits, and every schedule returns the same bits. Measured with a
probe counting live tables, caller's pool at 8 threads: running the N-at-once walks on the caller's
pool held **8 tables at N = 2**; a pool sized to the machine, 8; N = 2 run as N = 1, 1; removing
`self.prepared = None`, **2 at N = 1**; a walk splitting its tracts inside the N-thread pool, **6 at
N = 4**. The suite stayed green under all five. Fix: count live tables under `cfg(test)` and assert
the peak.

**M2: the checksum test cannot reach the several-at-once schedule.** **Categories:** extras,
reliability. Confidence: High. Its cohort holds one stratum (its progress line reads "fitting 1
strata … one stratum at a time"), and several at once needs more than one. So its unchanged
digests are evidence about the default schedule only. Fix: state the limit; the schedule's parity
is pinned by the six-stratum tests.

### Minor

- **Mi1** `estimate_parameters.rs` — nothing tests that the parsed value reaches `SsrFitConfig`;
  replacing it with 1 passes all 23 command and checksum tests.
- **Mi2** `neither_arm_moves_with_the_width_of_the_pool` cannot fail at N = 4: that schedule now
  builds its own pool and never sees the caller's width.
- **Mi3** `a_refreshed_scorer_holds_the_tables_a_fresh_one_builds` passes with the drop-before-build
  removed; its doc claims it checks that.
- **Mi4** the pool-build fallback is unreachable from any test.
- **Mi5** units: the help says 7.9 GB, the progress line prints the same table as 7.4 GiB; below
  about 50 MiB the line prints `0.0 GiB`, so a cohort of a few hundred samples gets no number to
  choose N from.
- **Mi6** the library's progress line names `--str-param-estimates-at-once`.
- **Mi7** help text: the `-h` summary says "each on one thread", wrong at the default; "stratum" and
  91 are not explained; the placeholder renders as `<STR_PARAM_ESTIMATES_AT_ONCE>`, not `N`.
- **Mi8** the default of 1 is written twice with nothing linking them.
- **Mi9** stale docs: `AcrossStrata` says "every stratum at once", and two docs bound memory by
  "threads × the largest stratum".
- **Mi10** a single stratum at N ≥ 2 prints "one stratum at a time" without saying why.
- **Mi11** `WhereTheThreadsGo` stays `pub` though no public item takes it.
- **Mi12** the pool is not capped at the work available; N above the walks idles threads.
- **Mi13** the `ng_joint_fit/strata` bench changed schedule under the same name.
- **Mi14** the report points at the commit message for the full-suite figures.

### Nits

- A `mut String` filled in one `match` and read in another.
- A test binding still named `every_class_at_once`.
- The test doc's "seven … leaves a thread idle" — 18 walks, so 7 threads are never idle.
- The example's module doc; a four-field tuple in a test; `..SsrFitConfig::default()` in a new
  literal; literal `**` in rendered help.

## 7. Out of scope observations

- Nobody has measured the wall-time cost of the new one-at-a-time default.

## 8. Missing tests to add now

- `n_strata_at_once_hold_at_most_n_tables` (M1, Mi3), `a_single_stratum_is_fitted_one_at_a_time_whatever_was_asked`,
  `a_pool_that_cannot_be_built_is_one_stratum_at_a_time_and_says_why` (Mi4), the flag reaching the
  config (Mi1), N above the walks available.

## 9. What's good

- `fit_strata` bounds memory by the pool it builds rather than by trusting the caller's.
- The refresh fix is one line and its comment says why it changes no number.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo test --lib -- ssr_fit estimate_parameters cross_platform census_fit`
