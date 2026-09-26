# Issue 2: the run says how many strata are fitted at once

*Implementation report, 2026-09-26. Plan:
[estimation_memory.md](../../implementation_plans/estimation_memory.md) §5. Branch
`strata-at-once`.*

## Plan

The repeat-tract fit has two schedules, which return the same bits:

- **one stratum at a time**, every thread on that stratum's tracts — memory: one likelihood table;
- **several strata at once**, each fitted on one thread — memory: one table a stratum running.

Until now the second was the default, spread over every thread of the machine, and commit
`b435adfe` switched to the first when the machine's `MemAvailable` said the tables would not fit.
The owner's rule is that a setting deciding whether a run fits in memory is stated by the user, not
guessed from the machine. So:

- `estimate-parameters --str-param-estimates-at-once N`, default 1;
- N = 1 is one stratum at a time; N ≥ 2 runs N strata at once on a pool of exactly N threads built
  for the step;
- the guess (`arm_that_fits`, `memory_available`) is deleted, not kept as a fallback;
- the stage's first progress line prints the largest stratum's table and what N of them cost;
- `Scorer::refresh` drops its old tables before building new ones, so a walk holds one table, not
  two, and the printed size is the real one.

## Assumptions and deviations

- **`SsrFitConfig::threads` is replaced, not extended.** The field is now
  `strata_at_once: NonZeroUsize` (default 1). The plan allowed either; one field keeps the two
  schedules from being described twice. `WhereTheThreadsGo` stays as the internal name for the two
  schedules, which the scorer and the quadrature still branch on.
- **A pool that cannot be built falls back to one stratum at a time**, and the progress line says
  so with rayon's error. The plan did not say. This keeps `fit_strata`'s signature, and the
  fallback is the schedule with the smallest memory, which returns the same bits.
- **N ≥ 2 with a single stratum runs one stratum at a time**, as before this change: with one
  stratum there is nothing to spread N threads over.
- **The table-size test from the deleted `the_arm_that_fits` module is kept**, in a module renamed
  `the_largest_table`, because `largest_table_bytes` stays for the progress line. The plan said to
  delete the module; its other two tests tested only the deleted guess.
- **Benchmarks and examples that take `SsrFitConfig::default()` now time one stratum at a time**,
  the new default. `examples/ng_ssr_strata_threads.rs` keeps its two arms: `across_classes` sets
  `strata_at_once` to `THREADS`, `across_tracts` to 1.
- The progress line names the command-line flag from inside the library, as the plan's example
  line does.

## Changes made

- [ssr_fit.rs](../../../../src/parameter_estimation/joint/ssr_fit.rs): `SsrFitConfig::strata_at_once`;
  `fit_strata` builds an N-thread pool and runs `every_stratum_from_every_start` inside it, or fits
  one stratum at a time; the new progress line; `arm_that_fits` deleted; `Scorer::refresh` sets
  `self.prepared = None` before building.
- [progress.rs](../../../../src/parameter_estimation/progress.rs): `memory_available` deleted
  (nothing else used it).
- [estimate_parameters.rs](../../../../src/cli/estimate_parameters.rs): the option, passed into the
  fit's config.
- `cross_platform_digests.rs` and `estimate_parameters/tests.rs` build the args struct literally,
  so each gains the new field at its default. The checksum constants are not touched.
- [ng_ssr_strata_threads.rs](../../../../examples/ng_ssr_strata_threads.rs): the two arms expressed
  as a count.

The first progress line has this shape. The figures here are illustrative — the kimura cohort's
141 strata, and a 5,000-tract stratum's table at 2,169 samples by the arithmetic above — not a
measurement. Sizes below 1 GiB are printed in whole MiB, above it in GiB to one decimal.

```text
estimating: repeat-tract fit: fitting 141 strata from 3 starting point(s) each, one stratum at a time; the largest stratum's table is 7.4 GiB, so --str-param-estimates-at-once N needs about N × 7.4 GiB
```

The command's help, as the review left it, is in the `--str-param-estimates-at-once` doc comment
in `estimate_parameters.rs`.

## Tests added or updated

- `ssr_fit::tests::any_number_of_strata_at_once_gives_the_same_bits` — six strata of unequal size
  fitted at 1, 2, 4 and 7 at once (7 is more than there are strata) give identical bits.
- `ssr_fit::tests::a_refreshed_scorer_holds_the_tables_a_fresh_one_builds` — a scorer that moved from
  slippage level 0.05 to 0.12 holds, bit for bit, the tables a new scorer builds at 0.12, and they
  differ from the 0.05 tables.
- `cli::estimate_parameters::tests::strata_are_fitted_one_at_a_time_unless_the_run_asks_for_more`
  — the default is 1, `4` is read as 4, `0` is refused when the command line is parsed.
- The existing parity tests now compare 1 against 4 at once, and 1 and 4 at pool widths one and
  eight.
- **The checksum test passes unchanged** at the default. **It reaches only one stratum at a
  time**: its cohort holds one stratum, and several at once needs more than one. Its unchanged
  digests are therefore evidence about the default schedule; the other schedule's bits are pinned
  by the six-stratum tests above.

After review ([review](../reviews/estimation_memory_issue2_2026-09-26.md),
[fixes](../reviews/fixes_applied_2026-09-26_v2.md)):

- `n_strata_at_once_hold_at_most_n_tables` — counts live likelihood tables (a test-only census,
  `table_census`) and requires at most N, and at least 2 when N ≥ 2, at N = 1 to 4 with the
  caller's pool at 8 threads. **It is the only test that sees what the option is for**, since every
  schedule returns the same bits. It fails, as it should, on the two defects review measured:
  without the drop-before-build ("1 at once held 2 tables together") and with the walks on the
  caller's pool ("2 at once held 8 tables together").
- `a_single_stratum_is_fitted_one_at_a_time_whatever_was_asked` — one stratum at N = 4 holds one
  table and gives the bits of N = 1.
- `a_pool_that_cannot_be_built_is_one_stratum_at_a_time_and_says_why` and
  `a_table_size_reads_in_mebibytes_below_a_gibibyte`.
- The CLI test also checks that the parsed count reaches the fit's config.
- `neither_arm_moves_with_the_width_of_the_pool` became
  `one_stratum_at_a_time_does_not_move_with_the_width_of_the_pool`: several at once no longer runs
  on the caller's pool, so varying that pool's width tests nothing there.

## Validation

In the dev container, after review:

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings` → clean;
- `cargo test --lib -- ssr_fit estimate_parameters cross_platform census_fit` → 65 passed, 0 failed;
- `cargo test --all-targets --all-features --no-fail-fast` → 4,871 passed, 3 failed; the three are the pre-existing failures in `examples/ng_generic_loci_dump.rs` and `examples/ng_ssr_loci_dump.rs` described in [issue 1's report](estimation_memory_issue1_2026-09-26.md).

## Tradeoffs and follow-ups

- At N below the core count, the other cores are idle during the repeat-tract fit, as the option's
  help says. Plan §5 accepts that as the cost of the simple scheme.
