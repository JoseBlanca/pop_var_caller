# Fix Application Report: estimation_memory_issue2_2026-09-26.md

**Date:** 2026-09-26
**Source review:** `doc/devel/reports/reviews/estimation_memory_issue2_2026-09-26.md`
**Source state reviewed against:** `4178f4b0`, branch `strata-at-once`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

- Review totals: 0 Blockers, 2 Majors, 14 Minors, 7 Nits.
- Applied: 2 Majors, 11 Minors, 3 Nits. Won't fix: 3 Minors (Mi6, Mi11, Mi12) and 4 Nits, reasons
  below.
- Unresolved high-priority findings: none.

## 2. Findings table

| ID | Final status | What was done |
|---|---|---|
| M1 | Applied | A test-only census of live tables (`table_census`, keyed by the evidence's address so concurrent tests do not move it) and `n_strata_at_once_hold_at_most_n_tables`: at N = 1 to 4 with the caller's pool at 8, at most N tables and, for N ≥ 2, at least 2 |
| M2 | Applied (documented) | The implementation report and the commit say the checksum test reaches only one stratum at a time; the six-stratum parity tests pin the other schedule |
| Mi1 | Applied | `repeat_tract_config(args)` extracted; the CLI test asserts 4 reaches the config and the default is `DEFAULT_STRATA_AT_ONCE` |
| Mi2 | Applied | Renamed `one_stratum_at_a_time_does_not_move_with_the_width_of_the_pool`, over N = 1 only, with the reason in its doc |
| Mi3 | Applied | Doc reworded to what it checks; the memory claim is M1's test |
| Mi4 | Applied | `a_pool_for(at_once, strata, build)` extracted; a test builds with a spawn handler that fails |
| Mi5 | Applied | Help and config doc in GiB (0.34 GiB and 7.4 GiB); the progress line prints whole MiB below 1 GiB (`table_size`, tested at 347 MiB and 7.4 GiB) |
| Mi6 | Won't fix | The plan's example progress line names the flag; the line exists to tell a user of that command what to type |
| Mi7 | Applied | Help rewritten: summary "1 uses the least memory", stratum and the 91 explained, `value_name = "N"` |
| Mi8 | Applied | `pub const DEFAULT_STRATA_AT_ONCE`, used by `SsrFitConfig::default()` and by clap's `default_value_t` |
| Mi9 | Applied | `AcrossStrata` and `every_stratum_from_every_start` docs rewritten; the second now states that the bound rests on a walk making no call into the pool |
| Mi10 | Applied | "one stratum at a time, since there is 1 to fit" |
| Mi11 | Won't fix | Making `WhereTheThreadsGo` crate-private breaks intra-doc links from the public config's docs; no public item takes it, so it costs nothing to leave |
| Mi12 | Won't fix | N is the run's statement; an idle rayon thread holds no table. The same-bits test covers N = 7 and 64 from one starting point (6 walks) |
| Mi13 | Applied | The bench's config doc says its results from before 2026-09-26 timed the other schedule |
| Mi14 | Applied | The report states the full-suite figures |
| Nits | 3 applied (the `mut String` is gone with `a_pool_for`; `four_at_once`; the "seven idles a thread" doc corrected and one-start cases added), 4 won't fix (the example doc already names both arms; a four-field tuple kept in one test; the `..default()` literal is the pattern throughout; `**` in help matches the file's other options) | |

## 3. Questions asked and answers

None.

## 4. Verification of the new tests against the defects they exist for

Run on this tree, each defect put in by hand and the file restored byte for byte afterwards
(`cmp` against a copy):

- `self.prepared = None` removed → `n_strata_at_once_hold_at_most_n_tables` fails: "1 at once held
  2 tables together".
- the several-at-once walks run on the caller's pool instead of the built one → fails: "2 at once
  held 8 tables together".

## 10. Commands run

- `./scripts/dev.sh cargo test --lib -- ssr_fit estimate_parameters cross_platform census_fit` →
  65 passed, 0 failed.
- `./scripts/dev.sh cargo test --lib -- n_strata_at_once_hold` under each defect above → 1 failed,
  as stated.
