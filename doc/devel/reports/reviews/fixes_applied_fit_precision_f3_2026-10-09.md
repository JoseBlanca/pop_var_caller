# Fixes applied: fit_precision_f3
**Date:** 2026-10-09 · **Review:** [fit_precision_f3_2026-10-09.md](fit_precision_f3_2026-10-09.md)

| finding | outcome |
|---|---|
| B1 — no test reaches the new branches | **fixed** for the best-point return: `every_start_at_the_pass_limit_returns_the_best_point_it_reached` checks every start of both fits on the 30-sample cohort; switching the return off fails it (start 2: −12,292.7 against a best of −12,267.3). The fallback has no test because it was removed. |
| M1 — the fallback replays one cycle | **resolved by removal (owner, 2026-10-09).** First re-judged against the cycle before rather than the best, with a plain cycle's landing left to stand: it still fired 19 times a start on that cohort and gained at most 0.4 units, because a loss from the alternation's own step recurs when the step is taken again. Recorded in the plan. |
| M2 — the go-back target can be a jumped-to point | **resolved by removal.** (Before removal it went back to the best cycle start only.) The best point a start returns still includes a jump that held. |
| M3 — the cap clamps each shape on its own | **fixed.** The cap is a box the Newton step is kept inside, projected as at a bound; tested exactly (`no_shape_moves_more_than_a_factor_of_two_in_a_pass`). |
| Stale docs | **fixed**: the module doc, `maximise`, `JUMP_SLACK`, `trace_the_returned_fit`, `fit_trace`. |
| "Converged" with a best point returned | **fixed**: only a start at the pass limit returns its best point. |
| Last pass's time covers two passes | **fixed**: the clock restarts before the pass at the best point. |
| Cap test asserts through the constant | **fixed**: from (1, 1), asked for (4, 1), the first shape stops at exactly 2.0. |
| `clamp` panics on a non-number shape | **fixed**: such a shape moves nothing (tested). |
| Tuple / names | **fixed**: `MeasuredPoint`, `best_above_the_last`; the flags went with the fallback. |
| Misplaced test doc | **fixed.** |
| Departures not recorded | **recorded** in the plan and the implementation report (the cap on the log scale with the step unchanged; the fallback's removal). |
| Nits | Units added to the log line; `MAX_SHAPE_FACTOR_A_PASS` gives its source; `maximise`'s length **kept** (the removal shortened it). |

Re-run after the fixes: fmt, clippy, the fit tests (119 passed), the cross-platform checksums (one fixture
re-recorded, explained in its doc), the oracle (unchanged), the full suite.
