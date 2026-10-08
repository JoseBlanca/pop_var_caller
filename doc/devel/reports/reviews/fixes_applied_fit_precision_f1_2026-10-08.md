# Fixes applied: fit_precision_f1
**Date:** 2026-10-08 · **Review:** [fit_precision_f1_2026-10-08.md](fit_precision_f1_2026-10-08.md)

| finding | outcome |
|---|---|
| B1 — blended fall-off untested | **fixed.** The writer test has a third slippage group: level from the curve whole, shorter share the stratum's, fall-off a blend. Measured: putting the old rule back for the fall-off alone now fails it (`fall_off: None` against `Some(0.09)`). |
| M1 — reader's share rule unpinned | **fixed.** The fixture's blended shorter share carries an error (0.013, both golden files); the validator test refuses an error beside a shorter share taken whole from the curve and accepts a blended fall-off with one. |
| M2 — four stale doc comments | **fixed**, each rewritten to the amended rule. |
| M3 — wrong figures and mechanism | **fixed.** (a) the plan names each start's furthest parameter and its distance (464.9, 44.3, 407.8); (b) the weight range is scoped to the 22 one-base strata in the plan and spec; (c) the weight is described as how much of the number is the curve's, in the file note, the spec and the plan; (d) the doc now gives both ends, measured: 2 × 10⁻⁵ to 1 × 10⁻⁴ at kimura's one-base strata, 0.05 to 0.9999 (median 0.6) over the four-accession oracle's 44 blends; (e) the note says an origin "can carry" an error, absent where the tracts do not determine the number. |
| Mi1 — deny-list rule | **fixed.** `LevelSource::holds_the_own_fit`, `ShareSource::holds_the_own_fit`, and the validator's `a_level_the_own_fit_went_into` / `a_share_the_own_fit_went_into` are exhaustive matches; the validator's parameter is now positive (`the_own_fit_went_into_the_number`). |
| Mi2 — bare row indices | **fixed**, named constants; the misdirected "Its fall-off" comment rewritten. |
| Mi3 — spec table row | **fixed**: the row states the amended rule; the amendment quotes the old wording. |
| Mi4 — level from the curve untested | **fixed** by B1's third group (`level: None`). |
| Mi5 — the log call untested | **won't fix.** Asserting a log line needs a test that captures the run's log, which this module has not got; the line was seen on the oracle cohort's real run instead: `strata: 68 in all; 15 fitted …, 22 furnished …, 31 refused — 19 …, 12 …, holding 1 to 5`. |
| Nits | **fixed**: the clause variable renamed, the line reworded with no verb agreement and one phrasing for the refused, "holding 3" for a single stratum (tested). **Kept:** the configured floor is printed; every refusal's own floor is that same value. |

Re-run after the fixes: fmt, clippy `--all-targets --all-features -D warnings`, the affected modules' tests, the
oracle and the cross-platform checksums (results in the implementation report and the commit).
