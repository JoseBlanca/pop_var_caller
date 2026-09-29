# Fit precision, step B1 — the settled test as specified stops too early; stopped for the owner

> **Superseded — kept as the record of why the rule changed.** The owner approved the Newton rule of §5
> on 2026-09-29; step B1 was then reviewed, bounded to the parameters' intervals, and committed. What was
> built is in [fit_precision_b1_2026-09-29.md](fit_precision_b1_2026-09-29.md); "not committed" below
> describes the tree when this was written.

**Date:** 2026-09-29. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step B1.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §2, §3.3, §3.6 item 4. **Branch:** `fit-precision`,
**not committed**: the implementation is in the working tree (and `tmp/fit_precision/b1/b1_as_specified.patch`).

## 1. What was built

The rule exactly as spec §2 and §3.3 state it: the log-likelihood rule (a cycle gaining less than 10⁻⁴ a
position) as the trigger; the next cycle's first pass sums the information, again every 10 cycles; each
parameter's projected distance `|Δ|·λ/(1−λ) + |Δ|`, `λ` capped at 0.95 and zero when two moves point opposite
ways; settled below 0.1 of its standard error, a parameter without an error settled by definition; converged at
the end of a cycle whose log-likelihood rule holds and in which every parameter is settled. The relative-move
rule is gone (`JointFitConfig::stillness` becomes `settled_fraction`).

- `fit/settled.rs` (new): `SETTLED_FRACTION`, `MAX_CONTRACTION`, `ERROR_REFRESH_CYCLES`, `projected_distance`,
  `SettledTest`, `CycleVerdict`, nine unit tests.
- `fit/standard_errors.rs`: `StandardErrors::by_coordinate`, the errors in the fit's vector order; `named` uses it.
- `fit.rs`: `maximise` runs the test; `Alternation::step` can keep the information; the progress line names the
  pass the errors come from; `largest_relative_move` and the scale floors removed.
- `fit/information.rs`: `a_settled_fit_lands_within_a_tenth_of_an_error_of_a_long_one`, which fails (§2).

## 2. What was measured — the validation fails

**The oracle cohort** (spec §3.6 item 4: four tomato accessions at about 3 reads; `tmp/fit_precision/oracle_b1/`,
against a probe build run to 1,000 passes, `tmp/fit_precision/b1/fit_1000/`; `compare_settled.py`):

| start | new rule: passes, log-likelihood | 1,000 passes | shortfall |
|---|---|---|---|
| 1 | 57, −799,202.4 | −799,128.8 | 73.7 units |
| 2 | 36, −799,205.8 | −799,128.8 | 77.0 |
| 3 | 42, −799,203.0 | −799,128.7 | 74.3 |

The old rule's 200-pass limit left the best start at −799,130.2. Under the new rule every one of the 19
parameters with an error is more than 0.1 of it from the 1,000-pass fit; the invariant share 28 errors
(0.99441 against 0.99834, ± 0.00014), one homozygote excess 2.5, the fixed non-reference share 2.2.

**Why.** On the 1,000-pass trace the density's first shape sits near 0.25 to pass 60 (0.227 at pass 20, 0.254 at
57), is 0.49 at pass 100 and 15.5 by pass 200; the log-likelihood gains 72 units on the way. At pass 57 the fit
was on that plateau and still gaining 3.0 units a cycle (passes 57 to 60) — below the trigger's 200 units a
cycle (10⁻⁴ × 2 million positions). Two things then combine:

1. **The errors judged against were the plateau's, not the maximum's**: the invariant share ± 0.0048 at the
   returned point, against ± 0.00014 at the 1,000-pass point — 34 times wider — and the first shape ± 0.32
   against ± 36.6. Moves of about 10⁻⁵ a cycle look settled against errors that wide.
2. **The projection under-states a slow approach.** Where a path turns, its moves change sign and `λ` is taken
   as zero; where it crawls, the true contraction (about 0.98 to 0.99 a cycle on the drawn cohorts below) is
   above the 0.95 cap.

**Drawn cohorts** (the new test, 15 cohorts: five regimes × three seeds, each against the same cohort fitted
600 passes with `settled_fraction` 0; `tmp/fit_precision/b1/sweep*.log`):

| regime | spec rule: passes / furthest, errors | cap 0.99 | cap 0.99, no sign rule |
|---|---|---|---|
| 20 samples × 3 reads | 18–21 / 0.004–0.062 | 21–39 / ≤ 0.010 | 24–39 / ≤ 0.010 |
| 20 × 8, duplicated class | 28–39 / 0.001–0.81 | 31–69 / ≤ 0.025 | 37–69 / ≤ 0.015 |
| 4 × 3 | 12–43 / 0.06–0.81 | 33–49 / 0.06–0.76 | 33–51 / 0.04–0.76 |
| 4 × 8, duplicated class | 30–42 / 0.23–0.49 | 39–60 / 0.23–0.43 | 60–70 / 0.17–0.32 |
| 2 × 3 | 27–56 / 0.25–1.90 | 27–74 / 0.25–1.77 | 51–145 / 0.11–1.39 |

Nine of the 15 stop with a parameter at or beyond 0.1 of an error under the spec's rule; the parameters left
behind are the invariant share and the density's shapes (log-likelihood 0.03 to 0.86 units short). A higher cap
repairs the 20-sample cohorts and not the small ones, and neither variant touches the oracle's cause 1.

## 3. What is left for the owner

A design decision (spec §2, §3.3). The recommendation is in the chat report of this date; nothing here is
committed, and B2, B3 wait on it. The fixture checksums were not re-recorded (the oracle's fit lines move to
`759a48db…` / `956c5e09…` under the rule as specified).

## 4. The amended rule, measured (owner's go-ahead, 2026-09-29)

After merging main (`c37e14e9`; both cross-platform checksums unchanged on the merge alone). The candidate: the
contraction cap at 0.99, and the fit may stop only in a cycle gaining less than 0.005 log-likelihood units
(half of 0.1², what one parameter a tenth of an error from its maximum still has to give). Probe knobs, not
committed; `tmp/fit_precision/b1/sweep_cap*_gain*.log`, `oracle_cap0.99_gain0.005/`.

**Oracle:** starts 1 and 2 reach the 200-pass limit, start 3 stops at pass 191; the winner is 1.6 units below
the 1,000-pass fit, 4 of 19 parameters beyond 0.1 of an error (the carrier shapes 0.26 and 0.25, the duplicated
share 0.11, one homozygote excess 0.11); the invariant share 0.015. The 1,000-pass fit is itself still gaining
(0.2 units from pass 700 to 999).

**Drawn cohorts**, furthest parameter from the 600-pass fit, in errors (passes in brackets):

| regime | spec rule | cap 0.99 + gain 0.005 | cap 0.95 + gain 0.005 | cap 0.99 + gain 0.05 |
|---|---|---|---|---|
| 20 × 3 | 0.004–0.062 (18–21) | ≤ 0.010 (21–39) | ≤ 0.056 (18–21) | ≤ 0.010 (21–39) |
| 20 × 8, duplicated | 0.001–0.81 (28–39) | ≤ 0.015 (31–69) | ≤ 0.033 (28–57) | ≤ 0.025 (31–69) |
| 4 × 3 | 0.06–0.81 (12–43) | 0.06–0.54 (39–63) | 0.06–0.56 (12–54) | 0.06–0.76 (33–49) |
| 4 × 8, duplicated | 0.23–0.49 (30–42) | 0.22–0.42 (45–66) | 0.22–0.43 (45–60) | 0.23–0.43 (39–60) |
| 2 × 3 | 0.25–1.90 (27–56) | 0.23–1.42 (27–84) | 0.23–1.44 (27–72) | 0.25–1.77 (27–74) |

At 20 samples the cap is what matters; the gain condition changes nothing there. **At 2 and 4 samples no rule
tried meets 0.1 of an error**: the fits stop while gaining under 0.005 units a cycle with 0.01 to 0.6 units
still to gain, a contraction near 0.99 or above. A rule that judges progress by its rate cannot certify the
distance left when the fit crawls that slowly.

## 5. The Newton distance, measured (owner's go-ahead, 2026-09-29)

Replacing the projection: the information pass also sums each parameter's slope, and the information
solved against it (`I⁻¹ g`, the whole matrix up to 20 samples, the arrow of blocks above) gives each
parameter's distance to the maximum; a parameter at a bound whose step points out of its interval is held
there and the rest solved again. Settled when every distance is below 0.1 of that parameter's error; the fit
stops at the end of that cycle. The information pass runs every cycle once the log-likelihood gate has held.
Probe-wired (`newton_step` in `fit/standard_errors.rs`; `tmp/fit_precision/b1/sweep_newton.log`,
`oracle_newton/`).

| regime | converged | passes | furthest from the 600-pass fit, errors |
|---|---|---|---|
| 20 × 3 | 3 of 3 | 18–33 | 0.004–0.043 |
| 20 × 8, duplicated | 3 of 3 | 25–69 | 0.022–0.037 |
| 4 × 3 | 3 of 3 | 12–180 | 0.047–0.088 |
| 4 × 8, duplicated | 1 of 3; 2 at the limit | 129–199 | 0.045–0.0995 |
| 2 × 3 | 0 of 3; all at the limit | 198–200 | 0.085–0.28 |

**Every fit it called converged is within 0.1 of an error** (the furthest 0.088); the slow ones run to the
limit and say so. **Oracle:** the Newton distance sees the plateau — the invariant share 14 to 18 errors from
the maximum at every check to pass 103 — and no start converges: all three stop at the 200-pass limit, at the
old rule's point (1.5 units below 1,000 passes; the carrier shapes 0.26 and 0.25 errors off, the invariant share
0.015). At pass 198 the winner's Newton distance is still 6.3 errors on the duplicated share, where the
1,000-pass fit differs from it by 0.11: far from the maximum the quadratic model over-states some distances
(a distance longer than the parameter's interval is not clipped). The cost of an information pass every cycle
is not yet measured: 186 of about 600 passes on the oracle, at about 1.4 times a plain pass each (step A5).
