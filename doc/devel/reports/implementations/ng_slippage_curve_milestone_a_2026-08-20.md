# ng — the slippage curve, Milestone A: the curve, fitted and proven on real cells

**2026-08-20**, branch `ng-str-slippage-curve`. Steps A1–A5 of
[`str_slippage_level_curve.md`](../../ng/impl_plan/str_slippage_level_curve.md), against the design
in [`../../ng/spec/str_slippage_level_curve.md`](../../ng/spec/str_slippage_level_curve.md).

**Nothing in the pipeline calls any of this yet.** Milestone A builds the curve as a pure function
and proves it; wiring it into the fit is Milestone B.

---

## 1. Plan

Five steps, one commit each: the types, evaluating a curve, fitting one line at a fixed shape
number, choosing the shape number across a period's slippage groups, and running the whole thing
against the per-cell tables two real cohorts produced.

---

## 2. Assumptions — choices the design left open

**One, and it changes a fitted number at one period.** Spec §4.3 says the shape number is chosen
by "leaving each contributing cell out in turn, refitting that group's line without it, and
predicting it", and does not say whether that prediction *continues* the line or *holds it flat*
as a deployed curve does (spec §6).

It matters because leaving out the lowest or highest cell puts that cell outside what remains.
Measured three ways over both cohorts:

| | tomato p1 | HG002 p1 | HG002 p2 | HG002 p3 | HG002 p4 |
|---|---|---|---|---|---|
| continuing the line | **0.00** | **1.00** | **0.80** | **0.00** | **1.00** |
| holding flat | 0.00 | 1.00 | **0.70** | 0.00 | 1.00 |
| scoring only cells that stay inside | **0.15** | 1.00 | 0.80 | 0.00 | 1.00 |

**Chosen: continue the line.** It is the only variant that keeps every cell in the score — holding
flat predicts an end cell with its neighbour's value, which says nothing about the shape, and
scoring inside-only leaves 3 of tomato's 5 cells scored. It also reproduces the table the report
already published. **This does not leak into what a caller reads:** `level_at` still holds flat,
and the score is discarded once the rung is picked.

**Raised for the owner at Checkpoint A** — it belongs in spec §4.3, and the plan-driven loop does
not edit the design.

---

## 3. Changes made

**New: [`src/ng/parameter_estimation/joint/slippage_curve.rs`](../../../../src/ng/parameter_estimation/joint/slippage_curve.rs)**,
declared in [`joint/mod.rs`](../../../../src/ng/parameter_estimation/joint/mod.rs).

| what | it is |
|---|---|
| `RiseShape` | the shape number, a validated newtype over `[0, 1]`; 0 means each extra repeat multiplies the level, 1 means it adds |
| `FittedCell` | one cell's own answer as the fit sees it — repeat count, level, slipped reads — plus `relative_standard_error()` |
| `SlippageCurve` | a fitted line with its shape number, the repeat counts it was fitted over, its held-out error and its cell count; `level_at` and `reach` |
| `NoCurve` | why a cell set produced none: too few cells, one repeat count only, or a level that would fall |
| `PeriodCurves` | one period's shared shape number and one line per slippage group |
| `fit_line` | weighted least squares of `level ^ rise_shape` on repeat count |
| `choose_rise_shape` | the grid search, scored by leaving each cell out in turn |
| `LevelSource`, `CurveReach` | where an emitted level came from, and whether the cell sat inside the curve's range |
| `LEVEL_FLOOR`, `LEVEL_CEILING`, `RISE_SHAPE_RUNGS`, `MIN_CELLS_FOR_A_CURVE`, `DISAGREEMENT_KNEE` | the constants, each with its source in its doc comment |

**New fixtures:** `tests/data/slippage_cells/tomato_63_accessions_8mb.csv` (71 cells, 6 fitted) and
`hg002_300x_tier.csv` (137 cells, 55 fitted) — the raw `SSR_CELL_TABLE` output of the two runs
reported in [`str_slippage_shape_2026-08-20.md`](../../ng/reports/str_slippage_shape_2026-08-20.md),
copied unchanged.

**Two departures from the plan, both recorded in their commits:**

1. **The types live in their own module beside `ssr_fit.rs` rather than inside it.** The curve is
   self-contained arithmetic with its own tests, and `ssr_fit.rs` is already 1,800 lines.
2. **`SlippageCurveConfig` carries no knob for the per-cell weight.** It is the cell's slipped
   reads; the measurement says the ranking of families is identical under all four weights tried,
   and a knob nobody sets is how a measured decision gets changed by accident. The config does
   carry `draw_curves`, the switch Milestone B's parity oracle needs.

---

## 4. Tests added

**29 unit tests** in the module and **5 integration tests** in
[`tests/slippage_curve_on_real_cells.rs`](../../../../tests/slippage_curve_on_real_cells.rs).

The ones that would catch a real defect:

- **The fit finds the shape its cells were drawn with**, at both ends of the grid, to a held-out
  error under one part in a million. Neither end has more freedom than the other, so only held-out
  prediction can separate them.
- **A positive slope never lets a longer tract come back less slippery** — asserted at every one
  of the 21 rungs, over repeat counts 1 to 60, which spans well outside every fitted range.
- **Holding flat, with the number that made it the rule:** a curve fitted over 8–12 repeats
  reaches 244.7 at 30 repeats if continued, and stays at its 12-repeat value of 0.1827 when held.
- **A falling fit is refused**, and the refusal names which of the two floors turned the cells
  away — two cells for arithmetic, four for meaning.
- **The two cohorts return opposite shape numbers from their real cells**, which is the finding
  the whole family exists for.

**The per-period numbers are asserted, not narrated**, so the report's table cannot go stale:

| cohort, period | contributing cells | shape number | held-out error |
|---|---:|---:|---:|
| tomato, homopolymer | 5 | 0.00 | 12.36% |
| HG002, homopolymer | 23 | 1.00 | 4.39% |
| HG002, dinucleotide | 20 | 0.80 | 3.79% |
| HG002, trinucleotide | 4 | 0.00 | 31.10% |
| HG002, tetranucleotide | 7 | 1.00 | 53.34% |

**Every one of the five agrees with the independent analysis the report was written from**, to the
precision the report prints — a Rust reimplementation landing on the same answers as the Python
that produced the design.

---

## 5. Validation

Run in the container (`./scripts/dev.sh`):

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --all-features -- -D warnings` | clean |
| `cargo test --lib slippage_curve` | 29 passed |
| `cargo test --test slippage_curve_on_real_cells` | 5 passed |

**⚠ The aggregate `cargo clippy --all-targets --all-features` gate is red and pre-dates this
branch.** Three files this work does not touch fail it: `src/ng/run/cohort_merge/observation_cache.rs`
(a test assertion), `examples/ng_joint_duplicated_in_fit.rs`, and `examples/shared/synthetic_alignment.rs`.
`PROJECT_STATUS.md` already records the red gate.

**⚠ One commit briefly carried a clippy error and was amended.** `A4` was committed while a
`manual_is_multiple_of` lint was outstanding, because the gate's output was masked by a shell
pipeline in the check. Caught, fixed, and the commit amended before the next step.

---

## 6. Tradeoffs and follow-ups

- **Two-stage, not one-stage.** The curve is fitted to the per-cell answers rather than to the
  reads. Spec §1.2 states it; the research plan's D1 owns measuring the gap.
- **`fit_line` recomputes the transform per cell per rung.** 21 rungs × (cells + 1) fits per
  period, each a two-pass weighted least squares over at most a few hundred cells — microseconds
  against a repeat-tract fit measured in minutes. Not optimised, deliberately.
- **The four-cell floor is arithmetic and the tests now say what it costs.** HG002's
  trinucleotides clear it with exactly four cells and predict a held-out cell to 31.10%, against
  3.79% at the twenty-cell dinucleotides — a test asserts that ratio, so spec §11's doubt about the
  floor is pinned rather than remembered.
- **Nothing calls this yet.** Milestone B wires it into `fit_strata`, and B3 is the step whose
  failure would be a wrong genotype rather than a panic.
