# ng — the slippage level as a curve across repeat count: implementation plan

**Status:** draft, 2026-08-20. The build order for the design settled in
[`../spec/str_slippage_level_curve.md`](../spec/str_slippage_level_curve.md): the slippage level
stops being an independent per-cell number that thin cells borrow, and becomes a curve in repeat
count with its two coefficients fitted per slippage group and its shape number shared across
groups at each motif period. **This is not a place for new design** — every open question is in
the spec §11, and none of them blocks a step here.

It builds the first piece of the research plan
[`str_slippage_across_repeat_count.md`](str_slippage_across_repeat_count.md); that plan's steps A
and B are done and reported in
[`../reports/str_slippage_shape_2026-08-20.md`](../reports/str_slippage_shape_2026-08-20.md).

---

## Scope

**In:** a `SlippageCurve` type and its fit inside
[`src/ng/parameter_estimation/joint/ssr_fit.rs`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs);
the per-period choice of the shape number; the blend between a cell's own level and the curve's;
holding flat outside the fitted range; the emitted outcome type gaining a cell that carries a
level with no shares; the provenance the spec §8 requires; narrowing the borrowing rule to the two
shares; the amendments to
[`../spec/parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) §4.3 and §4.5 and to its
architecture sibling.

**Out (with the plan that owns each):**

- **Smoothing the direction split, the fall-off and the substitution rate** — the research plan's
  C4. They keep today's per-cell fit and today's borrowing.
- **Fitting the curve directly against the reads** (one-stage) — the research plan's D1.
- **Partial pooling of the two coefficients across periods**, for periods 5 and 6 — the research
  plan's D2.
- **A read-group axis that carries library identity across a cohort** — the census specification
  [`../spec/parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md).
  Until it exists, a cohort of single-library samples fits one curve, and the spec §3 says so.
- **The monotonicity merge** — it does not exist on this route; `merge_until_monotone` is the
  per-sample route's
  ([`ssr/mod.rs:1479`](../../../../src/ng/parameter_estimation/ssr/mod.rs#L1479)).
- **Genotype concordance under the new parameters** — the research plan's F2, which needs the
  caller and not just the fit.

---

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The arithmetic before the plumbing.** The curve fit is a pure function from a list of
  (repeat count, level, weight) to three numbers. It is built and tested against drawn cells with
  a known curve — no census, no reads — before anything in `fit_strata` calls it. A wrong curve
  is a wrong genotype, not a panic, so it gets proven where a test can see it.
- **Verify against ground truth, and the ground truth already exists.** Both cohorts' per-cell
  tables were produced before a line of this was written
  ([the report](../reports/str_slippage_shape_2026-08-20.md)). With the curve switched off, every
  cell's own fitted level must come back unchanged — this plan touches stage one nowhere.
- **Isolate the step whose failure is silent.** B3 changes what number reaches a consumer. It
  lands as its own commit with the parity oracle green before and after, so a bisect can find it
  if a fitted level moves.
- **Incremental, with pauses.** One milestone, then stop.
- **Container builds.** `cargo` through `./scripts/dev.sh` (CLAUDE.md); this machine has Apple
  `container`, so `target-container/` is where the binaries land.

---

## Preconditions (already in place)

- **`StratumEvidence` carries the evidence the weight reads** — `spanning_reads()`, and from
  2026-08-20 `bases_compared`, `mismatching_bases`, `substitution_rate()` and
  `reads_off_reference_length()`
  ([`ssr_fit.rs:169`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L169)).
- **The walk writes a per-cell table** under `SSR_CELL_TABLE`, and a second one under
  `SSR_CELL_TABLE_BORROWED`
  ([`ng_joint_records_walk.rs`](../../../../examples/ng_joint_records_walk.rs)).
- **Both cohorts' independent-fit tables exist** and are the parity oracle:
  `tmp/slippage_curve/tomato_a1_noborrow.csv` (6 fitted cells) and
  `hg002_a1_noborrow.csv` (55 fitted cells). **A5 copies them under `tests/data/` so the oracle
  is not scratch.**
- **The drawn-stratum generator** exists for tests that need evidence with a known answer
  ([`bench_fixtures`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs), `cfg(test)` and
  the `bench-fixtures` feature).

---

## The steps

### Milestone A — the curve: its type, its fit, and what consumes it

**A1. The types, no logic.**  ✅
`SlippageCurve` (`rise_shape`, `intercept`, `slope`, the fitted repeat-count range, the held-out
error, the contributing-cell count) and `LevelSource` (`Cell` / `Curve` / `Blend { curve_weight }`).
Add a `SlippageCurveConfig` carrying the shape-number grid, the minimum contributing cells, and
which weight is used — all named `pub const`s with their source in the doc comment.
***Depends:*** none. ***Source:*** spec §2, §4, §7.

**A2. Evaluate a curve, and hold flat outside its range.**  ✅
`SlippageCurve::level_at(repeat_count) -> f64`: invert `level ^ rise_shape = intercept + slope · n`
at `rise_shape > 0` and `exp(intercept + slope · n)` at `rise_shape = 0`, clamp the result into
`(0, 1)`, and hold at the nearest fitted end outside the fitted range. Unit tests: a curve
evaluated at its own fitted cells returns them; below and above the range it returns the end
values; `rise_shape = 0` and `rise_shape = 1` each reproduce a hand-computed value.
***Depends:*** A1. ***Source:*** spec §2, §6.

**A3. Fit one slippage group's line at a fixed shape number.**  ✅
Weighted least squares of `level ^ rise_shape` on repeat count, with the weight the cell's slipped
reads (`level × spanning_reads`). Reject a fit whose `slope` is not positive, returning no curve
rather than a falling level. Unit tests: cells drawn exactly on a line at `rise_shape = 1` and
exactly on an exponential at `rise_shape = 0` are both recovered to within rounding; a
deliberately falling set returns no curve.
***Depends:*** A2. ***Source:*** spec §4.2, §9.

**A4. Choose the shape number for a period, across slippage groups.**  ✅
For each rung of the grid, fit every group's line (A3), then score the rung by leaving each
contributing cell out in turn, refitting that group's line without it, and predicting it; keep
the rung with the lowest median relative error, ties to the larger rung. Record the winning
rung's held-out error on the curve. Refuse to draw a curve for a period with fewer than the
minimum contributing cells.
***Depends:*** A3. ***Source:*** spec §4.1, §4.3, §4.4.

**A5. The oracle fixtures, and the end-to-end test of A4 on them.**  ✅
Copy the two cell tables under `tests/data/slippage_cells/` and add a test that reads them, runs
A4, and asserts the fitted shape number per period against the report's table — 0.00 at tomato's
period 1, 1.00 and 0.80 at HG002's periods 1 and 2 — and that HG002's period 1 curve predicts a
held-out cell to within 5%.
***Depends:*** A4. ***Source:*** [the report](../reports/str_slippage_shape_2026-08-20.md) §4.1.

> **Checkpoint A: the curve is fitted and proven against both cohorts' real cells, and nothing in
> the pipeline calls it yet. Pause for review.**

### Milestone B — wiring it into the fit, and what comes out

**B1. The outcome type gains a cell that has a level and no shares.**  ✅
Extend `StratumOutcome`
([`ssr_fit.rs:322`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L322)) so a cell
below the refusal floor can carry a curve-derived level, and add the provenance the spec §8
requires: the cell's tracts, reads crossing and slipped reads; the `LevelSource`; whether the cell
sat inside the curve's fitted range; and the curve's held-out error and cell count. Types and
their construction only — no caller changes yet.
***Depends:*** A1. ***Source:*** spec §8.

**B2. Blend the cell's own level with the curve's.**  ✅
Inverse-variance on the log scale: the cell's relative standard error is `1 / sqrt(slipped reads)`
and the curve's is its held-out error, with the curve's weight divided by `(gap / 2.5)²` where the
gap in combined errors exceeds 2.5. A cell with no fit takes the curve whole; a period with no
curve keeps the cell's own level; a cell whose fitted level is zero takes the curve whole and
never enters the logarithm. Unit tests pin the weights against the spec's two worked figures — at
a curve error of 4.4% the curve carries 93% of the weight at 40 slipped reads and 6% at 8,000 —
the knee standing the curve down at a gap of 9.3 combined errors, and the three degenerate cases.
***Depends:*** B1, A2. ***Source:*** spec §7, §7.1, §7.2.

**B3. `fit_strata` draws the curves and emits blended levels.**  ✅
After every cell is fitted on its own tracts, group the fitted cells by period, run A4, and emit
each cell's level through B2. **Own commit, do not bundle** — this is the step that changes the
number a consumer reads, and its failure is a wrong genotype rather than a panic. The oracle:
with the curve disabled by config, both cohorts' cell tables must come back identical to
`tests/data/slippage_cells/`.
***Depends:*** B2, A5. ***Source:*** spec §5.

**B4. Narrow borrowing to the two shares.**  ✅
`borrow_up_to_the_floor`
([`ssr_fit.rs:711`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L711)) stops
supplying the level and keeps supplying the direction split and the fall-off, as
[`../spec/parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) §4.5 already specifies.
**Leave the pooled-set dedup in place** — with borrowing off every set is a singleton and it looks
like dead code; it is what makes the borrowing arm affordable and that arm stays selectable.
***Depends:*** B3. ***Source:*** spec §5, §9.

> **Checkpoint B: the fit emits curve-derived levels, the independent-fit oracle is unmoved, and
> borrowing serves only the shares. Pause for review.**

### Milestone C — measure what it bought, and say so in the specifications

**C1. Re-run both cohorts and report the movement per cell.**  ☐
One walk each, writing the cell table with the curve on and with it off. Report how far each
cell's level moved, which cells gained a level they did not have, and the curve's held-out error
per period. Expect the tomato run to take about 70 minutes and HG002 about 60.
***Depends:*** B4. ***Source:*** spec §9.

**C2. Amend the specifications.**  ☐
[`../spec/parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) §4.3 (borrowing and merging
— the level leaves it), §4.4 (what the summary over strata must now say), §4.5 (the split floor,
now split differently), and the architecture sibling
[`../arch/parameter_prepass_ssr.md`](../arch/parameter_prepass_ssr.md). Point at the new spec
rather than restating it.
***Depends:*** C1. ***Source:*** spec §1.2, §5, §8.

**C3. Update `PROJECT_STATUS.md`** — the in-scope feature's block only, per the *Project
status protocol* in [`ai/skills/apply-code-review-fixes/SKILL.md`](../../../../ai/skills/apply-code-review-fixes/SKILL.md).  ☐
***Depends:*** C2.

> **Checkpoint C: measured, documented, and the status file says where this stands. Pause for
> review.**

### Milestone D — the shares' floor, and retiring the pooled refit

**Added 2026-08-20 (owner), after C1 showed the level's curve reaches no thin stratum until the
shares do.** Design: spec §5.1. This is what makes spec §1.1's first goal — every populated
stratum gets a level — actually arrive, and it removes the run's expensive arm on the way.

**D1. The floor, the source, and the outcome a stratum gets when nothing was fitted.**  ✅
`MIN_SLIPPED_READS_TO_FIT_SHARES = 4_000` with its derivation in the doc comment; a `SharesSource`
recording whether the two shares are the stratum's own or the repeat count they were copied from;
and a `StratumOutcome::Derived` for a stratum whose numbers were all supplied from elsewhere —
carrying the slippage numbers and the provenance, and **no length spectrum, concentration or
log-likelihood**, because nothing was fitted and a fitted-looking zero would be a lie. Types only.
***Depends:*** B1. ***Source:*** spec §5.1, §8.

**D2. Copy the shares, and stop pooling.**  ✅
A stratum whose own slipped reads reach the floor keeps its own shares; below it, the two are
copied from the nearest stratum at the same period that clears it, nearer repeat count first and
the shorter tract winning a tie. `borrow_up_to_the_floor` and the pooled refit go. **Own commit,
do not bundle** — this changes what the shares are at every thin stratum, and its failure is a
wrong genotype rather than a panic.
***Depends:*** D1. ***Source:*** spec §5.1.

**D3. Emit a stratum that was never fitted.**  ✅
A stratum with spanning reads but no fit of its own comes back `Derived` with a curve level and
copied shares. Only two refusals survive: no read crossed it, and its period has neither a curve
nor a stratum clearing the shares' floor.
***Depends:*** D2. ***Source:*** spec §5.1, §1.1.

**D4. Measure it on both cohorts.**  ✅
How many strata gained a complete parameter set, what the shares were copied across, and what the
run now costs against the pooled arm's 1,036.8 s against 155.5 s.
***Depends:*** D3. ***Source:*** spec §5.1, §9.

> **Checkpoint D: every populated stratum carries a level and a pair of shares, each saying where
> it came from. Pause for review.**

### Milestone E — one treatment for all three fitted numbers

**Added 2026-08-20 (owner), after D4 measured what the floor-and-copy rule delivers: 13 furnished
strata on HG002 and none at all on tomato, because only one motif period of six has a stratum
clearing 4,000 slipped reads.** Design: spec §5.1. This replaces that rule with the level's own
machinery applied to the two shares — a curve per period fitted from *every* stratum weighted by
its precision, and each stratum departing from it by inverse variance.

**E1. Choose each share's family, by measurement.**  ✅
Compare a constant, logit-linear and logit-quadratic in repeat count, per period and per
parameter, on the held-out-cell criterion, over both cohorts' ±8 cell tables. **Produces the
answer to the research plan's C4** — a period whose held-out error is lowest at the constant has
no trend to fit, and says so. **No code ships from this step**; its output is the report that fixes
E2's family list.
***Depends:*** D4's tables. ***Source:*** spec §5.1.
***Output:*** [`../reports/str_slippage_share_families_2026-08-20.md`](../reports/str_slippage_share_families_2026-08-20.md)
— all three families are needed and are chosen at run time; two questions go back to the owner
before E2 (which precision weights a stratum, and whether a curve that bends twice may be drawn
through four strata).

**E2. A share's own precision, and its curve.**  ✅
`sqrt((1 − p) / (p · S))` on `S` slipped reads — the model the 4,000 came from, which reproduces
the architecture's 1,400 and 4,000 to within 3%. The curve fit itself is the level's, generalised
over the families E1 chose, on the logit scale.
***Depends:*** E1. ***Source:*** spec §5.1.
***Two recorded departures, both ruled by the owner 2026-08-20.*** The weight is the inverse
variance of the **logit**, `S · p · (1 − p)`, not of the share: the two disagree by `1/(1 − p)²`
and the fit runs on the logit ([report](../reports/str_slippage_share_families_2026-08-20.md) §7.1
and §9). And **a curve always comes back** — a period too thin to choose a shape gets a flat mean,
one with nothing gets the run's other periods, and a run that fitted nothing anywhere gets a
built-in default, each recorded in the provenance. These numbers are a prior, so answering coarsely
beats refusing.

**E3. Blend each share, and retire the floor-and-copy rule.**  ☐
The blend is §7's with the logarithm replaced by the logit. `SharesSource` becomes the same
three-way provenance the level carries. `strata_lending_their_shares`, `copy_shares_from_a_neighbour`
and `nearest_lender` go. **Own commit, do not bundle** — it changes both shares at every stratum.
***Depends:*** E2. ***Source:*** spec §5.1, §7.

**E4. Lower the refusal floor so thin strata contribute.**  ☐
Fit strata far below 50 tracts so they feed their period's curves, weighted by their own
precision. **Gated on spec §11's open question** — fit drawn strata down to a handful of tracts
first and report how often the climb converges and how often the level returns exactly zero. If
that measurement says the thin fits are unusable, this step does not happen and the floor stays.
***Depends:*** E3. ***Source:*** spec §5.1, §11.

**E5. Measure it on both cohorts, and retire `Derived` if nothing uses it.**  ☐
How many strata carry each of the three numbers from their own fit, from a curve, and from a
blend; what the run costs; and whether any stratum is still furnished-from-nothing.
***Depends:*** E4. ***Source:*** spec §5.1, §8, §9.

> **Checkpoint E: all three fitted numbers smoothed the same way, with one mechanism rather than
> two. Pause for review.**

---

## Verification summary

| milestone | proven by |
|---|---|
| A | unit tests on drawn cells with a known curve (A2, A3), then the fitted shape number per period reproduced against both cohorts' real cell tables (A5) |
| B | the independent-fit parity oracle — every cell's own fitted level byte-identical with the curve disabled (B3) — plus unit tests pinning the blend crossover and its three degenerate cases (B2) |
| C | a re-run of both cohorts reporting per-cell movement and per-period held-out error (C1) |

**The parity oracle is the point of B and must not be weakened.** Nothing in this plan touches
stage one, so a cell's *own* fitted level moving is a defect in the plumbing, not a consequence of
the design.

---

## Out of scope (next plans)

- **The other three numbers** — the research plan's C4, which has measurements pointing both ways
  already ([the report](../reports/str_slippage_shape_2026-08-20.md) §5).
- **One-stage fitting** — the research plan's D1.
- **Whether the flattening is the polymerase or the census's ±4 recording window** — the spec §11's
  first open question. It changes what the fitted numbers *mean*, not what the code does, and the
  measurement is one constant and one re-walk.
- **Genotype concordance** — the research plan's F2.
