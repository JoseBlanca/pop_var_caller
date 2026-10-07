# Code Review: fit_precision_a7
**Date:** 2026-09-28
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step A7 of `fit_precision.md` — the allele-frequency shapes follow the likelihood's own slope
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** the diff `41e0f9bf..8be6683f` (`8be6683f`, "WIP: … for review", a
  review-only commit on no branch).
- **In-scope files:** [fit.rs](../../../../src/parameter_estimation/joint/fit.rs) (`step_beta_shapes`,
  which replaces `fit_beta_shapes`; the slope sums in `Statistics` and `one_position`;
  `PassModel`'s rules and slopes; `BetaQuadrature::prior_slopes`; the monotonicity test),
  [fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs) (`RuleSlopes`
  shared with the pass, the coverage test's re-recorded checks, a trace test's seed),
  [cli/cross_platform_digests.rs](../../../../src/cli/cross_platform_digests.rs) (both checksums
  re-recorded), the A7 implementation report.
- **Categories dispatched**, three reviewers, each in its own worktree at the review commit:
  - *design and claims* — naming, idiomatic, module structure, refactor safety, smells, defaults;
    conformance to plan step A7; every figure in the prose (11 mutations, a probe at a bound);
  - *reliability and errors* — tests and mutation (32 mutations, two review tests);
  - *correctness, float portability, cost* — an independent derivation of the slopes, bit
    identity at 1 and 4 threads, pass counts against the parent, a probe at a bound (8 mutations).

## 2. Verdict

**Request-changes**, for one defect in the update and for slopes nothing tested.

**The slopes are right.** The correctness reviewer derived the segregating branch's slope in a shape
independently — each node's posterior times the node weight's slope, plus how the node's product
over samples moves along the frequency axis times the node's own slope — and the carrier branch's,
coverage odds included; at 12 fitted cohorts the pass's summed slopes equal the scorer's to 6 or
more significant figures. **Away from the bounds the fit now stops at the likelihood's flat point**:
at 20 samples the full-matrix Newton step still left in each shape and in the invariant share is at
most 0.0024 errors. **The same bits at 1 and 4 threads**, with and without the duplicated class, on a
12-chunk cohort. Only `+ − × ÷` outside `crate::float`.

**The defect is at a bound**, where the oracle cohort and the checksum fixture both sit (the density's
second shape at its upper bound, 50): the step clamped each shape of a joint step whose cross term
assumes both move, so the free shape came to rest where its own slope was not zero — up to 34 errors
from its maximum with the other held, and up to 372 log-likelihood units below it on drawn cohorts.

## 3. Execution status

- Each reviewer: `cargo test --release --lib parameter_estimation::joint::fit` — 69 passed, 1 ignored
  (plus each reviewer's own ignored probes).
- Mutations: design 11 run (11 killed, 3 only by a test that asserts which start wins); reliability
  32 run (30 changed behaviour, all failed some test, 6 only the checksum and that same winner test; 2
  changed nothing); correctness 8 run (1 survived, latent: no production caller sets coverage odds).
- Pass counts against the parent `41e0f9bf`, 12 cohorts × 3 starts, with and without acceleration
  (correctness reviewer).
- Not re-run by reviewers: the oracle (its data was not mounted), the 35-minute coverage test, the
  full suite. The per-pass timing was checked against the author's log, not re-timed.

## 4. Open questions and assumptions

- **Whether "never lowering the log-likelihood" must be enforced or may be measured** is the
  owner's to rule on (design reviewer): the plan states it as a property of the update; the author
  reported it measured. The reviewer recommends ruling only after the bound fix, with a bound cohort
  in the measurement. Carried to checkpoint A′.
- **The slower convergence** (M2 below) matters for kimura's wall time and for the stopping rule of
  plan step B; carried to checkpoint A′.

## 5. Top 3 priorities

1. **B1** — the step at a bound.
2. **B2** — the pass's slope sums pinned by a test.
3. **M2** — the pass counts, measured and stated.

## 6. Findings

### Blocker

**B1: [fit.rs](../../../../src/parameter_estimation/joint/fit.rs) `step_beta_shapes` — at a bound, the
free shape comes to rest where its slope is not zero, and the log-likelihood can fall.**
**Categories:** design, correctness (convergent, each measured independently). The step solves the
two shapes' Newton step jointly, then clamps each. When one shape is held at its bound, the other
still takes the joint step, whose cross term (`h12 = −n·ψ′(a+b)`, negative) assumes the held shape
moved; it rests where `slope_a = slope_b · h12 / h22`, strictly negative when the held shape's slope
pushes out. Measured:

| cohort (truth) | clamped: first shape, its slope in errors | the projected step: first shape | log-likelihood gained |
|---|---|---|---|
| 20 samples × 3 reads (15, 200), 3 seeds | 6.87–6.95, −29.0 to −29.5 | 3.44–3.47 | 341 to 355 |
| 20 × 8 reads (2, 120), 3 seeds | 3.71–3.78, −31.6 to −33.8 | 0.46–0.47 | 331 to 372 |
| 4 × 3 reads (15, 200), 2 seeds at the bound | 6.2–6.3, −5.0 to −5.3 | 2.4–2.5 | 28 to 31 |
| 20 × 8 reads, duplicated class, carrier's second shape at 50 | carrier's first 28.2, 33.3 | 12.1, 13.4 | 9.0, 16.2 |

On the last cohort the error spread to the density: the step still left in the invariant share was
3.0 errors. The clamped fit ended **156 units below its parent** on one bound cohort, and plain passes
of it **lowered the log-likelihood by 36 units** between passes 100 and 200 on a 4-sample cohort
drawn with b = 150 (design) and by up to 0.94 units on the duplicated cohort (correctness). The flaw
predates A7 — the digamma update clamped inside the same joint step — but A7 made it worse there.
**Fix (both reviewers, the same):** project the step — hold the shape whose step leaves its interval
at the bound, and give the other the step the same curvature gives with that one held.

**B2: [fit.rs](../../../../src/parameter_estimation/joint/fit.rs) `one_position` — the slope sums the
update steps on are compared with nothing.** **Categories:** reliability (Blocker), design and
correctness (Major; convergent). The pass keeps its own copy of the slope formula beside the scorer's,
which is finite-differenced. Six wrong accumulations — the along-frequency term dropped, its sign
flipped, the other shape's node slope, for each branch — fail only the fixture checksum and
`the_trace_files_the_returned_fit_after_the_winning_starts_last_pass`, whose assertion that the middle
start wins moves with any change to the fit; the checksum is re-recorded whenever fitted numbers move,
as this step did. The carrier's coverage odds are exercised by no fixture. **Fix:** a test that the
pass's summed slopes equal the scorer's over a multi-chunk cohort, with and without the duplicated
class, and one with coverage odds set.

### Major

**M1: the monotonicity claim — "no plain pass lowers the log-likelihood" — rests only on cohorts with
interior shapes.** **Categories:** correctness, design. True on the author's eight cohorts (0 of 1,920
passes); false at a bound (B1). **Fix:** fix B1, re-measure with bound cohorts, put one in the test.

**M2: [fit.rs](../../../../src/parameter_estimation/joint/fit.rs) `step_beta_shapes`' doc — the update
is not "as fast near the answer as the expectation-maximisation update it replaces".**
**Category:** correctness (cost). Passes to convergence, parent against A7, three starts:

| cohort | accelerated | plain alternation |
|---|---|---|
| 20 samples × 3 reads, 3 seeds | 12–18 → 18–30 | 18–27 → 21–84 |
| 4 × 3 reads, seed 1 | 39–55 → 63–70 | 171–174 → 315–324 |
| 20 × 8 reads, duplicated, seed 0 | 33–39 → 39–57 | 120–123 → 135–138 |
| 4 × 3 reads, duplicated, seed 1 | converged at 165–178 → the 198-pass limit | 678–687 → the 1,500 limit |

With 4 to 7% more a pass, a converged 20-sample fit costs about 1.4 to 2.1 times what it did.
**Fix:** state the measured ratio; raise it at the checkpoint (kimura's 9 h 25 min, plan step B).

**M3: [fit.rs](../../../../src/parameter_estimation/joint/fit.rs)
`plain_passes_never_lower_the_log_likelihood` passes nine wrong updates.** **Category:** reliability. A
frozen carrier, a half-size step and wrong slopes lower nothing either; the one-sample cohort cannot
fail even under the digamma update (0 of 240 passes fell). **Fix:** say on the test what it cannot
see; replace the one-sample cohort with one that can fail.

**M4: [fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs) — the only
behavioural check that the density stops at its flat point is the ignored 35-minute coverage test, and
it never fits the duplicated class.** **Category:** reliability. **Fix:** B2's tests as the fast
guard; say in the coverage test's doc that the carrier is not covered there.

**M5: [cli/cross_platform_digests.rs](../../../../src/cli/cross_platform_digests.rs) — the re-recording
note says both shapes' "errors come out wider than their range"; the fit reports both "not
identified".** **Categories:** design, reliability (convergent). **Fix:** say what the fit reports.

### Minor

- **m1** (design, correctness) — the report explains the oracle's fit but not how its calls moved
  (records 6,693 → 6,706, 32 shared records changed genotype, all four inbreeding coefficients, the
  ordinary-site prior's concentration, the duplicated share tripled), nor that every start of the
  oracle's fit stops at the pass limit, before and after.
- **m2** (design, reliability) — `step_beta_shapes`' doc says the old update stopped "0.4 to 3.1
  standard errors beside the maximum"; that range is the distance from the truth. The distance from
  the flat point was 0.65 to 2.3.
- **m3** (design) — the module doc still says each iteration cannot lower the likelihood.
- **m4** (design) — spec §3.2's score table still names `fit_beta_shapes` and the digamma terms. The
  spec is the owner's; offer a one-line amendment.
- **m5** (design, reliability) — `PassModel` holds the carrier's rule and its slopes as two `Option`s
  that must agree; a mismatch would freeze the carrier's shapes silently. The carrier's slopes are
  sized by the density rule's node count.
- **m6** (design; reliability, cross-category) — the slope formula is written twice, in the pass and in
  the scorer, with different thresholds for skipping a node (`share ≤ 10⁻¹²` against `= 0`).
- **m7** (design) — the parameter `weight` does not say it is the posterior count of positions.
- **m8** (reliability) — a non-finite slope passes the determinant guard and becomes NaN shapes.
- **m9** (correctness) — "at rest only where the slope is zero" is stated without its bound exception.
- **m10** (design, cross-category) — `step_beta_shapes` has no limit on its step's length; with few
  segregating positions one step might reach a bound.

### Nit

- **n1** (correctness) — the determinant guard returns the shapes unmoved without a word; a
  `debug_assert!` would say it is unreachable. (The reliability reviewer shows it is reachable at a
  vanishing count of positions.)
- **n2** (design) — `RuleSlopes`' fields became `pub(super)` for the pass; an accessor would keep its
  "records the shapes it was computed at" check meaningful.
- **n3** (design) — the report's test table says the step's bound case clamps both shapes, which
  cannot see B1.
- **n4** (reliability) — the step's unit test builds its expectation from the same trigamma formula
  as the code; one value worked by hand would anchor it.
- **n5** (reliability) — the report's "of 24 seeds, 10 give start 1 …" has no log behind it.

## 7. Numbers checked

Every figure in the report's §2–§3 and in the coverage test's doc was checked against the author's
logs by the design and reliability reviewers and found correct, except the three already listed (M2's
"as fast", M5's "wider than their range", m2's "beside the maximum"). The oracle's checksums, its
29-unit gain and the density's a = 15.6 ± 36.9 are correct as figures but describe a fit at the bound
(B1).
