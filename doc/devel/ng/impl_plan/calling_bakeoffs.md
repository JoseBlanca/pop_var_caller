# ng calling bake-offs — implementation plan

**Status:** draft, 2026-08-21. The build order for **the measured alternatives to the default
loop**: spec Q1's arms B and C with the exhaustive scorer and the local search, spec Q2's
slippage re-fit at its three pull-back settings, and spec Q3's discovery at its three modes —
plus the harness that reports each comparison. Design is settled in
[`calling_em_loop.md`](../spec/calling_em_loop.md) (Q1–Q3, §4.1, §5.1) and
[`../arch/calling_em_loop.md`](../arch/calling_em_loop.md) §3, §6. This plan turns that design
into build order; it is **not** a place for new design — every arm, setting and reported number
below is the spec's.

**Where this sits.** Six plans build calling:
[`calling_prerequisites`](calling_prerequisites.md) ∥
[`calling_foundations`](calling_foundations.md) → [`calling_prior`](calling_prior.md) ∥
[`calling_read_likelihoods`](calling_read_likelihoods.md) → [`calling_loop`](calling_loop.md) →
`calling_bakeoffs`. It follows the loop plan, which ended with ng calling genotypes; everything
here is measurement and the machinery measurement needs.

**The settled ground rule (spec Q1, 2026-08-21): every Q1 arm runs at `F = 0`.** The joint
priors are written over allele counts and have no per-sample slot; the exact composition is a
`2^k`-term mixture at `k` homozygotes, so there is no cheap joint counterpart to arm A's
mixture. Consequently `JointAssignmentPrior` implementations **take no inbreeding argument**,
and the arm-B/C genotyper **rejects a non-zero `InbreedingF` loudly** rather than ignoring it —
a silently dropped `F` on a selfing panel is the failure this guards. Whether a whole-cohort
scorer can carry inbreeding at all is the spec's recorded follow-up, not this plan's.

---

## Scope

**In:** `calling/inference/assignment.rs` — `JointAssignmentPrior`,
`CohortDirichletMultinomialPrior` (arm B), `EwensArrangementPrior` (arm C),
`AssignmentGenotyper` with `Exhaustive` and `LocalSearch` enumeration;
`calling/inference/discovery.rs` and the slippage re-fit body in `summarise_condition.rs`
(lifting the loop plan's reject-non-default gates); the `bench/` harness runs Q1–Q3 specify.

**Out (with owners):**

- **Real-cohort end-to-end runs remain blocked on candidate selection (step 6, no spec)** — the
  same blocker [`calling_loop.md`](calling_loop.md) records. What runs today: everything
  fixture- and simulator-based, and any bench whose candidates can be supplied. The tomato-panel
  and GIAB/HG002 halves of Q1–Q3's reports run only once selection exists; each is flagged below
  where it appears.
- **The sibling specs' change measurements** (read_likelihoods §12 items 14–19; the priors' GIAB
  regression under both seam impls) — same blocker, same instrument; they are named here so they
  are not dropped, and they need an owner and a plan once step 6 exists.
- **A normalised Ewens variant** — a fourth arm, not a bug fix (spec Q1); built only if someone
  schedules it.
- **A six-number re-fit arm** — askable only once the part-repeat estimator has an owner
  ([`read_likelihoods.md`](../spec/read_likelihoods.md) §10; spec §5.1).

## Principles (how the order was chosen)

- **The spec's own build order for Q1** (spec Q1, "build it in four pieces"): priors as pure
  functions → exhaustive scorer → search → swap the second prior in. The exhaustive scorer comes
  before the search because it is the piece that separates "the model disagrees" from "the
  search is too narrow".
- **One code path per axis.** Frozen / pulled-back / free are three values of
  `SlippageRefitConfig`; off / frozen-frequencies / full-convergence are three values of
  `DiscoveryConfig`. No setting is a second implementation (spec §5.1, §4.1).
- **Verify against ground truth.** The exhaustive scorer is the search's oracle (spec §13 test
  9); Ewens' non-normalisation is pinned at its closed-form total, not "fixed" (test 10);
  discovery and the re-fit are pinned by plant/terminate/append tests against the loop plan's
  instrumented counters.
- **Build natively, not by porting freebayes' layout.** freebayes copies the whole assignment
  per neighbour — quadratic in cohort size on its exhaustive-local path; store **moves**, and if
  the measurement reports a quadratic curve, first check whether ng copied the copy (spec Q1's
  prediction).
- **Isolate the silent steps.** A re-fit that drifts on an empty profile, a rebuilt-instead-of-
  appended column, and the emission-call count at the re-fit setting all fail as quietly-wrong
  or quietly-slow; own commits, marked below.
- **Container builds.** All `cargo` via `./scripts/dev.sh`; bench runs on the machines the
  benchmarks already use.

## Preconditions (already in place)

- **The loop plan merged:** arm A end-to-end on fixtures, the instrumented emission-call and
  allocation counters, the config gates that reject non-default `SlippageRefitConfig` /
  `DiscoveryConfig` (this plan lifts them), and the two loop oracles.
- The ported Dirichlet-multinomial primitive ([`calling_prior.md`](calling_prior.md) B1) — arm
  B's prior is that primitive on cohort counts; one primitive, two callers.
- freebayes' Ewens implementation as the 51-line *reference for behaviour*, reimplemented
  natively ([`freebayes/src/Ewens.cpp`](../../../../freebayes/src/Ewens.cpp)); the discovery
  mechanism built **from spec §4.1's description, never from HipSTR's GPL tree** (the licence
  rule, [`read_likelihoods.md`](../spec/read_likelihoods.md) §4.2).
- Production's re-fit pieces to port: `refine_theta_locus`
  ([`ssr/cohort/stutter.rs:184`](../../../../src/ssr/cohort/stutter.rs)) and the attribution
  whose input changes to posteriors (`attribute_locus`,
  [`em.rs:1192`](../../../../src/ssr/cohort/em.rs)); the pull-back target is the fitted curve's
  value at the cell — `blend_level` via the `StratumFits` gather
  ([`calling_prerequisites.md`](calling_prerequisites.md) F), not a cell's own estimate (spec
  Q2 preamble).
- The STR simulator that produced the Model-A table
  ([`sim.rs`](../../../../src/ssr/cohort/sim.rs)) — the harness's exact-truth rung.

## Branch and merge (sequential — no worktree)

- **Branch** `ng-calling-bakeoffs`, from `main` **after `ng-calling-loop` has merged**. Nothing
  runs beside it, so no worktree — a plain branch in the primary checkout. (If two people pick
  it up, the natural split is A–C against D–E: the Q1 machinery and the STR mechanisms share no
  file until F.)
- Conflict surface: none — this branch is the sole calling editor while it lives.

---

## The steps

### Milestone A — the two joint priors (pure functions)

**A1. `JointAssignmentPrior` + arm B's prior.**  ☐
The trait — log-prior of one whole-cohort assignment from its per-allele chromosome counts plus
the per-sample genotypes (the arrangement term needs the multiplicities), **no inbreeding
argument** — and `CohortDirichletMultinomialPrior`: the ported primitive evaluated on the
cohort's own counts, nothing conditioned, no leave-one-out. Test: at one sample, arm B's prior
over that sample's genotypes is a hand-checkable Dirichlet-multinomial row. *Source:* spec Q1
(arm B); arch §3.2.

**A2. `EwensArrangementPrior` — arm C, native.**  ☐
freebayes' arrangement term × Ewens' sampling formula on the assignment's partition pattern,
reimplemented natively. Test (spec §13 test 10): summed over every assignment at a two-allele,
one-diploid-sample locus the exponentiated prior gives **`(2 + θ)/(1 + θ)` — about 2, not 1** —
asserted at the closed form, because a test asserting it sums to one fails against a correct
implementation. **Reported, never silently normalised.** *Depends:* A1. *Source:* spec Q1, §2.1;
arch §3.2.

> **Checkpoint A:** both priors are pure, pinned functions — including the pinned
> non-normalisation. Pause for review.

### Milestone B — the exhaustive scorer

**B1. `AssignmentGenotyper` + `Exhaustive`.**  ☐
Enumerate every assignment (bounded by `max_assignments`), sum joint posteriors, marginalise
each sample's genotype. Rejects non-zero `InbreedingF` loudly (the ground rule above). Tests: on
two or three samples with hand-enumerable assignments, the marginals match a by-hand sum; arm B
exhaustive vs arm A on a small cohort — **differences recorded, not asserted away** (the two
differ by design in shape; what the test pins is that the harness attributes them). *Depends:*
A1, A2. *Source:* spec Q1 piece 2; arch §3.1.

### Milestone C — the local search

**C1. `LocalSearch` — moves, not copies.**  ☐
Start at the likelihood-preferred assignment, climb by bounded moves, sum the neighbourhood —
storing **moves** rather than assignment copies, so retained memory stays linear in cohort size.
Tests: **exhaustive-vs-search agreement on enumerable cohorts** (spec §13 test 9 — a difference
is the search's narrowness and is reported as such, never as a model difference; and **no**
"exhaustive agrees with arm A at one sample" assertion, which is false by the priors differing);
a retained-memory check that grows linearly, not quadratically, with samples. *Depends:* B1.
*Source:* spec Q1 pieces 3–4; arch §3.1.

> **Checkpoint B/C:** arms B and C run under both enumerations; the oracle and the search agree
> where enumeration is possible. Pause for review.

### Milestone D — the slippage re-fit (Q2's one code path)

**D1. Posterior-weighted attribution + the re-fit.**  ☐
Port the attribution with its input changed — every read weighted by the **genotype
posteriors**, never the called genotype (HipSTR's choice; production's deferred refinement) —
and the shape re-fit (`refine_theta_locus`'s port) moving **three numbers**, part-repeat
placeholders held fixed. The level's pull-back target is **the curve's `blend_level` at the
cell** (spec Q2 preamble — the curve redefinition postdates the spec's Q2 wording, and "far from
its stratum" is now measured against the curve). *Source:* spec §5.1; arch §6.1.

**D2. The nested rounds, live.**  ☐
Lift the loop plan's gate: an outer round re-fits, rebuilds the `Lg` table, reruns the frequency
loop; rounds stop on the re-fitted numbers' movement or the cap, breaking **before** the rebuild
on convergence (production's order). Frozen / pulled-back (50, 20) / free (0, 0) are three
configurations of this one path. Tests: **the no-op collapse** — at a locus with no slips, the
re-fitted numbers return the per-stratum values and genotypes are bit-for-bit the frozen
setting's (spec §13 test 6); **the emission-call count at both settings** — `builds` equals the
instrumented rebuild count, one less than rounds run when the last round leaves before
rebuilding (test 5's re-fit half; run at both settings, because a per-pass rebuild hides behind
a per-round one when only the re-fit arm is tested). **Own commit, do not bundle** — a drifting
no-op re-fit and a double-paid table are both silent; the two instrumented tests are the oracle.
*Depends:* D1. *Source:* spec §5.1, §6, §13 tests 5–6; arch §6.1.

> **Checkpoint D:** three pull-back settings, one code path, collapse and cost both pinned.
> Pause for review.

### Milestone E — discovery (Q3's one mechanism)

> **2026-09-02 — ownership note.** The STR loop's spec now carries discovery's build and its
> measurement as part of the STR caller's own scope
> ([`calling_loop_ssr.md`](../spec/calling_loop_ssr.md) §3.5, owner's ruling), sequenced
> after tract selection and the driver branch. This milestone's design and F3's report shape
> are what that work builds to — one implementation, whichever plan schedules it first;
> coordinate before starting either.

**E1. `discovery.rs` — the mechanism.**  ☐
From spec §4.1's description (never the GPL tree): after convergence, retrace what the model
explains as slippage; admit a tract length clearing **both halves** of the bar (`min_reads = 2`,
`min_spanning_read_share = 0.15`, inherited and soft); **append** the new alleles' emission
columns, never rebuild; re-run per `DiscoveryMode` (`AgainstFrozenFrequencies` — scoring passes
against held frequencies, converge once at the end; `AgainstFullConvergence` — the outermost
repeat); stop when a round adds nothing or at the allele cap; then prune alleles no sample's
best genotype used and rerun the frequency loop. Extent never grows — alleles within the locus
only. *Source:* spec §4, §4.1; arch §6.2.

**E2. The three discovery pins.**  ☐
Tests: **plant/terminate/free-when-off** (spec §13 test 11) — a planted hidden allele is found
and the sample ends heterozygous where `Off` calls it homozygous; a locus supporting nothing
ends after one round and a runaway stops at the cap (assert the round count); with `Off`, the
emission count and genotypes are **bit-for-bit** the plain loop's. **Append-only columns**
(test 12): after a round admits `k` alleles the emission count has risen by exactly those
alleles' entries and nothing else. **Own commit, do not bundle** — a rebuilt table is quietly
slow and a leaky `Off` quietly changes the default; the instrumented counter is the oracle.
*Depends:* E1. *Source:* spec §13 tests 11–12; arch §6.2.

> **Checkpoint E:** discovery finds, terminates, appends, and costs nothing when off. Pause for
> review.

### Milestone F — the measurement harness and the runs

**F1. Q1's report.**  ☐
Per arm: peak memory and wall-clock per locus, at one sample, at 63 samples, and at a thousand
(synthetic cohorts — real-cohort candidates are blocked on step 6, and the report says so on
every affected row); **the count of genotypes that differ from arm A, reported as B−A and C−B
separately** — the shape difference and the prior difference; a single A−C number is the one
figure that cannot be acted on. Include the quadratic-curve check against C1's prediction.
*Depends:* B1, C1. *Source:* spec Q1 ("what to report, per arm").

**F2. Q2's report.**  ☐
The three pull-back settings on the STR simulator (exact truth, slippage settable) now; on the
HG002 tandem-repeat bundle and tomato's recurrence standard **once step 6 unblocks them** —
three numbers per setting: genotype accuracy given detection, the share of loci whose re-fitted
numbers land far from the curve's value, wall-clock (each round rebuilds the table). *Depends:*
D2. *Source:* spec Q2; §13's bench shape.

**F3. Q3's report.**  ☐
The three modes: **how often it fires at all, first** — if a handful of loci in ten thousand,
everything after is a rounding error; then, on firing loci: genotypes changed, alleles surviving
the prune, rounds per locus, wall-clock against `Off`. Sweep the bar and report **the depth at
which the fraction half starts to bind** (about 13 reads at the inherited values). Simulator
now; HG002 (should pay) and tomato at 3 reads (should be dangerous) apart, once unblocked.
*Depends:* E2. *Source:* spec Q3.

> **Checkpoint F:** the three questions have their reports, each number attributable to one
> axis; the real-data rows are either filled or explicitly blocked on step 6. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | hand-checkable one-sample rows; **the pinned Ewens total `(2 + θ)/(1 + θ)`** |
| B | by-hand marginal sums on enumerable cohorts; attributed A-vs-B differences |
| C | **exhaustive-as-oracle agreement**; linear retained-memory check |
| D | **no-op collapse bit-for-bit**; instrumented `builds` count at both settings |
| E | plant/terminate/free-when-off; **append-only column count** |
| F | the reports themselves, with B−A / C−B separated and blocked rows named |

## Out of scope (next plans)

- **Candidate selection (step 6)** — the spec that unblocks every real-data row above; it needs
  writing before it needs planning.
- **The sibling specs' remaining change measurements** — the dropped coefficient, STR
  contamination on/off, the zero-slippage path comparison, the per-locus re-fit distribution
  study, the GIAB regression under both prior impls — same instrument, same blocker; they need
  an owner once step 6 exists.
- **Whether a whole-cohort scorer can carry inbreeding** — the spec's recorded follow-up if a
  joint arm wins at `F = 0`; a spec question, not a plan step.
- **The tract QUAL calibration experiment** —
  [`calling_loop_ssr.md`](../spec/calling_loop_ssr.md) §3.3's obligation (owner, 2026-09-02):
  before the repeat-tract site quality is designed, observe the arms' behaviour — the
  inherited fold, a tract-specific correction, production's emission decision as comparator —
  on GIAB tract ground at 30×/50× and on the STR simulator, reporting calibration and a QUAL
  threshold sweep per arm, split by period and by parameter provenance. Runs once the STR
  loop's driver branch emits tract records; **belongs as a step of the STR loop's own
  implementation plan when that is written**, and is parked here so it is not dropped before
  then.
- **Q4–Q7's instrument-and-tune items** — the `passes` distribution, the cost crossover, the
  one-sample skip, the flat-first-pass count — cheap runs the loop's instruments already
  support, schedulable independently.
