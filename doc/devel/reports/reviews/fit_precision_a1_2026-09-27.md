# Code Review: fit_precision_a1
**Date:** 2026-09-27
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step A1 of `fit_precision.md` — each SNP/indel parameter's slope (score) at one position
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the diff `a81da5df..190c549e` on branch `fit-precision` (commit `190c549e`,
  "WIP: fit precision A1, for review", made only for the review and never referenced by a branch).
- **In-scope files:** [src/parameter_estimation/joint/fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs)
  (new), [src/parameter_estimation/joint/fit.rs](../../../../src/parameter_estimation/joint/fit.rs)
  (four `Scratch` fields filled in `one_position`, `mod information;`), the A1 implementation report
  and the PROJECT_STATUS entry.
- **Out of scope:** the rest of `fit.rs`, except the model functions the scores differentiate
  (`one_position`, `genotype_frequencies`, `ReadLogs::of`, `maximise_error_rate`, `fit_beta_shapes`),
  which were read as callees.
- **Categories dispatched**, three reviewers, each in its own worktree, each told to mutation-test:
  - *design conformance and claims* — naming, idiomatic, module structure, refactor safety, smells,
    defaults; conformance to plan step A1 and spec §3.2; the prose against the writing skill;
  - *reliability* — tests, and a 42-mutation run over `score_position` and `RuleSlopes`;
  - *correctness and numerics* — errors, float portability, extras; an independent re-derivation of
    every score, a check that no fitted bit moves, edge cases at the parameters' bounds.

  `unsafe_concurrency` was not dispatched: the step adds no threading (the scorer runs per position
  on the caller's thread); determinism across pool widths becomes a question at A2 and is recorded
  there.

## 2. Verdict

**Approve-with-changes.** No reviewer found an error in the mathematics: the correctness reviewer
re-derived every score from `one_position`, clamps included, and all agree with the code. The change
to `fit.rs` moves no fitted bit — a drawn 5-sample, 4,000-position cohort fitted at the review commit
and with the parent's `fit.rs` gave byte-identical `Debug` output (sha1 `9fc6dc7b…` both). What the
review found is test gaps on paths the drawn cohort never reaches, and prose whose numbers were
measured under a looser measure than it claimed.

## 3. Execution status

- `cargo test --release --lib parameter_estimation::joint::fit::information` in each reviewer's
  container: 6 passed, 0 failed, at `190c549e`.
- Reliability: 42 mutations run, 35 killed, 7 survived — 4 changed no behaviour on the fixture (two
  unreachable there, one a proven equivalence, one no-op), 3 moved the slopes by less than any test
  resolves (a deliberate tolerance probe and two posterior-skip mutations).
- Correctness: 14 of 36 planned mutations run before the loop was stopped, 13 killed, 1 survived
  (M9 below). The 22 not run are marked *Needs verification* in its file; the reliability run covers
  the same surface.
- Full-suite, clippy and fmt results are the author's (fixes-applied report).
- Findings labelled *Needs verification*: 1 (the slope at a homozygote excess of exactly 1, §7).

## 4. Open questions and assumptions

1. **A sample's second and later read groups.** The likelihood scores a sample's reads under its first
   read group's rates only, so the others carry no information and will have no standard error.
   Whether they are then written `defaulted` (the owner's standing rule for no-information
   parameters) or with their first group's value and error is an owner decision for checkpoint A.
   Affects Mi4.
2. **The maximisation's Beta-shape update** solves the digamma form, 2.9% off the log-likelihood's
   slope at the shipped 16 nodes. So the density-shape scores will not be centred on zero at the
   fitted shapes, which the outer-product estimator assumes. Owner decision at checkpoint A, with the
   size of that mean measured at A3. Cross-category note from correctness and reliability.

## 5. Top 3 priorities

1. **M1** — test the scorer on held-out reads and non-unit coverage odds; both terms are in its
   arithmetic and no fixture reaches them.
2. **M2** — test the scorer with the duplicated class switched off; a mutation that writes a slope
   into a class that does not exist survives.
3. **Mi2** — make "relative disagreement" relative and tighten the tolerance from 10⁻⁵ to 10⁻⁶; as
   written, a uniform error of 5 × 10⁻⁶ passes every test.

## 6. Findings

### Major

#### M1: src/parameter_estimation/joint/fit/information.rs:295, 445-451 — Held-out reads and coverage odds never reach the scorer in any test
- **Categories:** reliability, correctness (Minor there, coverage odds), design (cross-category note)
- **Confidence:** High
- **Problem:** `read_slope` subtracts `sample.on[4]` (reads on an indel, a spanning deletion or `N`),
  and the carrier branch weights by `odds`. The drawn cohort never emits a held-out read (0 of 2,400
  sample-positions at depth 8 and at 132), and every test passed `coverage_odds = &[]`. Counting
  held-out reads as ordinary errors, or dropping `odds` from the carrier's slope, left the suite green.
- **Suggested fix:** a test whose positions carry held-out reads and whose samples carry odds that are
  not 1. Run in review: passes on the tree at 1.38 × 10⁻⁷, fails the two mutants at 292 and 0.269.

#### M2: src/parameter_estimation/joint/fit/information.rs:329 — No test runs the scorer with the duplicated class off
- **Categories:** correctness; design and reliability (cross-category notes)
- **Confidence:** High
- **Problem:** every fixture fitted the duplicated class. Mutating `if carrier.is_some()` to `if true`
  (M9) writes about −1 a position into the duplicated-share slot of a run with no such class, and all
  six tests passed. The code itself is correct: a review-only test over all 17 coordinates of that
  layout agreed to 1.27 × 10⁻⁷, with the class's three slots exactly zero.
- **Suggested fix:** add that test, mapping coordinates by the layout without the class.

### Minor

#### Mi1: PROJECT_STATUS.md (Current focus) — "to 10⁻⁷ or better" overstates the precision
- **Categories:** design · **Confidence:** High
- **Problem:** the worst measured disagreement is 1.54 × 10⁻⁷ as printed, 1.7 × 10⁻⁷ once the carrier
  shape's is expressed relatively (Mi2).
- **Suggested fix:** "to within 2 parts in 10 million (1.7 × 10⁻⁷ at worst)".

#### Mi2: information.rs (test helper `largest_disagreement`) — "Relative disagreement" is absolute below a slope of 1, and the tolerance is 65 times the worst reach
- **Categories:** design, reliability (convergent) · **Confidence:** High
- **Problem:** the helper divided by `max(|numeric|, 1)`, so the carrier shapes (slopes −0.203 and
  0.063) were held to an absolute 10⁻⁵ — 5 × 10⁻⁵ and 1.6 × 10⁻⁴ relative — and the report's table
  called the result relative (carrier `a` printed 3.45 × 10⁻⁸, relatively 1.7 × 10⁻⁷). The bar of 10⁻⁵
  passed M39, a uniform 5 × 10⁻⁶ error in every excess slope.
- **Suggested fix:** divide by `max(|numeric|, 10⁻³)`; assert below 10⁻⁶; correct the report.

#### Mi3: report §2 item 2, information.rs module doc — the digamma comparison quoted at 12 nodes as though 12 were shipped
- **Categories:** design · **Confidence:** High
- **Problem:** the fit ships with 16 nodes, where the digamma form is 2.88% off, not 4.71%.
- **Suggested fix:** quote 2.9% at 16 beside 4.7% at 12, and pin both in the test.

#### Mi4: report §2 item 1; spec §3.2 — the later-read-groups deviation has no recommendation and no test
- **Categories:** design, correctness, reliability (convergent) · **Confidence:** High
- **Problem:** the report stated the consequence for the parameters file and left the decision bare;
  every fixture has one read group a sample, so nothing would notice a scorer reading a group other
  than the first. The spec's "k is 1 + 2 × its read groups" and the plan's A2 wording describe the
  wider block.
- **Suggested fix:** a recommendation in the report; a test that a second group's rate moves the
  log-likelihood by nothing; the spec note to follow the owner's ruling (Open question 1).

#### Mi5: information.rs `cohort` index module — its doc claimed `Parameters::coordinates`' order, false with the duplicated class off
- **Categories:** design · **Confidence:** High
- **Suggested fix:** say the orders match only with the class fitted; without it the last three slots
  are zero here and absent there.

#### Mi6: information.rs `cohort`/`sample` modules — bare `usize` constants where later steps need names
- **Categories:** design · **Confidence:** Medium
- **Problem:** A3, A5 and B3 print per-parameter names; nothing enumerates the slots or makes a match
  over them exhaustive.
- **Suggested fix:** enums `CohortParameter` and `SampleParameter` with `ALL` and `name()`.

#### Mi7: information.rs `score_position` — must be handed the same inputs `one_position` was, and nothing checks it
- **Categories:** design · **Confidence:** High
- **Suggested fix:** one value bundling the per-pass inputs, taken by both functions; `RuleSlopes`
  recording the shapes and node count it was built for.

#### Mi8: information.rs — the scorer re-derived parts of the model inline, apart from the functions they differentiate
- **Categories:** design; correctness (Nit) · **Confidence:** High
- **Suggested fix:** the derivatives beside `genotype_frequencies` and the read probabilities in
  `fit.rs`, with a unit test against a finite difference of each.

#### Mi9: information.rs `score_position` — 230 lines, nested five deep
- **Categories:** design · **Confidence:** High
- **Suggested fix:** one helper per branch. Optional if the parallel shape with `one_position` is
  preferred, but say so.

#### Mi10: report §5 — "see the commit message" for evidence the commit message did not hold
- **Categories:** design · **Confidence:** High
- **Suggested fix:** quote the suite's counts and the two checksum tests, with the log's path.

#### Mi11: PROJECT_STATUS.md — the status line read as if the whole programme were built
- **Categories:** design · **Confidence:** High
- **Suggested fix:** "through step A1 of 17", then the programme.

#### Mi12: report and doc comments — undefined terms for a geneticist reader
- **Categories:** design · **Confidence:** High
- **Problem:** "the E2 ruling", "SQUAREM", "the M-step", "the digamma form" and `ψ`, "the pass's
  attribution", "Gauss–Jacobi rule", each without a definition at first use.
- **Suggested fix:** define each once where first used.

#### Mi13: information.rs `score_position` doc — "No posterior is dropped for being small" is a rule no test can see
- **Categories:** reliability · **Confidence:** High
- **Problem:** adopting the pass's 10⁻¹² skip moves the summed slopes by about 10⁻¹¹ relative or not at
  all, because every share in the fixture is at least 0.006; at a share near its bound the total
  log-likelihood's finite difference drowns in rounding.
- **Suggested fix:** a per-position finite-difference test near a bound, or record the limitation on
  the shares test.

#### Mi14: information.rs test docs — two statements do not match the measured slopes
- **Categories:** reliability, design (Nit) · **Confidence:** High
- **Problem:** "thousands of log-likelihood units a unit of share" (measured −125 to 240, the
  duplicated share's −13); "every slope is large" (the carrier shapes' are −0.20 and 0.063).

### Nits

- "The three shares" followed by a list of four (module doc; the spec has the same miscount).
- The range test's "most sample-positions" was asserted as `> 1_000` of 2,400.
- Report §1's "4,880 lines" is the parent's length.
- One `///` over three `Scratch` fields; rustdoc attaches it to the first.
- `RuleSlopes.node` / `.ln_weight` hold derivatives; `node_slope` / `ln_weight_slope`.
- `super::super::` in the tests where crate paths read better.
- A literal `8` meaning `2 × samples` beside `COHORT_PARAMETERS`.
- `per_share` reads as "per unit of share"; it is the log posterior with the branch's own share
  divided out (correctness).
- The module name `information` holds only scores at A1 (fine if A2's blocks land there).

## 7. Out of scope observations

- **`fit.rs` `MAX_RECORDED_SPREAD` (pre-existing).** Its doc says the widest recorded range is 76–97,
  but the census ladder's rung 125–159 is 35 depths, so a census whose cap leaves that rung wider than
  32 panics in `next_position`'s assertion. The shipped cap of 124 never reaches it. Follow-up: owner,
  at checkpoint A.
- **`fit.rs` `fit_beta_shapes` (pre-existing).** Open question 2.
- **At a homozygote excess of exactly 1** (correctness, *Needs verification*): the analytic slope was
  −2.86 × 10⁶ against a one-sided difference of −6.44 × 10⁵; the reviewer judges the analytic value
  right and the difference unable to resolve a slope that steep. Either way, a parameter on that bound
  can carry per-position scores of order 10⁶ into the information. For A2/A3.
- **`read_slope` is recomputed for each genotype inside the node × sample loop**, though it depends
  only on (class, candidate, sample, genotype). Performance, once A2 runs the scorer over two million
  positions and thousands of samples.
- **Determinism across pool widths** does not arise at A1; at A2 the summed products must follow the
  pass's fixed chunk bounds and serial join, tested at one and several threads.

## 8. Missing tests to add now

- `the_slopes_hold_with_held_out_reads_and_coverage_odds` — M1; held-out reads at 583
  sample-positions and odds of 0.25, 3, 1 and 8.
- `the_slopes_hold_with_the_duplicated_class_off` — M2; every coordinate of the layout without the
  class, and the class's three slots exactly zero.
- `a_second_read_groups_rates_have_no_slope` — Mi4; a second read group's rates leave the
  log-likelihood bit-for-bit unchanged, and the first group's slope still matches.
- `the_models_slopes_are_the_derivatives_of_its_functions` — Mi8; each helper derivative against a
  finite difference of the function it differentiates.
- *Proposal, not run:* a per-position finite difference near a share's lower bound, to see the
  no-drop rule (Mi13).

## 9. What's good

- Every score is checked against a central difference of the pass's own log-likelihood, so the
  oracle is the function the fit maximises, not a restatement of the scorer
  ([information.rs](../../../../src/parameter_estimation/joint/fit/information.rs)).
- The Beta-shape slope is taken as the derivative of the quadrature rule actually used, with a test
  that pins why the textbook form would fail the oracle.
- The duplicated branch's stored value is computed beside its original sum, not reused in it, so the
  sum's association — and every fitted bit — is unchanged, and the correctness reviewer measured it
  byte-identical ([fit.rs](../../../../src/parameter_estimation/joint/fit.rs)).
- The depth-range test counts that its fixture really produces ranges (1,729 of 2,400) before trusting
  the comparison, and kills a mutation the depth-8 test cannot.

## 10. Commands to re-verify

- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit::information -- --nocapture`
- `scripts/dev.sh cargo fmt --check`
- `scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `scripts/dev.sh cargo test --all-targets --all-features --no-fail-fast`

Per-category findings and mutation evidence: `tmp/review_2026-09-27_fit_precision_a1/`
(`design.md`, `reliability.md`, `reliability_mut/`, `correctness.md`, `correctness_evidence/`).
