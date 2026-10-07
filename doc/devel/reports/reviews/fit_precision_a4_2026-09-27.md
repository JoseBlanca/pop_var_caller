# Code Review: fit_precision_a4
**Date:** 2026-09-27
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step A4 of `fit_precision.md` — the coverage of the SNP/indel fit's standard errors
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** the diff `121128f3..42400835` (`42400835`, a review-only commit on no branch).
- **In-scope files:** [src/parameter_estimation/joint/fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs)
  (tests only: `the_errors_mean_what_they_say`, `coverage_of`, `Coverage`, `label_of`, the true
  values, `full_matrix` lifted out of the step A3 comparison, `fitted_parameters_and_convergence`),
  the A4 implementation report, `PROJECT_STATUS.md`'s Step 4 block.
- **Out of scope:** the fit's Beta-shape update (`fit_beta_shapes` in `fit.rs`), which the review
  confirmed is the root cause of the density-shape finding; it moves fitted numbers and is for the
  owner.
- **Categories dispatched**, three reviewers, each in its own worktree:
  - *design and claims* — naming, idiomatic, module structure, refactor safety, smells, defaults;
    conformance to plan A4 and spec §3.6 item 3; every number in the prose (6 mutations; runs with
    the stop thresholds tightened and with the spread and a Newton-corrected coverage printed);
  - *reliability and errors* — tests and mutation (10 mutations);
  - *correctness and statistics* — what the measurement measures, the reading of it, float
    portability; direct probes of the flat point, the fixed point and the node count (6 mutations).

## 2. Verdict

**Request-changes.** Every number in the report and the doc comments was re-derived and is right
(all 32 table rows, every prose range, the probe's figures), and the mechanism the report proposes
for the density shapes was confirmed directly. What is wrong is three readings, each of which feeds a
checkpoint decision, and a test that cannot fail on anything it measures.

## 3. Execution status

- Each reviewer built its own worktree with its own `dev.sh`, release profile; each ran the coverage
  test with a shortened cohort count and its own instrumentation. The reliability reviewer ran
  `the_block_errors_against_the_full_matrix` on this commit: 1 passed, output identical to its doc.
- Mutations: reliability 10 run, 9 survived with a behaviour change, 1 killed; design 6 run, 5
  survived, 1 changed no behaviour; correctness 6 run, 5 survived, 1 killed. The one killed in both
  cases was including the duplicated class's slots, which panics in `label_of`.
- Findings labelled *Needs verification*: 0.

## 4. Open questions and assumptions

1. **Spec question 3** (is the outer-product estimator close enough?) is not answered by A4 for the
   kinds the spec doubted — the density shapes — because the fit does not stop at their maximum
   (M4, M5). Checkpoint A.
2. **The fix for the Beta-shape update** (out of scope) — the review measured that the gap between the
   fit and the likelihood's maximum grows, in errors, as the square root of the positions (M4).
   Checkpoint A, before plan step B1.

## 5. Top 3 priorities

1. **M1** — make the coverage test assert what it measured, so a later change cannot silently
   reverse the conclusions.
2. **M4** — correct the report's scaling inference: the gap grows with the positions, not with the
   samples; the kimura census has 10⁵–10⁶ positions.
3. **M3, M5** — correct the two readings checkpoint A decides from: the full matrix is too wide, not
   right; the shapes' errors hold at the maximum, the invariant share's at 4 samples do not.

## 6. Findings

### Major

#### M1: src/parameter_estimation/joint/fit/information.rs — the coverage test cannot fail on any defect in what it measures
- **Categories:** reliability, design, correctness (convergent)
- **Confidence:** High
- **Problem:** the only assertion is `coverage.estimates[0] > 0`. Surviving mutations, each with the
  printed numbers changed: a wrong truth index (clean rates −162 errors off), a wrong error slot, the
  Newton step's sign flipped, blocks and full matrix swapped, the full matrix reduced to the blocks,
  shapes `a` and `b` swapped, a within-one threshold of 1.5, every block error doubled. Three of these
  would reverse the report's conclusions with plausible-looking numbers.
- **Suggested fix:** at the default cohort count, where the draws are fixed by their seeds and the
  numbers change only with the code: coverage bands for the kinds that cover, the ordering of the two
  matrices where they differ, the Newton step's identity at 20 samples, a convergence count.

#### M2: information.rs — "every kind was estimated" is checked only for kinds already in the tally
- **Categories:** reliability
- **Confidence:** High
- **Problem:** dropping the homozygote excess from the tally passed; a cohort count of 0 printed
  "0 of 0 fits converged" four times and passed.
- **Suggested fix:** assert the tally holds exactly the eight expected kinds, each estimated by both
  matrices.

#### M3: report §2.2, §6; information.rs (the A3 comparison's doc); PROJECT_STATUS — "the full matrix's errors hold" overstates: they are too wide where they differ from the blocks
- **Categories:** reliability, design (convergent)
- **Confidence:** Medium (the pooled shares are not fully independent)
- **Problem:** at 4 samples and 3 reads the full matrix puts 0.743, 0.739, 0.780 within one error and
  0.973 within two — three to four sampling spreads wide; at 20 samples the blocks sit on nominal for
  the mismapped rates (0.677, 0.690) and the full matrix is the wide one (0.728, 0.722). A3's own
  mechanism predicts it: chance correlations between two samples' scores enter the full matrix as
  noise, which widens its inverse on average.
- **Suggested fix:** state the sizes and which side is off; in §6 state the trade (the full matrix is
  conservative at 4 samples and not needed at 20).

#### M4: report §2.3 — the scaling inference names the wrong axis; the gap grows with the positions
- **Categories:** correctness
- **Confidence:** High on the two measured directions; Medium on the extrapolation
- **Problem:** measured in `a`'s own units, the distance from the fit to the likelihood's maximum
  shrinks with the samples (mean 0.47 at 4 samples, 0.15 at 20) and does not shrink with the positions
  (0.15 at 5,000 positions, 0.14 at 20,000, 20 samples, 3 reads). In errors it grows as the square root
  of the positions: 1.30 to 2.47. The report said the opposite axis. The mechanism predicts it: the
  update's error and the information are both sums over positions.
- **Suggested fix:** replace the paragraph with the measurement; in §6, recommend fixing the shape
  update before B1; say that more nodes is not the fix (at 96 nodes the gap is still 0.16 to 0.88
  errors on four cohorts, while the likelihood itself is the same to four decimals from 16 nodes to
  192).

#### M5: report §2.3, §6; the test's doc — "the errors are not what is wrong" is supported for the shapes, not for the invariant share at 4 samples; the spread is never reported
- **Categories:** design, correctness (convergent)
- **Confidence:** High for the measurements, Medium for the reading
- **Problem:** spec §3.6 item 3 asks for the spread of the estimates, which the test does not print.
  The fitted shapes scatter 2 to 7 times less than their errors (spread of (estimate − truth)/error
  0.14 to 0.64). Moved by the Newton step to the likelihood's maximum, 70 to 85 in 100 shapes land
  within one error — their errors hold there — but the invariant share at 4 samples lands only 35 to
  45 in 100. The correctness reviewer found why: at 4 samples the likelihood's own maximum puts that
  share about 1.3 errors above the truth (six cohorts, a damped Newton search to the flat point), and
  the single Newton step lowers the likelihood there in all six. §6's answer to spec question 3 does
  not say that A4 cannot answer it for the shapes.
- **Suggested fix:** print the spread and the coverage at the Newton-corrected estimate; split the
  claim by kind; answer spec question 3 per kind.

#### M6: information.rs — `full_matrix` has no test of its own
- **Categories:** design
- **Confidence:** High
- **Problem:** made to drop the products between two samples — so that it equals the blocks — it
  survived both the A4 test and the A3 comparison, whose tolerances admit a zero gap. It is the
  reference behind A4's central finding.
- **Suggested fix:** a unit test on hand-built scores for two samples, checking a cross-sample entry.

### Minor

#### Mi1: information.rs — a bad cohort count is silently replaced by 200, and 0 passes
- **Categories:** reliability · **Confidence:** High
- **Suggested fix:** panic on a value that is not a positive whole number.

#### Mi2: information.rs, report — the cohort-count variable never reaches the container
- **Categories:** design · **Confidence:** High
- **Problem:** `dev.sh` forwards `CARGO_*`, `NG_*` and a few `RUST*` names only, so the documented
  shortcut silently runs all 800 cohorts.
- **Suggested fix:** rename it with the `NG_` prefix.

#### Mi3: information.rs — the printed line hides the count of full-matrix errors
- **Categories:** reliability · **Confidence:** High
- **Suggested fix:** print both counts.

#### Mi4: information.rs — every sample is drawn at the same values, so an estimate paired with the wrong sample cannot be detected
- **Categories:** reliability · **Confidence:** High (the present code is correct)
- **Suggested fix:** draw each sample at its own rates.

#### Mi5: information.rs — truths and estimates are read by position, labels by the layout's names
- **Categories:** design · **Confidence:** High
- **Suggested fix:** match on the named slots in both.

#### Mi6: report §2.3 — the 48-node probe is compared against another set of cohorts, and omits its most direct result
- **Categories:** reliability, design, correctness (convergent) · **Confidence:** High
- **Problem:** 40 cohorts at 48 nodes against 200 at 16; the estimates themselves moved toward the
  truth (shape `a` +1.29 to +0.51 at 20 samples and 3 reads), which is not quoted. The correctness
  reviewer measured why it only halves: the digamma form's error falls as about N^−2a with the node
  count N, 4.0 to 5.4-fold from 16 to 48 nodes at these cohorts' fitted `a`.
- **Suggested fix:** the same 40 seeds at 16 nodes; quote the estimates' shift and the node-count law.

#### Mi7: report §3 — the corner the A3 review asked for (two to four samples, 30,000 positions) is not the one measured
- **Categories:** design · **Confidence:** High
- **Suggested fix:** add a regime below 4 samples, or say it is unmeasured.

#### Mi8: report §2.3 — the stop rule is not ruled out as the cause, though the review measured it
- **Categories:** design · **Confidence:** High
- **Problem:** with the stop thresholds 10⁴ times tighter (375 to 987 passes a start instead of 48 to
  117), shape `a` sat at +1.558 errors against +1.557 and the flat point at −2.047 against −2.045
  (4 samples, 3 reads, 20 cohorts). The correctness reviewer found the same: 300 further plain passes
  move `a` by less than 10⁻⁵.
- **Suggested fix:** add it; it separates this defect from plan step B1's.

#### Mi9: report §2.3 — the invariant share's coupling to the shapes stated as fact, without its direction
- **Categories:** design · **Confidence:** Medium

#### Mi10: report §2.3, §6 — "1 to 2 errors" where the measured range is 1.05 to 2.30
- **Categories:** design · **Confidence:** High

#### Mi11: the test's doc — the digamma cause stated as fact where the report says "likeliest"
- **Categories:** design · **Confidence:** High
- **Note:** the correctness reviewer's probes now establish it (the fit is an exact fixed point of the
  update, where the digamma form is zero to 0.012 and the likelihood's slope in `a` is −8 to −20).

### Nits

- `Coverage`'s `[_; 2]` arrays indexed by the literals 0 and 1; `estimates` names a count; the
  flat-point sum divided by the full matrix's count without saying so.
- `label_of(i)`: `slot` matches the layout's vocabulary.
- "each fitted until it converges": one of 800 stopped at the pass limit.
- "a Beta cut off at 0.2": the draw floors at 0.2.
- The true values appear as literals in three tests.
- The seeds of the 3-read and 30-read regimes overlap (173 of 200); harmless, since the depth draw
  diverges the streams.
- "0.94 to 0.97" within two is 0.938 to 0.965.

## 7. Out of scope observations

- `fit.rs` `fit_beta_shapes` — the Beta-shape update solves the digamma form, so the fit's fixed point
  is not the likelihood's maximum (M4, M5). Moves fitted numbers; owner, checkpoint A.

## 8. Missing tests to add now

- `coverage_of_tallies_every_fitted_kind_in_every_regime` (M2).
- `the_full_matrix_sums_the_products_between_two_samples` (M6).
- The coverage assertions at the default count (M1): bands for the kinds that cover; the full
  matrix's within-one share above the blocks' at 4 samples and 3 reads; |off centre + flat point|
  below 0.5 for the shapes at 20 samples; a convergence count.

## 9. What's good

- The flat-point column separates "the errors are wrong" from "the fit stops beside the maximum", and
  the correctness reviewer's independent search to the flat point agreed with it to about 5% at 20
  samples.
- Every figure in the report was transcribed correctly from its log, and the log is kept.
- The measurement reports both matrices side by side, which is what exposed that each is off in a
  different direction.

## 10. Commands to re-verify

- `scripts/dev.sh env NG_FIT_PRECISION_COVERAGE_COHORTS=3 cargo test --release --lib the_errors_mean_what_they_say -- --ignored --nocapture`
- `scripts/dev.sh cargo test --release --lib the_errors_mean_what_they_say -- --ignored --nocapture`
- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit -- --nocapture`

Per-category files: `tmp/review_2026-09-27_fit_precision_a4/` (`design.md`, `reliability.md`,
`correctness.md`).
