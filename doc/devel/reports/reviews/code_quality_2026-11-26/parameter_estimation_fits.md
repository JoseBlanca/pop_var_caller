# Code-quality map: the model fits in `parameter_estimation/joint/` (2026-11-26)

Scope:

- the SNP/indel fit: `fit.rs` and `fit/` (`information.rs`, `standard_errors.rs`, `settled.rs`, `determined.rs`);
- the STR fit: `ssr_fit.rs` and `ssr_fit/` (`penalised_fit.rs`, `penalised_fit/run_log.rs`, `expected_counts.rs`, `read_group_stutter.rs`, `tract_exclusion.rs`).

This review is read only. Nothing was built or run, every finding is reasoned from the code, and every `path:line` was checked against the worktree at `f84129cf`.

**How big these files are without their tests.** Most of each file is tests:

| File | Production lines | Tests start at line |
|---|---|---|
| `fit.rs` | 5,517 | 5518 |
| `fit/standard_errors.rs` | 1,687 | — |
| `fit/information.rs` | 1,069 | — |
| `ssr_fit/penalised_fit.rs` | 4,097 | — |
| `ssr_fit/expected_counts.rs` | 1,746 | — |
| `ssr_fit.rs` | 1,646 | — |

The whole area holds about 18,800 production lines.

**Words used below.**

- **Pass**: one sweep over every census position at fixed parameters.
- **EM step**: a pass followed by re-fitting every parameter to the counts that pass collected.
- **SQUAREM**: an accelerator that takes two EM steps and then jumps further along the direction they took.
- **Newton step**: a jump to where a quadratic fitted to the log-likelihood's slope and curvature peaks.
- **Newton decrement**: the log-likelihood that Newton step would gain. The SNP/indel fit stops when it falls below 0.01.

## Summary

In priority order:

1. **The SNP/indel fit re-fits each sample's homozygote excess against quadrature points it did not count at.** The excess is re-fitted from genotype counts the pass collected at the allele-frequency points of the old Beta shapes. The fit reads the frequencies at those points from a rule rebuilt after the shapes have already moved (`fit.rs:5074-5128`). Each EM step's excess update is therefore slightly wrong. The fit still stops at the right place, because its stopping rule reads the true likelihood, but it may take extra passes. This may also be the source of the non-monotone EM that `JUMP_SLACK` was introduced to tolerate (a hypothesis, not tested). The fix is one line; it changes output in the last digits.
2. **Two STR-fit defects in the SQUAREM cycle can leave reported numbers that no round checked.**
   - After a jump, the next round reuses class counts gathered before the jump.
   - A fit that reaches its round limit on a jump returns the jumped numbers, while its curves and standard errors still describe the round before.
   - Separately, a fit that runs out of rounds is reported as an ordinary progress line, not as a warning.
3. **`fit.rs` should split, along the steps of the algorithm rather than by size.** Move out the Beta quadrature (200 lines, self-contained), the maximisation step (300 lines) and the per-position expectation pass (about 1,100 lines, the natural sibling of `information.rs`). The driver — starts, SQUAREM, Newton tries — stays in `fit.rs`. `fit_jointly` (610 lines) and `one_position` (580 lines) should also be cut into their phases.
4. **The SNP/indel fit and the STR fit share almost no numerics worth merging.** Their golden-section searches, SQUAREM steps and Newton steps each stop under different rules, and merging any pair would change output somewhere. The real duplication is elsewhere:
   - `fit.rs` against `contamination.rs`: the slice log-sum-exp, and the Jacobi eigen-solver;
   - inside `standard_errors.rs`: the arrow elimination, about 100 lines, written twice;
   - in a few small places inside each fit.
5. **Smaller items:**
   - The SNP/indel error-rate maximisation repeats the read-probability formulas inline, beside the functions that exist so it would not.
   - The STR side carries about 210 lines of test-only scoring code that its tests exercise in place of production.
   - The STR binomial tail sum can walk every read count when its first term underflows.

## Bugs

### 1. The homozygote excess is re-fitted at frequencies the pass did not count at

- **Where:** `fit.rs:5074-5086` moves the density's shapes `a`, `b`. `fit.rs:5114-5118` then builds `BetaQuadrature::new(parameters.density.a, parameters.density.b, …)` from the moved shapes. `fit.rs:5124-5127` maximises each excess against `quadrature.nodes[node]`. The counts it maximises over are `statistics.genotypes[s][node][j]`, which were accumulated at `fit.rs:4858` against `model.quadrature`, the rule at the shapes before this step.
- **What:** the M-step (the part of an EM step that re-fits a parameter) for the excess has to maximise `Σ count(node, genotype) · ln P(genotype | f_node, excess)` with each `f_node` the frequency at which that count was collected. The code pairs each count with the frequency of the same-numbered point of a different rule. While the shapes move, every pass's excess update is therefore a little off the true M-step.
  - **Scenario:** any cohort of two or more samples whose density shapes are still moving, which is every early pass.
  - **Why the answer survives:** when the shapes stop moving the two rules coincide, and the stopping rule (the Newton decrement on the full likelihood, `settled.rs:73`) does not depend on the M-step being exact. The converged answer is still the maximum. What the mismatch can cost is passes, and possibly the small non-monotone drift that `JUMP_SLACK` (`fit.rs:3849-3859`) tolerates. That comment blames "one of the maximisations is close to, but not exactly, the maximum of what it is handed". This mismatch would produce exactly that, but it is a hypothesis.
- **Proposal:** build the excess's rule from the shapes before `step_beta_shapes` changes them. Capture `(a, b)` at the top of `maximisation`, or move lines 5113-5129 above line 5074. Then check whether `JUMP_SLACK` can drop: rerun the ten-sample cohort described in its comment.
- **Downside:** the iteration path changes, so every fitted number moves within the stopping tolerance, and test fixtures that pin exact values may need updating. If the drift does not go away, the `JUMP_SLACK` question stays open.
- **Effect:** output changes in the last digits; the fit may get faster. Needs owner approval as an output change.
- **Cost:** S. **Confidence:** high on the mechanism (traced end to end); medium on the link to `JUMP_SLACK`. **Status:** likely.

### 2. STR fit: after a SQUAREM jump, the next round reuses counts gathered before the jump

- **Where:**
  - `ssr_fit/penalised_fit.rs:1981-1987`: the `reuse` closure keeps a class's `last_counts` when `class_change[point] < tolerance && multiplier_change < tolerance`;
  - `:1797-1813`: the jump calls `set_state(&jumped)`;
  - `:1817-1822`: the fallback resets `class_change` to infinity.
- **What:** the jump can move each number by up to `2·50·|r| + 50²·|v|`, where `r` and `v` are SQUAREM's first step and the change between its two steps, and 50 is its largest step multiple. The jump does not reset `class_change`. In the round after it, a class that had nearly settled reuses counts taken at the numbers from before the jump. Those counts then feed both the curve update and `composite_log_likelihood`, which decides whether the jump is kept (`:1817`). The fallback branch does reset the change, so the need for a reset was known and missed here.
- **Proposal:** in the `OneRound` arm that jumps, set every `class_change` to infinity, and `multiplier_change` too when multipliers are fitted.
- **Downside:** more class-rounds are counted near convergence, which is slightly slower.
- **Effect:** possible small output change; slightly slower.
- **Cost:** S. **Confidence:** medium. **Status:** suspected (path traced, effect not measured).

### 3. STR fit: reaching the round limit right after a jump returns numbers that no round checked

- **Where:** `ssr_fit/penalised_fit.rs:1958` (`while self.rounds < fit.max_rounds`) and `:2239` (`self.squarem = SquaremCycle::Idle`, which discards the fallback).
- **What:** if the last allowed round is the one that jumps, the slippage, share and concentration numbers are left at the jumped values. `last_counts` and the curve values are still the previous round's. The reported numbers come from the first set, while the curves and standard errors are built from the second, so they disagree. The limit is 2,000 rounds, shared across every call to `count`; how often production reaches it is unknown.
- **Proposal:** when the loop ends in `SquaremCycle::Jumped { fallback, .. }`, call `set_state(&fallback)` before leaving.
- **Downside:** none found.
- **Effect:** output changes only for fits that hit the limit.
- **Cost:** S. **Confidence:** medium. **Status:** suspected.

### 4. STR fit: a fit that runs out of rounds does not warn

- **Where:**
  - `ssr_fit/penalised_fit/run_log.rs:357-366` writes "used all N rounds" as an ordinary progress line;
  - `ssr_fit/penalised_fit.rs:1225-1229` sends only the heavy-exclusion lines as `Line::Warning`;
  - `calling/parameters_file/to_toml.rs:530` records `settled` only for the stutter-tail rule.
- **What:** a motif length whose fit did not converge produces numbers indistinguishable from a converged one, both in the log and in the parameters file.
- **Proposal:** send the "used all N rounds" line as a `Line::Warning`; that changes the log only. Writing a per-motif-length `settled` flag into the file would be an output change, so that is the owner's decision.
- **Effect:** log only (none for the file unless the owner approves the flag).
- **Cost:** S. **Confidence:** high. **Status:** likely.

### 5. STR fit: a NaN in a fitted number can neither stop the loop nor fail the jump test

- **Where:**
  - change measures fold with `f64::max`, which drops NaN (`ssr_fit/penalised_fit.rs:2202`, `:2204`, `:1807`, `:1850`, `:995-996`);
  - `log_likelihood < to_beat` is false for NaN (`:1817`);
  - `state()` takes `ln(p/(1−p))` with no clamp (`:1721`);
  - `set_state` normalises the shares by `exp` without taking out the largest first (`:1752-1757`);
  - the multipliers are `exp(next())` with no clamp (`:1762`).
- **What:** a logit at the ±40 clamp rounds `logistic` to exactly 1.0. The next `state()` then holds +∞, the jump arithmetic turns every number into NaN, and a NaN round reports zero change, so the loop "converges" on NaN. No input is known to reach this.
- **Proposal:**
  - clamp the logit in `state()`;
  - normalise the shares with log-sum-exp;
  - clamp the multipliers;
  - count a NaN change as infinite;
  - refuse a jump whose log-likelihood is not finite.
- **Downside:** a few more branches.
- **Effect:** none on ordinary runs.
- **Cost:** S. **Confidence:** medium on the mechanism, low that it is reachable. **Status:** suspected.

### 6. SNP/indel fit: the census length is read from the first sample's first section only

- **Where:** `fit.rs:1058-1063` (`EvidenceCursor::position_count`) is used by `fit_jointly` (`:2163`) and by every pass (`:4103`).
- **What:** a sample can lend no ordinary-position section. The code accepts this:
  - `fit.rs:2474-2477` and `information.rs:2410-2412` describe it, as a sample whose census holds repeat tracts only;
  - the test at `fit.rs:7951` puts such a sample last.
  If such a sample came **first**, `position_count` would return 0. Every pass would then visit no position, and every rate would come back `Defaulted`, silently. As a side effect, `fit.rs:2588-2590` would divide by zero positions and write NaN genotype rates.
  - **Reachability:** the production census builder gives every declared read group an ordinary-position section (`census.rs:2980-2990`), so this needs a sample with no declared read group, or a hand-built census.
- **Proposal:** take the length from the first sample that has a section, and assert that every section agrees. Use `statistics.positions.max(1.0)` at `:2588-2590`, as `maximisation` already does at `:5035`.
- **Effect:** none on current outputs.
- **Cost:** S. **Confidence:** high on the mechanism, low on reachability. **Status:** suspected.

### 7. STR exclusion test: the binomial tail sum can walk every read count when its first term underflows

- **Where:** `ssr_fit/tract_exclusion.rs:280-292`.
- **What:** the loop stops when `term < 1e-17 * total`. If the first term underflows to 0, then `total == 0`, the test `0 < 0` is false, and the loop runs to `trials`, calling `lgamma` three times a step. For example, at 60,000 cohort-summed reads, a slip rate of 0.005 and 10,000 off-length reads, the first term's logarithm is about −26,000. That is about 180,000 `lgamma` calls for one tract, repeated in each stutter-tail round (up to 10).
- **Proposal:** once past the mode (`k > n·p`), stop when a term is 0. Do not stop on the first term alone: below the mode, later terms can still be non-zero.
- **Effect:** none on output; faster on deep, large cohorts.
- **Cost:** S. **Confidence:** high on the mechanism; the cost was not measured. **Status:** likely.

## Architecture

### 8. Split `fit.rs` along the steps of the algorithm

- **Where:** `fit.rs`, 5,517 production lines plus a 3,100-line test module. Its sections:

  | Lines | What |
  |---|---|
  | 83-1000 | public results and configuration |
  | 1000-1290 | the evidence cursor |
  | 1290-1592 | the read model |
  | 1594-2028 | parameters and pass statistics |
  | 2029-3535 | the driver: `fit_jointly`, `maximise`, Newton tries, log lines |
  | 3535-3881 | the EM step, SQUAREM and the coordinate bounds |
  | 3883-5001 | the per-position expectation pass |
  | 5003-5306 | the maximisation step |
  | 5308-5515 | Beta quadrature, the eigen-solver, trigamma |

- **Proposal:** split in this order, each split its own reviewed change:
  1. `fit/beta_quadrature.rs`: `BetaQuadrature`, `gauss_jacobi`, `symmetric_eigen`, `trigamma` and the test-only `digamma` (`:5312-5515`). This code is pure, takes numbers and returns numbers, and is already used from `information.rs` (`RuleSlopes`).
  2. `fit/maximisation.rs`: `maximisation`, `maximise_error_rate`, `maximise_hom_excess`, `golden_section`, `step_beta_shapes`, `fits_homozygote_excess` (`:5003-5306`). It reads `Statistics` and writes `Parameters`, nothing else.
  3. `fit/expectation.rs`: `PassModel`, `BranchShares`, `CarrierRule`, `Scratch`, `one_position`, `expectation_pass_every`, the chunking constants, plus the read model of `:1290-1592`. This is the half `information.rs` already reads field by field, so the two become siblings.

  `fit.rs` keeps the public types, `Parameters`, and the driver (about 2,500 lines). Moving the public types out is possible but buys little.
- **Downside:**
  - `information.rs` reads about 20 private `Scratch` fields. After step 3 they need `pub(super)`, which is annotation churn, not new coupling.
  - The test module uses `super::*` and would need a few more imports.
  - Steps 1-2 have no real downside. In step 3, `one_position` must stay `#[inline]`-friendly; the move does not change that.
- **Effect:** none (code moves only).
- **Cost:** S + S + M. **Confidence:** high.

### 9. Cut `fit_jointly` and `one_position` into their phases

- **Where:** `fit_jointly` at `fit.rs:2029-2640` (610 lines); `one_position` at `fit.rs:4424-5001` (580 lines); `maximise` at `:2772-3113` (340 lines).
- **What and proposal:**
  - **`fit_jointly`** does four things in sequence:
    1. log what the evidence holds (`:2096-2156`);
    2. run the starts and choose the winner (`:2157-2388`);
    3. compute the standard errors and contamination (`:2397-2484`);
    4. turn `FittedCohort` into `JointFit` (`:2499-2639`).

    Each phase should be its own function. Phase 4 has no dependence on the evidence borrow, so it can be a `JointFit::from_fitted(...)`.
  - **`one_position`** already marks its phases with comment banners: read likelihoods (`:4475`), the frequency integral (`:4555`), the carrier integral (`:4585`), branches and classes (`:4629`), attribution (`:4734`). Each can become a function over `&mut Scratch`.
  - **`maximise`:** move the Newton-try block (`:2887-2984`) into a function returning "accepted / rejected"; that alone removes 100 lines from the loop.
- **Downside:** `one_position` is the innermost loop. Splitting it risks lost inlining, so run a speed check before and after (the "possible" below). The other two are not hot.
- **Effect:** output none if the arithmetic order is kept; speed possible for `one_position`.
- **Cost:** M. **Confidence:** medium.

### 10. No wrong-way dependencies; leave the module graph alone

- **Where:** imports at `fit.rs:45-81`, `fit/*.rs`, and the STR files.
- **What:** nothing in this area imports `crate::run` or `crate::cli`. The only mention of `crate::run` is a doc link (`fit.rs:3684`). The STR fit does not import `fit`: the homozygote excess reaches it as a plain `&[f64]` built in `run/census_fit.rs:348-362`, as `joint/mod.rs:32-34` claims.
- **Proposal:** none.

### 11. Narrow `invert_identified`'s visibility

- **Where:** `fit/standard_errors.rs:1594` and `:1607` are `pub(in crate::parameter_estimation::joint)`, but every caller is in `standard_errors.rs`.
- **Proposal:** make them private.
- **Effect:** none. **Cost:** S. **Confidence:** high.

## Duplication

**The SNP/indel fit against the STR fit: leave their optimisers apart.** Here is each pair:

| Numeric | SNP/indel fit | STR fit | Verdict |
|---|---|---|---|
| Golden section | `fit.rs:5214`, stops at a width of 1e-9, returns the bracket's middle (which the Newton step's box relies on, `fit.rs:3698-3705`) | `penalised_fit.rs:3590-3626`, exactly 12 steps over `CurveFit` values, after a grid | leave apart |
| SQUAREM | `fit.rs:3644` and `:3827`, step length clamped to `[1, reach]` with reach ×4 or ÷4, judged against the cycle's start with 1 unit of slack, re-tried shorter | `penalised_fit.rs:1770-1827`, clamped to `[1, 50]`, judged one round later, falls back once | the formula is the same 3 lines; leave apart |
| Newton | box-constrained active-set solve over the information (`standard_errors.rs:1067`), then up to 4 tries along it (`fit.rs:3192`) | banded penalised Newton with up to 40 halvings (`penalised_fit.rs:3718`); one-dimensional for the multipliers (`read_group_stutter.rs:100`) | leave apart |
| Convergence | Newton decrement below 0.01 | largest change below 1e-5 | different questions |

Sharing any of these couples two models that tune their stopping rules independently, saves 20-40 lines each, and changes output. One shared helper is still worth having (item 13).

### 12. One slice log-sum-exp, three identical copies

- **Where:** `fit.rs:1586`, `contamination.rs:1404`, `calling/quality/mod.rs:856`.
- **What:** three identical functions of 6-7 lines each, token for token. The STR fit has none: it takes the largest out of a vector and keeps the scaled vector, a different operation.
- **Proposal:** add `pub fn ln_sum_exp(values: &[f64]) -> f64` to `crate::float`, beside `ln` and `exp`, and call it from all three.
- **Downside:** `calling` and `parameter_estimation` then share a helper, which is the right home for it anyway.
- **Effect:** none (the arithmetic is identical).
- **Cost:** S. **Confidence:** high.

### 13. One dense inverse by Cholesky factorisation, two copies (SNP/indel and STR)

- **Where:** `fit/standard_errors.rs:1644` (`inverse_by_cholesky`) and `ssr_fit/penalised_fit.rs:3943` (`dense_inverse`), about 40 lines each.
- **What:** both use the same factorisation; they invert in a different summation order.
- **Proposal:** move `inverse_by_cholesky` to a small `joint/linear_algebra.rs` and call it from the STR sandwich estimate (`penalised_fit.rs:3925-3939`). To keep the file identical, either give the shared function the STR side's order as an option, or accept last-digit changes in the STR standard errors with owner approval.
- **Downside:** shared code between the two fits, but it is pure linear algebra with no model in it.
- **Effect:** none if the STR order is kept; otherwise the STR standard errors change in the last digits.
- **Cost:** S. **Confidence:** high.

### 14. The Jacobi eigen-solver, twice

- **Where:** `fit.rs:5429-5484` (`symmetric_eigen`) and `contamination.rs:1414-1467` (`leading_eigenvectors`).
- **What:** the same cyclic Jacobi rotation, 40 lines each. They differ in:
  - the stopping threshold: 1e-22 against 1e-18;
  - the skip threshold: 1e-18 against 1e-15;
  - what they return: all eigenpairs ascending against the leading vectors.
- **Proposal:** add one `jacobi_rotations(matrix, n, stop_below, skip_below) -> (diagonal, vectors)`, with each caller keeping its own thresholds and ordering. It belongs in the same `linear_algebra.rs` as item 13, or in `beta_quadrature.rs` if contamination is left out.
- **Downside:** the two thresholds become arguments, and a reader must check that the right pair goes to each caller.
- **Effect:** none (same operations with the same thresholds).
- **Cost:** S. **Confidence:** high.

### 15. The standard errors' block elimination is written twice

- **Where:**
  - `fit/standard_errors.rs:301-426` (`of_the_blocks`) against `:1329-1440` (`arrow_solve`);
  - `:238-298` (`of_the_whole_matrix`) against `:1287-1320` (`whole_matrix_solve`).
- **What:** the matrix has an "arrow" shape: the cohort's parameters along one edge and each sample's own parameters in a diagonal block. Both functions in each pair do the same things:
  - choose the parameters with information;
  - invert each sample's block;
  - form `A_s⁻¹ B_s`, where `A_s` is the sample's own block and `B_s` its pairing with the cohort's parameters;
  - subtract `Σ B_sᵀ A_s⁻¹ B_s` from the cohort block and invert what is left.

  The loops at `:319-374` and `:1342-1399` are the same code. The only differences are that one also carries a slope vector and one honours a `held` mask. That is about 100 duplicated lines. The two copies must also choose identical parameters, or the step and the errors disagree. Today they do: compare `:325-334` with `:1351-1357`.
- **Proposal:** add an `ArrowElimination::new(sums, held)` holding each sample's kept parameters, inverse and `A_s⁻¹ B_s`, and the reduced cohort inverse. The errors read its diagonals; a solve applies it to a slope. Do the same for the whole matrix.
  - **Side benefit:** `cohort_covariance` (`:877-899`) runs the elimination 8 times on the same matrix, and `newton_step` repeats what `StandardErrors::of` just did. Built once, it is reused.
- **Downside:** one more type to learn, and `bounded_step` changes `held` every round, so it still builds a new elimination each round. The speed gain is small next to a pass.
- **Effect:** none, provided the summation loops are kept in the same order. Verify with the output-identity check.
- **Cost:** M. **Confidence:** medium.

### 16. The error-rate maximisation repeats the read-probability formulas inline

- **Where:** `fit.rs:5148-5150` against `reference_read_probability` and `candidate_read_probability` at `:1458-1475`.
- **What:** the two functions exist, with their slopes beside them, "so the standard errors' slopes follow any change to the read model" (`:1463-1464`). `maximise_error_rate` writes the same two expressions out again. A change to the read model made in the functions would silently miss the maximisation.
- **Proposal:** call the functions. The expressions are the same token for token, so the bits are the same.
- **Effect:** none. **Cost:** S. **Confidence:** high.

### 17. Smaller copies inside each fit

All of these have no effect on output (the arithmetic is the same) and cost S each.

- **SNP/indel, the product over samples with rescaling, written twice:** `fit.rs:4562-4581` (frequency points) and `:4599-4624` (carrier points). A helper that multiplies a sequence of terms with the 1e150 rescale, `fn ln_of_product(terms) -> f64`, serves both.
- **SNP/indel, "credit a library's tally with expected reads", written three times:** `fit.rs:4795-4807`, `:4880-4892` and `:4965-4977`. Make it one helper.
- **STR, the stable binomial log-likelihood, written twice:** `penalised_fit.rs:3099-3111` (`binomial_at_logit`) and `:3164-3180`. The second should call the first.
- **STR, the "one tract of a class" `StratumEvidence` literal, written twice:** `penalised_fit.rs:2457-2473` duplicates `one_tract_of` (`:2965-2983`).
- **STR, the second-difference smoothness penalty, assembled three times:** dense at `penalised_fit.rs:2514-2529` and `:2658-2668`, banded at `:3827-3849`. A helper taking a closure that adds one entry serves all three. Its downside is the closure argument. Confidence: medium.

## Readability

### 18. STR: about 210 lines of test-only scoring code sit in the production file, and four tests exercise only that code

- **Where:** `ssr_fit.rs:858-1070`: `TractLikelihoods`, `add_ln_reads`, the test-only parts of `ReadDistributionsByLevel`, the test-only `ln_gamma` (`:1078-1103`) and `normalise` (`:832-839`).
- **What:** production scores reads through `expected_counts::ReadModel`. These four tests check code that no run executes:
  - `:1891`;
  - `:1928`, the test that the order a sample's libraries are listed in does not matter;
  - `:1960`;
  - `:2017`.

  The only production use of `ReadDistributionsByLevel::of` is the pooled arm of `ReadModel::of` (`expected_counts.rs:581-589`). It always runs at multiplier 1, so it can compute the distributions directly, as the per-read-group arm does at `:626-628`.
- **Proposal:** delete `TractLikelihoods` and `ReadDistributionsByLevel`, port the library-order test to `ReadModel` rather than dropping it, and move `ln_gamma` and `normalise` into `bench_fixtures`.
- **Effect:** none. **Cost:** S-M. **Confidence:** high.

### 19. STR: dead items and stale references

All of these have no effect on output, cost S and are high confidence.

- **Never or barely constructed:**
  - `StratumError::NotIdentified` (`ssr_fit.rs:599`) is never constructed;
  - `StratumError::NoShare` (`:610`) is constructed only in tests.
- **Public but used only by tests, yet compiled into the binary:** `read_group_stutter.rs:94` (`fit_multiplier`) and `:235` (`scale_to_a_typical_read_group`).
- **Visibility wider than needed:** `CurveState`, `CurveData`, `StrengthSearch` and `MultipliedReads` are `pub(crate)` but used only inside `penalised_fit.rs`.
- **Empty or broken:** `impl ParameterCurve {}` (`penalised_fit.rs:600`) is empty, and the doc at `:620` cites a `ParameterCurve::at` that does not exist.
- **Doc references to functions that no longer exist:**
  - `expected_counts.rs:24`, `:48` and `:54`;
  - `ssr_fit.rs:595` and `:860`;
  - `penalised_fit.rs:735`.
- **Wrong claim about the bench fixtures:** `ssr_fit.rs:1428-1434` says they are also built under a `bench-fixtures` feature and cites `examples/ng_str_fit_bench.rs`. The attribute at `:1435` is `#[cfg(test)]` only, and the example does not exist. The same missing example is cited at `run/census_fit.rs:345`, `cli/estimate_parameters.rs:706`, `benches/ng_joint_fit_perf.rs:15` and `Cargo.toml:257`.
- **Misplaced doc comment:** `penalised_fit.rs:1712-1716` documents the multipliers' step but sits on `state()`; it belongs on `update_multipliers` (`:1829`).

### 20. STR: one quantity under four names

- **Where:** `penalised_fit.rs`, throughout. The concentration `c` also appears as "spread" (`ln_spread`, `spread_held`, `spread_from_pairs`), as the between-locus share `1/(1+c)`, and as "F". `pair_tables` passes `-logit` into an argument named `ln_spread` (`:2646`, `:2783`).
- **Proposal:**
  - keep "concentration" and "between-locus share";
  - rename the `spread` identifiers;
  - add `share_of(concentration)` and `concentration_of(share)`, which also replace the duplicated clamps at `:2180-2184` and `:2860-2863`.
- **Downside:** the rename touches tests.
- **Effect:** none. **Cost:** S-M. **Confidence:** medium.

### 21. STR: name the clamp values

- **Where:** fall-off `0.99` is written four times (`penalised_fit.rs:1747`, `:2177`, `:2340`, `:2856`). The level's `0.95` is written once, in `set_state` only (`:1745`).
- **What:** an ordinary round (`count`, `:2175`) does not clamp the level, so a jump and a round obey different bounds.
- **Proposal:** make them constants beside `LARGEST_UNBOUNDED`. Whether `count` should also clamp the level would change output, so that question goes to the owner.
- **Effect:** none for the constants. **Cost:** S. **Confidence:** high.

### 22. STR: split `penalised_fit.rs` at the curve fitting only

- **Where:** `penalised_fit.rs`, 4,097 production lines. The curve-fitting code (`CurveData`, `CurveState`, `BendPrior`, the strength search, Newton) is at `:2939-3893`, and the dense and banded matrices at `:3895-4095`. Both meet the rest through `CurveState::refit` and `penalised_sandwich` only.
- **Proposal:** move them to `penalised_fit/curve.rs` (about 950 lines) and `penalised_fit/matrices.rs` (about 200 lines; or the shared `linear_algebra.rs` of item 13).
- **Leave in place:**
  - `Period`'s standard errors: they read about 15 of its private fields, and moving them would spread one type over two files;
  - `Period::count` (`:1928-2243`, 315 lines): extract its maximisation step (`:2050-2161`) and leave the loop where it is.
- **Effect:** none. **Cost:** M. **Confidence:** medium.

## Looked at, left alone

- **The SNP/indel fit's thread-count independence holds:** chunks are cut from the census length, never from the pool's width (`fit.rs:4110-4116`), and are joined in a fixed tree or in position order (`:4205-4252`). The STR side is collected in order as well.
- **`Coordinate::newton_range`'s claim that a golden-section parameter never lands exactly on an end holds** for every path that sets `parameters`: maximisation output, a Newton target clamped inside, a Newton try between the two. A SQUAREM jump's `bounded()` can round to an end, but its output is never kept without a maximisation after it.
- **`fit_jointly`'s `best.expect("the first start always competes")`** holds for every run order: the Continuing run that is set aside re-queues the starting points, and a Finish that does not converge competes and re-queues the others. It would panic only on an empty `starting_points`, which no command line can produce.
- **`Parameters::continuing` does not clamp an earlier run's numbers into bounds** (`fit.rs:1642-1686`). This is harmless: the first maximisation clamps them, and `Coordinate::free` clamps before taking a logit.
- **`count_times_ln`'s comment says a zero count needs no guard** (`fit.rs:1330-1333`): every logarithm is built clamped at `MIN_POSITIVE`. Correct.
- **The `!position_ln.is_finite()` branch in `one_position`** (`fit.rs:4692-4724`) is unreachable from well-formed evidence. Its comment says so and keeps the bookkeeping; leave it.
- **`step_beta_shapes`** refuses non-finite slopes and a non-positive determinant (`fit.rs:5270-5280`). `symmetric_eigen`'s `expect("no NaN")` cannot fire, because the shapes stay within `[0.02, 50]`.
- **`RuleSlopes::of`** finite-differences the quadrature rule at a relative step of 1e-4 (`information.rs:93`). It is documented and tested against a whole-likelihood difference.
- **The two golden sections in the SNP/indel area** (`fit.rs:5214` and `contamination.rs:1303`) return different points: the bracket's middle, against the better interior point. Merging them changes output, and the Newton box relies on the first's behaviour.
- **The STR side flooring before a logarithm at 1e-300** (`expected_counts.rs:425`, `:1463`) where the SNP/indel side uses `MIN_POSITIVE` is harmless.
- **Slippage groups threaded through the STR types while production fits one** (`penalised_fit.rs:1173-1180`): leave it, because the parameters file and the calling step read the per-group vectors.
- **One sample:**
  - SNP/indel: the excess is held and written `Defaulted` (`fit.rs:5119-5129`, `:2559-2565`);
  - STR: the between-locus share is reported absent;
  - nothing divides by a zero count on either side.
- **`fit_trace.rs`** holds process-wide state, but it is a diagnostic behind an environment variable and changes nothing it observes.
