# Fit precision, step A7 — the allele-frequency shapes follow the likelihood's own slope

**Date:** 2026-09-28. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step A7.
**Source:** checkpoint A's decision 2. **Branch:** `fit-precision`. **Review:**
[fit_precision_a7_2026-09-28.md](../reviews/fit_precision_a7_2026-09-28.md); **fixes:**
[fixes_applied_fit_precision_a7_2026-09-28.md](../reviews/fixes_applied_fit_precision_a7_2026-09-28.md).

## 1. What was built

The SNP/indel fit's update of the allele-frequency density's two Beta shapes, and of the duplicated
class's carrier Beta, no longer solves the digamma form. Checkpoint A measured why that mattered: the
form's integral on the quadrature rule misses the slope of the likelihood the rule computes, so the
fit came to rest beside the likelihood's maximum — 0.65 to 2.3 standard errors from its flat point,
and the shapes and the invariant share 0.4 to 3.1 errors above the truth, the gap growing as the
square root of the positions ([A4 report](fit_precision_a4_2026-09-27.md) §2.3).

**The update now takes one step along the log-likelihood's own slope in the two shapes**
(`step_beta_shapes`), so it is at rest only where that slope is zero:

- **the slope** is the likelihood's, taken on the rule the likelihood is computed on: each pass
  accumulates, at every node of the rule, the node's posterior times how its weight moves with the
  shape, plus how the node's product over samples moves along the frequency axis times how the node
  moves (`Statistics::density_shape_slopes`, `carrier_shape_slopes`) — the same slope step A1's scorer
  checks against finite differences, now summed on every pass;
- **the step** divides that slope by the curvature the Beta would have if each position's frequency
  were seen — the posterior count of positions the Beta covers times the trigamma matrix the old
  update used. That is one step of the EM gradient algorithm (Lange 1995);
- **at a bound the step is projected.** A shape whose step would leave its interval (0.02 to 50) is
  held at the bound it crosses, and the other takes the step the same curvature gives with that one
  held, so the free shape also comes to rest only where its own slope is zero. The first version
  clamped each shape of the joint step on its own, and the review found that it left the free shape
  where its slope was not zero (§3.2).

The digamma form's statistic, Σ posterior · ln f, is now kept in test builds only, for the test that
measures how far that form lands from the likelihood's slope; `fit_beta_shapes` is gone. The carrier
Beta's rule and its slopes travel together (`CarrierRule`), so a pass cannot hold one without the other.

## 2. What moves

**Fitted numbers move**, as the plan says this step must.

**The oracle** — four tomato accessions at about 3 reads, one library each
(`tmp/fit_precision/oracle_a7b/`, against the pre-step `oracle_a6b/`):

- **checksums:** `fitted.parameters.toml` `c590689d…` → `393e7cd356811ddab00e033a42903384`,
  `from_fit.comparable` `5cdeafd0…` → `1f16dec86c23eff1286c736d634f0a50`; the five other checksums
  unchanged;
- **the fit ends 30 log-likelihood units higher**: its best start at −7.991302 × 10⁵ against
  −7.991599 × 10⁵. **Every start stops at the 200-pass limit, before and after** (198 to 200 passes);
- **the density:** the first shape moved from 19.5 to 15.4 and now has an error (± 36.5); the second
  is at its bound, 50, both before and after. The carrier's shapes moved from 0.99 and 3.79 to 0.12
  and 1.92, the duplicated share from 6.8 × 10⁻⁴ to 2.4 × 10⁻³;
- **the parameters file:** 62 lines move — all four inbreeding coefficients, the largest from 0.952
  to 0.938; the ordinary-site prior's reference concentration from 1.662 to 1.850 and its alternative
  total from 7.56 × 10⁻⁴ to 8.07 × 10⁻⁴; the four read groups' error multipliers in their fourth
  figure; and repeat-tract rows, which read the inbreeding coefficients;
- **the calls:** 6,693 records become 6,706 (15 new, 2 gone); of the 6,691 shared records, 30
  change some sample's genotype.

Against the clamped version of this step: 0.4 units higher, the first shape 15.4 against 15.6, and
2 of the 6,706 records called differently. On the oracle the bound fix changes little
because the first shape's error, ± 36.5, is larger than the distance the clamp left.

**The cross-platform checksums** (`src/cli/cross_platform_digests.rs`), re-recorded in this commit with
the explanation in their doc: on the two-sample fixture two lines of the parameters file move — the
ordinary-site prior's reference concentration from 3.9591 to 39.089, its alternative total from
7.30 × 10⁻¹⁰ to 6.60 × 10⁻³ — read off a density whose shapes moved from a = 10.8, b = 47.9 to a = 0.288,
b = 50 (its bound); the fit's best log-likelihood rises from −44.919 to −44.704. The fit reports neither
shape as identified. In the calls, three records' QUAL rise — the SNPs at chrV:121 (480.2 to 541.9) and
chrV:456 (481.0 to 542.7), each with the other sample's GQ from 65 to 73, and the repeat tract at
chrV:201 (0.0 to 7.4). No genotype moves. Checksums `70e63484…` → `9d421df0bc755522f8a7b5617b4c30f8`
and `dad5ff61…` → `7a570c442a4060fcd62f4af85b50a4fa`.

## 3. What was measured

### 3.1 The fit stops at the likelihood's flat point

A4's coverage test, 200 cohorts a regime (`tmp/fit_precision/a7b_coverage_full.log`, 2,211 s). Every
line is the same as the clamped version's, since no drawn shape reaches a bound:

| regime | flat point, shapes and invariant share (errors) | before (A4) | within one / two errors, shapes and invariant share | before (A4) |
|---|---|---|---|---|
| 20 samples, 3 reads | −0.001 to −0.002 | −0.86 to −1.35 | 0.705 to 0.735 / 0.970 to 0.975 | 0.270 to 0.465 / 0.835 to 0.930 |
| 20 samples, 30 reads | −0.001 to −0.002 | −0.65 to −1.29 | 0.700 to 0.705 / 0.965 to 0.980 | 0.325 to 0.550 / 0.860 to 0.890 |
| 4 samples, 30 reads | −0.037 to −0.041 | −2.0 to −2.3 | 0.670 to 0.695 / 0.950 to 0.970 | 0.000 / 0.055 to 0.605 |
| 4 samples, 3 reads | −0.100 to −0.108 | −1.8 to −1.9 | 0.690 to 0.810 / 0.950 to 0.985 | 0.000 to 0.040 / 0.110 to 0.995 |
| 2 samples, 3 reads | −0.27 to −0.28 | −0.82 to −0.83 | 0.875 to 1.000 / 0.955 to 1.000 | 0.095 to 1.000 / 0.825 to 1.000 |

At 20 samples and at 4 samples and 30 reads **every kind** of parameter lands between 0.62 and 0.74
within one error and 0.925 and 0.98 within two. At 2 samples 146 of 200 fits converge within the pass
limit, against 166 before: the step converges more slowly (§3.3).

### 3.2 At a bound

The correctness reviewer's cohorts drawn with the density's second shape above its bound, re-run on
this code (`tmp/fit_precision/a7b_bounds.log`, three seeds each):

| cohort (truth) | first shape: clamped → projected | its slope at rest, in errors | log-likelihood gained over the clamped step | over the parent |
|---|---|---|---|---|
| 20 samples × 3 reads (15, 200) | 6.87–6.95 → 3.44–3.47 | +0.007 | 341 to 355 | 192 to 212 |
| 20 × 8 reads (2, 120) | 3.71–3.78 → 0.46–0.47 | −0.03 | 331 to 372 | 283 to 309 |
| 4 × 3 reads (15, 200) | 6.2–6.3 → 2.4–2.5 (two seeds; the third ends inside) | −0.05 to −0.07 | 28 to 31 | 35 to 36 |

The test `the_shapes_rest_where_the_likelihoods_slope_is_zero_at_a_bound_too` holds a 20-sample cohort
drawn at b = 150 to it: the step still left in the first shape is 2.7 × 10⁻⁵ at a = 0.1146.

### 3.3 Passes to convergence

Parent (`41e0f9bf`) against this step, the same cohorts and starts, three starts each
(`tmp/fit_precision/a7b_passes.log`, the correctness reviewer's probe on this code, against its
`passes_parent_full.log`):

| cohort | accelerated, as the fit runs | plain alternation |
|---|---|---|
| 20 samples × 3 reads, 3 seeds | 12–18 → 18–30 | 18–27 → 21–84 |
| 4 × 3 reads, seed 1 | 39–55 → 63–70 | 171–174 → 315–324 |
| 20 × 8 reads, duplicated class, seed 0 | 33–39 → 39–57 | 120–123 → 135–138 |
| 20 × 8, duplicated, seeds 1 and 2 (carrier's second shape at its bound) | 32–78 → 28–47 | 75–372 → 69–132 |
| 4 × 3, duplicated, seed 1 | converged at 165–178 → the 198-pass limit | 678–687 → the 1,500-pass limit |
| 4 × 3, duplicated, seeds 0 and 2 | the 200-pass limit → converged at 167–197 on 5 of 6 starts | the limit → 1,113–1,389 |

**Where the shapes are inside their bounds the step takes more passes than the update it replaced**
— 1.3 to 2 times with acceleration at 20 samples — because it takes one step a pass where the old
update solved its equations to convergence. Where a shape sits on its bound it takes fewer, because the
old update suffered the clamping defect too. With the 4 to 7% more a pass measured below, a converged
20-sample fit with interior shapes costs about 1.35 to 2.1 times what it did.

### 3.4 No plain pass lowers the log-likelihood

Measured with a probe (plain passes, no acceleration, from each of the three starts;
`tmp/fit_precision/a7b_mono.log`): on eight cohorts with interior shapes, 1 to 63 samples at 3 to 30
reads, with and without the duplicated class, 80 passes each, **none of 1,920 passes** lowers it (the
digamma update lowered it at 357, by up to 0.0095 units, `a7_mono_old.log`); on nine cohorts with a
shape at its bound, 300 passes each, 161 of 8,100 passes lower it, **by at most 2.9 × 10⁻¹¹ units** —
the rounding of a log-likelihood of 10⁴ to 10⁵ that has stopped moving. **This is measured, not
guaranteed**: the step is not a full maximisation of a surrogate, so no inequality forbids a fall. The
test `plain_passes_never_lower_the_log_likelihood` holds three cohorts to it, one at a bound.

### 3.5 The cost of a pass

4 to 7% more (single thread, median of 11 passes, two interleaved rounds,
`tmp/fit_precision/a7_timing.log`, measured before the bound fix, which adds no work to a pass):
135.0–135.7 → 144.7–145.3 ms at 8 samples × 40,000 positions with the duplicated class, 105.7–106.8 →
111.6 ms at 32 × 10,000, 85.9–86.0 → 89.3–90.4 ms at 8 × 40,000 without it.

## 4. Deviations and assumptions

1. **The curvature is the complete-data one** (the trigamma matrix times the posterior count), not the
   observed information: the plan says "Newton steps on the rule's slopes" and leaves the curvature
   open. It makes the update the EM gradient algorithm, whose fixed point is the likelihood's flat
   point, and it costs nothing a pass; its price is the pass counts of §3.3.
2. **"Never lowering the log-likelihood" is measured, not enforced** (§3.4). Enforcing it would need a
   pass at the candidate shapes before accepting them. The acceleration still refuses a jump that
   loses more than one unit. **For the owner to rule on at checkpoint A′.**
3. **One step a pass**, where the old update solved its equations to convergence each pass; the
   acceleration carries the rest.
4. **The projected step at a bound** (§1) is an addition the plan's "kept within the shapes' bounds"
   did not spell out; clamping each shape was the first reading, and the review showed it wrong.
5. **A test's fixture changed**: `the_trace_files_the_returned_fit_after_the_winning_starts_last_pass`
   now reads which start won from the trace's rows instead of asserting it, since the winner moves with
   any change to the fit.

## 5. Changes

- [fit.rs](../../../../src/parameter_estimation/joint/fit.rs): `step_beta_shapes` (projected at a
  bound) replaces `fit_beta_shapes`; `Statistics` gains `density_shape_slopes` and
  `carrier_shape_slopes` and keeps `sum_ln_f` in test builds only; `BetaQuadrature` gains
  `prior_slopes` (and keeps `ln_nodes` in test builds only); `PassModel` gains the density rule's
  slopes, and `CarrierRule` holds the carrier's rule with its slopes; `one_position` sums each node's
  share of the slopes; `digamma` is test-only; the module doc says the step is measured, not
  guaranteed, not to lower the likelihood.
- [fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs): `RuleSlopes`'
  two tables visible to the pass; `ScoringTables` reuses the pass's slopes; two tests of the pass's
  slopes; the coverage test's checks re-recorded (§3.1); the trace test reads its winner.
- [cli/cross_platform_digests.rs](../../../../src/cli/cross_platform_digests.rs): both checksums
  re-recorded, with the explanation.

## 6. Tests

| test | what it shows |
|---|---|
| `the_shapes_step_follows_the_slope_by_the_complete_data_curvature` (fit.rs), replacing `the_beta_shapes_come_back_from_their_own_log_means` | A slope built as the curvature times a known move lands on it to 10⁻¹²; one step worked by hand at (1, 1) to 10⁻⁶; a zero slope, a NaN slope and a vanishing count move nothing; a slope towards (100, 100) stops at both bounds; **with the second shape at its bound and pushing out, the first takes its own step**, `slope / curvature`. |
| `the_shapes_rest_where_the_likelihoods_slope_is_zero_at_a_bound_too` (fit.rs) | A 20-sample cohort drawn at b = 150: the fit returns b at 50, and the first shape where its slope asks for a step below a hundredth of it (measured 2.7 × 10⁻⁵ at 0.1146). |
| `plain_passes_never_lower_the_log_likelihood` (fit.rs) | Plain passes from each start on three drawn cohorts (20 samples at 3 reads; 8 at 8 reads with duplications; 4 at 5 reads drawn at b = 150, 150 passes): the log-likelihood never falls. The 20-sample cohort fails with the digamma update, the bound cohort with the clamped step. |
| `the_passs_shape_slopes_are_the_likelihoods` (information.rs) | The pass's summed density and carrier slopes equal the scorer's within 10⁻⁷ over a cohort of more than two chunks, with and without the duplicated class (measured at most 1.8 × 10⁻⁹). |
| `the_passs_carrier_slopes_follow_the_coverage_odds` (information.rs) | The same for the carrier with each sample's coverage odds set. |
| `the_errors_mean_what_they_say` (information.rs, ignored) | Re-recorded (§3.1): every kind covers at 20 samples and at 4 samples and 30 reads; flat points within 0.02 errors at 20 samples and 0.25 at 4; shapes within 0.25 errors of the truth at 20 samples. |

Going back to clamping each shape fails three of these (`tmp/fit_precision/a7fix_mutations.log`).
`cargo test --release --lib parameter_estimation::joint::fit`: 72 passed, 1 ignored.

## 7. Validation

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- `cargo test --all-targets --all-features --no-fail-fast` (`tmp/fit_precision/suite_a7b.log`):
  4,930 passed, 3 failed, 5 ignored; the three failures are the pre-existing ones
  (`examples/ng_generic_loci_dump.rs` × 2, `examples/ng_ssr_loci_dump.rs` × 1); both
  `cli::cross_platform_digests` tests pass with the re-recorded checksums.
- The coverage run at its full count passed every re-recorded check (§3.1).
- The oracle: §2. Its exit 1 is the stale baseline's two fit lines, as since step A1 — which now
  differ from main's by this step too.
