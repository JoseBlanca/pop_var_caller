# Fit precision, step C2 — whether a stratum's errors mean what they say

**Date:** 2026-09-30. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step C2.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §4.5 item 1.
**Branch:** `fit-precision`. **Review:** [fit_precision_c2_2026-09-30.md](../reviews/fit_precision_c2_2026-09-30.md);
**fixes:** [fixes_applied_fit_precision_c2_2026-09-30.md](../reviews/fixes_applied_fit_precision_c2_2026-09-30.md).

## 1. What was built

A test, ignored by default, that asks whether step C1's standard errors describe how far a stratum's fitted
numbers actually land from the truth: `a_stratums_errors_mean_what_they_say` in `ssr_fit.rs`.

- It draws strata from known slippage (level 0.05, shorter share 0.8, fall-off 0.3), a known length spectrum
  (`spectrum_of`) and concentration (0.5), fits each by the fit's own rule, and for every estimate with an error
  records its distance from the truth in that error.
- For each kind of number it prints how many estimates land within one error and within two — about 68 and 95
  in 100 for an error that means what it says — the mean distance (a bias, in errors), the spread (about one),
  and how many came back without an error, by why.
- Three regimes, 20 samples each: three allele classes over 300 tracts at 3 reads a sample and at 30 (100
  draws each by default), and the production span, thirteen classes, over 200 tracts at 3 reads (25 draws).
  `NG_FIT_PRECISION_STRATUM_DRAWS` sets the count; 20 to 40 minutes in the container at the default (three runs).
- The three-class regimes are held to bounds taken from the run below with a margin, at the default count; the
  thirteen-class regime is printed and not held to anything (§3).

No production code changes.

## 2. What was measured

`tmp/fit_precision/c2/coverage100_keep.log`, 100 draws a three-class regime, 25 at thirteen:

| regime | number | within one error | within two | mean distance | spread |
|---|---|---|---|---|---|
| 3 classes, 3 reads | slippage level | 0.70 | 0.98 | −0.11 | 1.00 |
| | shorter share | 0.67 | 0.94 | +0.22 | 1.06 |
| | fall-off | 0.67 | 0.96 | +0.08 | 0.95 |
| | concentration | 0.67 | 0.95 | −0.03 | 1.06 |
| | reference-length share | 0.59 | 0.91 | +0.11 | 1.16 |
| | other length shares | 0.61 | 0.90 | −0.13 | 1.34 |
| 3 classes, 30 reads | slippage level | 0.69 | 0.99 | +0.01 | 0.94 |
| | shorter share | 0.61 | 0.91 | +0.23 | 1.10 |
| | fall-off | 0.76 | 0.94 | +0.10 | 0.96 |
| | concentration | 0.65 | 0.94 | −0.02 | 1.05 |
| | reference-length share | 0.47 | 0.79 | +0.24 | 1.61 |
| | other length shares | 0.53 | 0.85 | −0.23 | 1.61 |
| 13 classes, 3 reads | slippage level | 0.32 | 0.72 | +1.42 | 1.17 |
| | shorter share | 0.64 | 0.80 | −0.63 | 1.21 |
| | fall-off | 0.64 | 0.88 | −0.25 | 1.42 |
| | concentration | 0.00 | 0.04 | +8.45 | 3.04 |
| | reference-length share | 0.04 | 0.04 | −18.41 | 10.18 |
| | other length shares | 0.12 | 0.22 | +1.24 | 11.73 |

- **At three classes the errors mean what they say** for the slippage numbers and the concentration: 61 to 76
  in 100 within one error, 91 to 99 within two, each within a quarter of an error of the truth on average. With
  100 draws a share within one error is itself uncertain by about ±0.05.
- **The class shares' errors are too small at 30 reads**: 47 and 53 in 100 within one, the distances spread 1.6
  times wider than the errors say. At 3 reads they are close (59 and 61).
- No number came back without an error at three classes; 6 of the 300 other-class shares at thirteen were not
  identified.
- Of the climbs, 55 in 100 (3 reads) and 60 in 100 (30 reads) settled within their five rounds at three classes,
  2 in 25 at thirteen — the rule step C3 replaces.

## 3. At thirteen classes the fit's answer is off, because its integral is too coarse

**The thirteen-class numbers are not, in the main, a failure of the errors: the fitted numbers are far from the
truth.** The reference-length share sits a mean of 18 errors from it, the concentration 8, the slippage level 1.4
(C1 saw the level 3.0 to 4.5 errors off on three strata).

**What the fit integrates, in words.** At each tract the fit does not know how the stratum's alleles are shared
out among its chromosomes; it averages the tract's likelihood over every way they could be, weighted by the
stratum's length spectrum and concentration. That average is computed on a fixed set of 256 sample points
(`QUADRATURE_POINTS`, `parameter_prepass_joint_fit.md` §4.2), checked against an exact grid at three allele classes
only. Production fits thirteen.

**It is not the climb stopping short.** Refitting two drawn strata at 20 and 60 rounds, where the default is 5,
both settled by 20 and moved no number at 60; the reference-length share stayed 0.175 and 0.171 against 0.296
(`tmp/fit_precision/c2/rounds.log`).

**It is the 256-point average.** The correctness reviewer evaluated one drawn 13-class stratum's log-likelihood a
tract at the truth and at the 256-point fit's answer, with the average taken over more and more points
(`tmp/review_2026-09-30_fit_precision_c2/evidence_correctness/review/probeA.log`):

| points in the average | at the truth | at the 256-point fit's answer | fit's answer minus truth |
|---|---|---|---|
| 256 | −40.811 | −38.469 | +2.342 |
| 1,024 | −37.841 | −37.911 | −0.071 |
| 4,096 | −37.571 | −37.568 | +0.003 |
| 16,384 | −37.250 | −37.444 | −0.194 |
| 65,536 | −37.122 | −37.401 | −0.280 |

- At 256 points the fit's answer looks 2.34 units a tract better than the truth (468 over the stratum's 200
  tracts); with 65,536 points it is 0.28 a tract **worse** than the truth. The fitted maximum is made by the coarse
  average, not by the reads.
- The 256-point average is 3.7 units a tract too low at the truth. Four other placements of the 256 points gave
  the same pattern.
- **No point count tried is known to be enough**: 4,096 points still leave the average 0.45 a tract low at the
  truth (90 units over the stratum). Refitting at 1,024 and 4,096 points moved the reference share to 0.229 and
  0.309 (truth 0.296), but the 1,024-point climb did not settle in 20 rounds, and the 4,096-point one landing near
  the truth is partly luck. A fit at 4,096 points took 468 s against 48 s at 256.
- **It starts by five classes**: over 12 draws at five classes the other length shares landed within one error 12
  times in 100 and the concentration sat 1.95 errors off, the slippage numbers unaffected.
- **It reaches the three-class shares too, in part**: refit at 4,096 points, three drawn three-class strata's
  reference-share errors grew 1.29 to 1.56 times — the under-coverage of §2 at 30 reads.

## 4. For the owner

The coarse average is the fit's model, which this plan does not change (spec §1.3) — and raising the point count
moves every stratum and reverses the decision recorded in `parameter_prepass_joint_fit.md` §4.2. But it is not
outside this plan in effect: the stratum errors C1 computes, C3's check against a 20-round climb, D's precision
target on the level, and the errors E2 writes all stand on this likelihood.

**Recommendation: a separate investigation of the average before step D, and continue with C3 now.**

- **C3 is independent of it.** C3 changes when the climb stops, and is validated against a longer climb of the
  same objective, whatever that objective's error.
- **D and E2 are not.** D grows each stratum's sample subset until the level is known to 2%, and the level's error
  at thirteen classes is not yet trustworthy (32 in 100 within one error); E2 would write these errors into the
  parameters file.
- **The investigation**: measure the average's error at the truth against the point count, the number of classes
  and the concentration, and weigh more points (4,096 cost about ten times the time) against placing them where a
  tract's reads say the frequencies are.

## 5. Tests

- `a_stratums_errors_mean_what_they_say` (ignored): the table of §2. At the default 100 draws the three-class
  regimes are held to bounds from these runs with a margin (`hold_to_the_three_class_bounds`): every draw fitted,
  every estimate with an error; the slippage numbers and the concentration within one error 55 to 82 times in 100,
  within two at least 88, spread 0.85 to 1.2, mean distance under 0.4 errors; the class shares within one at least
  40 in 100 and within two 75. Those bounds cannot see every error scaled by 0.88 to 1.18 (review measurement), and
  the shares' accept their under-coverage. Run at 100 draws with them in place (`coverage100_final.log`, 1,211 s):
  passed, every figure the same as §2's.
- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- Full suite (`cargo test --all-targets --all-features --no-fail-fast`): 4,992 passed, 3 failed (the pre-existing
  `examples/ng_generic_loci_dump.rs` × 2 and `examples/ng_ssr_loci_dump.rs` × 1), 6 ignored (the new test is the
  sixth).
