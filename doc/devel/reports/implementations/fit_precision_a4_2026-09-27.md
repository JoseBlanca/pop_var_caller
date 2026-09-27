# Fit precision, step A4 — the errors mean what they say

**Date:** 2026-09-27. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step A4.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §3.6 item 3. **Branch:** `fit-precision`.
**Review:** [fit_precision_a4](../reviews/fit_precision_a4_2026-09-27.md); this report is the
version after its fixes.

## 1. What was built

A test that checks the SNP/indel fit's standard errors against what a standard error promises: if
the same kind of data were drawn again and again, about 68 estimates in 100 would land within one
error of the truth, 95 within two, and the distances from the truth, in errors, would scatter with a
standard deviation of about one. The test draws 200 cohorts from known values, fits each by the
fit's own rule, and tallies how far each estimate lands from the value it was drawn from — per kind
of parameter, for the errors the fit reports and for those of the full matrix.

- **The errors the fit reports** come from the information summed in blocks (step A2): each
  sample's parameters are paired with its own and with the cohort's, never with another sample's.
- **The full matrix's errors** pair every parameter with every other, the products between two
  samples' scores included. It is affordable only on small cohorts, which is why the fit reports
  the blocks.

Five regimes: 4 samples over 20,000 positions and 20 samples over 5,000, each at 3 and at 30 reads a
position, and 2 samples over 30,000 positions at 3 reads. The test is `the_errors_mean_what_they_say`
in [fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs), ignored by
default because it runs for 20 minutes in the container (1,184 s,
`tmp/fit_precision/a4_coverage4.log`). No code outside the tests changed; no fitted number moved.

Beside the coverage, the test prints per kind, all in units of the error:

- **off centre** — the mean of (estimate − truth) / error;
- **spread** — their standard deviation, about one for an error that is the right size;
- **the flat point** — how far the log-likelihood's flat point lies from where the fit stopped, by
  one Newton step (the full matrix's inverse times the log-likelihood's slopes summed over the
  positions). A fit that stops at the likelihood's maximum gives zero;
- **coverage at the flat point** — the estimate moved by that step, against the truth. It tells
  whether the errors describe the likelihood's maximum when the fit does not stop there.

**The test holds each regime to what it measured** (`check_the_coverage`), so a change to what it
measures fails it: every kind is tallied with an error from both matrices; at the full cohort count,
the kinds that cover stay in their bands, the full matrix covers more than the blocks at 4 samples
and 3 reads, and the Newton step carries the density shapes to the truth at 20 samples.
`NG_FIT_PRECISION_COVERAGE_COHORTS=N` shortens a run (it reaches the container by its `NG_` prefix).

## 2. Results

Converged: 166 of 200 fits at 2 samples, 199 of 200 at 4 samples and 3 reads, 200 of 200 in the
other three regimes. Each cohort-level parameter has 200 estimates a regime; each sample-level kind
pools every sample (400, 800 or 4,000). At 2 samples a few parameters have no error (14 mismapped
shares, 2 mismapped rates, one each of three others); everywhere else none is missing.

"Within one", "within two" and "spread" give the blocks' value first and the full matrix's second.
"Off centre" is in the blocks' error; the flat point and the columns after it in the full matrix's.

| samples | reads | parameter | within one | within two | off centre | spread | flat point | at the flat point: within one · two | its off centre |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 3 | clean error rates | 0.537 · 0.800 | 0.848 · 0.917 | +0.193 | 1.398 · 0.958 | −0.371 | 0.640 · 0.880 | −0.007 |
| 2 | 3 | mismapped error rates | 0.643 · 0.872 | 0.869 · 0.992 | −0.271 | 1.453 · 0.674 | −0.344 | 0.704 · 0.940 | −0.225 |
| 2 | 3 | homozygote excess | 0.627 · 0.625 | 0.932 · 0.927 | +0.181 | 1.098 · 1.110 | −0.028 | 0.615 · 0.925 | +0.147 |
| 2 | 3 | mismapped share | 0.543 · 0.817 | 0.833 · 0.914 | −0.201 | 1.772 · 1.476 | +0.385 | 0.683 · 0.876 | −0.124 |
| 2 | 3 | fixed share | 0.593 · 0.593 | 0.975 · 0.975 | +0.849 | 0.560 · 0.560 | −0.728 | 0.578 · 0.935 | +0.121 |
| 2 | 3 | invariant share | 0.095 · 0.090 | 0.825 · 0.825 | +1.620 | 0.716 · 0.986 | −0.815 | 0.530 · 0.805 | +0.854 |
| 2 | 3 | density shape a | 1.000 · 1.000 | 1.000 · 1.000 | +0.421 | 0.037 · 0.037 | −0.830 | 0.709 · 0.940 | −0.409 |
| 2 | 3 | density shape b | 1.000 · 1.000 | 1.000 · 1.000 | +0.413 | 0.043 · 0.043 | −0.828 | 0.698 · 0.945 | −0.415 |
| 4 | 3 | clean error rates | 0.595 · 0.743 | 0.895 · 0.973 | −0.068 | 1.234 · 0.893 | −0.098 | 0.721 · 0.978 | −0.055 |
| 4 | 3 | mismapped error rates | 0.589 · 0.739 | 0.874 · 0.958 | −0.386 | 1.464 · 0.938 | −0.018 | 0.729 · 0.954 | −0.188 |
| 4 | 3 | homozygote excess | 0.676 · 0.677 | 0.950 · 0.950 | +0.120 | 1.018 · 1.017 | −0.142 | 0.679 · 0.951 | −0.023 |
| 4 | 3 | mismapped share | 0.460 · 0.780 | 0.775 · 0.950 | +0.114 | 1.667 · 1.022 | +0.012 | 0.745 · 0.950 | −0.151 |
| 4 | 3 | fixed share | 0.605 · 0.605 | 0.905 · 0.905 | +0.732 | 0.973 · 0.973 | −0.606 | 0.645 · 0.925 | +0.127 |
| 4 | 3 | invariant share | 0.000 · 0.000 | 0.110 · 0.110 | +2.869 | 0.668 · 0.670 | −1.765 | 0.410 · 0.780 | +1.112 |
| 4 | 3 | density shape a | 0.000 · 0.000 | 0.995 · 0.995 | +1.590 | 0.161 · 0.162 | −1.939 | 0.760 · 0.975 | −0.346 |
| 4 | 3 | density shape b | 0.040 · 0.040 | 0.955 · 0.955 | +1.533 | 0.289 · 0.290 | −1.875 | 0.775 · 0.970 | −0.339 |
| 4 | 30 | clean error rates | 0.690 · 0.690 | 0.951 · 0.951 | −0.062 | 0.992 · 0.991 | +0.001 | 0.691 · 0.951 | −0.062 |
| 4 | 30 | mismapped error rates | 0.655 · 0.661 | 0.944 · 0.948 | +0.018 | 1.047 · 1.037 | −0.000 | 0.656 · 0.949 | +0.017 |
| 4 | 30 | homozygote excess | 0.681 · 0.682 | 0.963 · 0.963 | +0.180 | 0.956 · 0.954 | −0.127 | 0.682 · 0.966 | +0.053 |
| 4 | 30 | mismapped share | 0.665 · 0.665 | 0.950 · 0.950 | −0.116 | 1.007 · 1.004 | +0.000 | 0.655 · 0.950 | −0.116 |
| 4 | 30 | fixed share | 0.650 · 0.650 | 0.945 · 0.945 | +0.338 | 0.982 · 0.982 | −0.443 | 0.700 · 0.930 | −0.105 |
| 4 | 30 | invariant share | 0.000 · 0.000 | 0.055 · 0.055 | +3.085 | 0.738 · 0.738 | −2.010 | 0.460 · 0.740 | +1.075 |
| 4 | 30 | density shape a | 0.000 · 0.000 | 0.530 · 0.530 | +2.001 | 0.235 · 0.235 | −2.303 | 0.705 · 0.970 | −0.302 |
| 4 | 30 | density shape b | 0.000 · 0.000 | 0.605 · 0.605 | +1.910 | 0.337 · 0.337 | −2.192 | 0.740 · 0.975 | −0.283 |
| 20 | 3 | clean error rates | 0.683 · 0.686 | 0.945 · 0.945 | −0.061 | 1.038 · 1.030 | −0.019 | 0.688 · 0.945 | −0.079 |
| 20 | 3 | mismapped error rates | 0.677 · 0.728 | 0.938 · 0.954 | −0.151 | 1.067 · 0.965 | +0.006 | 0.727 · 0.954 | −0.134 |
| 20 | 3 | homozygote excess | 0.697 · 0.714 | 0.951 · 0.958 | +0.077 | 1.002 · 0.973 | −0.028 | 0.712 · 0.956 | +0.047 |
| 20 | 3 | mismapped share | 0.610 · 0.655 | 0.950 · 0.965 | +0.017 | 1.099 · 1.001 | −0.035 | 0.655 · 0.970 | −0.025 |
| 20 | 3 | fixed share | 0.655 · 0.655 | 0.960 · 0.960 | +0.083 | 0.998 · 0.998 | −0.028 | 0.660 · 0.955 | +0.055 |
| 20 | 3 | invariant share | 0.465 · 0.465 | 0.835 · 0.835 | +1.126 | 0.883 · 0.883 | −0.861 | 0.645 · 0.910 | +0.266 |
| 20 | 3 | density shape a | 0.270 · 0.270 | 0.930 · 0.930 | +1.290 | 0.485 · 0.485 | −1.346 | 0.740 · 0.970 | −0.056 |
| 20 | 3 | density shape b | 0.430 · 0.430 | 0.915 · 0.920 | +1.074 | 0.667 · 0.667 | −1.130 | 0.715 · 0.985 | −0.056 |
| 20 | 30 | clean error rates | 0.681 · 0.683 | 0.951 · 0.952 | −0.026 | 1.009 · 1.006 | +0.000 | 0.684 · 0.952 | −0.026 |
| 20 | 30 | mismapped error rates | 0.690 · 0.722 | 0.955 · 0.969 | −0.039 | 0.996 · 0.929 | +0.000 | 0.721 · 0.968 | −0.036 |
| 20 | 30 | homozygote excess | 0.682 · 0.695 | 0.955 · 0.962 | +0.031 | 1.006 · 0.978 | −0.018 | 0.697 · 0.961 | +0.012 |
| 20 | 30 | mismapped share | 0.670 · 0.670 | 0.965 · 0.965 | −0.130 | 0.999 · 0.999 | −0.003 | 0.670 · 0.965 | −0.133 |
| 20 | 30 | fixed share | 0.665 · 0.665 | 0.950 · 0.950 | −0.008 | 1.030 · 1.030 | −0.008 | 0.665 · 0.950 | −0.016 |
| 20 | 30 | invariant share | 0.550 · 0.550 | 0.880 · 0.880 | +0.815 | 0.960 · 0.960 | −0.653 | 0.645 · 0.930 | +0.161 |
| 20 | 30 | density shape a | 0.325 · 0.325 | 0.890 · 0.890 | +1.280 | 0.527 · 0.527 | −1.290 | 0.745 · 0.985 | −0.011 |
| 20 | 30 | density shape b | 0.500 · 0.500 | 0.860 · 0.860 | +1.054 | 0.741 · 0.741 | −1.053 | 0.725 · 0.980 | +0.001 |

If the estimates were independent, a share of 0.68 over 200 of them would have a sampling spread of
about ±0.033, and over 800 about ±0.016; at 0.95 about ±0.015 and ±0.008. The samples of one cohort
share its fit, so the pooled kinds are not quite independent.

### 2.1 Where the errors mean what they say

**The error rates, the homozygote excesses and the mismapped share, and the fixed share at 20
samples**, at 20 samples at both depths and at 4 samples at 30 reads: between 0.61 and 0.70 land
within one of the blocks' errors and between 0.938 and 0.965 within two, they scatter with a spread
of 0.96 to 1.10, and they sit within 0.18 of an error of the truth on average. The homozygote excess
covers in every regime of 4 samples or more (spread 0.96 to 1.02).

### 2.2 Where the blocks' errors are too small: few samples at low depth

**At 4 samples and 3 reads the errors the fit reports are too small** for the error rates and the
mismapped share: their estimates scatter 1.23 to 1.67 times as widely as the blocks' errors say, and
0.46 to 0.60 of them land within one error. **The full matrix's errors are close to the right size
there, and on the wide side**: a spread of 0.89 for the clean rates (errors about 12% too wide,
which is what puts 0.74 of them within one error), 0.94 for the mismapped rates and 1.02 for the
mismapped share, and 0.95 to 0.97 within two. The mismapped share's 0.78 within one error with a
spread of 1.02 means its distances are more peaked than a normal distribution's, not that its
errors are wide.

At 2 samples and 3 reads it is the same, further: the blocks' spreads are 1.40 to 1.77 for these
kinds, the full matrix's 0.67 to 1.48, and 34 of the 200 fits stopped at the pass limit. **At 4
samples and 30 reads both matrices hold.** At 20 samples the blocks are a little small for the
mismapped kinds (spread 1.07 and 1.10 at 3 reads) where the full matrix is a little wide (0.97 and
1.00), and both are close.

This is what step A3's review predicted: two samples at one position share its unknown allele
frequency and class, and with few samples that shared part is a real share of each sample's
information, which the blocks leave out (A3 report §4). The homozygote excess is not affected (spread
1.02 against 1.02 at 4 samples and 3 reads).

### 2.3 Where the fit stops beside the likelihood's maximum: the density shapes

**The density's two shapes and the invariant share sit above the truth in every regime**, by 0.4 to
3.1 errors on average; none of the 200 shape `a`s lands within one error of the truth at 4 samples,
at either depth. At 2 and 4 samples the fixed share sits above it too, by 0.34 to 0.85. **The fitted
shapes also scatter far less than their errors say**: a spread of 0.04 at 2 samples, 0.16 to 0.34 at
4, 0.49 to 0.74 at 20.

**The fit's update for the shapes stops where the likelihood still slopes.** The log-likelihood's
flat point lies below where the fit stopped by about as much as the estimates sit above the truth:
0.65 to 2.3 errors for the shapes and the invariant share, 0.44 to 0.73 for the fixed share. The
review established the cause directly, on six cohorts of 20 samples at 3 reads:

- **the fit is a fixed point of its own update, not a fit that stopped early**: 300 further plain
  passes from the fit move shape `a` by less than 10⁻⁵, and with the fit's stop thresholds made 10⁴
  times tighter (375 to 987 passes a start instead of 48 to 117; 4 samples, 3 reads, 20 cohorts)
  shape `a` sits at +1.558 errors against +1.557;
- **at that point only one equation is out of balance**: the update solves the digamma form — the
  slope of the log-likelihood as if each position's frequency were known, averaged on the quadrature
  rule — which is zero there to within 0.012, while the likelihood's own slope in `a` is −8 to −20.
  The slopes in `b` and in the invariant share are zero to within 0.002. The digamma form uses
  `ln f`, which the rule integrates badly where a first shape below one piles mass near zero (step
  A1's finding 3);
- **the flat point really is higher**: a damped Newton search reaches it, gaining in log-likelihood
  what the outer-product matrix predicts to within about 5%.

**At the flat point the shapes' errors hold**: moved there by the Newton step, 0.70 to 0.78 of them
land within one error and 0.94 to 0.985 within two, in every regime. **The invariant share's do not
at 4 samples** (0.41 and 0.46 within one at the flat point, +1.1 errors off centre): there the
likelihood's own maximum sits above the truth — the review's search to the flat point found it +1.26
errors above on average over six cohorts, and the single Newton step lowered the likelihood in all
six, so the one-step columns are least reliable at 4 samples. Fixing the shapes' update would not
bring that share's coverage to 68 in 100 at few samples.

**The gap grows with the positions.** Measured by the review in `a`'s own units, the distance from
the fit to the flat point shrinks as samples are added (0.47 at 4 samples, 0.15 at 20, both at 3
reads) and does not shrink as positions are added (0.15 at 5,000 positions, 0.14 at 20,000, 20
samples). The error shrinks as one over the square root of the positions, so in errors the gap grows
as its square root: 1.30 errors at 5,000 positions, 2.47 at 20,000. Both the update's error and the
information are sums over positions, which is why. A real census has hundreds of thousands of
positions or more; on this law the shapes would sit tens of errors from the likelihood's maximum
there. That is an extrapolation from two sizes, not a measurement.

**More quadrature nodes are not the fix.** The likelihood itself is already exact at 16 nodes: the
review found it the same to four decimals from 16 to 192 nodes. Only the digamma form's error shrinks
with more nodes, as about N^−2a with N the node count — 4.0- to 5.4-fold from 16 to 48 nodes at
these cohorts' fitted `a`, and sevenfold at step A1's `a = 0.9`. Measured on the same 40 cohorts a
regime (`a4_probe16_40.log` at 16 nodes, `a4_probe48.log` at 48), the flat point in `a` moved from
−1.39 to −0.62 errors at 20 samples and 3 reads, from −1.29 to −0.51 at 30 reads, from −2.05 to −1.19
and from −2.25 to −1.16 at 4 samples, and the estimates themselves moved toward the truth (shape `a`
+1.26 to +0.51 and +1.57 to +0.94 errors off centre at 20 and 4 samples, 3 reads). At 96 nodes the
review still measured 0.16 to 0.88 errors on four cohorts, and that gap would grow with the positions
the same way. The fix is to maximise the rule's own likelihood in the two shapes — for example Newton
steps on the slopes step A1 already computes (`RuleSlopes`).

The invariant and fixed shares move with the shapes: plausibly, the fitted first shape is too large,
which moves mass away from frequency zero, and the invariant share rises to take it. Not measured.

**What it means for the stopping rule (plan step B1).** The settled test stops a fit when its moves
are small against its errors, and it will stop here: the update has stopped moving. It cannot see
that the point is not the maximum. B1 would then report a fit as settled at a point its own errors
place far from the answer.

## 3. Assumptions and deviations

1. **The duplicated class is off in the draw and in the fit.** The drawing generator floors each
   carrier frequency at 0.2 (every draw below it becomes exactly 0.2), which is outside the model's
   family, so no carrier shape would be true. Its errors were checked in step A3 against the full
   matrix instead.
2. **"Converged" is the fit's own rule, unchanged.** Fits that stopped at the pass limit are kept in
   the tally, since a run keeps what it has: 34 of 200 at 2 samples, 1 at 4 samples and 3 reads.
3. **The truth is the value drawn from, not the likelihood's maximum on each cohort.** So a displaced
   fit (§2.3) counts against coverage, which is what the owner reads the errors for: how far the
   written value is from the real one. The at-the-flat-point columns give the other reading.
4. **The Newton-step columns use the full matrix** (its inverse and its error), because the blocks
   have no products between two samples to step with.
5. **The regimes are not the A3 review's exactly.** It asked for two to four samples at 3 reads with
   30,000 positions or more; 2 samples run at 30,000 positions, 4 samples at 20,000, to keep the test
   near twenty minutes. 2 samples run at 3 reads only: at 4 samples both matrices already hold at 30
   reads, and the 2-sample regime alone costs about 11 minutes. One sample needs no measurement of the
   block approximation: with one sample the blocks are the full matrix.
6. **Every sample is drawn at the same values**, so a slip that read one sample's estimate against
   another's error would barely show. The review confirmed the indexing is right; drawing each sample
   at its own rates is deferred (it would redraw every cohort).

## 4. Changes

[fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs), tests only:

- `the_errors_mean_what_they_say` (ignored), `check_the_coverage`, `coverage_of`, `Coverage`,
  `Tally`, `TrueValues` and `TRUTH`, `value_at`, `label_of`, `COVERAGE_KINDS`, `COVERAGE_COHORTS`.
- `full_matrix`, lifted out of the step A3 comparison so both tests build the full matrix the same
  way, and `the_full_matrix_sums_the_products_between_two_samples`, its own test;
  `fitted_parameters_and_convergence`, which also says whether the fit converged
  (`fitted_parameters` now calls it, and is private to the tests again — nothing outside them used
  it).
- The step A3 comparison's documentation now states A4's answer.

## 5. Validation

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- `cargo test --all-targets --all-features --no-fail-fast` (`tmp/fit_precision/suite_a4b.log`): counts
  in the commit message; the three pre-existing failures only; both cross-platform checksum tests pass
  unchanged (nothing outside the tests changed).
- `the_errors_mean_what_they_say` at the full count: passed, every check held
  (`tmp/fit_precision/a4_coverage4.log`, 1,184 s).

## 6. For checkpoint A

- **Spec question 3 — is the outer-product estimator close enough?** Per kind:
  - *error rates, homozygote excesses, mismapped share, fixed share:* yes, for the full outer-product
    matrix, in every regime measured (spread 0.89 to 1.04 at 4 samples and more). **The block
    approximation is not, at few samples and low depth**: estimates scatter 1.23 to 1.67 times its
    errors at 4 samples and 3 reads, 1.40 to 1.77 at 2. The choice is between the full matrix for
    small cohorts, where it is cheap and at worst about 12% too wide, and the blocks' errors, 19 to
    40% too small at 4 samples and 3 reads and up to 44% at 2. At 20 samples the two spreads are
    about 10% apart.
  - *density shapes:* yes, at the likelihood's maximum (0.70 to 0.78 within one error there). But the
    fit does not stop at the maximum, so the errors it writes describe a point that is not the
    estimate's.
  - *invariant share at 4 samples:* no — the likelihood's own maximum sits about 1.3 errors above the
    truth there; this is not the estimator's fault and a fix to the shapes would not remove it.
  - Louis's method, the alternative spec question 3 names, changes the estimator and not which
    products are summed, so it does not bear on the block question.
- **The Beta-shape update** (§2.3): the fit stops 0.8 to 2.3 errors from the likelihood's maximum on
  these cohorts, and the gap grows as the square root of the positions. Recommended before plan
  step B1: its stopping rule would otherwise settle there. The fix moves fitted numbers and is outside
  this plan.
