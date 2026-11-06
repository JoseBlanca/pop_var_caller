# Does the SNP/indel fit's curvature estimate cost it passes? (2026-11-05/06)

**Outcome: the fit and its standard errors stay as they are.** The fit estimates how sharply its
log-likelihood curves from the outer product of each position's slopes, the slopes of every
parameter multiplied pairwise and summed over the positions. That estimate is wrong in both
directions, by up to seven times at the answer and by far more in the middle of the climb. But
stepping with the true curvature makes the 59-sample fit slower, not faster. Reporting standard
errors from the true curvature works and is exact, but costs about 7% of the fit's time at 59
samples. It also changes which parameters get an error at all. The owner kept today's errors to
avoid the cost.

Everything below was measured on the tomato cohort (about three reads a position): one sample, ten
and 59. The code is not on `main`; it is kept under the tag `fit-curvature-experiment`, every part
behind an environment switch that is off by default.

## 1. What was asked

The fit stops when the Newton step still to take is worth less than 0.01 log-likelihood units (the
*decrement*). Near the end it shrinks only by a factor of 0.4 to 0.7 a pass, where an exact Newton
step would shrink it much faster. The question
([`fit_precision.md`](../../ng/spec/fit_precision.md) §8, question 3) was whether the outer product
is the cause, and whether to replace it.

## 2. The true curvature, against the outer product

The true curvature was measured by central differences of the summed slope, every parameter in
turn. The analytic slope agrees with the log-likelihood's own difference to a thousandth of a
standard error, so the comparison is sound.

Each figure below is how much sharper (above 1) or flatter (below 1) the true curvature is than the
outer product, along the directions where the two differ most:

| where | flattest | sharpest | the directions |
|---|---|---|---|
| 1 sample, at the answer | 1.0 | 2.0 | mismapped share with both error rates |
| 10 samples, at the answer | 0.14 | 2.25 | 0.14: mismapped share with every sample's mismapped-class error rate; 2.25: two samples' homozygote excess with the frequency shapes |
| 59 samples, at the answer (every 8th position) | 0.28 | 1.13 | mismapped share with every sample's mismapped-class error rate |
| 10 and 59 samples, mid-climb | — | 30 to 2,190 | duplicated share with the carrier-frequency shapes |

So along the mismapped share the outer product is up to seven times too sharp, and the reported
error for it is too small: by a factor of 1.4 to 1.75 at 10 and 59 samples.

## 3. Stepping with the true curvature makes the 59-sample fit slower

| fit | finish passes, outer product | finish passes, true curvature |
|---|---|---|
| 1 sample | 8 | 7 |
| 10 samples | 20 | 10 and 16 (two runs) |
| 59 samples (curvature from every 8th position) | 15 | 70 (blocks only), 81 (whole matrix) |

At 59 samples the true curvature had a direction where the likelihood curves upward, along
allele-frequency shape `a` resting on its lower bound of 0.02. Every step needed damping, and then
the steps overshot, because curvature from every 8th position did not match the slope from every
position. Taking the curvature from every position would cost eight times as much a step. Answers
stayed within 0.055 standard errors of today's on all but 4 of 6,476 numbers.

**A second cause of slow passes, untouched by the curvature:** in the 10-sample run that took 16
passes, passes 5 to 15 walked a ridge where the duplicated share rose from 0.010 to 0.018 while
carrier shape `a` fell towards its bound of 0.02. More duplicated positions each carried by fewer
samples fit almost equally well. The run converged in one step once the duplicated share left the
set of parameters being solved for.

## 4. Standard errors from the true curvature: exact, but not free

Built behind `PVC_OBSERVED_INFORMATION=1`. The final pass also sums each position's likelihood
second derivatives (divided by the likelihood) in the blocks the errors already use. The true
curvature is then the outer product less that sum. The fit's steps, its stopping rule and the
starts' agreement still read the outer product, so the fitted values do not move; at 1, 10 and 59
samples they were identical. Every entry was checked against central differences, agreeing to 1
part in 10⁸ for each sample's own parameters (including two libraries, held-out reads and ranged
depths) and to 4 parts in 10⁵ for the frequency shapes.

**Cost:** the final pass took 3.2 plain passes at 10 samples and 5.1 at 59, against 1.9 today.
That is about 2% more fit time at 10 samples and 7% at 59. Nearly all of it is work done for each
sample at each of the 16 frequency nodes, at every position.

**Effect on the parameters file** (error with the switch on over today's):

| | 1 sample | 10 samples | 59 samples |
|---|---|---|---|
| invariant share | – | 0.66 | – |
| allele-frequency shape b | – | 0.88 | 1.05 |
| carrier-frequency shape b | – | 0.76 | 1.20 |
| expected heterozygosity | – | 0.18 | error lost |
| fixed non-reference share | 4.13 | 0.93 | 0.98 |
| per-sample error rates and inbreeding, median | 0.79 (one sample) | 1.00 | 1.00 |

The mismapped share, the parameter the true curvature corrects most, is not written to the
parameters file.

**Two side effects:**

- **A parameter on its bound has to be left out.** At 59 samples both shapes `a` sit on 0.02, and
  with them in, the true curvature cannot be inverted. Leaving them out works, but then every
  quantity derived from them loses its error (expected heterozygosity, the shares above the stated
  frequency), and so do inbreeding coefficients on their bound: 4 at 59 samples, 1 at 10. Leaving
  carrier shape `a` out is also what takes expected heterozygosity's error to 0.18 at 10 samples;
  with it in, the ratio was 0.63.
- **Which barely determined parameters get an error changes.** At one sample the outer product
  drops the frequency shapes as undetermined; the true curvature keeps them, with errors wider than
  their range, and that widens the fixed share's error four times. Over the parameters today keeps,
  the true curvature gives the fixed share the same error as today.

## 5. If this is reopened

- The code is under the tag `fit-curvature-experiment`: `information::SecondDerivatives` and its
  tests, `standard_errors::holding_parameters_on_their_bounds`, and the measurement switches
  `PVC_CURV*` and `PVC_DECREMENT` in `fit.rs`.
- To adopt the errors, two things are missing: the reports would have to treat a parameter on its
  bound as known, so derived quantities keep their errors; and the cost would need to come down.
  Skipping frequency nodes of negligible posterior weight is the obvious lever, but the terms they
  multiply can be as large as 10⁹, so the cut-off would have to be around 10⁻¹² and may skip
  little. That was not measured.
- The ridge in §3 is a separate lead on fit speed.
