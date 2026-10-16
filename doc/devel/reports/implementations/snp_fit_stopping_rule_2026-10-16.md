# Stopping the SNP/indel fit on the Newton decrement: measured at 1, 4 and 59 samples

**Date:** 2026-10-16
**Branch:** `snp-fit-convergence` (from `main` 51b7be99). Nothing is committed; the prototype was
built from a scratch copy of the tree (`tmp/snp_stop/proto/`), and `src/` is untouched.
**Question:** the lead proposed stopping each start of the SNP/indel fit when the **Newton
decrement** — the log-likelihood the remaining Newton step is worth, `½·gᵀI⁻¹g` within the
parameters' box — falls below a threshold, instead of today's rule that every parameter be within a
tenth of its standard error of the maximum. What does it cost and what does it leave, at thresholds
1, 0.1 and 0.01 (and 0.001, added when 0.01 fell short), and which threshold should the spec adopt?

The spec text this feeds is [`fit_precision.md`](../../ng/spec/fit_precision.md) §2 (the amendment
of 2026-10-16) and [`parameter_prepass_joint_fit.md`](../../ng/spec/parameter_prepass_joint_fit.md)
§2.1.3 and §3.3.1.

---

## 1. Summary

**Recommendation: replace today's rule with the decrement, at a threshold of 0.001.** It is the only
threshold measured that brings every determined quantity — the mismapped share, each read group's
error multiplier, each sample's homozygote excess, the share segregating above 1%, expected
heterozygosity, the mean alternative frequency — within a tenth of its standard error of the fully
converged value on every input: the worst was always the mismapped share, at 0.021 to 0.085 of its
error. It costs the first start three to six passes beyond 0.01, and the whole fit two to three
passes beyond today's rule (§4.1). **0.01 is the cheaper alternative** if a quarter to a third of an
error on the mismapped share is acceptable: it stops where today's rule does on the oracle inputs.

**0.1 and 1 are too loose.** At 0.1 the mismapped share was left 0.55 to 1.17 of its error away on
the 4- and 59-sample inputs, and on 59 samples read one chunk in four it left expected
heterozygosity 0.21 and the mean alternative frequency 0.29 of their errors away. At 1 the 4-sample fit read one chunk in eight stopped 4 log-likelihood units short
with the mismapped share 11 errors away.

**Today's rule does not deliver what it says either.** Where it first held, the mismapped share was
still 0.13 to 0.47 of its error from the converged value, though the rule promises 0.1. Both rules
read the same Newton distances, and those distances are short by a factor of about two near the
maximum, because the information they are solved with is three to six times too sharp along the
slowest direction (§3.2). The threshold above is chosen from the measurement, not from the rule's
own promise.

**The decrement does not stop the fit on the 59-sample ridge.** Along the ridge (§3.4) the decrement
stayed between 26 and 45 units, and the fit gained 68 units there; it fell below 1 only once `a`
reached its lower bound. On these inputs the decrement's advantage over today's rule is that one
number bounds every parameter and every derived quantity at once, not that it stops sooner.

**One trap to carry into the code**: a negative decrement means the box solve failed, and must not
stop the fit (§3.3).

---

## 2. Inputs and how the fits were run

| name below | what it is | threads |
|---|---|---|
| **1 sample** | `SRS1839219.psp` from the identity oracle, all positions | 4 |
| **4 samples, 1 chunk in 8** | the oracle's four psps, the fit reading every 8th chunk of positions; contamination skipped | 4 |
| **4 samples** | the oracle's four psps, all positions | 4 |
| **59 samples, 1 chunk in 4** | tomato2's 59 psps, every 4th chunk; contamination skipped | 18 |
| **59 samples** | tomato2's 59 psps, all positions — the reference run only (§4.2) | 18 |

The oracle psps are `pop_var_caller/tmp/snp_fit_baseline/oracle_psps/`, tomato2's
`pop_var_caller/tmp/tomato2/psps/`. Every run is `pvc estimate-parameters` from the diagnosis's
prototype (`tmp/snp_conv/prototype_knobs.diff`) with the Newton step of decision 1 on (four tries,
the remembered fraction), extended with:

- `PVC_SNP_DECREMENT=<threshold>`: a start stops at a judging pass whose decrement is in
  `[0, threshold)`, at the point that pass judged; today's per-parameter rule is then not used.
- `PVC_SNP_JUDGE_TSV`: every judging pass appends the log-likelihood, the decrement, the count of
  parameters not settled under today's rule, every parameter's value, error and Newton distance, and
  the cohort's 8 × 8 covariance; the returned fit is appended too, after the final pass.

**The reference** for each input is a run at threshold 10⁻⁴ (pass limit 400). Its last judging pass
is within about 5 × 10⁻⁴ log-likelihood units of the maximum (its decrement, under 10⁻⁴, times the
understatement of §3.2), far below any difference reported here.

**Distances are in the reference's standard errors.** For a parameter, its own error. For the share
segregating above 1%, expected heterozygosity and the mean alternative frequency, the delta method
over the reference's covariance of the four density numbers. "Error multipliers" are the read
groups' error rates at ordinary positions; the multiplier is that rate over a constant, so its
distance in errors is the same. Each row quotes the **largest** distance over the read groups or the
samples.

**Two ways of reading a threshold, both reported.** *Post hoc* (§3): the first start's trajectory in
the reference run is the same as in a run at any larger threshold up to the pass where that run stops,
so the reference's judging passes say where each threshold would have stopped the first start and
what it would have returned; the runs at 0.001 reproduced these numbers exactly. *Actual runs* (§4): one run per threshold, all three starts, which also
shows how the later starts' agreement test reacts.

**Wall times are noisy.** The machine was shared with other jobs during the oracle and 59-sample runs
(load average up to 29 on 18 cores); passes are the measure to compare.

---

## 3. Where each threshold stops the first start

### 3.1 The tables

Pass at which the first start stops (its final pass follows), log-likelihood short of the reference,
and the largest distance in errors of each determined quantity. "Today" is the first judging pass at
which every parameter was within a tenth of its error by today's rule.

**1 sample** (reference: 26 passes, log-likelihood −111,299.006). The density's shapes are not
identified at one sample, so the derived quantities have no error, and the homozygote excess is held.

| threshold | stops at pass | log-likelihood short | mismapped share | error multipliers |
|---|---|---|---|---|
| 1 | 14 | 0.27 | 1.13 | 0.74 |
| 0.1 | 14 | 0.27 | 1.13 | 0.74 |
| 0.01 | 18 | 0.027 | 0.32 | 0.22 |
| 0.001 | 22 | 0.0023 | 0.075 | 0.051 |
| today | 19 | 0.015 | 0.23 | 0.16 |

**4 samples, 1 chunk in 8** (reference: 48 passes, −102,443.820).

| threshold | stops at pass | short | mismapped share | error multipliers | homozygote excesses | share above 1% | `Hexp` | mean alt. frequency |
|---|---|---|---|---|---|---|---|---|
| 1 | 30 | 4.0 | 11 | 3.5 | 0.24 | 0.048 | 0.027 | 0.012 |
| 0.1 | 40 | 0.24 | 0.55 | 0.29 | 0.054 | 0.009 | 0.005 | 0.002 |
| 0.01 | 43 | 0.022 | 0.13 | 0.074 | 0.021 | 0.003 | 0.002 | < 0.001 |
| 0.001 | 46 | 0.001 | 0.021 | 0.012 | 0.004 | < 0.001 | < 0.001 | < 0.001 |
| today | 43 | 0.022 | 0.13 | 0.074 | 0.021 | 0.003 | 0.002 | < 0.001 |

**4 samples** (reference: 57 passes, −796,940.963).

| threshold | stops at pass | short | mismapped share | error multipliers | homozygote excesses | share above 1% | `Hexp` | mean alt. frequency |
|---|---|---|---|---|---|---|---|---|
| 1 | 41 | 2.25 | 2.73 | 1.01 | 0.022 | 0.005 | 0.002 | < 0.001 |
| 0.1 | 44 | 0.42 | 1.17 | 0.45 | 0.008 | < 0.001 | < 0.001 | < 0.001 |
| 0.01 | 49 | 0.028 | 0.27 | 0.11 | 0.002 | < 0.001 | < 0.001 | < 0.001 |
| 0.001 | 53 | 0.003 | 0.069 | 0.027 | < 0.001 | < 0.001 | < 0.001 | < 0.001 |
| today | 49 | 0.028 | 0.27 | 0.11 | 0.002 | < 0.001 | < 0.001 | < 0.001 |

**59 samples, 1 chunk in 4** (reference: 52 passes, −2,020,462.594).

| threshold | stops at pass | short | mismapped share | error multipliers | homozygote excesses | share above 1% | `Hexp` | mean alt. frequency |
|---|---|---|---|---|---|---|---|---|
| 1 | 32 | 2.54 | 1.20 | 0.23 | 0.22 | 0.64 | 1.33 | 1.61 |
| 0.1 | 35 | 0.31 | 0.70 | 0.094 | 0.050 | < 0.001 | 0.21 | 0.29 |
| 0.01 | 39 | 0.049 | 0.32 | 0.042 | 0.007 | < 0.001 | 0.031 | 0.044 |
| 0.001 | 45 | 0.005 | 0.085 | 0.011 | 0.002 | < 0.001 | 0.003 | 0.003 |
| today | 37 | 0.11 | 0.47 | 0.063 | 0.018 | < 0.001 | 0.080 | 0.11 |

On this input the threshold of 1 also left the Beta's second shape `b` 1.8 of its error away; at 0.1,
0.43.

**59 samples, all positions** (reference: 48 passes, −7,985,014.368).

| threshold | stops at pass | short | mismapped share | error multipliers | homozygote excesses | share above 1% | `Hexp` | mean alt. frequency |
|---|---|---|---|---|---|---|---|---|
| 1 | 35 | 2.39 | 1.93 | 0.30 | 0.10 | 0.016 | 0.30 | 0.39 |
| 0.1 | 38 | 0.24 | 0.72 | 0.11 | 0.017 | 0.003 | 0.052 | 0.070 |
| 0.01 | 41 | 0.031 | 0.25 | 0.037 | 0.003 | < 0.001 | 0.010 | 0.013 |
| 0.001 | 45 | 0.002 | 0.047 | 0.007 | < 0.001 | < 0.001 | 0.001 | 0.001 |
| today | 41 | 0.031 | 0.25 | 0.037 | 0.003 | < 0.001 | 0.010 | 0.013 |

### 3.2 Why the thresholds land where they do: the decrement understates what is left

The decrement predicts the gain still to come. Measured against the reference, the log-likelihood
actually still to gain was, once the decrement was below 0.1:

| input | log-likelihood still to gain ÷ decrement |
|---|---|
| 1 sample | 3.2 to 4.0 |
| 4 samples, 1 chunk in 8 | 2.5 to 3.3 |
| 4 samples | 3.4 to 4.2 |
| 59 samples, 1 chunk in 4 | 4.2 to 5.9 |
| 59 samples | 3.0 to 3.5 |

(the last two or three passes before the reference's own stop are left out: there the reference's
own remaining 10⁻⁴ dominates.) In the middle of the climb on the subsampled oracle it reached 10.

Near the maximum every whole Newton step is accepted, and each pass's decrement is 0.4 to 0.6 times
the last on 1 and 4 samples, 0.68 on 59 — linear convergence, not the quadratic convergence an exact
Newton step gives. So the outer-product information the step is solved with is three to six times as
sharp as the likelihood along the slowest direction, and the step covers only a fifth to a third of
the remaining way. **Today's rule has the same blind spot**: its distances are the same Newton step,
which is why it left the mismapped share 0.13 to 0.47 of its error away while saying 0.1.

So the threshold that meets "a tenth of an error" is about `0.1² / (2 × 6)` ≈ 0.001, and that is what
the tables show. *This is also evidence on [`fit_precision.md`](../../ng/spec/fit_precision.md) §8
question 3 (is the outer-product estimator close enough to the observed information?): along the
slowest direction, on these inputs, it is three to six times sharper, which would make the errors
along it a half to two-fifths of what the observed information gives. Not measured directly here.*

### 3.3 A negative decrement is a failed box solve, not convergence

On the 4-sample, one-chunk-in-8 input, starts 2 and 3 both reached pass 30 with a decrement of −2,959
and −2,804. The quadratic model's maximum within the box is never below zero (staying put is in the
box), so a negative value means the active-set solve ended at its round limit without finding it;
at that pass the step was sending several error rates at mismapped positions and the carrier Beta's
shapes to their bounds. The first build of the prototype read it as converged and stopped both starts
13 and 14 log-likelihood units below start 1. Guarded (`0 ≤ decrement` required), both went on and
agreed with start 1. Today's rule did not stop there either: it counted 16 parameters unsettled.

### 3.4 What the decrement does on the ridge

On 59 samples (one chunk in 4) the first start walked the ridge between the invariant share and the
Beta's first shape `a` from pass 17 to pass 31: the invariant share fell from 0.969 to 0.858 and `a`
from 0.109 to its bound of 0.02. The decrement over those passes was 26 to 45 units, and the
log-likelihood still to gain at pass 17 was 68 units, so the ridge was not flat and neither rule
should have stopped on it. Once `a` sat on its bound the decrement fell to 0.97 in one pass. **No
input here tests a ridge that is flat without a bound**; on these inputs the bound is what ends it.

---

## 4. Actual runs

### 4.1 All three starts

Passes per start (c converged, a agreed with start 1), SNP/indel stage time, CPU seconds of the whole
command. "Today" is the diagnosis's run with the Newton step and today's rule
(`tmp/snp_conv/runs/{s1_newton,s8_newtonA_all,o4_newton}`) or, on 59 samples, this measurement's.

| input | rule | passes by start | all | stage | CPU s |
|---|---|---|---|---|---|
| 1 sample | 10⁻⁴ (reference) | c26 c24 c25 | 75 | 1 m 14 s | 301 |
| | 1 | c14 c16 c15 | 45 | 42 s | 175 |
| | 0.1 | c14 c16 c15 | 45 | 42 s | 177 |
| | 0.01 | c18 c17 c17 | 52 | 50 s | 207 |
| | today | c21 c20 c20 | 61 | — | 218 |
| 4 samples, 1 in 8 | 10⁻⁴ | c48 a46 a44 | 138 | 29 s | 126 |
| | 1 | c30 c39 c39 | 108 | 23 s | 103 |
| | 0.1 | c40 a43 a41 | 124 | 27 s | 116 |
| | 0.01 | c43 a45 a43 | 131 | 29 s | 124 |
| | today | c45 a45 a43 | 133 | — | 118 |
| 4 samples | 10⁻⁴ | c57 a46 a39 | 142 | 4 m 24 s | 1,055 |
| | 1 | c41 a37 c43 | 121 | 3 m 55 s | 941 |
| | 0.1 | c44 a42 c44 | 130 | 5 m 53 s | 1,367 |
| | 0.01 | c49 a45 a37 | 131 | 5 m 03 s | 1,188 |
| | today | c51 a45 a37 | 133 | 3 m 48 s | 916 |
| 59 samples, 1 in 4 | 10⁻⁴ | c52 a46 a37 | 135 | 5 m 33 s | 4,224 |
| | 1 | c32 c44 a34 | 110 | 3 m 57 s | 3,377 |
| | 0.1 | c35 c46 a35 | 116 | 4 m 07 s | 3,566 |
| | 0.01 | c39 a46 c40 | 125 | 4 m 46 s | 3,910 |
| | today | c39 c50 a37 | 126 | 5 m 37 s | 3,907 |
| **1 sample** | **0.001** | c22 c20 c21 | 63 | 57 s | 236 |
| **4 samples, 1 in 8** | **0.001** | c46 a45 a44 | 135 | 29 s | 124 |
| **4 samples** | **0.001** | c53 a46 a37 | 136 | 3 m 59 s | 959 |
| **59 samples, 1 in 4** | **0.001** | c45 a46 a37 | 128 | 3 m 54 s | 3,608 |

At one sample the later starts never stop by agreement under any rule; each converges on its own.

**What the returned fit was** — the best-scoring start's answer, after its final pass — against the
reference: log-likelihood short, and the largest distance in errors over the determined quantities
(the mismapped share unless named).

| input | 1 | 0.1 | 0.01 | 0.001 | today |
|---|---|---|---|---|---|
| 1 sample | 0.058; 0.41 | 0.058; 0.41 | 0.022; 0.29 | 0.002; 0.075 | 0.01; — |
| 4 samples, 1 in 8 | 1.82; 2.0 | 0.24; 0.55 | 0.022; 0.13 | 0.001; 0.021 | 0.02; — |
| 4 samples | 0.30; 0.31 | 0.017; 0.20 | 0.028; 0.27 | 0.003; 0.069 | 0.02; — |
| 59 samples, 1 in 4 | 1.12; 1.23 (mean alt. frequency) | 0.15; 0.36 (mean alt. frequency) | 0.016; 0.19 | 0.005; 0.085 | 0.007; 0.13 |

("today" on the three smaller inputs is the diagnosis's runs, which kept no judge table: only their
log-likelihood, as the diagnosis quotes it, is compared.)

At 0.001 the returned fit was the first start's on every input, with the distances §3 predicts; in
all, 0.001 took two to three passes more than today's rule (63 against 61, 135 against 133, 136
against 133, 128 against 126).

At looser thresholds the returned fit is often better than the first start's (§3), because a later
start that climbs on its own may end higher, and the best start wins. That is luck of the starts, not a property of the
threshold: at 0.01 on 59 samples start 3 converged on its own and its answer was returned.

**A looser threshold can cost passes.** A later start stops when its projected endpoint agrees with
the first start's answer within half that answer's error ([`fit_precision.md`](../../ng/spec/fit_precision.md)
§3.4). The less settled that answer, the likelier a later start misses it and climbs to its own end:
on the 4 oracle samples at 0.1, start 3 converged on its own after 44 passes where at 0.01 it agreed
after 37. At threshold 1, on the subsampled oracle, starts 2 and 3 converged on their own to a higher
log-likelihood than start 1, and the fit returned start 3's answer, 1.8 units short.

### 4.2 59 samples, all positions

One run, at 10⁻⁴, the only run of this measurement at the time, though other jobs shared the
machine: start 1 converged after 48 passes, and starts 2 and 3 agreed with it after 43 each — 134
passes, against 133 for the diagnosis's run with today's rule. Its SNP/indel stage took 24 min 05 s
and 16,756 CPU-seconds, against 15 min 51 s and 15,355 for the diagnosis's run; the difference is the
shared machine and the extra passes, not a measured cost of the rule. Its first start's trajectory
gives §3.1's last table; the thresholds were not run separately on all positions.

---

## 5. The test for "not separately determined"

The same runs give the inputs of decision 2's test
([`parameter_prepass_joint_fit.md`](../../ng/spec/parameter_prepass_joint_fit.md) §2.1.3): hold the
share segregating above 1% fixed and ask what halving or doubling the share of segregating positions
costs, from the covariance of the density's four numbers at the reference fit.

| input | segregating share | share above 1% ÷ segregating | its relative error | halving it costs, all else refitted | … along the trade, share above 1% held |
|---|---|---|---|---|---|
| 1 sample | 0.0057 | 0.41 | 0.07, shapes held | — | shapes not identified |
| 4 samples, 1 in 8 | 0.0014 | 1.00 | 0.40 | 1.5 units | no trade |
| 4 samples | 0.0016 | 1.00 | 0.20 | 5.8 units | no trade |
| 59 samples, 1 in 4 | 0.155 | 0.059 | 0.84 | 0.34 units | 0.36 units |
| 59 samples | 0.158 | 0.053 | 0.43 | 1.31 units | 1.35 units |

The middle column of costs is why the test is taken along the trade: on the subsampled oracle the
segregating share itself is poorly measured (halving it costs 1.5 units), but every segregating
position there is above 1% — `a` is 17.9 — so there is nothing below 1% for it to trade with, and
the invariant share and `a` are not confounded. On 59 samples, with the determined share held, the
segregating share can be halved for 1.35 log-likelihood units (0.36 when read one chunk in four),
under the two-unit line: the two numbers are not separately determined. There the invariant share
and `a` are correlated at 0.998 (1.00 on the subsample). The full cohort's 1.35 is not far below the
line. If the outer-product information is as much too sharp along the trade as along the slowest
direction of §3.2, the true cost is lower than 1.35, not higher.

---

## 6. Files

All under `/Users/jose/devel/pop_var_caller-snp-fit-convergence/tmp/snp_stop/`:

- `proto/` — the scratch copy the prototype was built from: the diagnosis's patch plus
  `newton_decrement` and `cohort_covariance` (`src/parameter_estimation/joint/fit/standard_errors.rs`)
  and the stopping knob and judge table (`src/parameter_estimation/joint/fit.rs`, search `PROTOTYPE`).
- `fit.sh`, `queue.sh`, `queue2.sh`, `queue3.sh` — how each run was made (`queue.sh` ran the
  one-sample runs; `queue2.sh` the rest after the negative-decrement guard; `queue3.sh` the 0.001 runs).
- `prototype_stopping.diff` — the prototype against `main` 51b7be99 (the diagnosis's patch included).
- `runs/` — each run's log, parameters file, judge table and `/usr/bin/time` output, named
  `<input>_d<threshold>` (`s1`, `s8`, `o4`, `t59s4`, `t59`) and `t59s4_settled` for today's rule.
- `analyse.py` (post hoc and actual comparisons, the test of §5), `passes.sh`, `summary.sh`.
