# Fit precision, step C3 — stopped: the climb's projected gain stops it short

**Date:** 2026-10-01. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step C3.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §4.3, §4.5 item 2. **Branch:** `fit-precision`, at
`189a4350`; the step's code is **uncommitted** in the worktree (`tmp/fit_precision/c3/c3_uncommitted.patch`).

**Why this report and not a commit.** C3 was built as spec §4.3 states it, and its own validation fails it in
the way the handoff said to stop on: on drawn strata and on the four-accession oracle cohort, the rule stops
climbs more than a tenth of an error short of a climb run to twenty rounds. At the production span of thirteen
allele classes it stops some of them several log-likelihood units short. What follows is the measurement and a
recommendation, for the owner's decision.

## 1. What was built

The repeat-tract fit climbs each stratum's numbers by coordinate ascent: one **round** moves every number once
(the three slippage numbers, every allele class's share, the concentration), each by a golden-section search,
and a **walk** is the rounds from one of three starting points. A walk used to stop when a round raised the
log-likelihood **a tract, on average** by less than 10⁻⁶, or after five rounds.

As built (spec §4.3):

- **The total, not the mean.** The rule works on the stratum's total log-likelihood, so that it means the same
  at 10 samples and at 3,000.
- **A target from the errors.** Near the answer, the total still to gain is about half the sum of each number's
  squared distance in its own errors. So "every number within a tenth of its error" leaves at most
  ½ · p · 0.1² to gain. For p numbers that is 0.08 at one slippage group and thirteen classes (p = 16).
- **A projection from the last two rounds.** A walk stops once the gain still to come, projected from its last
  two rounds' gains, is below that target. The projection assumes the gains shrink by a steady factor λ, so
  `gain · λ / (1 − λ)` is still to come. λ is capped at 0.95, and a first round is projected at that cap.
- **A losing round.** A round that ends lower than it began has its moves undone; the walk stops at the best
  point it stood at and records that it did, rather than calling it convergence.
- **The record.** `StratumFit` now says how its winning walk ended — settled, stopped at a round that lost, or
  out of rounds — and records every walk's ending and rounds. The log gains one line with the counts.

Unit tests of the rule itself pass (scripted gains; the total's units; a losing round; tracts without reads
leaving every walk the same). The comparison below is an ignored test in the same diff,
`the_climbs_stop_loses_nothing_against_twenty_rounds`.

## 2. What was measured

**The comparison.** Each stratum is fitted by the rule. The same three walks are then climbed twenty rounds with
no rule, and the best point any of them reached is taken as the answer to compare against. Every number's
distance is measured in that twenty-round answer's own standard error (step C1's), as the plan asks. "Short"
means how many total log-likelihood units below the twenty-round answer the fit stopped. "The old rule" is the
same walks stopped where the 10⁻⁶-a-tract rule would have stopped them, read off the same twenty-round paths.

### 2.1 Drawn strata, five rounds at most (`tmp/fit_precision/c3/twenty1.log`)

| stratum drawn | rounds, three walks a draw, all draws | numbers more than a tenth of an error short | furthest | most short |
|---|---|---|---|---|
| 3 classes, 30 tracts × 8 samples × 3 reads, 20 draws | 192 (old rule 295) | none of 136 | 0.071 errors | 0.003 units |
| 3 classes, 300 × 20 × 3, 20 draws | 187 (old 288) | 8 of 60 class shares; none of the others | 0.283 | 0.039 |
| 3 classes, 300 × 20 × 30, 20 draws | 150 (old 287) | 11 of 60 class shares | 0.501 | 0.149 |
| 3 classes, 400 × 63 × 3, 5 draws | 54 (old 72) | 2 of 15 class shares, 1 of 5 shorter shares | 0.366 | 0.092 |
| 13 classes, 200 × 20 × 3, 5 draws | 72 (old 75) | 58 of 84 numbers | 10.5 | 17.6 |

- **At three classes the rule saves a quarter to a half of the rounds**, and at 20 samples the slippage numbers
  and the concentration stay within 0.075 errors of twenty rounds (at 63 samples one shorter share is 0.112
  short). The class shares do not: 21 of the 195 compared end 0.1 to 0.5 errors short. The old rule left every
  number within 0.060 errors.
- **At thirteen classes five rounds are too few under either rule.** 9 of the 15 walks ran out of rounds, and
  three of the five draws stopped 3.5 to 17.6 units short.

### 2.2 Two reasons, each measured

**The target bounds a sum, not each number.** 0.08 is what is left when every number is a tenth of an error
away. But the same 0.08 can sit entirely in one number, which is then √(2 × 0.08) = 0.4 errors away. A target
of 0.005 (a tenth of an error for one number alone) was run on the same three-class draws
(`variant_3class_target005.log`): 186 to 209 rounds a regime against the old rule's 287 to 295, every number
within 0.086 errors at 3 reads, but **at 30 reads 5 of 60 class shares were still 0.1 to 0.30 errors short**.
A tighter target helps and does not suffice.

**The projection cannot see a climb that crosses a plateau.** At thirteen classes with up to twenty rounds
(`variant_13class_rounds20.log`), each walk's projection when it stopped is set against what its own path went
on to gain:

| that walk's gains, round by round (total log-likelihood units) | projected still to gain | its path still gained |
|---|---|---|
| …, 106.6, 13.7, 0.98 | 0.076 | 5.54 |
| …, 226.1, 15.0, 1.50, 0.25 | 0.052 | 15.2 |
| …, 4.78, 2.67, 0.81, 0.23, 0.072 | 0.034 | 12.8 |

The gains collapse for a round or two and then resume, and a projection from the last two rounds cannot tell
that from the end. This is how the same projection failed the SNP/indel fit at step B1 (spec §2's amendment: "the
projection cannot see a slow approach"). With twenty rounds allowed, 12 of 15 walks still stopped by the rule,
and two of the five draws stopped 5.5 and 17.4 units short.

**Stopping at a losing round costs as much.** In the draw that stopped 17.4 short, all three walks lost a round
(at rounds 4, 6 and 8) and stopped there. Refusing a coordinate move that scores worse, so that no round can lose,
was also run (`variant_13class_rounds20_guard.log`): no round lost, and the same draws stopped 36.3 and 5.5
units short. The projection is the defect; the losing-round handling is secondary.

### 2.3 The four-accession oracle cohort (`tmp/fit_precision/oracle_c3p/fit.log`)

Fifteen strata are fitted on their own tracts (four samples at about three reads; thirteen classes). Each was
also climbed twenty rounds by a probe in the same run, which changed no output.

- **Rounds: 166 against 217 under the old rule**, 23% fewer. Of the 45 walks, 36 settled, 6 stopped at a round
  that lost and 3 ran out of rounds; the winning walk settled in 12 of the 15 strata.
- **Neither rule reaches the twenty-round answer on several strata**: 4 of the 15 stop more than 0.08 units short
  under both rules, and 10 of the 15 (the old rule: 8) have a number more than a tenth of an error away. The worst
  stratum, motif length 1 at 8 repeats, stops 4.5 units short under this rule and 17.8 under the old one.
- **Its checksums moved**: the fitted parameters file and the calls made with it differ from the baseline. They
  were not re-recorded, since the step is not committed.

## 3. For the owner

**Recommendation: do not ship the projected gain.** Replace it with the rule step B1 adopted for the SNP/indel
fit, built from step C1's errors: a walk has settled when the Newton step from its current point — its
curvature, which C1 computes, solved against its slope — is within a tenth of each number's error.

- **Why that rule.** It measures each number's distance directly, in that number's own error, so it is neither
  fooled by a plateau nor able to let one number carry the whole target. On the SNP/indel fit it replaced the
  same projection for the same reason, and every fit it called converged was within 0.088 errors of 600 passes.
- **What it costs.** A curvature is 513 evaluations at thirteen classes, about 1.7 rounds of the climb, plus 32
  for the slope. It need be taken only when a round's gain falls below a trigger (the projection as built makes a
  reasonable one), and the winning walk's last one is the curvature C1 already computes for its errors. At three
  classes that is close to the rounds the projection saves.
- **The round limit must rise with it.** At thirteen classes five rounds are too few under any rule: 9 of 15
  walks ran out (§2.1). How many are enough is not measured — the twenty-round climbs are the yardstick here, not
  shown to have converged — and a limit of 20 would make the climb up to four times as long at the production
  span. Step D's subsets cut the cost of a round.
- **A losing round** should be undone and the walk left to the settled test, not stopped (§2.2) — not yet
  measured under the Newton rule.
- **What I would do next**, if you agree: build the Newton test on C1's `curvature_of` and `invert_identified`
  behind the same trigger, rerun this comparison and the oracle, and amend spec §4.3 with you as §2 was amended.
  The spec and the plan are not edited.

**What stands from this step** if it is redone that way: the total-log-likelihood units, the record of how each
walk ended, the log line, and the comparison test.
