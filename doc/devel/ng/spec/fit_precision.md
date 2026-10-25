# Fitting to the precision the data support

**Status:** design, 2026-09-27; amended at checkpoints A and A′ (§1.3, §3.2; 2026-09-29), at plan step B1 (§2, §3.3, §3.4; 2026-09-29) at plan step C3 (§4.3; 2026-10-01; §4.1, §4.2, §4.5 and §6 brought into line at checkpoint C, 2026-10-02) and at plan step D2 (§4.4; 2026-10-02); on 2026-10-16 (owner) §2's per-parameter stopping test is replaced by one on the Newton decrement, §3.3 gains the Newton step the SNP/indel fit takes ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.3.1), and §3.5, §3.6 and §8 follow; on 2026-10-17 (owner, after the review of that build's first steps) §2 gains a box solve that always finishes and a decrement measured over the parameters that have a standard error. Not built yet; build order [`../impl_plan/snp_fit_convergence.md`](../impl_plan/snp_fit_convergence.md).
Build order:
[`../../implementation_plans/fit_precision.md`](../../implementation_plans/fit_precision.md). It
amends the SNP/indel fit ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.3,
which says only "repeat until the fitted values stop moving"), the repeat-tract fit
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.2) and the parameters file
([`parameters_file.md`](parameters_file.md) §2, §3). There is no separate architecture document:
the types it touches are few and are given here.

Decided with the owner on 2026-09-27, after `estimate-parameters` on 2,169 tomato accessions on
kimura took 9 h 25 min in the SNP/indel fit and about two hours a stratum in the repeat-tract fit
(the owner's run log; its figures are quoted in §1.1, the log itself is not in the repository):

> "What we need is not for the parameters that can be precisely estimated to be estimated with a
> high degree of precision, what we need is an approximate estimate for those and, most
> importantly, error estimates. Once we can calculate both we're fine and we can stop."

> **Part B was replaced, 2026-10-16; Parts A and C stand.** Part B (§4, and goals 3 and 4 of §1.2)
> is written for the climb: each repeat-tract stratum fitted on its own by coordinate ascent, its
> errors taken from that climb's curvature, its stopping rule, and the subset of samples a stratum
> was read from. Penalised expectation maximisation (EM),
> [`str_penalised_fit.md`](str_penalised_fit.md), replaced that fit; it fits every stratum of one
> motif length together, reads every sample, stops when no number moves by more than a fixed
> tolerance (its §2.7), and gives every number a standard error by its own method (its §2.6). Commit
> `3198d55d` deleted the climb, its errors, its stopping rule and the subsets; the links below to
> that code now name the commit. Also no longer applicable: §1.4's *round*, Part B's lines in §6, the
> first two items of §7 (the curves are now drawn inside the fit, not through its answers), and
> questions 4 and 5 of §8.
>
> **Still true:** the principle (§2); Part A, the SNP/indel fit (§3); Part C, standard errors carried
> in `Estimate<T>` and the parameters file (§5) — the repeat-tract rows still carry errors, now the
> new fit's; and §1.1's kimura measurements, which are why both fits were changed.

---

## 1. What this is

**Every parameter the fits produce gets a standard error, and the fits stop when the answer is
settled to a small fraction of that error — not to a fixed fraction of the value.** A parameter the
data determine well is fitted closely; one the data barely determine is fitted roughly, and its
large error says so. Separately, **the repeat-tract fit reads each stratum from as many samples as
its precision needs**, not from every sample in the cohort.

Three parts, each a module someone can build on its own:

- **Part A — the SNP/indel fit** (`src/parameter_estimation/joint/fit.rs`): standard errors, a
  stopping rule measured in them, and stopping a starting point early when it has reached the same
  answer as an earlier one.
- **Part B — the repeat-tract fit** (`src/parameter_estimation/joint/ssr_fit.rs`): standard errors
  for each stratum's fit, a stopping rule for its climb that does not depend on the cohort's size,
  and fitting each stratum on a reproducible subset of samples, grown until the stratum's slippage
  level is known well enough.
- **Part C — the parameters file** (`src/calling/parameters_file/`): writing the standard errors.

### 1.1 Why now: what the kimura run showed

| | measured on kimura, 2,169 samples, 2,651 read groups, median depth 10.8 |
|---|---|
| SNP/indel fit | three starting points, **each stopped at the 200-pass limit**, 9 h 25 min; all three ended at log-likelihood −1.293107 × 10⁹ |
| why it never stopped | the log-likelihood rule was met from about pass 100 (gains of 5 to 30 a cycle against a threshold of 200); **the parameter rule never was** — some parameter kept moving by 1% to 10% of itself a cycle |
| repeat-tract fit | one stratum at a time with every core busy: 7 strata of 141 in 11 h 20 min |

**The parameter rule measures each move relative to the value itself**, with a floor of 10⁻⁶
(`largest_relative_move`, [`fit.rs:2017`](../../../../src/parameter_estimation/joint/fit.rs#L2017);
floors at [`fit.rs:2052-2054`](../../../../src/parameter_estimation/joint/fit.rs#L2052)), over
**every** coordinate, including each sample's homozygote excess. Every start sets every sample's
homozygote excess to 0 ([`fit.rs:1833`](../../../../src/parameter_estimation/joint/fit.rs#L1833)). So
a sample whose excess stays near 0 turns an absolute move of 10⁻⁹ into a "relative move" of 10⁻³,
and with 2,169 samples there is always such a sample. **Which parameter was moving on kimura is not
known** — the log prints only the size of the largest move — so this is the likely cause, not a
measured one. A rule in units of each parameter's own error removes it either way: a move of 10⁻⁹ in
a quantity known to ±0.05 is nothing.

**The repeat-tract fit's cost grows with the number of samples** (each score sums over every sample
with reads at every tract, `ssr_fit.rs`
`ln_tract`, deleted in commit `3198d55d`), while a stratum's slippage is estimated from thousands of tracts and does not need
thousands of samples to be known well.

### 1.2 Goals

1. **Every fitted parameter has a standard error**, computed from the data the fit read, or is said
   to have none because the data carry no information about it. Never a number presented as
   precise when nothing determined it.
2. **The SNP/indel fit stops when every parameter is within a small fraction of its own standard
   error of where it is heading**, and stops a starting point early when it is heading where an
   earlier start already arrived.
3. **The repeat-tract climb stops by a rule that means the same thing at 10 samples and at 3,000.**
4. **The repeat-tract fit reads a stratum from a subset of samples**, chosen without looking at the
   data, grown until the stratum's slippage level has a relative standard error below a target.
   **Below the subset's starting size nothing changes**, so a small cohort is fitted exactly as
   today.
5. **The parameters file carries the standard errors** beside the values they belong to.
6. **Everything stays deterministic**: the same input gives the same bytes on every platform and at
   every thread count ([`src/cli/cross_platform_digests.rs`](../../../../src/cli/cross_platform_digests.rs)).

### 1.3 Non-goals, and what this does not do

- **It does not change the models** — **amended at checkpoint A (owner, 2026-09-28; recorded
  2026-09-29): two corrections do.** Each library's reads are now scored under that library's own
  error rates, where a sample's libraries shared the first one's (plan step A6), and the
  allele-frequency density's and the carrier Beta's shapes now step along the likelihood's own
  slope, where they solved the digamma form and stopped beside the maximum (step A7). Both move
  fitted numbers. Otherwise the same likelihoods, parameters and starting points; what changes is
  when the fits stop, which samples a stratum is read from, and what is reported.
- **It does not make calling read the standard errors.** Consumers "combine warrants; they do not
  branch on them" ([`parameters_file.md`](parameters_file.md) §2); whether calling should weight by
  a standard error is its own design (§7).
- **It does not re-weight the repeat-tract curves.** They weight each stratum by a count-based
  error, `1/√(slipped reads)` ([`slippage_curve.rs:133`](../../../../src/parameter_estimation/joint/slippage_curve.rs#L133)),
  which its own documentation says understates the error. Replacing it with Part B's error is a
  natural next step and is deferred (§7), because it moves every curve.
- **It does not subsample samples in the SNP/indel fit.** Every sample needs its own homozygote
  excess and every read group its own error rate; the subsample proposal of
  [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §7 stays open there.
- **It does not touch the repeat-tract guard**, which dropped 45,432 of 86,688 tracts on kimura (§7).
- **It does not run the starting points in parallel.** Deferred until the passes a start takes
  under the new rule are measured (§7).

### 1.4 Vocabulary

- **Standard error (SE).** How far the estimate would typically move if the same kind of data were
  drawn again: the spread of an estimate, in the estimate's own units. Computed here from how sharply
  the log-likelihood falls away from its peak — flat means a large error, sharp a small one.
- **Information.** The sharpness itself: the curvature of the log-likelihood at its peak. A standard
  error is one over the square root of it. A parameter the data say nothing about has zero
  information and no standard error.
- **Score.** The slope of the log-likelihood with respect to one parameter. Zero at the peak.
- **Cycle.** One accelerated step of the SNP/indel fit: three passes over the data or more
  (SQUAREM, Varadhan & Roland 2008; [`fit.rs:1778-1924`](../../../../src/parameter_estimation/joint/fit.rs#L1778)).
- **Round.** One sweep of the repeat-tract climb over all of a stratum's parameters, about 306
  evaluations of its likelihood at one slippage group
  (`ssr_fit.rs:1048-1100`, deleted in commit `3198d55d`).

---

## 2. The principle, and where it comes from

**Stop when the distance still to travel is small compared with the uncertainty of where the answer
is.** Two pieces, both with long statistical precedent:

1. **The distance still to travel is projected, not read off the last step.** These fits approach
   their answer in steps that shrink by a roughly constant factor λ each cycle (EM converges linearly,
   at a rate set by the fraction of missing information: Dempster, Laird & Rubin 1977). If the last
   step moved a parameter by Δ, the rest of the way is about Δ·λ/(1 − λ) — Aitken's extrapolation
   (Aitken 1926), used as an EM stopping rule by Böhning, Dietz, Schaub, Schlattmann & Lindsay (1994)
   and standard in mixture-model software. A small step in a slow fit is not closeness (Lindstrom &
   Bates 1988); the projection is what corrects for that.
2. **It is compared with the standard error, not with the value.** Finishing a computation more
   precisely than the data determine the answer buys nothing (Bates & Watts 1981, the "relative
   offset" criterion). The standard errors come from the same pass structure the fit already has
   (Louis 1982; the link between convergence rate and missing information is Meng & Rubin 1991).

*The citations are from memory and are to be checked before this document is cited elsewhere.*

**The rule, precisely.** For each parameter *j*, after each cycle:

- Δⱼ is its move over the cycle, on its natural scale;
- λⱼ = |Δⱼ| / |Δⱼ from the previous cycle|, **capped at `MAX_CONTRACTION` = 0.95** (soft — a
  starting value), and taken as 0 when the two moves have opposite signs (an oscillation is not a
  slow approach);
- the projected distance is dⱼ = |Δⱼ| · λⱼ / (1 − λⱼ) + |Δⱼ| — the step just taken counts as well,
  since the SE is computed at the cycle's start;
- **the parameter is settled when dⱼ < `SETTLED_FRACTION` × SEⱼ**, with `SETTLED_FRACTION` = 0.1
  (soft — to be chosen from the measurement in the plan's first checkpoint).

**The fit has converged when every parameter is settled.** A parameter with no standard error
because it has no information is settled by definition: nothing the fit does to it matters to the
likelihood.

**Why the natural scale and not logit.** The accelerated step works in logit coordinates
([`fit.rs:1993-2012`](../../../../src/parameter_estimation/joint/fit.rs#L1993)), and there a
homozygote excess heading to 0 moves forever (logit −10, −11, …) while its standard error in logit
units grows without bound too. On the natural scale the move is 10⁻⁹ and the error is a few
hundredths; the comparison is well behaved at the boundary, which is exactly where the current rule
fails.

**Why the cap on λ.** As λ approaches 1 the projection explodes, and those are the parameters with
the most missing information — the poorly determined ones, whose standard errors are large anyway.
The cap keeps one near-flat direction from holding the fit forever; with it, the projection is at
most twenty times the last step.

**Amended at plan step B1 (owner, 2026-09-29): for the SNP/indel fit the distance still to travel is
Newton's, not the projection above.** Built as written, the projection stopped the fit far too early
(the plan's B1 report, `fit_precision_b1_stopped_2026-09-29.md`): on the 4-accession oracle cohort
every start stopped at 36 to 57 passes, 74 log-likelihood units below a fit run to 1,000 passes, the
invariant share 28 of its errors away; and on 15 drawn cohorts nine stopped with a parameter more than a
tenth of an error short, up to 1.9 errors at two samples. Two causes. The errors in force had been
computed on a plateau the fit crosses slowly, where the invariant share's was 34 times wider than at the maximum. And
the projection cannot see a slow approach: where a path turns, the moves change sign and λ is taken as
zero; where it crawls, λ is 0.98 to 0.99, above the cap. A higher cap and a log-likelihood condition
repaired the 20-sample cohorts and not the 2- and 4-sample ones.

**The rule, as amended.** On a pass that sums the information (§3.3) the pass also sums each
parameter's score — the slope of the whole log-likelihood, `g`. The information solved against it,
`dⱼ = (I⁻¹ g)ⱼ`, is each parameter's distance to the likelihood's maximum as a Newton step estimates it,
on its natural scale; it does not depend on how fast the fit happens to be moving.

- The matrix is the one the errors come from (§3.2): the whole one for a small cohort, the blocks
  above, solved as an arrow — each sample's own parameters given the cohort's, the cohort's from what
  the samples leave of its slope.
- **The step stays inside the parameters' bounds** — the maximum of the quadratic model
  `g·d − ½ dᵀ I d` within the box the fit keeps the parameters in: a parameter whose step would carry
  it past an end of its interval stops at that end, its distance the way there, and the others are
  solved again with its move taken out of their scores; a stopped parameter whose slope at the
  solution points back into its interval is freed (the active-set rule for a box). *Found in B1's
  review:* holding only a parameter exactly on an end left a clean error rate the golden section rests
  4.4 × 10⁻¹⁰ above its floor, and a density shape the fit walks towards its lower bound, unsettled for
  1,000 passes; with the box both converge. *The way this rule finds the box's maximum was replaced
  on 2026-10-17* (the amendment at the end of this section): it often failed to finish on small
  cohorts.
- **A parameter is settled when |dⱼ| < `SETTLED_FRACTION` × SEⱼ**, both from the same pass; one with
  no error, or no distance (not solved for), is settled by definition. **The fit has converged when
  every parameter is settled**, and stops at the end of that cycle.

Measured on the same 15 cohorts, before the box: every fit this rule calls converged is within 0.088 of
an error of the 600-pass fit; the ones still crawling run to the pass limit and say so (all three two-sample
cohorts and two of the three 4-sample ones with the duplicated class). On the oracle no start converges within 200 passes, and none
claims to. Part B's climb (§4.3) keeps its projection of the log-likelihood gain.

**Amended 2026-10-16 (owner): for the SNP/indel fit, a start stops when the remaining Newton step is
worth less than `DECREMENT_THRESHOLD` = 0.01 log-likelihood units** (0.001 until 2026-10-25; see
below). Proposed by the lead after the
owner asked what the fit should do on a plateau, measured in
[`../../reports/implementations/snp_fit_stopping_rule_2026-10-16.md`](../../reports/implementations/snp_fit_stopping_rule_2026-10-16.md),
and decided by the owner the same day at the threshold that measurement supports. It replaces the
per-parameter test above ("a parameter is settled when …"); everything else in the amendment of plan
step B1 stays — the judging pass, its information and scores, the bounded step `d`, the trigger.

- **The quantity.** The **Newton decrement** is the log-likelihood the quadratic model says the
  bounded step would gain: `λ = g·d − ½ dᵀ I d`, with `g`, `I` and `d` exactly as above. Without
  bounds it is `½ gᵀ I⁻¹ g`.
- **Where it is computed.** By `newton_step`
  ([`fit/standard_errors.rs:583`](../../../../src/parameter_estimation/joint/fit/standard_errors.rs#L583)),
  which already holds `g`, `I` and `d`: it returns the decrement beside the distances, from the flat
  step `bounded_step` gives with each parameter not solved for moving 0, as
  `Σⱼ gⱼ dⱼ − ½ Σⱼ dⱼ (I d)ⱼ`, the product `I d` from `information_times` (the same file). One box
  solve a judging pass; the prototype solved twice, once for the distances and once for `λ`.
- **The rule, and the order of a judging pass.** At each judging pass, at the parameters `θ` it
  judged:
  1. if `0 ≤ λ < DECREMENT_THRESHOLD` the start has **converged**. It stops **at `θ`**: the rest of
     the cycle, and the Newton step of
     [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.3.1, are not run. The
     final pass follows at `θ` as today; it is still needed, because only the final pass keeps each
     position's posteriors. *Until now a converged start finished its cycle first.*
  2. otherwise, for a later start, the agreement test of §3.4, unchanged, on the same `d`. **A start
     that has converged is never stopped by agreement**, as before: it competes to be the answer.
  3. otherwise the Newton step (§3.3.1 of the joint-fit spec), and, if none of its tries rises, the
     usual SQUAREM cycle from `θ`.
- **What it guarantees, in the model.** For any quantity `h` computed from the parameters — one
  parameter, or a derived number such as the share segregating above 1% or `Hexp` — the step moves
  it by at most `√(2λ)` of its own standard error (Cauchy–Schwarz on `I`). So one number bounds every
  parameter and every derived quantity at once: `λ` below 0.005 means all of them within 0.1 of their
  errors, which implies the per-parameter test; below 0.1, within 0.45. The per-parameter test,
  by contrast, does not bound a combination, and with 185 parameters each a tenth of an error from
  the maximum the remaining step can still be worth almost a unit.
- **Why it should stop on a plateau** (the lead's rationale, accepted by the owner). Along a direction
  the data barely determine, the information is small, so a long step there is worth little; `λ`
  weighs each direction by what moving along it is worth in likelihood rather than by how far it moves
  a parameter. *Not tested by the measurement:* on 59 samples the ridge between the invariant share
  and the Beta's shape `a` was not flat — `λ` was 26 to 45 units along it and 68 units were still to
  gain — and the fit stopped only after `a` reached its lower bound. No input measured has a ridge
  that is flat without a bound.

**Traps.**

- **A negative `λ` is not convergence.** The quadratic model's maximum within the box is never below
  zero, since staying put is in the box. The box solve as first built (the active-set rule of plan
  step B1 above, `3 × side + 3` rounds) could end at its round limit and return a step outside the
  box, and then `λ` had any sign. It happened on the 4 oracle samples read one chunk in eight: starts
  2 and 3 reached pass 30 with `λ` = −2,959 and −2,804, while several error rates and the carrier
  shapes were being sent to their bounds. Read as below the threshold, it stopped both 13 to 14
  log-likelihood units below start 1. **A non-negative `λ` did not prove the solve had finished**:
  on the two-sample fixtures of the cross-platform checks a failed solve gave 290. The solve of the
  2026-10-17 amendment below cannot return a negative `λ`, and reports one it did not finish as not
  a number. The test stays written `λ >= 0.0 && λ < threshold`: a NaN fails it, and so would a
  negative value from any later change to the solve. A start that fails it goes on to steps 2 and 3
  above, like any start not converged.
- **`λ` understates what is left, by a factor of two and a half to six near the maximum.** Measured
  on the converged references, the log-likelihood still to gain was 2.5 to 4.2 times `λ` once `λ` was
  below 0.1 on 1 and 4 samples, 3.0 to 5.9 times on 59, and up to 10 times in the middle of the climb
  on the subsampled oracle. Each full step is accepted and the decrement shrinks by 0.4 to 0.7 a
  pass, not quadratically: along the slowest direction the outer-product information is that much
  sharper than the likelihood (§8 question 3).
  **The per-parameter test has the same blind spot**, because its distances are the same `d`: where
  it first held, the mismapped share was still 0.23 (1 sample), 0.13 (4, subsampled), 0.27
  (4 samples), 0.47 (59, subsampled) and 0.25 (59 samples) of its error from the converged value. So
  the threshold is chosen from the measurement, not from `√(2λ)`.
- **A looser threshold can cost time.** A later start stops by agreeing with the first start's answer
  (§3.4); the less settled that answer, the likelier a later start is to miss it and climb to its own
  end. On the 4 oracle samples at 0.1, start 3 missed it and climbed 44 passes on its own, where at
  0.01 it agreed after 37: the looser threshold saved start 1 five passes and cost start 3 seven.
- **The constant is a config field, as `SETTLED_FRACTION` is.** `DECREMENT_THRESHOLD` lives beside
  `SETTLED_FRACTION` in [`fit/settled.rs`](../../../../src/parameter_estimation/joint/fit/settled.rs)
  and reaches `maximise` as `JointFitConfig::decrement_threshold`, defaulted to it, so a test can run
  a fit to `10⁻⁴` the way the measurement's references were.

**`DECREMENT_THRESHOLD` = 0.01 since 2026-10-25 (owner), re-measured after the starts were thinned.**
Since `parameter_prepass_joint_fit.md` §3.3.2 the starts explore on a seventh of the census and the
winner finishes on all of it, from close to the answer. Measured on the returned fit, against the
same fit run to 10⁻⁴, the largest distance of any number the parameters file gives an error, in
that error, and the fit's instructions retired against 0.001:

| input | 0.001 | **0.01** | 0.1 | instructions, 0.01 against 0.001 |
|---|---|---|---|---|
| 1 tomato sample | 0.035 | **0.012** | — | −11% |
| 8 tomato samples | 0.011 | **0.055** | 0.72 (the duplicated share above 1%) | −30% |
| 62 tomato samples | 0.003 | **0.049** | 0.13 (the second frequency shape) | −24% |

Calling from the 0.01 fits changed 4 of 495,712 genotypes at 8 samples and 11 of 10,491,082 at 62.
**0.01 now meets the tenth of an error that 0.001 was chosen for**; the measurement below, made on
the first start's own path before the starts were thinned, is why 0.01 was rejected then. Inputs are
still tomato only, at about 3 reads a position.

*The 2026-10-16 decision, kept for its measurement:* **`DECREMENT_THRESHOLD` = 0.001 — decided
2026-10-16 (owner), measured.** Where each threshold
would have stopped the first start, against a run to 10⁻⁴, largest distance in errors over the
determined quantities (the mismapped share, except at 1 on 59 samples read one chunk in four: the
mean alternative frequency) and the log-likelihood short:

| input | 1 | 0.1 | 0.01 | **0.001** | the per-parameter test |
|---|---|---|---|---|---|
| 1 sample | 1.13; 0.27 | 1.13; 0.27 | 0.32; 0.027 | **0.075; 0.002** | 0.23; 0.015 |
| 4 oracle samples, one chunk in 8 | 11; 4.0 | 0.55; 0.24 | 0.13; 0.022 | **0.021; 0.001** | 0.13; 0.022 |
| 4 oracle samples | 2.7; 2.3 | 1.17; 0.42 | 0.27; 0.028 | **0.069; 0.003** | 0.27; 0.028 |
| 59 tomato2 samples, one chunk in 4 | 1.6; 2.5 | 0.70; 0.31 | 0.32; 0.049 | **0.085; 0.005** | 0.47; 0.11 |
| 59 tomato2 samples | 1.9; 2.4 | 0.72; 0.24 | 0.25; 0.031 | **0.047; 0.002** | 0.25; 0.031 |

0.001 is the only threshold measured that meets the tenth of an error the per-parameter test was
meant to give. It costs the first start three to six passes beyond 0.01 and three to eight beyond the
per-parameter test (19 → 22, 43 → 46, 49 → 53, 37 → 45 and 41 → 45 on the five inputs); in whole
fits, two to three passes. *Rejected:* **0.01**, which left the mismapped share a quarter to a third
of an error away (the numbers calling reads were within 0.22 of their errors at one sample and 0.11
elsewhere); **0.1, first proposed**, which left the mismapped share 0.55 to 1.17 of its error away on
the 4- and 59-sample inputs, and on 59 samples read one chunk in four `Hexp` 0.21 and the mean
alternative frequency 0.29. The threshold is in log-likelihood units and the errors scale with the
data, so nothing in it depends on the cohort's size. Every input measured was tomato accessions: no
cohort above 59 samples, and no other depth, was measured.

The per-parameter test survives in two places: the progress line and a non-converged start's report
still name the parameter furthest from settled, in errors (§3.5), and `SETTLED_FRACTION` keeps that
meaning there. **It no longer stops anything.**

**Amended 2026-10-17 (owner): the box solve always finishes, and the decrement counts only the
directions the data determine.** Decided after the review of the build's first three steps
([`../../reports/reviews/snp_fit_convergence_AB_2026-10-16.md`](../../reports/reviews/snp_fit_convergence_AB_2026-10-16.md),
findings B1 and M1). That review traced why the cross-platform checks' small cohorts — two samples
without sequencing errors, two with, and 21 — stopped converging once the decrement rule was built:
every start of all three ran to the 200-pass limit. Two causes stack, and each has its own decision.

**1. The box solve is the primal active-set method.** The rule of plan step B1 above fixes every
parameter that crosses an end of its interval in one round and frees every one whose slope points
back in, for at most `3 × side + 3` rounds. Nothing makes that sequence finish, and on information
that is nearly singular it cycles. On the two-sample fixture with sequencing errors 95 of its 100
solves ended at the round limit (75 of 102 on the fit module's own 3-sample test cohort), and the step
it then returned lay outside the box, so the decrement read from it meant nothing: 290, −36 and
−1,100 were all seen. It is replaced, with no fallback to it, by the textbook method (Nocedal & Wright,
*Numerical Optimization*, §16.5):

1. start from staying put, `d = 0`, which is inside the box;
2. solve the free parameters' Newton step from the current point, the fixed ones held where they are;
3. move along it until the first free parameter meets an end of its interval; fix that parameter
   there — and any other that meets an end at the same point — and go back to 2;
4. after a whole step, free the one fixed parameter whose slope points furthest back into its
   interval and go back to 2; when none does, that point is the box's maximum.

A slope counts as pointing back in only when it exceeds 10⁻⁹ of the size of the two terms it is the
difference of, so that rounding cannot free and re-fix the same parameter forever. Every point the
solve passes through is inside the box, and the quadratic model's value never falls along the way —
each move is a Newton step on a subspace, or part of one — so **the decrement is never negative**.
When no bound binds, the first solve is the whole step and the result is the old one, bit for bit.
The round limit stays, at ten rounds a parameter plus ten; a solve that reaches it still returns a
step inside the box, which the Newton step of §3.3.1 may try, but **its decrement is reported as not
a number**, so the start cannot stop on it and runs on.

*Measured* (review, B1): no solve unfinished on any digest fixture, any fit test of the module, or
the 1- and 4-sample oracle; the 1-sample parameters file byte for byte as before; on 4 samples the
same passes (53 converged; 46 and 37 agreed), every SNP/indel number moved by at most 10⁻¹¹ of
itself, the repeat-tract slip rates by at most 2 × 10⁻⁵.

*Three departures from the textbook method, from the review of the built solve*
([`../../reports/reviews/snp_fit_convergence_B_2026-10-17.md`](../../reports/reviews/snp_fit_convergence_B_2026-10-17.md),
findings 1, 2 and 4; built 2026-10-17):

- **No parameter is freed twice in one solve.** The method assumes a freed parameter moves back into
  its interval and that the free set's maximum is unique; on nearly singular information neither
  holds. Of 6,000 random cases of 1 to 8 samples, 28 did not finish in ten rounds a parameter: one
  parameter freed, its next step pointing out, fixed again and freed again from round 21 to 199; or
  two parameters taking turns, each turn gaining 10⁻⁹ of a step worth 10⁻⁶. With each parameter
  freed at most once the solve takes at most three rounds a parameter plus one, so it always
  finishes, and the round limit is a guard. The price is that on such information the point it
  finishes at need not be the box's maximum. On information well conditioned no parameter needed a
  second freeing — 2,400 of the 6,000 cases and 134 of the property test's 400 — and the 1,200 best
  conditioned of the 6,000 end at the box's maximum, checked by the slope left at each parameter.
- **A parameter beside the end its step heads to is on it**, while the solve has not moved it:
  within 10⁻⁹ of one plus the end's size, the golden section's stopping width, which leaves a
  homozygote excess up to half that above 0. Fixing only the parameter meeting an end first spent a
  round on each such excess: at 2,000 outbred samples at their floor the solve took 1.83 s, against
  1.3 ms with the excesses inside; it now takes 1.9 ms.
- **The determined decrement is raised to zero when rounding makes it negative.** It is a small
  difference of large terms, and on exactly singular information it came out negative on 125 of
  3,546 random cases (as low as −1.4 × 10⁻¹⁰ against terms of 2.3 × 10⁴), which `has_converged`
  refuses. A decrement that is not a number stays one.

None of the three moves a number on the 1- and 4-sample oracle: both parameters files are byte for
byte those of the build before. On the 21-sample digest fixture the second moves the fitted numbers
by at most 6 × 10⁻⁶ of themselves, the same starts taking the same passes.

*A known limit (code review E6, 2026-10-17; not fixed).* **The solve fixes one parameter a round,
so when many parameters cross an end of their interval at different points along the step, its cost
grows with the square of the samples**: each round is one solve over every parameter, and there is
a round for each parameter that meets an end. Measured on a contrived case — one library a sample,
every sample's homozygote excess at 0.01 to 0.02 with a slope pushing it past 0, each at a different
fraction of the step: 10 ms a solve at 59 samples, 5.4 to 5.9 s at 2,000. With no parameter near an
end a solve at 2,000 samples costs 1.3 to 1.9 ms. A fit takes one solve a judging pass, so at 2,000
samples in that case each judging pass would gain several seconds; whether real cohorts come near it
is not measured. The fixes, in order of effort:

- **measure first**: a judging pass's share of wall time on a real cohort of a few thousand samples;
- **fix every parameter whose room is within a tolerance of the shortest in one round.** This does not
  help the case above, where the crossings are spread over the whole step, and moves parameters onto
  ends a little before the step reaches them, so it is not done;
- **a gradient-projection solve** (Moré & Toraldo, 1991; Nocedal & Wright §16.7), which fixes every
  parameter a projected step carries past an end at once, in a handful of rounds whatever the
  samples. It replaces this solve, and would need the oracle and the digests re-measured.

**2. The decrement counts only the parameters that have a standard error; the others adjust
freely.** With the solve fixed, the two-sample fixtures still ran to the limit. What was left to gain
lay along parameters the data do not determine — on the fixture with sequencing errors 12 of its 16
parameters have no standard error (the density's shapes against their bound, the duplicated class);
on the one without, all 16 — and there the decrement fell only to about 0.02 by pass 200. This
section already calls such a parameter "settled by definition"; the decrement did not.

- **The quantity, in words.** `½ dᵀ I d` is the step's length measured in standard errors: for any
  quantity `h` built from the parameters, the step moves `h` by at most `√(dᵀ I d)` of `h`'s own
  error. Split the parameters the step moves into those with an error, `K`, and those without, `U`.
  The decrement is now **the smallest that length can be when the determined parameters move as the
  step says and the undetermined ones may move anywhere** — the undetermined ones are profiled out.
  (Where the solve over the undetermined ones cannot tell one apart from the rest, it leaves that one
  still, and the length is a minimum over fewer directions: larger than `λ_K` below, on 12 of 3,546
  random cases by up to 3.3 times, never above 2.5 × 10⁻⁵.)
  That is the largest `½ (Δh / SE_h)²` over every quantity `h` built from determined parameters
  alone: the guarantee of `√(2λ)` above, restricted to the quantities it promises anything about.
- **Precisely.**

  ```text
  λ_K = min over e_U of ½ (d_K, e_U)ᵀ I (d_K, e_U)
      = ½ d_Kᵀ (I_KK − I_KU I_UU⁻¹ I_UK) d_K
      = ½ d_Kᵀ (Σ_KK)⁻¹ d_K
  ```

  `d_K` is the determined parameters' part of the whole bounded step; `I_KK`, `I_KU`, `I_UU` are the
  blocks of the information; `Σ_KK` is the determined parameters' block of its inverse, their joint
  covariance. The minimum is at `e_U = −I_UU⁻¹ I_UK d_K`. A parameter with an error that the step does
  not solve for is in `U` too; one with neither an error nor a distance (no information, held fixed,
  not told apart from the others) is in neither and stays put, as it does in the step.
- **When every parameter the step moves has an error, the decrement is the one above**,
  `g·d − ½ dᵀ I d`, computed exactly as before, so a cohort where everything is determined does not
  move by a bit. (Without a bound binding the two formulas give the same number, `½ dᵀ I d`; with one
  they differ, and the rule keeps the old formula there for that reason.)
- **It scales with the cohort as the step does.** The information is an arrow — each sample's own
  parameters paired with each other and with the cohort's eight, two samples' never paired (§3.2).
  `e_U` is one solve restricted to `U`, by the same block elimination the step's solve uses: each
  sample's own block, then the cohort's eight from what the samples leave. `λ_K` is then one product
  of the information with `(d_K, e_U)`. No dense inverse is formed; the cost is linear in the number of
  samples. *Measured* (plan step B4, on drawn arrows of one library a sample, the density's two
  shapes and one sample's rate in thirty without an error): 7 µs a judging pass at 59 samples (185
  parameters) and 160 µs at 2,000 (6,008), against 0.12 ms and 3.0 ms for the box solve and 0.08 ms
  and 1.8 ms for the standard errors the pass already computes. A dense inverse at 2,000 samples
  would be some 7 × 10¹⁰ operations a judging pass. On the 1 and 4 oracle samples every parameter
  the step moves has an error, so the old decrement is computed and the new one costs nothing.
- *Measured* (review, M1; together with 1 and with the error-rate swap deleted,
  [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.3): the two-sample fixture
  without errors converges in 16 passes a start where it ran to the limit; the one with errors in 22,
  13 and 13, best log-likelihood −482.27 where 200 passes of the corrected fit reach −482.16; the
  21-sample fixture in 84, 38 and 44, −5,045.75 against −5,045.64. On the 1 and 4 oracle samples,
  where every parameter the step moves has an error, nothing changes.
- *Rejected, measured:* **leaving the parameters without an error out of the Newton system.** It
  converges the fixtures too (16; 43, 21, 15; 95, 62, 59 passes), but it holds those parameters still
  while expectation–maximisation keeps moving them, and that drags the determined ones: on the
  4-sample test cohort with the duplicated class a parameter stopped 0.41 of its error from a long
  fit, and the 4 oracle samples took 56, 49 and 37 passes in place of 53, 46 and 37. **Treating an
  unfinished solve as "cannot judge"** without replacing it: no start can stop, and every digest
  start still ran to the limit.

**Trap: an undetermined parameter can still be moving when the start stops.** That is the point of
the decision, not a defect: the parameters file reports such a parameter without a standard error,
and the run log says why it has none (§3.5).

---

## 3. Part A — the SNP/indel fit

### 3.1 What is fitted, for reference

`Parameters` ([`fit.rs:1222-1233`](../../../../src/parameter_estimation/joint/fit.rs#L1222)), in the
order `Parameters::coordinates` fixes ([`fit.rs:2059-2089`](../../../../src/parameter_estimation/joint/fit.rs#L2059)):

| parameter | per | count |
|---|---|---|
| mismapped share | cohort | 1 |
| allele-frequency density: share invariant, share fixed non-reference, Beta shapes a, b | cohort | 4 |
| duplicated class: share, carrier Beta shapes (when on, the default) | cohort | 3 |
| error rate at ordinary positions, error rate at mismapped ones | read group | 2G |
| homozygote excess | sample | S |

8 + 2G + S in all: **7,478 on kimura** (G = 2,651, S = 2,169). Contamination is fitted afterwards
and is not in this vector.

### 3.2 Standard errors: the information from the data the fit already reads

**Decided: the information is estimated by summing, over positions, the outer product of each
position's score** — the "outer product of gradients" estimator. A position's score is the
expectation of its complete-data score under the posterior the E-step already computes (Fisher's
identity), so it needs no new pass structure, only new sums in an E-pass. The fit's log-likelihood
is a sum of one term per position, each computed on its own in the E-step (the pass counts them at
[`fit.rs:2497`](../../../../src/parameter_estimation/joint/fit.rs#L2497)), which is what makes the
sum of the positions' score products an estimate of the information.

**Alternatives considered.**

- *Numerical second derivatives of the total log-likelihood.* Two E-passes per parameter, 7,478
  parameters, about a minute a pass on kimura: rejected as unaffordable. Kept as the **oracle** on
  small cohorts (§3.6).
- *Louis's method* (complete-data information minus missing information): exact observed
  information, but needs every complete-data second derivative as well. More code for the same
  purpose; rejected for now. The outer-product estimator and the observed information agree at the
  maximum when the model is right and differ when it is not; the oracle comparison measures how much
  they differ here.

**Decided: a block approximation, not the full 7,478 × 7,478 matrix.** A sample's own parameters —
its homozygote excess, and its read groups' two error rates each — interact with each other and with
the eight cohort-level ones, and with other samples' parameters only through those eight. So for
each sample the fit accumulates the (k + 8) × (k + 8) block of its k own parameters and the cohort's
eight, and the cohort's own 8 × 8 block. A sample's standard errors come from its block with the
cohort's parameters accounted for (the Schur complement: the sample's information less what the
cohort-level parameters explain); the cohort's from its own block.

- **Cost.** Per position, per sample, a (k + 8)² update: about 2 million positions × 2,169 samples
  × 121 multiply-adds (k = 3), roughly 5 × 10¹¹ on kimura — of the order of one pass. **Memory:**
  (k + 8)² numbers a sample, a few kilobytes; nothing like the full matrix's 450 MB.
- **What it gives up:** the correlation between two samples' parameters through the cohort's. The
  oracle (§3.6) measures it on cohorts small enough to hold the full matrix.
- **Where a sample carries several read groups** k is 1 + 2 × its read groups; the block grows
  with it.

**Amended at checkpoints A and A′ (owner, 2026-09-28 and 2026-09-29) — what was built:**

- **The blocks are inverted exactly, as an arrow**, not by the Schur-complement sketch above: the
  cohort's errors come from `C − Σ_s B_sᵀ A_s⁻¹ B_s` inverted (the cohort's information less what
  each sample's own parameters explain), and a sample's from `A_s⁻¹ + A_s⁻¹ B_s V B_sᵀ A_s⁻¹` with
  `V` that inverse — the sample's own uncertainty plus the cohort's carried through the parameters
  they share. The sketch was measured no closer on any kind and further on the cohort's (plan step
  A3). A sample of k libraries carries 1 + 2k own parameters: its excess and each library's two
  rates (step A6).
- **Small cohorts use the whole matrix** (checkpoint A's decision 3, step A8): for a cohort of at most
  `FULL_MATRIX_SAMPLES` = 20 samples, and at most 188 parameters (checkpoint A′), every pairing is
  summed, two samples' included, and the matrix is inverted at once — the parameters taken in the
  blocks' order, each sample's own then the cohort's, so the same one is dropped when others mimic
  it. At 4 samples and 3 reads the blocks' errors were too small (the error rates and the mismapped
  share scattered 1.23 to 1.67 times them) and the whole matrix's about right (0.89 to 1.02). The
  parameter limit keeps the final pass's matrices to about 18 MB (20 samples of four libraries); a
  cohort of samples with more libraries takes the blocks.

**The scores to derive**, one per parameter kind, each new code:

| parameter | the complete-data score is … |
|---|---|
| error rates | the read tallies' derivative in the rate — `maximise_error_rate` maximises the same function by golden section ([`fit.rs:3137-3163`](../../../../src/parameter_estimation/joint/fit.rs#L3137)) |
| homozygote excess | the genotype prior's derivative in the excess, weighted by the posterior genotype counts `maximise_hom_excess` reads ([`fit.rs:3166-3182`](../../../../src/parameter_estimation/joint/fit.rs#L3166)) |
| density and carrier Beta shapes | the slope of the likelihood the quadrature rule computes, as the rule's nodes and weights move with the shapes (`RuleSlopes`, step A1) — **amended 2026-09-29**: not the digamma terms, whose integral on the rule misses that slope. Since step A7 the fit's update of the shapes steps along this same slope (`step_beta_shapes`), and `fit_beta_shapes` is gone |
| the three shares | the class posteriors over the share, as in their closed-form M-steps ([`fit.rs:3028-3048`](../../../../src/parameter_estimation/joint/fit.rs#L3028)) |

**Trap: the per-position score must be the observed-data score at that position**, the posterior
expectation of the complete-data score — not the complete-data score at the posterior mean. For
the Beta shapes the difference is the gap between E[ln f] and ln E[f], which is not small at three
reads.

**Trap: bounds.** Several parameters sit on or near a bound (`hom_excess` at 0, the shares near 10⁻⁹;
bounds at [`fit.rs:2160-2168`](../../../../src/parameter_estimation/joint/fit.rs#L2160)). There the
curvature-based error is not the spread of a normal distribution, because the estimate cannot move
past the bound. The error is still computed and written (it is still the right scale for the
stopping rule), and the file's reader is told what it is (§5).

**No information.** A sample with no reads, or a read group with none, has a zero score at every
position, so its row of the block is zero and its error does not exist. It is reported as absent,
never as zero and never as a large number (`Option::None`). **At one sample** the homozygote excess
is carried but never moved ([`fit.rs:3117-3127`](../../../../src/parameter_estimation/joint/fit.rs#L3117))
and has no error either; the error rates and the cohort's parameters still get theirs.

**Amended 2026-09-29 — two more reasons an error is absent**, found on drawn cohorts (plan steps A3
and A5), so a parameter without an error says which of four it is:

- **no information** — the case above;
- **held fixed** — the homozygote excess at one sample;
- **not identified** — other parameters mimic its effect at every position, so once they are
  accounted for less than 10⁻⁸ of its own curvature is left (at one sample the four density
  parameters, at two samples the duplicated class's share and carrier shapes). Only that parameter is
  dropped; the others are inverted without it;
- **wider than its range** — its error came out wider than the whole interval the fit keeps it in,
  so the data do not place it anywhere in that interval.

### 3.3 When the standard errors are computed

**Decided: not every cycle.** The fit keeps its current log-likelihood rule — the gain over a cycle
below 10⁻⁴ per position ([`fit.rs:1886`](../../../../src/parameter_estimation/joint/fit.rs#L1886)) —
as the **trigger**. When it first holds, the next E-pass also accumulates the information, and from
then on the settled test of §2 is applied each cycle with those errors. They are recomputed every
`ERROR_REFRESH_CYCLES` = 10 cycles (soft) while the fit continues. On kimura the trigger held from
about pass 100, where the current rule ran to 200 three times.

**The largest-relative-move rule is removed.** `stillness` no longer stops the fit; the two rules
are the log-likelihood trigger and the settled test. `max_passes` stays as a backstop, and a fit
that reaches it still returns what it has, as today.

**The final pass** ([`fit.rs:1909-1921`](../../../../src/parameter_estimation/joint/fit.rs#L1909))
also accumulates the information, so the errors reported are those at the returned parameters.

**Amended at plan step B1 (owner, 2026-09-29): once the trigger has held, every cycle's first pass
accumulates the information and the scores**, at the cycle's starting parameters, and the Newton
distances of §2's amendment are judged there; `ERROR_REFRESH_CYCLES` is gone. The fit can stop only
where fresh errors and distances say it is settled — errors ten cycles old were how the projection
came to judge the oracle's plateau by the plateau's much wider errors. The cost is one information
pass a cycle, at 1.48 to 1.80 times a plain pass (4 to 64 samples, one thread, measured in B1's
review), in cycles of three or more passes: about 11 to 13% more time on a small cohort's fit.

**Amended 2026-10-16 (owner): a judging pass's Newton step is also taken**, kept only if the
log-likelihood rises and halved otherwise
([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.3.1). An accepted step's pass
sums the information, so it is the next cycle's judging pass: from the first acceptance on, most
cycles are one judging pass long, and the information is summed on nearly every pass.

### 3.4 Stopping a starting point that is heading where another arrived

The three starting points exist because a start that puts the ordinary and the mismapped class close
together can collapse them and report convergence
([`fit.rs:496-503`](../../../../src/parameter_estimation/joint/fit.rs#L496)). They are run one after
another and the best log-likelihood wins
([`fit.rs:1617-1619`](../../../../src/parameter_estimation/joint/fit.rs#L1617)). On kimura all three
reached the same answer.

**Decided: a later start stops as soon as its projected endpoint is within `AGREEMENT_FRACTION` = 0.5
(soft) of the best earlier start's standard errors, on every parameter.** It then contributes
nothing new and is recorded as having agreed. A start that ends anywhere else runs to its own
convergence, as today, and the best log-likelihood still wins. So the saving is real only when the
starts agree, which is exactly when running them all to the end was wasted.

- **Trap: the test compares a start's projected endpoint** (its current value plus the projected
  distance, §2) **with the earlier answer**, not its current value — early in a start, the current
  value is far from anything. *Amended at plan step B1 (owner, 2026-09-29): the projected endpoint is
  the current value plus the Newton step of §2's amendment, from the cycle's information pass.*
- **Trap: the first start's errors are the yardstick**, so the first start must itself have
  converged under §2; if it hit `max_passes`, the others run to their own end.

### 3.5 What is reported

- **Per start:** converged, stopped at the limit, or agreed with an earlier start; passes; the
  log-likelihood; and, when it did not converge, **which parameter was furthest from settled, and by
  how many of its standard errors** — which spec §3.2 already asks for ("Reporting *not converged*
  without saying which parameter is still moving invites a correct number to be discarded",
  [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) lines 619-624).
- **The progress line** each cycle gains, once errors exist, the count of parameters not yet settled
  and the worst one's name and distance in errors.
- **`JointFit`** ([`fit.rs:264-361`](../../../../src/parameter_estimation/joint/fit.rs#L264)) gains
  each estimate's standard error through `Estimate` (§5.1), and a record per start. `converged` has no
  consumer outside `fit.rs` today; the parameters file records each start's outcome (§5.2).

*Amended 2026-10-16, with §2's decrement rule:*

- **The progress line** of a cycle that judged also prints the decrement and the threshold, before
  the count of parameters not yet settled: `the remaining Newton step is worth 3.1e-2
  log-likelihood units (stops below 1e-3); 4 of 185 parameter(s) still more than 0.1 standard
  errors …`. The per-pass trace (`PassSummary`,
  [`fit.rs:1631`](../../../../src/parameter_estimation/joint/fit.rs#L1631)) gains `decrement:
  Option<f64>` beside `not_settled`, `None` on a pass that did not judge.
- **What *converged* means changes, and three texts say the old meaning**: the parameters file's
  note on `snp_indel_fit_starts`
  ([`to_toml.rs:159`](../../../../src/calling/parameters_file/to_toml.rs#L159), and the golden file
  `testdata/every_shape_as_written.toml:106` that pins it), the doc comment of the file's
  `converged` outcome ([`parameters_file/mod.rs:714`](../../../../src/calling/parameters_file/mod.rs#L714)),
  and the module doc of [`fit/settled.rs`](../../../../src/parameter_estimation/joint/fit/settled.rs).
  The note says instead: *`converged` means the step still left to the best value the data allow was
  worth less than a thousandth of a log-likelihood unit, which on every cohort measured put every
  number within about a tenth of its standard error of that value.* The key and its values do not
  change, so this is not a format change.

### 3.6 How we know Part A works

1. **Scores against finite differences.** On a small drawn cohort, each parameter's summed score
   equals the numerical derivative of the total log-likelihood, to the precision of the difference.
2. **The block approximation against the full matrix.** On cohorts of 4, 20 and 63 samples the full
   outer-product matrix fits in memory; the block errors are compared with the full-matrix errors,
   and the difference is reported per parameter kind.
3. **The errors mean what they say.** Draw many cohorts from known parameters (the fit's own test
   generators), fit each, and compare the spread of the estimates with the reported errors: about
   68 in 100 estimates within one error of the truth, 95 within two. At 3 reads and at 30.
4. **The stopping rule loses nothing.** Fit the 4-accession oracle cohort
   (`scripts/promote_ng_oracle.sh`) and, where available, the 63-accession cohort, once with the new
   rule and once run to 1,000 passes; every parameter differs by less than `SETTLED_FRACTION` of its
   error, and the passes saved are reported.
   *Amended 2026-10-16:* under the decrement rule (§2), every determined quantity — each parameter
   with an error, and the share segregating above 1%, `Hexp` and the mean alternative frequency —
   within `SETTLED_FRACTION` of its error of a run to a decrement of 10⁻⁴, on 1 sample, the 4
   oracle samples and the 59 tomato2 samples. Unit tests: on a small information matrix with no
   bound binding, the decrement equals `½ gᵀ I⁻¹ g`; with one bound binding it equals
   `g·d − ½ dᵀ I d` for the box step `d`, and is smaller than the unbounded value. The verdict is a
   function of the decrement and the threshold alone, in `fit/settled.rs`, tested at −2,959, at NaN
   (neither converged), at 0 and just below the threshold (converged) and at the threshold (not).
5. **Determinism.** The same bits at 1, 4 and 8 threads; the new sums follow the fit's fixed chunk
   order ([`fit.rs:2235-2352`](../../../../src/parameter_estimation/joint/fit.rs#L2235)).

---

## 4. Part B — the repeat-tract fit

### 4.1 What one stratum's fit estimates, for reference

Per stratum (`ssr_fit.rs:803-809`, deleted in commit `3198d55d`):
three slippage numbers per slippage group (level, shorter share, fall-off), a length spectrum over 13
allele classes, and a concentration. Production pools every read group into one slippage group
(`every_read_group_pooled`, [`census_fit.rs:227-234`](../../../../src/run/census_fit.rs#L227)), so
17 numbers a stratum. The score is the **mean log-likelihood a tract**
(`ssr_fit.rs:2457`, deleted in commit `3198d55d`); the climb is
coordinate ascent by golden section, 18 evaluations a coordinate, up to 5 rounds from each of 3
starts, stopping when a round gains less than 10⁻⁶ on that mean
(`ssr_fit.rs:895`, deleted in commit `3198d55d`). **That was the climb
before plan step C3**; a walk now stops as §4.3's amendment says, within 40 rounds.

**Trap: the mean is over every tract in the stratum, including those with no reads**, and a tract's
log-likelihood sums over its samples. So the same 10⁻⁶ is a looser test in a stratum where most
tracts are empty and a stricter one as samples are added. That is why the climb's rule changes too.

### 4.2 Standard errors for one stratum

**Decided: the numerical curvature of the stratum's total log-likelihood at its answer**, by central
differences over all its numbers on the scales the climb uses (logit for the slippage numbers and
shares, log for the concentration, log-ratios for the spectrum), inverted to give the errors, then
carried to the natural scale.

- **Cost:** `1 + 2p²` evaluations over the p numbers the errors are taken over — 513 at p = 16 (one
  slippage group and thirteen classes; the shares sum to one, so one fewer than the 17 numbers fitted),
  against about 4,600 for the climb from three starts under the rule before plan step C3: roughly a
  ninth more, once per stratum at its final answer (not per start). Since C3 this curvature is usually
  the winning walk's last judgement (§4.3), and is not paid again.
- **Why all the numbers and not only the slippage three:** a read off the reference length is either
  a slip or a real allele, so the level and the spectrum trade against each other
  (`ssr_fit.rs` `PeriodLengthSpectrum`
  documentation, deleted in commit `3198d55d`). Holding the spectrum fixed would understate the level's error.
- **The step for each difference** is set relative to the coordinate's own scale, and the result is
  checked for a negative-definite curvature; where it is not (a flat or saddle direction), that
  number's error is reported as absent.
- **Trap: the climb's score is the mean a tract; the curvature must be of the total** (the mean
  times the number of tracts), or every error comes out too large by the square root of the tract
  count.

### 4.3 When the climb stops

**Decided: on the projected remaining gain of the stratum's total log-likelihood**, not on a fixed
gain of the mean. Near the answer the total log-likelihood still to be gained is about ½ Σ (dⱼ/SEⱼ)²
over the stratum's p numbers, so "every number within `SETTLED_FRACTION` of its error" is, in total,
a remaining gain of at most ½ · p · `SETTLED_FRACTION`² — 0.085 at p = 17. The climb projects the
remaining gain from successive rounds' gains (the same Aitken projection as §2, applied to the
log-likelihood: Böhning et al. 1994) and stops when it is below that.

- **Why the total gain and not the per-number test:** the per-number test needs the errors every
  round, which cost about two rounds each (§4.2). The total gain needs nothing new and means the same
  thing whatever the number of samples, because the total log-likelihood's curvature, not its mean,
  is what the errors are made of.
- **Trap:** a round that makes the score *worse* is counted as converged today
  (`ssr_fit.rs:895`, deleted in commit `3198d55d`); golden section
  never re-evaluates its starting value, so it can move a coordinate to a worse one. Under the new
  rule a round that loses is **not** convergence: the round's worse moves are undone and the walk
  stops at the best point it held, recorded as such.
- `max_rounds` stays as a backstop.

**Amended at plan step C3 (owner, 2026-10-01): a walk stops when a Newton step says every number is
within `SETTLED_FRACTION` of its error, not on the projected gain**, as §2 was amended for Part A. Built as
written, the projection stopped climbs short of a longer climb of the same objective
(`fit_precision_c3_stopped_2026-10-01.md`), for two reasons. The target bounds a sum: the 0.08 left when every
number is a tenth of an error away can sit in one number, which is then 0.4 errors away — at three allele
classes 21 of 195 class shares stopped 0.1 to 0.5 errors short. And the projection cannot see a climb crossing a
plateau: at thirteen classes, gains of 13.7 and then 0.98 were projected to leave 0.076, and the walk went on to
gain 5.5. Five rounds were also too few at thirteen classes: 9 of 15 walks ran out.

**The rule, as amended.**

- **The judgement.** At a point the walk has reached, the curvature and the slope of the stratum's total
  log-likelihood (§4.2's central differences; the slope comes from the same evaluations) give each number's
  distance to the maximum as a Newton step, `I⁻¹ g` over the numbers the curvature identifies, carried to each
  number's own scale as its error is. **The walk has settled when every number with an error is within
  `SETTLED_FRACTION` of it**; a number without one is settled by definition — not placed (its error on the climb's
  scale is wider than three units of logit or log), not identified (§3.2), or a class with no share. A judgement
  costs one curvature, `1 + 2p²` evaluations: 513 at thirteen classes, about 1.7 rounds, a round being 307
  evaluations (17 golden sections of 18, and the score where it ends).
- **When a walk is judged.** Once the projected remaining gain above is below ½ · p · `SETTLED_FRACTION`² — p is
  the numbers the errors are taken over, 16 at one slippage group and thirteen classes, so 0.08 — the projection
  being a trigger, not the test; and after a round that loses. A judgement that finds the walk unsettled puts
  the next off by one round more each time (one, then two, then three rounds later).
- **A round that loses** is undone back to the best point it stood at: the round's start, or a point part-way
  through it. A point part-way through is judged unless a judgement is put off, and unless it has settled the
  walk goes on from it. The round's start is always judged, put off or not, since a next round would repeat
  this one exactly: settled, the walk ends settled; not, it stops there, recorded as having lost a round, not as
  converged.
- **`max_rounds` is 40**, from 5. Measured with the limit at 60 on drawn strata: walks that settled took 2 to 6
  rounds at three classes and 6 to 36 at thirteen; at thirteen classes one walk of fifteen had not settled by
  60, and none settled between 37 and 60.
- **The winning walk's last judgement gives the stratum's errors** (§4.2) whenever its last round was judged at
  the point it returns, settled or not, so they cost no second curvature.
- §4.1's description of the climb — up to 5 rounds, stopping on a gain of the mean below 10⁻⁶ — is the rule
  before this amendment.

Measured with the limit at 60, against the best point the same walks reach in 60 rounds with no rule: at three
classes every number of every fit is within 0.058 errors; at thirteen classes four of five drawn strata are within
0.044 errors, and the fifth's winning walk stopped at a round that lost, 17.4 log-likelihood units short, and says
it has not settled.

### 4.4 A stratum read from a subset of samples

**Decided: each stratum is fitted on the samples that come first in a fixed, data-blind order,
starting with `FIRST_SUBSET` = 256 samples (soft), doubling until the stratum's slippage level has a
relative standard error below `LEVEL_RELATIVE_ERROR_TARGET` = 0.02 (soft), or every sample is in.**
Each larger subset starts its climb from the previous subset's answer.

- **The order.** Samples are ranked by a hash of the sample's **name** and a fixed seed — the same
  kind of fixed hash the census uses to choose positions
  ([`loci.rs:274-300`](../../../../src/parameter_estimation/joint/loci.rs#L274)) — so the subset is
  the same on every machine, does not depend on the order the psps were listed in, and is chosen
  without looking at depth or genotypes. **Decided: one order for all strata**, so that a sample is in
  every stratum's first 256 or in none: simpler to report and to reason about, at the cost that a
  sample with unusual chemistry is either in every stratum's fit or in none.
- **Every slippage group stays represented.** If a slippage group has reads in the stratum but none
  in the subset, that group's samples are added in hash order until it has `MIN_SAMPLES_A_GROUP` = 8
  (soft) or all of its own.
- **Small cohorts are untouched.** A cohort of at most 256 samples is fitted on all of them, exactly
  as today — including the cross-platform checksum's two-sample fixture.
- **What changes for the rest of the fit.** The stratum's evidence counts — reads crossing, slipped
  reads, tracts with reads — are those of the subset, so the curves across strata weight the stratum
  by the precision it was actually fitted to. The homozygote excess each sample carries is looked up
  by the sample's index as today (`ssr_fit.rs:2161-2165`, deleted in commit `3198d55d`);
  dropping rows from a stratum's evidence does not disturb it. The substitution counts are still
  summed over every sample when the evidence is gathered, since they cost nothing to keep.
- **Memory falls too:** a stratum's likelihood table has one row per sample with reads a tract, so
  it shrinks in proportion (5.9 GiB for the largest stratum on kimura at 2,169 samples).
- **Trap: the refusal floor is judged on the whole stratum**, not the subset: a stratum with reads in
  only the samples outside the subset must grow the subset rather than be refused. Refusal for too few
  tracts (`ssr_fit.rs:1135-1159`, deleted in commit `3198d55d`) looks
  at every sample's tracts before any subset is drawn.
- **The target is on the level, not on all three numbers.** One stratum can measure one of its
  numbers well and another to about 45% of itself ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md)
  §4.5, line 738), so a 2% target on every number would read every sample in every stratum.
  **Which number carries the target is soft** — decided here, and question 5 below.

**This is the §7 proposal of [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md)**
("bound the parameters fit's resident set by a subsample of *samples*, chosen the same data-blind way
the loci are", lines 1627-1634) applied to the repeat-tract half, with its question 8 — how large the
subsample must be — answered per stratum by the precision target rather than once for the cohort.

**Amended at plan step D2 (owner, 2026-10-02)**, after the step's review measured the rule above as first
built (`fit_precision_d2_2026-10-02.md`):

- **Each larger subset holds the smaller one's samples.** Built literally, a slippage group's added readers
  could leave the subset when it doubled.
- **Every slippage group with reads in the stratum is topped up to `MIN_SAMPLES_A_GROUP` of its readers**,
  not only one the first samples hold none of: one reader in 256 is barely better than none.
- **The subset takes every sample once doubling would hold more than three quarters of them.** At 2,169
  samples the subsets are 256, 512, 1,024, then every sample rather than 2,048.
- **A subset holding fewer tracts with reads than the refusal floor is not fitted; it grows.** The floor is
  still judged first on the whole stratum, so a stratum is never refused for its subset's thinness.
- **The target is judged only on the slippage groups that still have readers outside the subset.** A group
  whose every reader is in cannot be measured better by growing, and would otherwise hold every stratum to
  the whole cohort.
- **A larger subset is one walk from the last answer, which decides only whether to grow; the subset the
  answer is taken from is fitted from every starting point and from the last answer, the best winning.** On
  one drawn stratum of thirteen allele classes, one walk from the last answer ended 1.3 to 41
  log-likelihood units below three fresh starts on the same samples, in all four comparisons, once with its
  level 0.79 of an error away.
- **What else moves.** A subset-fitted stratum's evidence counts are the subset's, so the curves across
  strata weigh it by the subset's slipped reads, roughly 2,500 at thirteen classes once the target is met:
  against a curve with HG002's homopolymer held-out error (7.7%), the curve's share of a well-read stratum's
  blended level rises from under 1 in 100 to about 6 in 100. With several strata fitted at once, each stratum's
  subsets run one after another on one thread.

### 4.5 How we know Part B works

1. **The curvature errors mean what they say.** Draw strata at known slippage (the
   `bench_fixtures::draw_stratum` generator), fit each many times, and compare the spread of the
   fitted numbers with the reported errors, at 3 and at 30 reads.
2. **The new stop loses nothing.** On drawn strata and on the 4-accession oracle cohort, the fitted
   numbers under the new rule are within `SETTLED_FRACTION` of their errors of a longer climb of the
   same walks with no stopping rule (20 rounds as first written; 40 by the comparison test's default
   since plan step C3, when the round limit itself became 40).
3. **The subset loses little, and says how little.** On kimura's cohort (the owner's run), a handful
   of strata fitted on every sample and on the grown subset: the two levels differ by less than
   their errors, and the time saved is reported. Below 256 samples, byte-identical results to today.
4. **Determinism:** the subset is the same whatever order the psps are given in.

---

## 5. Part C — carrying the errors

### 5.1 `Estimate<T>` gains a standard error — reversing a recorded decision

`Estimate` ([`mod.rs:124-131`](../../../../src/parameter_estimation/mod.rs#L124)) is
`{ value, provenance, observations }`, and its documentation rules this out:

> "**No uncertainty interval.** These are priors; a caller mixes them into a genotype prior rather
> than reporting them, and an interval on a prior is not a quantity any consumer in the design
> reads." ([`mod.rs:120-122`](../../../../src/parameter_estimation/mod.rs#L120))

**Reversed, owner, 2026-09-27.** The error has two consumers now: the stopping rules of Parts A and
B, and the person judging a fit — a number with no error cannot be told apart from a well-measured
one. `Estimate` gains `standard_error: Option<f64>`: `None` where nothing determined the value (no
information, a defaulted or supplied value, or a quantity whose error this design does not compute).
The doc comment is rewritten to say so.

### 5.2 The parameters file, version 2

**Decided: `standard_error` beside `value` wherever the file writes an estimate this design gives an
error to**, and nowhere else:

| entry ([`parameters_file.md`](parameters_file.md)) | error from |
|---|---|
| inbreeding coefficient, per sample (§3.5) | Part A — the homozygote excess |
| base-quality calibration multiplier, per read group (§3.3) | Part A's error on the ordinary-position error rate, carried through the multiplier (the rate over the minted rate) |
| substitution rate, per read group and stratum (§3.7) | the binomial error of a count, √(p(1 − p)/n) over the bases compared |
| slippage numbers, per stratum and slippage group (§3.7) | Part B — **only where the number is the stratum's own fit**, in its `*_origin` block as `own_fit_standard_error`; a number drawn from or blended with a curve gets none, because the blend's own error is not computed here |

- **The file's format version goes from 1 to 2** ([`parameters_file/mod.rs:188`](../../../../src/calling/parameters_file/mod.rs#L188)).
  Every table refuses unknown keys (`deny_unknown_fields`, [`mod.rs:406`](../../../../src/calling/parameters_file/mod.rs#L406)),
  and the value tables accept only `value`, `warrant`, `observations`
  ([`mod.rs:2173`](../../../../src/calling/parameters_file/mod.rs#L2173)). A version-2 build reads
  version-1 files (the key is optional and absent); a version-1 build refuses version-2 files with its
  existing message ("written by a newer build", [`validate.rs:129`](../../../../src/calling/parameters_file/validate.rs#L129)).
  That is acceptable: files are refitted, not archived.
- **The error is written on the value's own scale**, with the documentation saying it is a
  curvature-based error and, near a bound, not the spread of a normal distribution.
- **Size.** About 30 bytes a row where it appears. §9 of the parameters-file spec sizes the largest
  table at up to 62 MB at 3,000 samples; the substitution-rate rows grow by a fifth.
- **The file records each start's outcome** — converged, stopped at the limit, or agreed with an
  earlier start, with its passes — in the section that already says what the fit read
  (`fitted_from`, [`mod.rs:493`](../../../../src/calling/parameters_file/mod.rs#L493)).

### 5.3 How we know Part C works

The golden files `testdata/every_shape.toml` and `every_shape_as_written.toml`
([`mod.rs:1857`](../../../../src/calling/parameters_file/mod.rs#L1857),
[`to_toml.rs:1465-1474`](../../../../src/calling/parameters_file/to_toml.rs#L1465)) gain the new
keys; a version-1 file still reads; a version-2 file round-trips; an absent error stays absent through
a write and a read.

---

## 6. What moves, and what the checksums do

- **Parts A and B change fitted numbers**, by less than a tenth of their errors by construction on
  converged fits, and more where the old rule stopped at its limit. **The cross-platform checksums
  move** and are re-recorded under their own policy: only after the change has been measured on the
  real cohort with `scripts/promote_ng_oracle.sh` and the difference explained
  ([`cross_platform_digests.rs`](../../../../src/cli/cross_platform_digests.rs) module doc).
- **Part B's subset does not move the checksums**: the fixture has two samples.
- **Part C moves the parameters file's bytes** (the new key and the version) without moving a fitted
  number; the calls do not move from Part C alone.
- **Memory.** Part A adds a few kilobytes a sample; Part B's subset lowers the repeat-tract fit's peak.
- **Time.** On kimura: Part A should cut the SNP/indel fit from 9 h 25 min to roughly a third, if the
  first start converges near pass 100 and the others stop on agreement; Part B's subset should cut a
  2,169-sample stratum by roughly the ratio of samples read, 256/2,169 at the first size. Both are
  estimates from the kimura log, not measurements. **Part B's stopping rule (§4.3) costs more than the
  rule it replaced**, measured on the four-accession oracle cohort: 394 rounds against 217, and the
  repeat-tract fit 5 min 13 s against 1 min 45 s, a ratio the host's load blurs (the unchanged SNP/indel
  half took 8 min 24 s against 14 min 54 s in the same two runs).

---

## 7. Deferred, with a recommended home

- **Weighting the repeat-tract curves by Part B's errors** instead of `1/√(slipped reads)`. Home: an
  amendment to [`str_slippage_level_curve.md`](str_slippage_level_curve.md) §4.2, after Part B's
  errors are validated; it moves every curve and every blended number.
- **An error for a blended or curve-derived slippage number.** Home: the same amendment.
- **Calling reading the errors** — for instance, widening a prior whose value is poorly known. Home:
  [`read_likelihoods.md`](read_likelihoods.md) / [`calling_priors.md`](calling_priors.md), as its own
  design.
- **Running the three SNP/indel starts in parallel.** Measured cost: one start holds about 1 GiB of
  its own beside the shared evidence (7.1 GiB resident when the evidence was read, 7.9–8.2 GiB during
  the starts, kimura log), so three at once add about 2 GiB. It helps only if all three settle early,
  since each pass already spreads over the cores. Home: this document's successor, once §3.6 item 4
  has measured the passes a start takes.
- **The repeat-tract guard** drops a tract when any one sample's read-group section has one read in
  ten off by a partial repeat unit (`guard_is_over_threshold`, [`census.rs:781-794`](../../../../src/parameter_estimation/joint/census.rs#L781),
  applied at [`ssr_fit.rs:2937-2940`](../../../../src/parameter_estimation/joint/ssr_fit.rs#L2937)).
  A section with one such read and no other off-length read trips it, so the share of tracts dropped
  grows with the number of samples: 45,432 of 86,688 on kimura. The records spec does not say whether
  the threshold is per section or pooled ([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md)
  §3.3). Home: its own investigation and an amendment to that section.

---

## 8. Open questions

1. **A parameter with no information is written with warrant `defaulted` and the default value**,
   not `fitted_here` and its starting value — **decided, owner, 2026-09-27.** A value the data never
   touched is a default, and calling it fitted is the confident-and-wrong case the project rules out.
   It changes what the file says for read groups and samples with no reads, not what any other gets.
   Built in the plan's step E2.
2. **The soft constants** — `SETTLED_FRACTION` (0.1), `MAX_CONTRACTION` (0.95),
   `AGREEMENT_FRACTION` (0.5), `ERROR_REFRESH_CYCLES` (10), `FIRST_SUBSET` (256),
   `LEVEL_RELATIVE_ERROR_TARGET` (0.02), `MIN_SAMPLES_A_GROUP` (8). All starting values, none
   measured. **Settled by:** the plan's checkpoints A and D.
   *Amended at plan step B1 (2026-09-29): `SETTLED_FRACTION` was kept at 0.1 at checkpoint A;
   `MAX_CONTRACTION` and `ERROR_REFRESH_CYCLES` belonged to the projection §2's amendment replaced and
   are no longer used.*
   *Decided 2026-10-16 (owner), changed to 0.01 on 2026-10-25 (§2):* `DECREMENT_THRESHOLD` (§2) = 0.001, measured rather than a
   starting value: the only threshold measured that leaves every determined quantity within
   `SETTLED_FRACTION` of its error on every input tried (1, 4 and 59 tomato samples). It replaces
   `SETTLED_FRACTION` as what stops the SNP/indel fit; `SETTLED_FRACTION` keeps only its reporting
   role (§3.5).
3. **Is the outer-product estimator close enough to the observed information here?** — OPEN.
   *Leaning:* yes for the well-determined parameters, possibly not for the Beta shapes at three reads.
   **Settled by:** §3.6 items 2 and 3; if it is not, Louis's method replaces it for the cohort-level
   block only (eight parameters, affordable).
   *Evidence, 2026-10-16 — a finding for later, not acted on by the 2026-10-16 changes*
   ([`../../reports/implementations/snp_fit_stopping_rule_2026-10-16.md`](../../reports/implementations/snp_fit_stopping_rule_2026-10-16.md)
   §3.2): near the maximum every whole Newton step solved with the outer-product information is
   accepted, yet the log-likelihood still to gain shrinks only by a factor of 0.4 to 0.6 a pass on 1
   and 4 samples and 0.68 on 59, where an exact Newton step would shrink it quadratically. So along
   the slowest direction the outer-product information looks three to six times sharper than the
   likelihood's own curvature. If the observed information is the right one there, the reported
   errors along that direction are a half to two-fifths of what they should be. Two things already
   rest on it: the decrement understates what is left (§2's traps; the threshold was chosen from the
   measurement for that reason), and the trade cost of
   [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §2.1.3 would be lower, not
   higher, so a pair that test marks stays marked. Not measured against the observed information
   itself; §3.6 items 2 and 3, or Louis's method on the cohort's eight, would.
4. **One subset order for all strata, or one per stratum?** — decided one for all (§4.4); revisit if
   Part B's validation finds a stratum whose error the fixed subset leaves much larger than a
   per-stratum draw would.
5. **Which of a stratum's numbers the subset's precision target is set on** — the slippage level
   (§4.4), **confirmed by the owner for now, 2026-09-27**, without a measurement behind the choice. *Alternative:* the largest relative error of
   the three slippage numbers, capped so a barely determined one cannot force every sample in.
   **Settled by:** checkpoint D's comparison, which reports all three numbers' errors at each subset
   size.
