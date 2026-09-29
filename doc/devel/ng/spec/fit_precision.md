# Fitting to the precision the data support

**Status:** design, 2026-09-27; amended at checkpoints A and A′ (§1.3, §3.2; 2026-09-29) and at plan step B1 (§2, §3.3, §3.4; 2026-09-29). Build order:
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
with reads at every tract, [`ssr_fit.rs`](../../../../src/parameter_estimation/joint/ssr_fit.rs)
`ln_tract`), while a stratum's slippage is estimated from thousands of tracts and does not need
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
  ([`ssr_fit.rs:1048-1100`](../../../../src/parameter_estimation/joint/ssr_fit.rs#L1048)).

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
  1,000 passes; with the box both converge.
- **A parameter is settled when |dⱼ| < `SETTLED_FRACTION` × SEⱼ**, both from the same pass; one with
  no error, or no distance (not solved for), is settled by definition. **The fit has converged when
  every parameter is settled**, and stops at the end of that cycle.

Measured on the same 15 cohorts, before the box: every fit this rule calls converged is within 0.088 of
an error of the 600-pass fit; the ones still crawling run to the pass limit and say so (all three two-sample
cohorts and two of the three 4-sample ones with the duplicated class). On the oracle no start converges within 200 passes, and none
claims to. Part B's climb (§4.3) keeps its projection of the log-likelihood gain.

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
5. **Determinism.** The same bits at 1, 4 and 8 threads; the new sums follow the fit's fixed chunk
   order ([`fit.rs:2235-2352`](../../../../src/parameter_estimation/joint/fit.rs#L2235)).

---

## 4. Part B — the repeat-tract fit

### 4.1 What one stratum's fit estimates, for reference

Per stratum ([`ssr_fit.rs:803-809`](../../../../src/parameter_estimation/joint/ssr_fit.rs#L803)):
three slippage numbers per slippage group (level, shorter share, fall-off), a length spectrum over 13
allele classes, and a concentration. Production pools every read group into one slippage group
(`every_read_group_pooled`, [`census_fit.rs:227-234`](../../../../src/run/census_fit.rs#L227)), so
17 numbers a stratum. The score is the **mean log-likelihood a tract**
([`ssr_fit.rs:2457`](../../../../src/parameter_estimation/joint/ssr_fit.rs#L2457)); the climb is
coordinate ascent by golden section, 18 evaluations a coordinate, up to 5 rounds from each of 3
starts, stopping when a round gains less than 10⁻⁶ on that mean
([`ssr_fit.rs:895`](../../../../src/parameter_estimation/joint/ssr_fit.rs#L895)).

**Trap: the mean is over every tract in the stratum, including those with no reads**, and a tract's
log-likelihood sums over its samples. So the same 10⁻⁶ is a looser test in a stratum where most
tracts are empty and a stricter one as samples are added. That is why the climb's rule changes too.

### 4.2 Standard errors for one stratum

**Decided: the numerical curvature of the stratum's total log-likelihood at its answer**, by central
differences over all its numbers on the scales the climb uses (logit for the slippage numbers and
shares, log for the concentration, log-ratios for the spectrum), inverted to give the errors, then
carried to the natural scale.

- **Cost:** about 2p² evaluations for p numbers — 578 at p = 17, against about 4,600 for the climb
  from three starts: roughly an eighth more, once per stratum at its final answer (not per start).
- **Why all the numbers and not only the slippage three:** a read off the reference length is either
  a slip or a real allele, so the level and the spectrum trade against each other
  ([`ssr_fit.rs`](../../../../src/parameter_estimation/joint/ssr_fit.rs) `PeriodLengthSpectrum`
  documentation). Holding the spectrum fixed would understate the level's error.
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
  ([`ssr_fit.rs:895`](../../../../src/parameter_estimation/joint/ssr_fit.rs#L895)); golden section
  never re-evaluates its starting value, so it can move a coordinate to a worse one. Under the new
  rule a round that loses is **not** convergence: the round's worse moves are undone and the walk
  stops at the best point it held, recorded as such.
- `max_rounds` stays as a backstop.

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
  by the sample's index as today ([`ssr_fit.rs:2161-2165`](../../../../src/parameter_estimation/joint/ssr_fit.rs#L2161));
  dropping rows from a stratum's evidence does not disturb it. The substitution counts are still
  summed over every sample when the evidence is gathered, since they cost nothing to keep.
- **Memory falls too:** a stratum's likelihood table has one row per sample with reads a tract, so
  it shrinks in proportion (5.9 GiB for the largest stratum on kimura at 2,169 samples).
- **Trap: the refusal floor is judged on the whole stratum**, not the subset: a stratum with reads in
  only the samples outside the subset must grow the subset rather than be refused. Refusal for too few
  tracts ([`ssr_fit.rs:1135-1159`](../../../../src/parameter_estimation/joint/ssr_fit.rs#L1135)) looks
  at every sample's tracts before any subset is drawn.
- **The target is on the level, not on all three numbers.** One stratum can measure one of its
  numbers well and another to about 45% of itself ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md)
  §4.5, line 738), so a 2% target on every number would read every sample in every stratum.
  **Which number carries the target is soft** — decided here, and question 5 below.

**This is the §7 proposal of [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md)**
("bound the parameters fit's resident set by a subsample of *samples*, chosen the same data-blind way
the loci are", lines 1627-1634) applied to the repeat-tract half, with its question 8 — how large the
subsample must be — answered per stratum by the precision target rather than once for the cohort.

### 4.5 How we know Part B works

1. **The curvature errors mean what they say.** Draw strata at known slippage (the
   `bench_fixtures::draw_stratum` generator), fit each many times, and compare the spread of the
   fitted numbers with the reported errors, at 3 and at 30 reads.
2. **The new stop loses nothing.** On drawn strata and on the 4-accession oracle cohort, the fitted
   numbers under the new rule are within `SETTLED_FRACTION` of their errors of a climb run to 20
   rounds.
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
  estimates from the kimura log, not measurements.

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
3. **Is the outer-product estimator close enough to the observed information here?** — OPEN.
   *Leaning:* yes for the well-determined parameters, possibly not for the Beta shapes at three reads.
   **Settled by:** §3.6 items 2 and 3; if it is not, Louis's method replaces it for the cohort-level
   block only (eight parameters, affordable).
4. **One subset order for all strata, or one per stratum?** — decided one for all (§4.4); revisit if
   Part B's validation finds a stratum whose error the fixed subset leaves much larger than a
   per-stratum draw would.
5. **Which of a stratum's numbers the subset's precision target is set on** — the slippage level
   (§4.4), **confirmed by the owner for now, 2026-09-27**, without a measurement behind the choice. *Alternative:* the largest relative error of
   the three slippage numbers, capped so a barely determined one cannot force every sample in.
   **Settled by:** checkpoint D's comparison, which reports all three numbers' errors at each subset
   size.
