# Why the SNP/indel fit of `estimate-parameters` never converges, and what fixes it

**Date:** 2026-10-16
**Branch:** `snp-fit-convergence` (from `main` 51b7be99). Nothing is committed; the prototype code is
reverted from `src/` and kept as a patch under `tmp/snp_conv/` (§8).
**Question:** on 59 tomato2 samples the SNP/indel half of `estimate-parameters` took 60 of the
command's 61 minutes, because each of its three starting points ran to the 200-pass limit without
converging (performance review of 2026-10-15, finding L5). Why, and which of the proposed remedies
works?

---

## 1. Summary

**The cause.** The fit climbs the likelihood by alternating a pass over the data with an update of
every parameter (an expectation–maximisation loop, accelerated by SQUAREM, a method that extrapolates
along the direction the last two updates took). That loop is slow wherever the data barely tell two
parameters apart, and this fit has such a pair at its centre: **the share of invariant positions and
the first shape of the allele-frequency Beta**. A position segregating at a frequency of one in a
million looks exactly like an invariant one in any cohort we can sequence, so the fit can trade
"more invariant positions" for "more rare-allele mass in the Beta" at almost no change in likelihood.
On 59 samples it spends its passes walking along that trade: the invariant share falls from 0.98
towards 0.84 while the shape falls towards its lower bound of 0.02, and the likelihood rises by
fractions of a unit per pass. The stop rule asks every parameter to be within a tenth of its
standard error of the maximum, so it never fires before pass 200. At 4 samples a second pair of the
same kind holds it back too (§3.2).

**The fix that works.** Take the Newton step the fit already computes to judge convergence, keep it
only when it raises the likelihood, and halve it until it does (finding L5b, with a step-halving guard
the review did not have). Measured:

| input | passes, current fit | passes, Newton step | wall time, current → Newton | log-likelihood gained |
|---|---|---|---|---|
| 59 tomato2 samples, all positions, 18 threads | 597, no start converged | 133, two starts converged, one agreed | 60 min → 15 min 51 s (SNP/indel stage) | +37 units |
| 4 oracle samples, all positions, 4 threads | 599, no start converged | 133, one converged, two agreed | 18 min 55 s → 3 min 48 s (SNP/indel stage) | +2.3 units |
| 1 oracle sample, 4 threads | 229, two converged, one agreed | 61, all three converged | 192 s → 58 s (whole command) | none (−0.01) |

**What it changes in the calls** on the 4 oracle samples: 2 records added (QUAL 9.4) to 6,680, and 5
sample genotypes changed out of 26,720, all at one locus (§5). The calling prior moves by 2% (§5.2).

**What it exposes.** Run to convergence, the fit puts the Beta's first shape **on its lower bound**,
0.02, and the invariant share wherever that bound leaves it: 0.84 on 59 samples, against 0.94 that the
current fit reported after 200 passes. The number the data determine is the share of positions
segregating at a frequency above about 1%, and that does not move: 0.0084 in both fits. So the
invariant share and the shape `a` are not separately measured on this cohort, and reporting them as
fitted numbers with standard errors is misleading (§6). That is a modelling decision, separate from
the speed fix.

**Recommendation.** Adopt the Newton step with step halving (§4.1). Do not build the other two
remedies the review proposed (§4.2, §4.3): once starts converge, the existing agreement rule already
stops the later ones, and judging less often conflicts with the Newton step, which needs the
information every cycle. Decide separately what to do about the invariant-share ridge (§6); my
recommendation there is to keep the model and report the identified quantity beside the two
unidentified ones.

---

## 2. Inputs and how the fits were run

| name used below | what it is | positions | where |
|---|---|---|---|
| **oracle, 4 samples** | the identity oracle's four tomato2 psps | 1,998,212 | `pop_var_caller/tmp/snp_fit_baseline/oracle_psps/` |
| **oracle, 4 samples, one chunk in 8** | the same psps with the fit's passes reading every 8th chunk of positions (16 of 128 chunks, contiguous stretches of genome) | about 250,000 | same, prototype setting |
| **tomato2, 59 samples** | the full cohort | about 2,000,000 | `pop_var_caller/tmp/tomato2/psps/` |
| **tomato2, 59 samples, one chunk in 4** | the same, every 4th chunk | about 500,000 | same, prototype setting |
| **1 sample** | `SRS1839219.psp` from the oracle, alone | 1,998,212 | oracle directory |

Every run is `pvc estimate-parameters` built from this branch with the prototype patch, whose
settings are all off by default; with none set the binary computes what `main` computes. The
per-pass parameter trace (`PVC_JOINT_FIT_TRACE`) was on in every run. Contamination was skipped on
the reduced inputs only, because it reads every position. The 59-sample baseline is the existing run
of 2026-10-15 (`t2_59.err`), not a rerun.

**How representative the small inputs are.** The 4-sample input with one chunk in 8 shows both slow
directions the full oracle shows and ends at nearly the same values (shape `a` 17.8 against 15.4, `b`
at its bound of 50 in both); a fit on it takes about 50 seconds a start, against 6 minutes. **For 59
samples, one chunk in 16 is not representative**: the density's three numbers come out "not
identified" there, the stop rule ignores parameters without an error, and the fit reports convergence
after 75 passes with the invariant share nowhere near where the full data put it. One chunk in 4 shows
the full cohort's behaviour.

---

## 3. What happens over the passes

### 3.1 On 59 samples: a walk along the invariant-share ridge

Start 1 on 59 samples, one chunk in 4, current fit, 400 passes (`runs/t59s4_base_start1_p400.trace.tsv`):

| pass | log-likelihood | invariant share | Beta shape `a` | Beta shape `b` |
|---|---|---|---|---|
| 3 | −2,024,892 | 0.983 | 0.278 | 4.28 |
| 26 | −2,020,514 | 0.973 | 0.149 | 3.95 |
| 51 | −2,020,491 | 0.965 | 0.104 | 3.70 |
| 101 | −2,020,475 | 0.946 | 0.063 | 3.46 |
| 126 | −2,020,474 | 0.943 | 0.058 | 3.44 |
| 400 | −2,020,467 | 0.916 | 0.039 | 3.33 |
| Newton step, converged at pass 48 | −2,020,463 | 0.845 | 0.020 (bound) | 3.23 |

From pass 26 on the log-likelihood gains 47 units over 374 passes, while the invariant share and `a`
keep moving together; at pass 400 the run's own judgement names the invariant share as the parameter
furthest from settled, 1.05 of its errors away, with 8 of 185 parameters not yet within a tenth of an
error. The full cohort's log (2026-10-15) says the same thing in its own terms: from
pass 73 to pass 199 of each start the parameter furthest from settled is the invariant share, its
distance falling from 8.5 to 4.1 of its standard errors.

**What stays fixed along the walk** (computed from the fitted values, `tmp/snp_conv/ridge.py`):

| fit | invariant share | `a` | share of positions segregating at all | share segregating above 1% | expected heterozygosity |
|---|---|---|---|---|---|
| 59 samples, current fit, pass 199 | 0.939 | 0.056 | 0.061 | 0.0084 | 0.00119 |
| 59 samples, Newton step, converged | 0.842 | 0.020 | 0.158 | 0.0084 | 0.00119 |

The share of positions "segregating" goes from 6% to 16%, but the additional 10% are segregating at
frequencies no read set can see. What the reads constrain — how many positions carry an allele at
more than 1%, and the heterozygosity — is the same to three figures.

### 3.2 On 4 samples: two slow directions, and a third parameter at a bound

Start 1, oracle with one chunk in 8 (`runs/s8_base_start1.trace.tsv`):

- **The Beta concentrates.** `a` and `b` climb from (0.2, 2.8) at pass 25 to (18, 50) at pass 170,
  where `b` meets its upper bound of 50 and stays. A Beta(18, 50) is a narrow peak near frequency
  0.26. Four inbred accessions show a segregating allele in one, two or three of four homozygous
  samples, so the fit is approximating a spike at a quarter. It gained 16 log-likelihood units on this
  walk, at about a tenth of a unit a pass.
- **The mismapped class shrinks and gets noisier.** Its share falls from 0.0108 (pass 25) to 0.0021
  (pass 200) while its error rate rises from 0.017 to 0.048 and the duplicated share rises from
  0.0006 to 0.0023. Still moving at pass 200; the same run continued converges at pass 215 with the
  share at 0.0018.
- After pass 109 the progress line names the mismapped share as furthest from settled, at 0.3 to 1.5
  of its errors. `a` and `b` are never named: from pass 170 on, both come out "no standard error"
  (wider than their range, or at a bound), and a parameter without an error is not judged.

### 3.3 SQUAREM is working, not failing

The accelerated jump is accepted in almost every cycle and often reaches far: on the 4-sample input
it reaches 16 to 119 times the length of one plain step in a third of the cycles after pass 37
(`runs/s8_base_start1.log`, lines `PROTO … squarem`). Long jumps that hold are the signature of a
direction the plain loop barely moves along. It is still slow because each cycle costs four passes
(two plain steps, the jump, one steadying step, and from about pass 15 the first pass also sums the
information), and each cycle gains a tenth of a unit to one unit.

Its guard accepts a jump that loses up to one log-likelihood unit (`JUMP_SLACK`). That happens in 5 of
the 64 cycles of 200 passes on the 4-sample input (losses of 0.06 to 0.69 units), so the fit is not strictly
climbing, but it is not oscillating: the next cycle recovers the loss.

### 3.4 Do the three starts reach the same maximum?

Yes, up to the precision the pass limit leaves. On 59 samples the three ended at −7,985,053,
−7,985,051 and −7,985,053; on 4 samples at −796,943.3, −796,943.3 and −796,944.5. With the Newton
step, start 1 converges and the later starts either stop by agreeing with it (the existing rule) or
converge themselves to the same log-likelihood within 1 unit (59 samples: −7,985,014, −7,985,015,
−7,985,014).

---

## 4. The remedies measured

### 4.1 L5b: take the Newton step the judging pass already computes

**What it does.** Once the fit judges convergence, the first pass of every cycle sums, over positions,
each parameter's slope and the products of slopes — the information the standard errors come from.
Solved against each other they give a step to where a quadratic approximation of the likelihood
peaks, kept inside each parameter's interval (`newton_step`, `fit/standard_errors.rs:583`). Today the
step is only used to measure distance. The prototype:

1. moves the parameters by the whole step and runs the next pass there, with the information;
2. keeps the new point if its log-likelihood is higher than the old one's, and that pass becomes the
   next cycle's judging pass, so an accepted step costs one pass;
3. otherwise halves the step and tries again, up to four times (the fraction that last worked is
   remembered: the next cycle starts at twice it, capped at the whole step);
4. if no fraction raises the likelihood, runs the usual SQUAREM cycle from the old point.

**Why the halving is needed.** The whole step overshoots badly while the fit is on the ridge: on 59
samples the first full steps lost 3,556, 2,901 and 2,299 log-likelihood units; a quarter of the step
gained 29, 11 and 4 units. The ridge is curved, and the information summed from slopes (the "outer
product of gradients") describes it as straight. Without halving the prototype rejected every step
on the 4-sample input and gained nothing. Close to the maximum the whole step is accepted every time
and the likelihood's gains roughly halve each pass.

**Measured.**

| input | fit | passes (by start) | how each start ended | best log-likelihood | time |
|---|---|---|---|---|---|
| 4 samples, one chunk in 8, start 1 only | current | 215 to converge (2,000 allowed) | converged | −102,444.19 | 73 s at 3 threads (machine shared) |
| same | Newton step | 51 | converged | −102,443.84 | 34 s at 4 threads |
| 59 samples, one chunk in 4, start 1 only | current | 400 | at the limit, invariant share 1.05 errors away | −2,020,467.4 | 8,580 CPU-seconds |
| same | Newton step | 48 | converged | −2,020,462.6 | 1,566 CPU-seconds |
| oracle, 4 samples, 3 starts, 4 threads | current | 200 + 199 + 200 | all at the limit | −796,943.28 | SNP/indel stage 18 min 55 s |
| same | Newton step | 51 + 45 + 37 | converged, agreed, agreed | −796,940.98 | SNP/indel stage 3 min 48 s |
| tomato2, 59 samples, 3 starts, 18 threads | current | 199 + 199 + 199 | all at the limit | −7,985,051 | SNP/indel stage 60 min |
| same | Newton step | 43 + 43 + 47 | converged, agreed, converged | −7,985,014 | SNP/indel stage 15 min 51 s |
| 1 sample, 3 starts, 4 threads | current | 144 + 24 + 61 | converged, converged, agreed | −111,299.01 | 192 s (whole command) |
| same | Newton step | 21 + 20 + 20 | converged ×3 | −111,299.02 | 58 s (whole command) |

The two reduced 59-sample runs overlapped on the machine, so their wall times are not comparable;
CPU-seconds are given instead. The rows for whole inputs (oracle, 59 samples, 1 sample) ran alone.

**How far the parameters move** (oracle, 4 samples; values the fits returned, in units of the current
fit's standard error, `tmp/snp_conv/compare_fits.py`): of 20 parameters, 9 move by more than a tenth
of an error, and two by more than half of one:

| parameter | current fit | Newton step | move |
|---|---|---|---|
| duplicated share | 0.0021 | 0.0139 | +3.3 errors |
| duplicated class's carrier shape `a` | 0.158 | 0.020 (bound) | −0.4 errors |
| mismapped share | 0.00515 | 0.00524 | +0.3 errors |
| invariant share | 0.998352 | 0.998338 | −0.1 errors |
| Beta shapes `a`, `b` | 15.4, 50 (bound) | 15.4, 50 (bound) | no error for either |

The duplicated class is the second weakly identified corner: with four samples its share and its
carrier shape trade against each other the way the invariant share and `a` do, and the converged fit
reports the duplicated share's error as wider than its whole range. On 59 samples the duplicated share
moves from 0.0085 to 0.0154 (0.3 of its error of 0.023) with its carrier shape on the same lower bound.

On 59 samples, each sample's homozygote excess moves by at most 0.05 of its error (median 0.01), and
each read group's base-quality multiplier by at most 0.07%.

**At one sample** both fits converge to the same log-likelihood (within 0.01 units) and the frequency
density is reported as not identified, as the spec says it must be. The invariant share differs —
0.9935 against 0.9954 — which reads as 4.7 of its reported error. It is not a real disagreement: that
error is computed with `a` and `b` held where they are, because both were dropped as not identified,
so it is the error of the invariant share given the shapes, not the error of the invariant share.
§6 says what to do about that.

**Code.** About 60 lines in `maximise` (`fit.rs`), around the existing cycle: compute the step from
the information the judging pass already has, try it with halving, carry an accepted pass into the
next cycle. No new pass type, no new sums, no change to the information or standard-error code. Every
accepted point raises the log-likelihood, so the fit's no-downhill guarantee holds without the
one-unit slack SQUAREM needs. Results stay independent of the thread count, because the passes are
the existing ones.

### 4.2 L5a: let a start agree with an earlier start that hit the limit

Measured on the oracle, 4 samples: start 1 runs to 200 passes, start 2 stops after 154, start 3 after
170 (524 passes against 599), and the SNP/indel stage takes 12 min 26 s against 18 min 55 s. The
answer is start 1's, unchanged, and still not converged. **The parameters file refuses the result**:
its reader rejects a start recorded as having agreed with one that did not converge, so building this
also means changing the file's rule and the spec's §3.4 trap that forbids exactly this. With the
Newton step it is not needed: start 1 converges and the existing rule stops the others (§3.4).

### 4.3 L5c: judge convergence only every third cycle

Measured on the oracle, 4 samples: the returned parameters and the VCF called with them are identical
to the current fit's, as the review predicted for a fit that runs to the limit, and the SNP/indel stage
takes 12 min 52 s against 18 min 55 s. It removes the cost of the information pass from two cycles in
three and does nothing about the passes themselves. It conflicts with the Newton step, which needs the
information every cycle.

---

## 5. Effect on the calls

### 5.1 The VCF

`pvc call-from-psps --parameters <file>` on the 4 oracle psps (as `scripts/promote_ng_oracle.sh`
calls them), compared record by record (`tmp/snp_conv/compare_vcfs.py`, `qualdiff.py`):

| fit the calls used | records | only in this call set | sample genotypes changed | QUAL change |
|---|---|---|---|---|
| current fit | 6,680 | — | — | — |
| Newton step | 6,682 | 2 (SL4.0ch01:85,900,664 and :85,900,669, QUAL 9.4 each) | 5 of 26,720, one sample at five adjacent positions (SL4.0ch01:30,161,572–577), 0/1 → 1/1 | median 0.1, 99th percentile 0.2, largest 8.1 |
| judge every third cycle (L5c) | 6,680 | 0 | 0 | none |

The L5a fit has no calls: its parameters file was refused (§4.2).

Calls on 59 samples were not run; that comparison would take the 59-sample calling run twice.

### 5.2 The calling prior on 59 samples

The SNP/indel prior's starting point, read off the fitted density, moves by about 2%: the reference
allele's concentration from 3.555 to 3.481, the alternative alleles' total from 0.00271 to 0.00267.
That is what the identified part of the density (§3.1) predicts: the prior depends on how many
positions carry an allele at a visible frequency, which the ridge does not change.

---

## 6. The invariant share and shape `a` are not separately measured

Run to convergence, the fit sets `a` to its lower bound on 59 samples and on the 59-sample subsample,
and the invariant share to whatever that bound implies. Move the bound from 0.02 to 0.01 and the
invariant share would move again with almost no change in likelihood. Today's file reports 0.939 ±
0.009 for the invariant share; the converged fit reports 0.842 ± 0.068. Neither number describes the
tomato panel: both describe where the fit stopped. At 4 samples the duplicated class's share and its
carrier shape behave the same way, and at one sample the spec already marks all four density numbers
as not identified.

This does not need to be solved for the speed fix, and the speed fix does not cause it; it makes it
visible, because the current fit stops part-way along the ridge at a value that looks like an answer.

Options, for a decision separate from §4.1:

1. **Keep the model and report what is identified.** Write, beside the four density numbers, the
   share of positions carrying an allele above a stated frequency (for example 1/(2N), or 1%) and the
   expected heterozygosity, and say in the file that the invariant share and `a` are not separately
   determined on this cohort. A few lines in the writer; no change to any fitted value.
2. **Reparameterise the density** so that one number is the identified share and the other is a shape
   that the data may leave undetermined. The calling prior is built from the identified part, so calls
   should not move; the fit and its standard errors change shape.
3. **Raise the lower bound of `a`.** Simple, and arbitrary: the reported invariant share would then be
   set by a different arbitrary number.

**Recommendation: option 1 now, option 2 only if a consumer needs the density's numbers themselves.**
The calls depend on the identified quantity, and §5 shows them nearly unchanged by a 10-point move in
the invariant share. A second, smaller defect belongs with this decision: a parameter whose partners
were dropped as not identified gets an error computed with those partners held fixed (the one-sample
invariant share in §4.1), which understates it. It should be reported as not identified too.

---

## 7. Recommendation

1. **Adopt the Newton step with step halving and the remembered fraction (§4.1).** Measured: 3.8 to
   4.5 times fewer passes, the SNP/indel stage from 60 to 16 minutes on 59 samples and from 19 to 4 on
   the 4-sample oracle, every start converging or agreeing, a log-likelihood never lower (37 units
   higher on 59 samples), and 2 records plus one locus's genotypes changed in 6,680 calls. It works at
   one sample (3.8 times fewer passes, same likelihood). The review's merge criterion — passes down by
   a third or more, final log-likelihood no lower — is met on every input measured.
2. **Do not build L5a or L5c.** L5a is made redundant by (1) and needs a file-format change; L5c
   conflicts with (1).
3. **Decide §6 separately**; my recommendation is to report the identified share beside the density's
   four numbers.

Before merging (1), the spec's §3.3 ("How the maximum is found") and the `maximise` doc comment need
the new step described, and the stop rule's tests need a case where the step is rejected and halved.
The pass budget of 200 can stay as a backstop: no measured start needed more than 51 passes.

---

## 8. Files

All under `/Users/jose/devel/pop_var_caller-snp-fit-convergence/tmp/snp_conv/`:

- `prototype_knobs.diff` — the prototype against `main` 51b7be99, every behaviour behind an
  environment variable, off by default: `PVC_SNP_NEWTON`, `PVC_SNP_NEWTON_TRIES`,
  `PVC_SNP_NEWTON_ADAPT` (the Newton step and its halving), `PVC_SNP_AGREE_LIMIT` (L5a),
  `PVC_SNP_JUDGE_EVERY` (L5c), `PVC_SNP_SUBSAMPLE` (one chunk in k), `PVC_SNP_STARTS`,
  `PVC_SNP_MAX_PASSES`, `PVC_SNP_SLACK`, `PVC_SNP_CYCLE_LOG` (one line a cycle on stderr).
- `fit.sh`, `call.sh`, `queue1.sh`, `queue2.sh` — how each run was made.
- `runs/` — every run's log, trace, parameters file and `/usr/bin/time` output; the names in §4 and §5
  are the file stems (`o4_base`, `o4_newton`, `o4_L5a`, `o4_L5c`, `s1_base`, `s1_newton`,
  `t59_newton`, `t59s4_*`, `s8_*`) and `call_o4_*.vcf`.
- `compare_fits.py`, `ridge.py`, `compare_vcfs.py`, `qualdiff.py`, `cmp_toml.py`, `trace_table.py` —
  the comparisons quoted above.
- `pvc_proto_v4` — the prototype binary the final measurements used.
