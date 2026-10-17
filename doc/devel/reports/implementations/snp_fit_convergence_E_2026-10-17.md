# The converging SNP/indel fit, checked: calls, precision and time (plan steps E2–E5)

**Date:** 2026-10-17
**Branch:** `snp-fit-convergence` at 2627af72 (the new fit), against `main` at 354754b8 (the current
fit). Plan: [`../../ng/impl_plan/snp_fit_convergence.md`](../../ng/impl_plan/snp_fit_convergence.md),
milestone E. Step E1 (the identity oracle and the cross-platform digests) is in commit 32961786.
**Machine:** macOS, Apple M5 Pro, 18 cores, 64 GB; every run alone, started with the one-minute
load below 3.

---

## 1. Summary

**The new fit is four times faster at 4 and 59 samples and lands where a much longer fit lands, but
it has one defect that changes thousands of calls at 59 samples, and that defect must be fixed before
the branch merges.**

- **The defect.** On the 59 tomato2 samples, two samples whose inbreeding coefficient the current fit
  puts at 0.999999999 (SRS1839222 and SRS1839240) are written to the parameters file as
  `inbreeding_coefficient = 0.0, warrant = "defaulted"` — scored as outcrossers. The Newton step
  lands their homozygote excess exactly on its upper bound of 1, and the code that turns an excess
  into a coefficient drops an excess of exactly 1 because a coefficient must be below 1; the sample
  then falls back to the default of 0. That code's own comment says the fit "cannot return exactly
  one" — true of the golden-section search the fit used alone, not of a Newton step clamped into
  `[0, 1]`. (§3.2)
- **What the calls look like without the defect.** With those two coefficients put back near 1 by
  hand, the 59-sample calls move by **7 records of 291,145 (3 lost, 4 added) and 52 sample genotypes
  of 17.2 million**; 51 of the 52 were called at GQ 5 or less before, the other had no GQ. QUAL moves by more than 1
  in 3,686 records. (§3.3)
- **What they look like with it**, as the branch stands: 94 records lost, 1,331 added, and 4,766
  sample genotypes changed, 4,314 of them in the two samples (mostly 1/1 → 0/1). (§3.1)
- **At 4 samples** the calls move exactly as B1 measured: 2 records added (QUAL 9.4) and 5 of 26,720
  genotypes changed, all at one locus. (§2)
- **Precision (E4).** At 59 samples every number the data determine lies within **0.047** of its
  standard error of a fit run on to a ten-times smaller stopping threshold; the plan's bar was 0.1,
  and the prediction 0.047. (§4)
- **Time (E5).** `estimate-parameters` at 59 samples: **15 min 49 s against 55 min 41 s**; 4 samples
  87 s against 366 s; 1 sample 24 s against 72 s. Peak memory is the same (410 MB against 417 MB at
  59 samples). (§5)

**Recommendation for the defect:** keep the Newton target for a homozygote excess inside what the
golden-section search can return, so the fit never hands back an excess of exactly 1. That restores
the assumption the coefficient code was written under, and on the 59 samples gives the two samples
a coefficient of 0.999999999 with a standard error, as the current fit does. The alternative —
letting the coefficient code turn an excess of 1 into the largest value below 1 — fixes this one
reader but leaves an exact bound in the fit's output for any other. Either is small; it is a change
to the fit's numbers, so it needs the oracle re-run (its 4 samples are not affected: none reaches the
bound there).

---

## 2. Calls on the 4 oracle samples (E2)

**What was compared.** The 4 oracle psps (main checkout `tmp/snp_fit_baseline/oracle_psps/`, the 20
oracle regions of 4 tomato2 accessions, about 3 to 13 reads a position) fitted twice — by `main`'s
`pvc` and by this branch's — then both parameters files called by this branch's `pvc call-from-psps`
with the psps in the same order, so only the parameters file differs.

`main`'s 4-sample file is byte-for-byte the file the branch started from (5ea720b6, the oracle's old
`fitted.parameters.toml`), and the branch's is the oracle's new one (3e60785a): the comparison is the
oracle's.

| | current fit | new fit |
|---|---|---|
| starts | 200, 199, 200 passes, all at the pass limit | 53 converged; 46 and 37 agreed with start 1 |
| records | 6,680 | 6,682 |

- **Records:** none lost; two SNVs added at SL4.0ch01:85900664 and 85900669, QUAL 9.4 each.
- **Genotypes:** 5 of 26,720 sample genotypes changed, all one sample at SL4.0ch01:30161572–30161577,
  0/1 → 1/1, GQ 3 on both sides. Those five records' QUAL fall by 6.6 to 8.2.
- **QUAL:** moves in 3,892 of the 6,680 shared records; apart from those five, by +0.1 or +0.2 (or
  −0.1 in 23 repeat-tract records).
- **What moved the calls:** the ordinary-site prior's two concentrations (reference 1.8730 → 1.8681,
  alternative total 8.096 × 10⁻⁴ → 8.167 × 10⁻⁴), the 4 inbreeding coefficients (SRS1839219
  0.9662 → 0.9681, a sixth of its standard error; the others by a twentieth of theirs or less) and the 4 error
  multipliers (at most 0.04%).

This is the size the diagnosis predicted (its §5.1: 2 records added, 5 genotypes changed).

---

## 3. Calls on the 59 tomato2 samples (E3)

**What was compared.** The 59 tomato2 psps (main checkout `tmp/tomato2/psps/`, whole-genome walks of
59 accessions at a median 5.5 reads a position) fitted by `main` and by this branch, both called by
this branch's binary, psps in the same (name) order. `main` was refitted rather than reusing
`tmp/snp_fit_baseline/t2_59.parameters.toml` of 2026-10-15: the refit gives **every number in that
file unchanged**, differing only in repeat-tract standard errors and notes (the later repeat-tract
commits D7 and D8), which calling does not read — so either file would have given the same calls.

### 3.1 As the branch stands

| | current fit | new fit |
|---|---|---|
| starts | 199, 199, 199 passes, all at the pass limit | 45 converged; 43 and 43 agreed with start 1 |
| records | 291,145 | 292,382 |

- **Records:** 94 lost (40 SNV, 40 indel, 14 repeat tract; median QUAL 289), 1,331 added (1,222 SNV,
  78 indel, 31 repeat tract; median QUAL 31, 198 below 20); 52 positions change their alleles.
- **Genotypes:** 4,766 of 17.2 million sample genotypes changed, in 4,193 records. **4,314 are in
  two samples**: SRS1839222 2,352 and SRS1839240 1,962, mostly homozygous to heterozygous (1/1 → 0/1
  2,518 times; 0/0 → 0/1 1,145 times). No other sample has more than 46.
- **Added records:** 809 of the 1,331 have no non-reference genotype outside those two samples.

### 3.2 Why: an excess on its bound becomes an outcrosser

The fit estimates, per sample, a **homozygote excess**: how much less heterozygous the sample is than
the cohort's allele frequencies predict, a number from 0 to 1 that becomes the sample's inbreeding
coefficient F. For SRS1839222 and SRS1839240 the data put it at the top: the current fit returns
0.9999999992 and 0.9999999975 (±0.010 and ±0.005).

- The current fit only ever reaches the excess through a golden-section search, which returns the
  middle of a bracket strictly inside `[0, 1]`, so it never returns exactly 1.
- The new fit also takes Newton steps, whose target is clamped into each parameter's interval
  (`Parameters::newton_target`, `fit.rs`). At 59 samples start 1's accepted Newton step puts both
  excesses at exactly 1.0, and the start converges there (the trace's returned row: `hom_excess_8`
  and `hom_excess_25`, 1.0; every pass before it shows 0.999999999).
- `fitted_inbreeding_of` (`run/census_fit.rs`) turns an excess into a coefficient, and a coefficient
  must lie in `[0, 1)`. An excess of exactly 1 has none, so the sample is left out and falls to the
  rung below: **F = 0, warrant `defaulted`**. Its comment states the assumption the Newton step
  broke: *"The shipped maximisation cannot return exactly one … so this guards the type's contract
  rather than an observed case."*

A sample that is almost fully inbred is then scored as an outcrosser: the genotype prior stops
favouring homozygotes, so its low-depth homozygous calls turn heterozygous and new records appear
where its reads alone suggest a second allele. The file says `defaulted` and the note above it
reads "nobody said how inbred this sample is", which is false.

The same happens in the run to a decrement of 10⁻⁴ (§4) and in checkpoint D's 59-sample run, so it
dates from B1, the step that added the Newton step. The 4 oracle samples and the one-sample run are
unaffected: no excess reaches 1 there.

### 3.3 With the defect removed by hand

To see what the new fit itself does to the calls, the branch's 59-sample file was edited to give the
two samples F = 0.999999999, warrant `fitted_here` (everything else untouched) and called again.

| | current fit | new fit, two coefficients restored |
|---|---|---|
| records | 291,145 | 291,146 |

- **Records:** 3 lost (SNVs at QUAL 15.3, 22.9 and 376.4) and 4 added (two indels at QUAL 245 and
  440, an SNV at 64, a repeat tract at 111); no position changes its alleles.
- **Genotypes:** 52 of 17.2 million sample genotypes changed, one each in 52 records (32 SNV, 12
  indel, 8 repeat tract); 51 of the 52 were called at GQ 5 or less before (median 3), the other
  with no GQ. No sample has
  more than 13.
- **QUAL:** moves in 132,476 of 291,142 shared records, by more than 1 in 3,686 and by more than 5
  in one (an indel, +11.6). Among those that move, half by 0.1 or less; 98 in 100 by between −0.1 and
  +1.4.
- **What moved the calls:** the ordinary-site prior's concentrations (reference 3.5545 → 3.4810,
  −2.1%; alternative total 2.711 × 10⁻³ → 2.669 × 10⁻³, −1.6%), which make a non-reference allele
  slightly less expected a priori — most QUAL changes are this, and in SNVs QUAL falls in 78,210 of
  the 116,563 that move; the 59 error multipliers by at most 0.05% (largest move 0.04 of a standard
  error); the 57 other coefficients by at most 0.05 of a standard error. Two coefficients go from
  3 × 10⁻¹⁰ to exactly 0 at the lower bound, which is allowed and changes nothing.

---

## 4. The stopping rule loses nothing at 59 samples (E4)

**What was compared.** The branch's 59-sample fit (stops when the remaining Newton step is worth less
than 0.001 log-likelihood units) against the same build with that threshold set to 10⁻⁴ (a scratch
build, not committed), on the same psps. The yardstick is each number's standard error in the 10⁻⁴
run. The 10⁻⁴ run's start 1 converges at 48 passes against 45, and ends 0.0019 log-likelihood units
higher.

**Largest distance, in standard errors of the 10⁻⁴ run** (from the trace's returned rows and the
parameters files):

| what | numbers | largest distance |
|---|---|---|
| mismapped share | 1 | **0.047** |
| error rates at mismapped positions | 59 | 0.018 |
| error multipliers (ordinary positions) | 59 | 0.007 |
| allele-frequency shape b | 1 | 0.004 |
| share segregating above 1%, expected heterozygosity, mean alternative frequency | 3 | 0.004 |
| share of positions duplicated in more than 1% of the panel; carrier shape b | 2 | 0.003 |
| fixed non-reference share | 1 | 0.0005 |
| inbreeding coefficients | 57 | 0.0004 |
| repeat-tract substitution rates | 5,428 | 0 |

Every determined number lies within 0.1 of its error, as the plan requires; the largest, 0.047,
is the figure the stopping-rule study predicted. The four numbers the data do not separately
determine (the invariant share and the allele-frequency shape a; the duplicated share and the
carrier shape a) carry no error and are not in the table; both runs mark the same pairs.

Against the prototype's reference run to 10⁻⁴ (`tmp/snp_stop/runs/t59_d1e-4`, made before the box
solve and error-rate fixes of B2–B4), the largest distance is the same 0.047, on the same number.

**The verdicts and tables at 59 samples match §2.1.3's example** of
`parameter_prepass_joint_fit.md`: both pairs marked not separately determined; the share segregating
above 1% 8.447 × 10⁻³ ± 2.04 × 10⁻⁴, expected heterozygosity 1.1894 × 10⁻³ ± 3.03 × 10⁻⁵, mean
alternative frequency 7.660 × 10⁻⁴ ± 2.02 × 10⁻⁵, the share duplicated in more than 1% of the panel
1.212 × 10⁻³ ± 3.38 × 10⁻⁵; halving or doubling the segregating share costs 1.35 log-likelihood
units, the duplicated share 0.116.

---

## 5. Time and memory (E5)

Each run alone, `/usr/bin/time -l`, 18 threads (`--threads 18` for `main`, all 18 cores by default on
the branch). The branch's runs also wrote the fit's trace (one row a parameter a pass; a few hundred
kilobytes). A macOS background process held about one core throughout.

| samples | build | passes by start | SNP/indel stage | `estimate-parameters` wall time | peak resident memory |
|---|---|---|---|---|---|
| 1 (SRS1839219) | main | 144 converged, 24 converged, 61 agreed | 1 min 05 s | 72 s | 131 MB |
| 1 | branch | 22, 20, 21, all converged | 17 s | 24 s | 127 MB |
| 4 (oracle) | main | 200, 199, 200, all at the limit | 5 min 59 s | 366 s | 203 MB |
| 4 | branch | 53 converged; 46, 37 agreed | 1 min 20 s | 87 s | 203 MB |
| 59 (tomato2) | main | 199, 199, 199, all at the limit | 55 min 17 s | 55 min 41 s | 417 MB |
| 59 | branch | 45 converged; 43, 43 agreed | 15 min 27 s | 15 min 49 s | 410 MB |
| 59 | branch, threshold 10⁻⁴ | 48 converged; 43, 43 agreed | 15 min 58 s | 16 min 21 s | 452 MB |

`main`'s 59-sample run took 55 min 41 s against 61 min on 2026-10-15; the branch's matches B5's
measurement (942 s, here 949 s, the same passes). Calling the 59 samples takes 33 s with either
file.

---

## 6. What was run, for reproduction

Scripts and outputs under the worktree's `tmp/e/`: `fit.sh` and `call.sh` (one run each),
`queue.sh` (the runs in order, each waiting for a quiet machine), `compare_calls.py` (records,
genotypes and QUAL between two VCFs), `gt_by_sample.py`, `e4.py` (distances in standard errors),
`filediff.py`; the runs in `tmp/e/runs/`. `main`'s binary was built from `git archive main` in
`tmp/e/main_src`; the 10⁻⁴ build from `git archive HEAD` with `DECREMENT_THRESHOLD` edited, in
`tmp/e/thr_src`.
