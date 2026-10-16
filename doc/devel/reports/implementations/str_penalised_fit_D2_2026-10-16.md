# D2: the production fit against the counting fit on tomato2's real reads, held out by whole locus (2026-10-16)

Plan step D2 of [`str_penalised_fit.md`](../../ng/impl_plan/str_penalised_fit.md) (spec §4 item 4):
the first test of the new repeat-tract fit on real reads. For the owner's decision at checkpoint D,
whether the new fit passes on real data before the old code is deleted.

## 1. The answer

**On 59 tomato accessions at about three reads a position, the production fit predicts held-out loci
better than the counting fit**: by 0.136 log-likelihood units a locus over the 2,349 loci both fits
scored (± 0.053, counting loci of one class as one unit), and by 0.213 (± 0.062) against the counting
fit without calling's stutter curves, over 2,231 loci. Both margins are 2.6 and 3.5 of their
standard errors above zero.

The margin comes from the two- and three-base classes. At one base, where classes are large and each
fit has plenty of reads, the two fits predict about equally (+0.029 ± 0.10). The production fit
also scores loci the counting fit leaves unfitted: 24 against the counting fit with curves and 142
against the one without, in classes with fewer than 8 loci with reads.

| motif length | loci both scored | classes | production fit, a locus | counting fit with curves, a locus | difference | ± by class |
|---|---|---|---|---|---|---|
| 1 | 986 | 19 | −77.552 | −77.581 | +0.029 | 0.102 |
| 2 | 1,023 | 35 | −46.160 | −46.351 | +0.191 | 0.049 |
| 3 | 304 | 18 | −18.490 | −18.792 | +0.302 | 0.197 |
| 4 | 15 | 1 | −9.319 | −9.676 | +0.357 | — |
| 6 | 21 | 1 | −22.212 | −22.186 | −0.026 | — |
| **all** | **2,349** | **74** | **−55.307** | **−55.443** | **+0.136** | **0.053** |

Against the counting fit without curves the per-motif margins are +0.074 (one base), +0.309 (two),
+0.381 (three); full table in `tmp/d2a/tomato2_59.paired.md`.

## 2. What the numbers are

- **Held-out score**: each class's loci are split into five folds of whole loci; each fold is
  scored by a fit made on the other four, as the log-likelihood of its loci's reads under the fitted
  numbers. A real cohort has no truth, so unlike the bench's drawn cells these are raw
  log-likelihoods, comparable only between fits on the same loci; hence the pairing by locus. Larger
  (less negative) is better.
- **"± by class"** is the standard error of the mean paired difference with loci of one class
  grouped, because they share their fitted numbers; treating loci as independent gives about 30%
  smaller errors (0.041 and 0.049).
- 2,392 loci were picked for scoring (at most 100 a class); 19 were scored by neither fit because
  their class had no reads in the other folds.

## 3. How it was run

- The bench's psp mode (step D2a, commit `830049c1`): the evidence comes from `estimate-parameters`'
  own reading path, the repeat-tract guard and the cohort's fixed sample order; on the oracle's four
  psps its production fit reproduces the oracle's fitted slip numbers bit for bit (D2a's review).
- Data: the 59 tomato2 psps of plan's preconditions (`../pop_var_caller-str-guard/tmp/tomato2/psps`,
  format 2.0), with each sample's homozygote excess read from an `estimate-parameters` fit of the same
  psps made after C10 (`--parameters tmp/perf_ep/t2_59.parameters.toml`; the SNP/indel fit takes an
  hour on these 59, and the file's coefficients equal the fit's on the 4-psp check).
- `--methods B,D_tract_prior`, five folds, at most 100 held-out loci a class, the default fold seed;
  8 threads: 98 s, peak 416 MB.
  Outputs: `tmp/d2a/tomato2_59.{tsv,loci.tsv,log}`; pairing: `tmp/d2a/heldout_real.py`.

## 4. Limits

- One cohort, one depth (about three reads a position), one fold seed.
- The homozygote excess is computed on the whole cohort before folding, as production does; both fits
  get the same values, so it favours neither (D2a's review).
