# Genotypes without reads: the `GQ` floor lowered from 10 to 9

*2026-11-08. Branches `fill-gq8`, then `fill-gq9`. Measured on rick, refine_sim cohorts, pvc aa39bbbf with refinement on.*

## The result

The refinement writes a genotype for a sample with no reads at a site only if the call's `GQ` is at
least a floor; below it the sample stays `./.`. The floor was 10. **It is now 9** — first 8, then
9 the same day, the conservative choice (see *8 or 9* below). In every cohort of 100 samples or more
measured, a floor of 8 cuts the SNP genotypes left missing in the shallowest samples by a fifth to a
third, and adds 1.5% to 4% to their SNP errors; 9 keeps about half of that gain.

The owner's question was whether pvc could leave fewer genotypes missing — Beagle leaves none, at
many more errors — without many more false positives and negatives. The owner set the target range:
the refinement is meant for cohorts of 100 samples and more; below that, turning it on is the
user's call (2026-11-08).

## How it was measured

A withheld genotype is written `./.` but keeps the pass's posterior as `GP`. So no rerun is needed
to see what a lower floor would do: every withheld genotype was filled with its most probable
genotype from `GP`, the filled calls were scored with pvc_eval exactly as the suite scores them, and
each filled genotype's outcome was set beside the `GQ` it would have been written with
(`-10 log10(1 − largest GP)`). One scoring run per cohort gives the outcome at every floor.

The cohorts, SNPs only (as in the user accuracy page's imputation section):

- the 400-sample refine_sim populations (100 samples each at 30×, 10×, 5×, 2×), outbred and selfing
  to F = 0.9 — the counts below are their 2× group;
- 100-sample cohorts cut from them — all at 2×, or all at 5× — each fitted on its own samples, as a
  user with only those samples would run pvc;
- a new F2 of 100 plants at 1× with its two parents at 30×, on the same genome as refine_sim's F2.

Three further cohorts — selfing with 100 samples all at 2×, and 200 samples at 2× and 5× in each
population — were stopped before their parameter fits finished, to free rick for benchmarks.

## What a floor of 8 adds

| cohort, group | SNP errors now | SNP genotypes missing now | filled at `GQ` 8–9 | of them wrong |
|---|---|---|---|---|
| outbred, 100 samples all at 2× | 7,420 | 9,734 | 2,107 | 129 |
| outbred, 100 samples all at 5× | 1,714 | 282 | 81 | 4 |
| F2, 100 plants at 1× | 1,473 | 7,324 | 1,550 | 61 |
| outbred 400, its 2× group | 1,921 | 3,554 | 1,098 | 28 |
| selfing 400, its 2× group | 643 | 769 | 93 | 19 |

*Errors are missed variants plus false variants, with and without reads summed.*

Fills at `GQ` 8–9 were wrong 3 to 6 times in 100 in the outbred cohorts and the F2. The selfing
population is the exception, 20 in 100, but it has few such fills: 19 more errors on 643.

**Why not lower.** At `GQ` 6–7 the fills were wrong 6 to 9 times in 100 in the outbred cohorts and
the F2, and 25 in 100 in the selfing population (27 of 109). At `GQ` 5, 11 to 18 in 100; below 5,
about 30 in 100.

The same `GQ` means a different error rate in different cohorts: a genotype written with `GQ` 8
should be wrong about 16 times in 100, and the outbred 400 cohort's were wrong about 3 in 100. So the
floor is set from these measurements, not from what the `GQ` says.

## What the change does not touch, and what it suggests next

**Two other rules also hold genotypes back**: the sample must have reads on both sides of the site
within one segment length, and must probably be copying some sample's haplotype rather than the
cohort's frequencies (§4.2 of `neighbour_evidence_pass.md`). In the 100-plant F2 they hold back 1,761
genotypes that already reach `GQ` 10, of which only 6 would be wrong; in the selfing 400 cohort's 2×
group, 371 with 23 wrong; in the outbred cohorts, under 100. Which of the two rules holds back the
F2's is not known — the VCF does not record it. Finding out needs a diagnostic build, and could cut
the F2's missing genotypes by about another quarter for almost no errors.

**Unchanged by the floor**: no phased allele moved. In the unit test of a phased window, four more
genotypes were filled, joining blocks, so 46 phased columns took a different phase-set id; that
test's recorded hash was re-recorded.

**Not yet rerun**: the refine_sim benchmark itself and the accuracy page's tables, which still show
the floor at 10.

## Files

Scripts: `tmp/fill_cutoff/` (not committed) — `fill.py` fills the `./.` genotypes from `GP`, `join.py`
sets each filled genotype's outcome beside its `GQ`, `sub.sh` builds and calls the 100-sample
cohorts, `summary.py` prints the tables above. Results on rick under
`~/pop_var_caller/tmp/fill_cutoff/`.

## 8 or 9

8 was merged first; the owner then chose 9, the more conservative floor (2026-11-08). The two bands
apart, SNP genotypes of the shallowest group:

| cohort, group | filled at `GQ` 8 | wrong | filled at `GQ` 9 | wrong |
|---|---|---|---|---|
| F2, 100 plants at 1× | 1,041 | 28 | 509 | 33 |
| outbred, 100 samples all at 2× | 1,088 | 83 | 1,019 | 46 |
| outbred 400, its 2× group | 676 | 20 | 422 | 8 |
| selfing 400, its 2× group | 49 | 15 | 44 | 4 |
| the four 5× groups together | 73 | 3 | 65 | 4 |

A floor of 9 gives up about 2,900 fills, 19 in 20 of them right, to avoid about 150 errors — chiefly
the selfing population's, where a fill at `GQ` 8 was wrong 31 times in 100.

## Which of the other two rules holds genotypes back

A diagnostic build (branch `fill-reason-diag`, `PVC_DIAG_FILL_REASON`) printed each rule's verdict
apart for every genotype without reads, on the 100-plant F2 (SNPs, floor 8):

| held back by | genotypes | of them wrong |
|---|---|---|
| the background rule alone | 2,419 | 34 |
| `GQ` below the floor alone | 1,822 | 313 |
| both | 2,178 | 395 |
| (filled) | 32,946 | 146 |

The rule wanting reads on both sides holds back none. In the F2 the background value does not
predict a wrong fill: of those held back by it alone, 388 between 0.1 and 0.2 had 22 wrong, 1,791
between 0.2 and 0.5 had 12, and 240 above 0.5 none. Not yet measured in the outbred and selfing
populations, so the limit is unchanged.
