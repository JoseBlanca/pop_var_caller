# Would phased templates help the neighbour-evidence pass? (2026-11-02)

**Question (owner, 2026-11-02).** The pass copies other samples' *genotypes*. A heterozygous
template only tells a copy "either allele, half and half", so in partly outbred material — the
owner's *S. pimpinellifolium* and *cerasiforme* accessions — its guide is blurred. Would templates
that are haplotypes do better, and can phasing at our depths produce them?

**Answer.** It depends on how many haplotypes the cohort shares, which we do not know for the
owner's material:

- From few founder haplotypes (16), with 10 of 50 samples at 30× among 3× ones, Beagle's
  haplotypes cut the pass's remaining errors by a quarter to a third, about half of what perfect
  phasing would.
- From many (64), even perfect phasing gains 15%, and Beagle gains 2%.
- With every sample at 3×, no phasing tried here helps much: Beagle phases little better than
  chance in the outbred cohorts.

The next measurement is the owner's own material, not more simulation (see the end).

## How it was measured

`phasing_ceiling.py` re-calls founder-simulator cohorts (`founder_grid.py`, 50 samples, 2 Mb) with
the Python reference walk (`hmm.py`), each target sample's 16 templates the same in every
arm, and only what the templates carry changed:

| arm | what each template sample contributes |
|---|---|
| as now | its allele shares from the `GP` pvc wrote — what the pass reads |
| perfect genotypes | its true genotype, unphased |
| perfect haplotypes | its two true haplotypes, as two templates |
| subtraction | its two haplotypes phased by the owner's subtraction idea: one copy follows the sample's closest relatives, the other is what its genotype leaves; two rounds |
| Beagle | its two haplotypes as Beagle 5.5 phased calling's `GT` (`beagle_phase.py`) |

The target is left out of every step that phases its templates — one Beagle run per target with
`excludesamples=` — so its own reads never return to it. The walk's lengths are fixed at the
simulator's (300 kb, 2 Mb) and the pass's rules for leaving a sample `./.` are not applied, so the
numbers compare the arms and are not the pass's own. On the one cohort checked, the "as now" arm
matches the Rust pass: 2,158 wrong against 2,155.

Scored: biallelic SNP records matching a truth record, one sample at one record a unit.

## Results

Wrong genotypes where samples have reads (calling's count first):

| cohort | calling | as now | perfect genotypes | perfect haplotypes | subtraction | Beagle |
|---|---|---|---|---|---|---|
| **all 3×** | | | | | | |
| outbred, 16 founders | 4,934 | 3,971 | 3,380 | 591 | 3,331 | 3,755 |
| outbred, 64 founders | 5,308 | 4,934 | 4,492 | 3,003 | 4,813 | 5,263 |
| inbreeding 0.5, 16 founders | 2,595 | 1,339 | 1,089 | 472 | 1,050 | 1,154 |
| inbreeding 0.8, 16 founders | 1,245 | 515 | 451 | 320 | 443 | 458 |
| inbreeding 0.95, 16 founders | 281 | 121 | 109 | 109 | 117 | 115 |
| **10 at 30×, 40 at 3×** | | | | | | |
| inbreeding 0.5, 16 founders | 2,187 | 1,106 | 918 | 392 | 850 | 787 |
| inbreeding 0.8, 16 founders | 1,053 | 409 | 357 | 251 | 335 | 310 |
| inbreeding 0.8, 64 founders | 1,080 | 729 | 688 | 617 | 728 | 715 |
| inbreeding 0.95, 64 founders | 235 | 159 | 156 | 149 | 162 | 161 |

How well each phased (`phasing_diagnostic.py`, `beagle_switch_errors.py`): how often the phase
flips wrongly between consecutive heterozygous sites it phased; a random phase flips half the time.

| cohort | subtraction | Beagle |
|---|---|---|
| all 3×, outbred, 16 founders | 21% (29% of heterozygous sites phased) | 21% |
| all 3×, outbred, 64 founders | — | 40% |
| all 3×, inbreeding 0.5 | 11% (44% phased) | 8% |
| mixed depth, inbreeding 0.5, 16 founders | — | 6% |
| mixed depth, inbreeding 0.8, 16 founders | — | 5% |
| mixed depth, inbreeding 0.8, 64 founders | — | 22% |

Missing genotypes filled, and wrong, mixed depth, inbreeding 0.5, 16 founders: as now 1,260 (1.5%
wrong), Beagle 1,636 (1.9%), perfect haplotypes 1,839 (0.8%). Beagle's haplotypes fill more but
less reliably: they come as hard alleles, so a phasing error is copied with full confidence.

## What the numbers say

1. **The ceiling is large where haplotypes are shared and heterozygosity is common** — outbred or
   inbreeding 0.5, 16 founders: perfect phasing leaves a sixth to a third of the errors the pass
   leaves now. It is nil for inbreeding 0.95, where genotypes already are haplotypes.
2. **Phasing 3× data alone does not reach it.** About 3 heterozygotes in 10 are called
   homozygous at 3×; Beagle reads hard calls only and phases them little better than chance
   (21–40% wrong flips) in the outbred cohorts. Giving it only the calls with `GQ` ≥ 20 was worse
   (41–49%): at 3× too few remain. The subtraction does about as well as unphased perfect
   genotypes and no better.
3. **Deep samples change that.** With 10 samples at 30×, Beagle's wrong flips fall to 5–6% and its
   haplotypes recover 45% (inbreeding 0.5) and 63% (0.8) of the gap between the pass as now and
   perfect phasing.
4. **Diversity sets the ceiling.** From 64 founders the ceiling itself is 15% at inbreeding 0.8,
   and Beagle phases with 22% wrong flips.

## What is not known, and the next measurement

How many haplotypes the owner's *pimpinellifolium* and *cerasiforme* accessions share over a few
hundred kilobases decides between the 16-founder and 64-founder answers, and the simulator cannot
say it. The measurement on the real material: the tomato hiding test (`tomato_fill_test.sh`) on a
cohort of those accessions with deep and shallow ones together, run twice — the pass as now, and
with templates' haplotypes from Beagle (phasing the deep accessions first). If Beagle's haplotypes
cut the hidden-genotype errors by a quarter or more there, a phaser inside the pass — one that
keeps each template's uncertainty, which Beagle's hard haplotypes lose — is worth building.
