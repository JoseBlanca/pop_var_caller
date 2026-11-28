# Refinement: drawn alleles for shallow samples only (2026-11-27 and 28)

Branch `no-drawn-alleles`, from main 859033f1. Specs changed: `neighbour_evidence_phasing.md` §2.2
(the decision), §2.3 and §5; `neighbour_evidence_pass.md` (memory per thread). User docs:
`neighbour_evidence.md`, `recipes.md`, `methods.md` §11.

## 1. What was wrong

Rerunning the founder-simulator cohorts of 2026-11-02 on current main, to check the user docs'
10-sample warning, showed refinement making genotypes worse than calling alone in four of seven
cohorts, two of them inbred cohorts where `--refine-genotypes auto` switches it on. In 10 samples
from 16 founders at F 0.95 and 3×, the wrong genotypes where samples had reads went from 117 with
refinement off to 320 with it on; on 2026-11-02 the same cohort went from 88 to 51. Nearly all the
added errors were true homozygous-reference genotypes called non-reference: 2 → 143.

## 2. Where it came from

A walk of main's merge history on that cohort (simulation and mapping done once, each commit's
`pvc` building its own psps, fit and calls):

| commit | wrong, off → on | homozygous reference called non-reference, off → on |
|---|---|---|
| f43519fc (2026-11-02) | 134 → 115 | 2 → 1 |
| 0019f331, 88f93350, ff4ea293, 1edce6ca | 134 → 108–109 | 2 → 3 |
| **11a97683 `refine-sources`** | **134 → 321** | **2 → 123** |
| c5aaa3c8 … aa39bbbf, 95569d47 (main) | 117–134 → 284–320 | 2 → 113–143 |

`refine-sources` (`refine_sources_2026-11-05.md`) made three changes. Switching two of them off in
current main, on the same cohort:

| setting | wrong, off → on | homozygous reference called non-reference |
|---|---|---|
| main | 117 → 320 | 2 → 143 |
| background landing at 1/K, as before | 117 → 243 | 2 → 90 |
| background landing at 1 in 100 | 117 → 267 | 2 → 110 |
| **emissions read the shares, not drawn alleles** | **117 → 111** | **2 → 10** |
| both | 117 → 105 | 2 → 3 |

The drawn alleles do most of the damage. `refine-sources` measured them only on 1× to 30× msprime
cohorts of 40 samples and an F2; none was a cohort where every sample is at 3×.

## 3. What was tried

Removing the draw (commit 54f35c4e) fixed every 3× cohort but cost 1× and 2× samples 3 to 9
genotypes right in 100, mostly left missing. Three ways to keep part of the draw were then built
behind a temporary switch (branch `exp-draw-quality`, never merged) and measured on 13 cohorts:

- **draw only from templates sure at the site**, a share of at least 90% or 99% for the surest
  allele (the owner's suggestion: let the templates' genotype quality choose where to copy from);
- **draw only for a shallow sample**: one with reads at fewer than 9 of the window's walked sites
  in 10 — at 1× about 63 sites in 100 have a read, 2× about 86, 3× about 95.

Drawing only from sure templates kept little of the 1× gain — F2 heterozygous offspring 87.1 in
100 right, against 87.6 with no draw and 90.9 drawing for shallow samples; outbred 1×
homozygous-reference genotypes 66.3, against 64.0 and 73.5 — because the gain comes from drawing
the *unsure* templates, which is what makes them plainly right or wrong. Drawing for shallow samples kept the gain and the protection, and is what was adopted.

## 4. The change

From round 2, a sample's walk reads alleles drawn from its templates only in a window where it is
shallow (`refine::DRAW_BELOW_READ_SHARE` = 0.9, through `refine::is_shallow`); otherwise its
emissions read the shares, as the prediction and the drawn haplotypes always do. The code of the
draw is `refine-sources`' own, restored unchanged; only the condition is new. The background's
landing share stays at 10⁻⁵.

## 5. Measurements

The final build on rick (`f65b459d`, the branch's three commits on 859033f1) reproduces the
experiment build's numbers exactly in all 13 cohorts. Off and on use the same psps and parameters.

**Founder cohorts.** `founder_grid.py`'s simulation, 2 Mb over two contigs, seed 20261102; one
unit is one sample at one matched SNP record where calling gave it a genotype. "No draw" is
54f35c4e.

| founders | F | samples | depth | off | main | no draw | sure ≥ 90% | **shallow only** |
|---|---|---|---|---|---|---|---|---|
| 16 | 0.95 | 10 | 3× | 117 | 320 | 111 | 108 | **123** |
| 16 | 0.95 | 50 | 3× | 295 | 167 | 69 | 67 | **69** |
| 16 | 0.95 | 200 | 3× | 1,359 | 597 | 354 | 358 | **359** |
| 64 | 0.95 | 10 | 3× | 118 | 415 | 124 | 125 | **143** |
| 64 | 0.95 | 50 | 3× | 359 | 447 | 174 | 173 | **174** |
| 64 | 0 | 10 | 3× | 973 | 1,088 | 935 | 943 | **935** |
| 16 | 0 | 50 | 3× | 4,959 | 2,188 | 1,780 | 1,823 | **1,814** |
| 64 | 0 | 50 | 3× | 5,031 | 5,127 | 4,149 | 4,211 | **4,149** |
| 16 | 0.95 | 50 | 1× | 231 | 499 | 395 | 408 | **499** |
| 64 | 0 | 50 | 1× | 3,404 | 3,953 | 4,959 | 4,822 | **3,953** |

In the two 10-sample inbred cohorts a few 3× samples fall below the line in some windows and draw;
the one from 64 founders ends 25 above calling (143 against 118, of about 6,800 genotypes).

**The `refine-sources` cohorts**, rebuilt on rick as its §7 builds them: an F2 (two parents at
30×, 20 offspring at 1×, 4 Mb), and 40-sample msprime cohorts, 10 samples each at 30×, 5×, 2×, 1×.
Genotypes right, in 100:

| cohort, samples, true genotype | off | main | no draw | sure ≥ 90% | **shallow only** |
|---|---|---|---|---|---|
| F2 offspring 1×, heterozygous (7,284) | 26.3 | 90.8 | 87.6 | 87.1 | **90.9** |
| outbred 1×, heterozygous (5,257) | 33.1 | 60.9 | 56.2 | 56.3 | **62.2** |
| outbred 1×, homozygous alternative (2,315) | 52.5 | 83.2 | 75.4 | 77.4 | **82.8** |
| outbred 1×, homozygous reference (18,688) | 49.6 | 73.2 | 64.0 | 66.3 | **73.5** |
| outbred 2×, heterozygous (5,215) | 54.0 | 80.4 | 76.6 | 78.3 | **80.2** |
| outbred 5×, heterozygous (5,268) | 85.5 | 94.5 | 93.3 | 93.8 | **93.4** |
| selfing 1×, heterozygous (838) | 12.1 | 81.0 | 75.7 | 79.4 | **80.7** |
| selfing 2×, heterozygous (202) | 29.7 | 78.2 | 77.2 | 77.7 | **77.7** |

Sure ≥ 99%, not in the table, kept no more of the gain: F2 heterozygous offspring 88.7, outbred
1× homozygous reference 64.0, selfing 1× heterozygous 76.8. Samples at 30× do not move under any
setting.

## 6. What it leaves

- **Where every sample is at 1×, refinement raises the wrong genotypes where samples have reads**,
  under every setting: 50 inbred samples from 16 founders, 231 with it off and 499 with it on.
  It also fills about 4,000 of their 8,175 genotypes without reads, 91 of them wrong, so a run is
  better on the whole, but a 1×-only cohort is a weak case. `auto` would switch it on there.
- **Small inbred cohorts from many founders**: 10 samples from 64 founders, 118 → 143.
- **The line is chosen from simulations**, and only one value was tried; 0.8 would exclude 2×
  samples and lose their gain (80.2 against 76.6 heterozygotes right).
- **The published refinement figures predate this.** `accuracy.md`'s imputation section was
  measured on 909ea135, with the draw for every sample; the `refine_sim` benchmark should be rerun
  on the merged branch.

## 7. Tests

- New: `refine::tests::a_sample_is_shallow_below_nine_sites_in_ten_with_reads`.
- Re-recorded: the phased-window phase fingerprint — that window's samples are not shallow, so
  its phase is the no-draw one (`0xbd48_9385_1932_57b9`, 2,894 columns). The posterior
  fingerprint's window does draw and keeps its value.
- `cargo test --release --lib neighbour_evidence`: 167 pass on macOS arm64 and Linux x86-64.

## 8. Reproducing

On rick, in `~/pop_var_caller/.claude/worktrees/founder-rerun/tmp/founder_bisect/`: `bisect.sh`
(the history walk), `exp.sh` (one founder cohort, `CELL=` chooses it, `TREE=` and `TAG=` point it
at a build), `suite.sh` (the `refine-sources` cohorts, `SUITE=` for a separate output tree), and
`exp.tsv` / `suite_fix/*/score_*.txt` with every number above. `exp.patch` and `draw.patch` are the
temporary switches. Copies of the scripts under the main checkout's `tmp/founder_bisect/`.
