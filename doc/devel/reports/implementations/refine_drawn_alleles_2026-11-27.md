# Refinement: the walk reads the template shares again, not alleles drawn from them (2026-11-27)

Branch `no-drawn-alleles`, from main 859033f1. Specs changed: `neighbour_evidence_phasing.md` §2.2
(the decision) and §2.3, `neighbour_evidence_pass.md` (memory per thread). User docs:
`neighbour_evidence.md`, `recipes.md`, `methods.md` §11.

## 1. What was wrong

Rerunning the founder-simulator cohorts of 2026-11-02 on current main, to check the user docs'
10-sample warning, showed refinement making genotypes worse than calling alone in four of seven
cohorts, two of them inbred cohorts where `--refine-genotypes auto` switches it on. In 10 samples
from 16 founders at F 0.95 and 3×, the wrong genotypes where samples had reads went from 117 with
refinement off to 320 with it on; on 2026-11-02 the same cohort went from 88 to 51.

Nearly all the added errors were true homozygous-reference genotypes called non-reference: 2 → 143
in that cohort.

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

The drawn alleles do most of the damage. Why they hurt here was not measured. These samples are
at 3×, and in a small cohort the 16 templates are themselves mostly 3× samples; a drawn allele
then makes an unsure template look plainly right or wrong at random, at sites where the sample's
own reads do count. The `refine-sources` measurements were all 1× to 30× msprime cohorts of 40
samples and an F2; none was a small cohort at 3×.

## 3. The change

From round 2 the walk's emissions read each template's shares, as the prediction and the drawn
haplotypes already did. Removed: the draw (`refine::solve_round`), `WalkSite::emission_shares`,
the allele-pair fast path that served definite alleles (`model::definite_emissions_into`) with its
scratch, and `StreamPurpose::DefiniteAlleles`. The background's landing share stays at 10⁻⁵: giving
it 1/K again cost the F2 a third of its heterozygotes (below).

## 4. Measurements

Built from the branch (`1e089142` on rick, the same tree as 54f35c4e) and from current main, the
same psps and parameters for off and on. Founder cohorts: `founder_grid.py`'s simulation, 2 Mb over
two contigs, 3×, seed 20261102; one unit is one sample at one matched SNP record where calling
gave it a genotype.

| founders | F | samples | off | main | branch |
|---|---|---|---|---|---|
| 16 | 0.95 | 10 | 117 | 320 | **111** |
| 16 | 0.95 | 50 | 295 | 167 | **69** |
| 64 | 0.95 | 10 | 118 | 415 | **124** |
| 64 | 0.95 | 50 | 359 | 447 | **174** |
| 64 | 0 | 10 | 973 | 1,088 | **935** |
| 16 | 0 | 50 | 4,959 | 2,188 | **1,780** |
| 64 | 0 | 50 | 5,031 | 5,127 | **4,149** |

At sites where a sample had no reads, the branch fills fewer and gets almost none wrong: in 64
founders, F 0.95, 50 samples, wrong fills 161 → 6 of 2,890, left missing 792 → 1,518.

The `refine-sources` cohorts, rebuilt on rick as its §7 builds them: an F2 (two parents at 30×,
20 offspring at 1×, 4 Mb), and 40-sample msprime cohorts, 10 samples each at 30×, 5×, 2×, 1×.
Genotypes right, in 100:

| cohort, samples, true genotype | off | main | branch | both changes reverted |
|---|---|---|---|---|
| F2 offspring 1×, heterozygous (7,284) | 26.3 | 90.8 | 87.6 | 55.2 |
| F2 offspring 1×, homozygous (7,196) | 54.0 | 93.6 | 93.5 | 92.9 |
| outbred 1×, heterozygous (5,257) | 33.1 | 60.9 | 56.2 | 52.4 |
| outbred 1×, homozygous alternative (2,315) | 52.5 | 83.2 | 75.4 | 70.9 |
| outbred 1×, homozygous reference (18,688) | 49.6 | 73.2 | 64.0 | 56.9 |
| outbred 2×, heterozygous (5,215) | 54.0 | 80.4 | 76.6 | 75.7 |
| outbred 5×, heterozygous (5,268) | 85.5 | 94.5 | 93.3 | 93.0 |
| selfing 1×, heterozygous (838) | 12.1 | 81.0 | 75.7 | 74.1 |
| selfing 2×, heterozygous (202) | 29.7 | 78.2 | 77.2 | 77.2 |

What the branch loses at 1× it mostly leaves missing, not wrong: outbred 1× homozygous reference
wrong 0.9 → 0.5 in 100, missing 6.4 → 15.9; heterozygous wrong 11.3 → 10.5, missing 23.2 → 28.8.
The F2's heterozygous wrong is 1.4 against 1.5 in 100. Samples at 5× and above move by at most 1.2
in 100; at 30× by nothing.

## 5. What it leaves

- **Small cohorts from many founders are still the weak case**: 10 samples from 64 founders,
  inbred, 118 → 124. Refinement there neither helps nor harms much.
- **The published refinement figures are from the defective code.** `accuracy.md`'s
  imputation section was measured on 909ea135, after `refine-sources`; the `refine_sim` benchmark
  should be rerun on the merged branch. Its selfing cohort shows the same pattern (false variant
  genotypes 73 → 178 with refinement on).
- A middle way — drawn alleles only for samples below some depth — was not tried: it needs a
  cut-off chosen from these simulations, and its own measurement round.

## 6. Tests

- Removed: `the_emissions_and_the_prediction_read_their_own_shares` and
  `definite_emissions_have_the_general_bits`, which tested only what is removed;
  `the_log_likelihood_alone_has_the_walks_bits` keeps its long window without definite sites.
- Re-recorded: the phased-window fingerprints (`refine::tests`), the arithmetic changing on
  purpose. They pass on macOS arm64 and Linux x86-64.
- `cargo test --release --lib neighbour_evidence`: 164 pass.

## 7. Reproducing

On rick, in `~/pop_var_caller/.claude/worktrees/founder-rerun/tmp/founder_bisect/`: `bisect.sh`
(the history walk), `exp.sh` (one founder cohort under the experiment switches, `CELL=` chooses it;
`TREE=` and `TAG=` point it at another build), `suite.sh` (the `refine-sources` cohorts, `SUITE=`
for a separate output tree), `exp.patch` (the two temporary switches, never merged). Copies under
the main checkout's `tmp/founder_bisect/`.
