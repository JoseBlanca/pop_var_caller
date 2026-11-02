# The neighbour-evidence pass, phase four: first measurements on the founder simulator

*Measurement report, 2026-11-02, branch `neighbour-pass` (code at `f43519fc`). Spec:
[`doc/devel/ng/spec/neighbour_evidence_pass.md`](../../ng/spec/neighbour_evidence_pass.md) §8.
Driver `doc/devel/ng/research/neighbour_evidence/founder_grid.py`, table
`grid_table.py`; cells in the worktree's `tmp/founder_grid/`. Ten cells of the planned seventeen;
the rest are listed in §5 and still running.*

## 1. What it says, in four lines

- **Where samples have reads, the pass removes errors in every cell at 3× and 10×** — from a
  quarter (an outbred cohort) to seven in eight (four founders).
- **Where samples have no reads, it calls nearly every site, and at 1× a third of those calls are
  wrong** — almost all homozygous for a common alternative, at `GQ` below 10.
- **Dropping calls without reads below `GQ` 10 removes most of that damage** (ruled, built): the
  pass then fills a quarter to nine tenths of the missing genotypes, wrong at most 3.5 times in 100
  — except with 10 samples at 1×, where one fill in six is wrong (decision 2).
- **The pass helps outbred and half-inbred cohorts that descend from 16 founders** — 4,934 wrong
  genotypes become 3,773 at inbreeding 0 — so the `auto` rule's median-inbreeding condition turns it
  off where it works. Decision 3, pending the 64-founder cells.

## 2. What was measured

Each cell is one simulated cohort: `founders` founder haplotypes, each sample's genome a mosaic of
them in stretches of 300 kb on average, a share `F` of it in autozygous runs (mean 2 Mb), sequenced
at `depth` with PCR-free libraries over 2 Mb in two contigs; mapped with minimap2, its parameters
fitted by `pvc estimate-parameters`, then called twice by `pvc call-from-psps`, with
`--neighbour-evidence off` and `on`. One unit is one sample at one SNP record whose position and
alleles match a truth record.

- **With reads**: calling gave the sample a genotype. Wrong genotypes counted before and after.
- **Without reads**: calling wrote `./.` (flat `PL`). The pass gives a genotype there; how many are
  wrong, and how many of those are true-reference genotypes called non-reference.
- **The §11 Q5 balance**: errors removed where samples have reads, against false non-reference calls
  added anywhere (a true-reference genotype called non-reference by the pass and not by calling).

## 3. The cells

| founders | F | samples | depth | with reads: wrong, calling → pass | without reads: sites; pass wrong (true reference called non-reference) |
|---|---|---|---|---|---|
| 4 | 0.95 | 50 | 3× | 228 → 32 of 35,048 | 1,952; 17 (6) |
| 16 | 0 | 50 | 3× | 4,934 → 3,773 of 45,556 | 2,444; 509 (157) |
| 16 | 0.5 | 50 | 3× | 2,595 → 1,173 of 47,107 | 2,643; 251 (67) |
| 16 | 0.8 | 50 | 3× | 1,245 → 374 of 46,862 | 2,538; 137 (40) |
| 16 | 0.95 | 10 | 1× | 46 → 44 of 1,635 | 955; 314 (216) |
| 16 | 0.95 | 10 | 3× | 88 → 51 of 8,068 | 402; 81 (31) |
| 16 | 0.95 | 10 | 10× | 19 → 6 of 8,949 | 1; 0 (0) |
| 16 | 0.95 | 50 | 1× | 57 → 67 of 4,180 | 2,570; 475 (374) |
| 16 | 0.95 | 50 | 3× | 281 → 90 of 46,976 | 2,674; 81 (29) |
| 16 | 0.95 | 50 | 10× | 29 → 6 of 50,045 | 5; 0 (0) |

At 1× pvc writes far fewer records — 4,180 sample-sites with reads against 46,976 at 3× for the
same cohort; why is calling's question, not measured here. The pass cannot add records (spec
§1.2).

**The no-read errors at 1× are one kind.** At 50 samples, all 374 true-reference genotypes called
non-reference are called homozygous for an alternative whose cohort frequency is 0.3 or more, and
360 of them have `GQ` below 10. At 1× a sample's neighbours carry too few reads to say which samples
it copies, so its posterior at a site without reads falls back towards the cohort's frequency, and
where that leans to the alternative the most likely genotype is the homozygous alternative. §4.2's
rule does not stop these: there *is* evidence within reach, it just does not identify the templates.

**`GQ` means what it says**, over all ten cells, for the pass's calls:

| GQ | with reads: wrong of called | without reads: wrong of called |
|---|---|---|
| 0–9 | 4,501 of 17,319 (26.0%) | 1,720 of 6,217 (27.7%) |
| 10–19 | 977 of 29,766 (3.3%) | 120 of 4,588 (2.6%) |
| 20–29 | 122 of 31,158 (0.4%) | 25 of 5,350 (0.5%) |
| 30 and above | 16 of 216,183 | none called |

A `GQ` of 10 to 19 promises 1.3 to 10 wrong in 100, and 20 to 29 promises 0.1 to 1; both sit inside.
So neighbouring records sharing reads (§11 Q7) does not make `GQ` overconfident here, and the
thinning Q7 holds in reserve is not needed on this evidence. Calls without reads are as well
calibrated as calls with them.

**Cost**: the pass added at most 1.5 seconds to a 50-sample run over 2 Mb (the column `pass s` of
`grid_table.py`'s output); the parameter fit dominates every cell, 2 to 13 minutes.

## 4. Decisions

### Decision 1 — no genotype at a site without reads below `GQ` 10. *Ruled yes by the owner, 2026-11-02; built (`NO_READ_MIN_GENOTYPE_QUALITY` in `refine.rs`).*

Today a sample with flat `PL` is written with the pass's genotype whenever it has an informative
record within `S` on both sides (§4.2, your ruling on §11 Q1). The proposal adds one condition: the
call's `GQ` is at least 10, a 1-in-10 chance of being wrong; below it the sample stays `./.` with the
pass's posterior as `GP`, as a withheld sample is today. Calls where the sample has reads are
untouched — `vcf_output.md` §7.1 forbids a quality floor there, and nothing here argues for one.

| cell | false non-reference added, no floor → floor 10 | no-read true non-reference recovered, no floor → floor 10 |
|---|---|---|
| 16 founders, 50 samples, 1× | 382 → 22 | 1,749 → 605 |
| 16 founders, 10 samples, 1× | 218 → 44 | 433 → 213 |
| 16 founders, 50 samples, 3× | 23 → −5 (fewer than calling) | 953 → 877 |
| 16 founders, F 0, 50 samples, 3× | 111 → −46 | 928 → 148 |
| 4 founders, 50 samples, 3× | −3 → −9 | 942 → 883 |

The cost is the recoveries the floor gives up, half to two thirds of them at 1× and under a tenth at
3× with 16 inbred founders; the gain is false calls cut by a factor of 5 to 17 at 1× and brought near
or below calling's own count at 3×. A floor of 20 cuts false calls little further and
gives up most recoveries (at 50 samples and 1×, 605 → 3). Either way it is one comparison in
`refine.rs`, made where the no-read rule is.

### Decision 2 — how the balance is judged. *Restated by the owner, 2026-11-02.*

The draft of this section asked whether a `./.` at a true non-reference genotype counts as an
error. The owner's answer: **a `./.` is not an error, and it is not a good outcome either.** The aim
is as little missing data as possible without filling it with wrong genotypes. So each depth is
judged by how much missing data the pass fills and how many of the fills are wrong, beside what it
does where samples have reads. With the `GQ` 10 floor:

| cell | `./.` before | filled | of them wrong | where samples have reads: wrong, calling → pass |
|---|---|---|---|---|
| 4 founders, 50 samples, 3× | 1,952 | 1,810 (93%) | 0 | 228 → 32 |
| 16 founders, F 0, 50 samples, 3× | 2,444 | 694 (28%) | 24 (3.5%) | 4,934 → 3,773 |
| 16 founders, F 0.5, 50 samples, 3× | 2,643 | 1,816 (69%) | 33 (1.8%) | 2,595 → 1,173 |
| 16 founders, F 0.8, 50 samples, 3× | 2,538 | 2,117 (83%) | 14 (0.7%) | 1,245 → 374 |
| 16 founders, 10 samples, 1× | 955 | 301 (32%) | 47 (16%) | 46 → 44 |
| 16 founders, 10 samples, 3× | 402 | 154 (38%) | 5 (3%) | 88 → 51 |
| 16 founders, 50 samples, 1× | 2,570 | 632 (25%) | 16 (2.5%) | 57 → 67 |
| 16 founders, 50 samples, 3× | 2,674 | 2,409 (90%) | 6 (0.25%) | 281 → 90 |
| 10× cells | 1 and 5 | 0 and 5 | 0 | 19 → 6, 29 → 6 |

**One corner fills with wrong data: 10 samples at 1×**, where one fill in six is wrong although each
claims `GQ` 10 or more — the floor of 20 does no better there (120 filled, 19 wrong). Everywhere else
a fill is wrong at most 3.5 times in 100, and in the inbred cohorts of 50 under 1 in 100. That
corner's `GQ` overstates its certainty, so no `GQ` floor fixes it; what does is the open question.
The 64-founder cell at 10 samples, still running, shows whether it is the small cohort or the low
depth; a second seed of the 10-sample cells shows whether it is this draw. *Counts above are the
written `GQ` after allele subsetting; the built floor reads the posterior before it, which moves a
handful of calls.*

### Decision 3 — the `auto` rule's inbreeding condition. I recommend deciding after the 64-founder cells.

`auto` runs the pass only where the cohort's median inbreeding is at least 0.5 (§5, starting value).
With 16 founders and 50 samples at 3× it helps at every inbreeding measured:

| F | wrong where samples have reads, calling → pass | false non-reference added, with the floor |
|---|---|---|
| 0 | 4,934 → 3,773 | −46 |
| 0.5 | 2,595 → 1,173 | 48 |
| 0.8 | 1,245 → 374 | 7 |
| 0.95 | 281 → 90 | −5 |

The evaluation report's outbred simulation gained nothing, and the difference is the founders: there
every site was drawn independently, so no two samples shared a stretch; here 16 founders give every
sample many close relatives, outbred or not. What decides the threshold is therefore how many
haplotypes the cohort's samples share, not their inbreeding — and the cells that test that, 64
founders at 10 and 50 samples, are still running. If the pass still helps there at F 0, I will
recommend dropping the inbreeding condition and keeping only `N_min`; if it does not, a condition on
something closer to haplotype sharing. Until then `auto` stays as specified.

## 5. Not yet measured

- **Still running**: 64 founders at 50 and at 10 samples (3×), 30× at 10 and 50 samples, 200
  samples at 1× and 3×. This report gains their rows when they land.
- **Not started**, from §8: the mixed inbred-and-outbred cohort (the simulator draws one `F` for the
  whole cohort); the autozygous posterior against `ancestry.tsv` (not in the VCF); template bias
  (§11 Q4); the circularity (§8.5); tomato1 held out and the two deep-and-shallow plants; cost at 500
  samples; accuracy at core edges against core centres. Thread-count invariance is covered by tests
  (§8.7), not by these runs.
- **What these numbers are**: one seed a cell, a 2-Mb genome, PCR-free reads. How far a cell's
  counts move between seeds is not measured; the differences the decisions rest on are hundreds of
  sample-sites, the smallest (10 samples at 1×) tens.
