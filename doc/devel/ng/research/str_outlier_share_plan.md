# Fitting the repeat-tract outlier share — research plan

**Started:** 2026-10-18. **Owner's request:** make the repeat-tract outlier share a number
`estimate-parameters` fits, instead of a stated constant. Simulations are allowed.
**Worktree:** `../pop_var_caller-str-group-b`, branch `str-group-b` (experiment branch, never
merged; this plan moves to a research branch when the program settles).

## 1. The question

At a repeat tract, the read likelihood mixes three sources for each read: the sample's own alleles
with stutter, a contaminating individual, and an **outlier term** for reads no allele explains. The
outlier term takes a share of the reads, today a stated **0.20**, and spreads it evenly over every
tract length the stutter model can reach — about 22 lengths at a homopolymer. So every read's
probability has a floor of roughly `0.20 / 22 ≈ 0.009` under it, whatever genotype is scored.

**The number does two jobs, and only one of them is a share of reads.**

- **Stray reads.** Some reads really did come from somewhere no allele explains: a paralogous
  tract, a chimera, a mis-anchored read. Their share is a property of the library and can be
  counted. Counted literally at HG002's homozygous-reference tracts at 30x it is **1 read in
  2,300 at homopolymers and 1 in 209 at longer motifs**
  (`doc/devel/reports/tract_genotype_investigation/readmodel.md` §2.3).
- **A cap on one read's pull.** A floor under every read's probability bounds how much any one
  read can count against a genotype. GATK does the same job with a Phred-45 cap. At 0.20 the floor
  sits between a one-repeat slip and a two-repeat one, so a read two or more repeats from an allele
  says almost nothing (`readmodel.md` §3).

**The best value depends on the rest of the model, which a parameter's true value would not.**
Measured on the HG002 tandem-repeat benchmark (one sample, 5x/10x/30x), counting correct tract
genotypes over the tracts every arm built:

| arm | best outlier share | what lowering it does |
|---|---|---|
| `--defaults` (shipped stutter, nearly flat prior) | 0.20 | +4 and +6 tracts at 5x and 10x, **−37 at 30x** |
| fitted stutter and prior | **0.03** | +65, +93 and +12 tracts at 5x, 10x and 30x from 0.20 |

A stated value tuned against one arm is compensating for something that arm does not model. **The
program's first job is to find what that is**; fitting the literal stray-read share alone would
ship a value near 0.001 and lose calls at depth.

## 2. What exists

- **The STR stutter fit has no outlier term.** Every whole-repeat read gets a stutter
  probability; a far length lands in the end bucket of the ±8 recorded offsets; a tract where more
  than 1 in 10 of the differing reads are a non-whole number of repeats off is dropped from the fit
  (`parameter_prepass_ssr.md` §5; `ssr_fit.rs`, `expected_counts.rs`). The fit already integrates
  over each sample's genotype under the fitted prior, so a new read-level term would sit in that
  same sum.
- **The generic path solved the analogous problem in its fit, and refused it at calling.** Its
  parameter fit carries two site classes: *clean* sites, where reads disagree with the reference at
  0.2–0.4%, and *noisy* sites, 0.4–1.4% of sites, where they disagree at 4–10%
  (`parameter_prepass_generic.md` §2.1; `research/noise_model_overdispersion_2026-08-10.md`). It
  exists so the fitted heterozygosity is not inflated (1.41 → 1.09 times the benchmark on HG002).
  **At calling the read likelihood carries no noisy class** (owner, 2026-08-19,
  `read_likelihoods.md`): from one locus's reads, a site where half the reads disagree because of a
  duplication and a heterozygous site are the same distribution. Depth separates them, and that is
  the duplicated-region filter's job.
- **The specs.** `read_likelihoods.md` §4.5 declares the outlier share inherited and unfitted;
  `parameters_file.md` §3.8 lists it among the constants no fit produces. Both still say 0.01 where
  the code ships 0.20 — a stale statement to correct whatever this program decides.

## 3. What the swept value may be compensating for

Each is tested by changing one thing at a time, in simulation and on real data.

| | candidate | what it predicts | first evidence (HG002 TR benchmark) |
|---|---|---|---|
| A | **Stray reads** — the literal share | best share near the counted rate, the same at every depth | counted 1/2,300 and 1/209; the best share is 10–100 times that |
| B | **Multi-repeat slips under-modelled** by the stutter shape | the floor matters for reads two repeats off; a better-fitted fall-off removes the need | two-repeat reads ≈ 1 in 100 at truth homozygotes |
| C | **A tract-level noisy class**: a few tracts whose reads show a second length for reasons other than the genotype (paralog, truth error) | over-confidence that grows with depth; a read-level floor can only cap it read by read | ~5 in 100 truth-homozygous tracts carry > 35% one-repeat-off reads; the likelihood is over-confident towards heterozygotes at 30x |
| D | **Reads not independent** at a tract | allele balance at heterozygotes wider than binomial | read-to-read correlation 0.01–0.03 at 30x, none measurable at 5x/10x |
| E | **The prior** — `--defaults`' near-flat prior needs a weaker likelihood | the best share falls when the prior is fitted | 0.20 with `--defaults`, 0.03 with the fit |

C is the generic path's noisy class seen at repeat tracts. If it is the main cause, the
2026-08-19 decision applies: the likelihood cannot tell such a tract from a heterozygote, and the
answer belongs to depth or to the cohort (a noisy tract is noisy in every sample), not to a
per-read share. The program has to say how far that carries over, since the STR path already has a
read-level class at calling and the generic path does not.

## 4. Phases

Each phase ends with a short report under `doc/devel/ng/research/` and a pause for the owner.

1. **Characterise the unexplained reads on real data, against truth.** At HG002 truth homozygotes
   (5x/10x/30x), the distribution of reads at each offset from the true length against what the
   fitted stutter predicts — the repeat-tract counterpart of the generic study's tail table. Fit
   three descriptions to it: one stray share; a beta-binomial (reads correlated within a tract); two
   tract classes. Split by period and depth. Repeat on the tomato cohort (tomato2, 3x), where an
   inbred line's heterozygote rate stands in for truth. **Output:** which of A–E carries the mass,
   and how big each is.
2. **Simulator.** Draws a cohort's tract genotypes from a prior, reads with stutter, and each of
   A–D as a separate knob, at any sample count and depth. It scores reads with **pvc's own
   repeat-tract likelihood** (a Rust example calling the library row), so calling arms differ only
   in parameters. Validated against pvc's `--explain-loci` likelihoods on real tracts.
3. **Truth-free estimators,** fitted the way `estimate-parameters` would fit them — from the
   reads, integrating over genotypes, never from truth:
   - the stray share inside the stutter fit's genotype sum;
   - the same plus a two-class tract model (share of noisy tracts, their off-length rate);
   - the same plus a within-tract correlation.
   Measured for bias and spread across the range grid: **1, 10, 100, 1,000 samples × 3, 10, 30,
   100, 300 reads**. Each estimator must say what it does where it is not identified (one sample at
   3x): fall back to a stated value, and report that it did.
4. **Calling in simulation.** Each estimator's output used at calling, against the stated 0.20 and
   against the true generating values, over the same grid. Scored as tract-genotype accuracy and as
   the calibration of heterozygote-against-homozygote odds.
5. **Real data.** The chosen estimator against `--defaults` and the current fit on the HG002 TR
   benchmark (tract genotypes over common tracts; RTG vcfeval; heterozygotes lost with the
   likelihood preferring the homozygote, and with the prior overruling it) and on tomato2 (the
   identity oracle's inputs; the heterozygote rate).
6. **Spec amendment and implementation plan**, if a fitted number wins: where it is fitted, the
   parameters-file format change, its warrant and fallback, and what `--defaults` keeps.

## 5. What would end the program early

- Phase 1 finds C carries nearly all the mass and A is negligible. Then a fitted stray share would
  sit near 0.001 and the real decision is about C — the 2026-08-19 decision's territory — and
  that decision comes back to the owner before any estimator is built.
- Phase 3 finds no truth-free estimator lands within a factor of two of the value calling needs,
  at any corner of the grid. Then the honest outcome is a stated value per arm, documented as a
  tuned cap rather than a share.

## 6. Already measured (inputs to phase 1)

From the calibration study of 2026-10-18 (HG002 TR benchmark; scripts in the worktree's `tmp/`):

- With `--defaults`, heterozygote-against-homozygote odds are calibrated near even odds at 5x and
  10x, and over-confident towards heterozygotes at 30x.
- Per read at a homopolymer, a read on the second allele adds 2.1 (natural log) to the
  heterozygote's evidence under `--defaults` and 2.5–3.4 under fitted stutter; truth implies
  2.4–3.0 or more.
- The fitted prior matches truth for homozygous reference and reference/one-repeat heterozygotes
  to within 2%, and under-predicts heterozygotes whose two alleles are both non-reference and one
  repeat apart about three-fold (167 against 50 at homopolymers).
- The held edge-deletion fix (a deletion that touches the tract edge counted as crossing it) is
  part of every arm from here on; without it a homozygous deletion's reads count as reference
  length and look like C.
