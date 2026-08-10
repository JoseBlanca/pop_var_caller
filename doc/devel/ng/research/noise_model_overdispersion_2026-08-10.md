# The generic path's error model has too thin a tail, and heterozygosity pays for it

**Date:** 2026-08-10. **Sample:** HG002, the 100 GIAB confident regions
(`benchmarks/giab/per_sample/bed/HG002_bench_azar_merged_100.bed`, 571,984 bases), at 30× and
at 300×. **Truth:** the whole-genome v4.2.1 benchmark VCF, restricted to those regions.
**Found while measuring G1's model-free anchors.**

## What was found

`spec/parameter_prepass_generic.md` §2 states the generic noise model as *"a per-base
substitution rate and nothing else"* — one `ε` per read group, and a read over a reference copy
shows another base with that probability. **The body of that distribution is right and its tail
is wrong by five orders of magnitude**, and the three-genotype mixture has only one class that
can absorb the surplus.

At the **550,976** loci where the benchmark records no variant of any kind, against a binomial
at the fitted rate:

| alternative reads | observed | binomial expects | ratio |
|---|---|---|---|
| 0 | 516,889 | 515,081 | 1.0 |
| 1 | 30,906 | 34,664 | 0.9 |
| 2 | 2,363 | 1,202 | 2.0 |
| 3 | 505 | 28.5 | 18 |
| 4 | 174 | 0.5 | 337 |
| 5 | 78 | 0.008 | 10,000 |
| ≥ 6 | 61 | 0.0001 | 626,000 |

**818 non-variant loci carry three or more alternative reads where the model predicts 29.**

## What it costs

Fitted heterozygosity comes out **1.41 times** the benchmark's count over the identical loci —
1.407 × 10⁻³ against 9.967 × 10⁻⁴, 776 heterozygous sites where the truth has 550. Splitting
the fit's heterozygous mass by alternative-read count shows where it comes from: sites with 3,
4 or 5 alternative reads hold **167 units of heterozygous mass against 15 real
heterozygotes**, and that band is exactly where the tail is fattest.

By depth, the excess per 10,000 loci falls monotonically — 27.7 at depth 11–15, 12.5 at 16–20,
7.6 at 21–25, 3.6 at 26–30, 1.5 at 31–40, 1.0 above 40. Deeper sites separate error from
variation; a 30× sample has a long shoulder below 25 where they cannot be told apart.

**The error rate barely notices.** 818 sites in 550,976 is 1.5 in a thousand, and the fitted
rate still lands **1.1% from a model-free count** — 2.239 × 10⁻³ against 2.263 × 10⁻³, one
fifth of a ladder rung, where the model-free number is 38,450 mismatching bases in 16,992,201
read observations at benchmark homozygous-reference positions. The parameter with the smallest
true value absorbs the misspecification; the parameter with the largest does not feel it.

## Three candidate models, fitted over the same 551,843 loci

Truth heterozygosity 9.9666 × 10⁻⁴, truth every-copy-non-reference 5.7444 × 10⁻⁴.

| emission | extra parameters | het / truth | hom-alt / truth | ln L |
|---|---|---|---|---|
| one error rate (today) | — | 1.417 | 0.919 | −153,066.3 |
| beta-binomial | 1 (`ρ`) | 1.190 | 0.926 | −150,409.8 |
| **two site classes** | 2 (`w`, `ε_noisy`) | **1.091** | 0.930 | **−149,984.9** |

**The comparison is against a faithful control.** The one-error-rate row is a re-implementation
of today's model and returns heterozygosity 1.4126 × 10⁻³ where production returns
1.407 × 10⁻³ — 0.4% apart — so what the other two rows are being compared against is the
shipping estimator and not a straw man.

Two site classes wins on both axes, and by 425 nats over the beta-binomial for one further
parameter, which no information criterion would refuse.

### What the winning model says

A site is **clean** with probability `1 − w` and **noisy** with probability `w`; the genotype
emission uses that site's error rate. Fitted independently at the two depths:

| | `ε_clean` | `ε_noisy` | `w` | het / truth |
|---|---|---|---|---|
| HG002 30× | 1.895 × 10⁻³ | 5.29 × 10⁻² | 0.88% | 1.091 |
| HG002 300× | 1.952 × 10⁻³ | 4.24 × 10⁻² | 1.28% | 1.141 |

**About one site in 110 disagrees with the reference at 5% rather than 0.19%.** That is what
mismapped reads and error-prone sequence contexts look like, and the current model has nowhere
to put them. The two depths were fitted separately and found the same population, which is the
main reason to believe it is a property of the data and not of the fitting.

### A second benefit, not looked for

**The clean rate is more depth-invariant than the single rate**, which is the property
`arch/parameter_prepass_generic.md` §9's coverage sweep exists to test. Across the same
four-fold change in the depth actually used:

| | 30× | 300× | drift |
|---|---|---|---|
| single error rate | 2.268 × 10⁻³ | 2.407 × 10⁻³ | +6.1% |
| `ε_clean` | 1.895 × 10⁻³ | 1.952 × 10⁻³ | **+3.0%** |
| heterozygosity, one rate | 1.4126 × 10⁻³ | 1.5401 × 10⁻³ | +9.0% |
| heterozygosity, two classes | 1.0877 × 10⁻³ | 1.1455 × 10⁻³ | **+5.3%** |

## The two controls, computed exactly

Both are run on cell probabilities computed in closed form under a known truth over the **real
depth distribution** of the HG002 30× walk — 69 distinct depths, 2,464 cells, probability mass
summing to 1.00000000 — so there is no sampling noise in either and any departure is bias.

**Control: a world that does not need the extension.** Genotype frequencies
(0.9885, 0.0105, 0.0010), one error rate 1.0 × 10⁻³.

| | error rate | heterozygosity | hom-non-ref |
|---|---|---|---|
| one error rate | −0.0009% | **+0.0000%** | −0.0000% |
| two site classes | **−1.10%** | **−0.0006%** | +0.0000% |

**The genotype frequencies survive intact and the error rate does not.** Given a world with no
noisy sites the extension still splits off a spurious class — `w` = 0.48% — and the clean rate
it reports is **1.10% below the truth**, about one fifth of a ladder rung. That is the price of
carrying the extension where it is not needed, and it is paid by the parameter that today has
the best anchor.

**Recovery: a world that really has one.** Same frequencies, `w` = 0.0100,
`ε_noisy` = 5.0 × 10⁻².

| | error rate | heterozygosity | hom-non-ref |
|---|---|---|---|
| one error rate | 1.428 × 10⁻³ for a truth of 1.0 × 10⁻³ | +4.09% | −0.03% |
| two site classes | 1.0001 × 10⁻³, `ε_noisy` 5.0012 × 10⁻², `w` 0.0100 | **−0.000%** | **+0.000%** |

## The bias is an absolute number of sites, not a percentage

The recovery world's one-rate heterozygosity is only **+4.09%** where the real data's is +42%,
and the reason is not that the synthetic tail is milder. **+4.09% of that world's 5,794
heterozygotes is 237 sites; +42% of HG002's 550 is 231 sites.** The same tail, the same depth
distribution, the same ~230 extra sites — expressed against a heterozygosity ten times larger.

**So the parameter is worst for the samples this caller is aimed at.** A 1% heterozygous
outbred sample would see the excess as 2%; HG002 at 0.1% sees 42%; a selfing tomato landrace,
where heterozygosity is lower still, would see more. The bias grows exactly as the quantity
being estimated shrinks.

## What it does not fix

**A residual of 9% at 30× and 14% at 300×.** Better than 42% and 54%, and not zero. The
homozygous-non-reference rate barely moves either — 0.919 → 0.930 — so whatever holds it 7%
below the benchmark is a different thing and is still unexplained.

## Provenance of every number here

The alignment side is ng's own walk: `count_whole_site` over the loci
`SampleLocusObservationsIterator` produced under `ReadFilterConfig::default()` (MAPQ ≥ 20, read
length ≥ 30, no BAQ), so the reads counted are the reads the estimator is fitted on. The
model-free error rate is the one exception and was counted with `samtools mpileup` under the
same policy; it agrees with the fit to 1.1%, which is itself the cross-check that the two
read populations are the same. The truth side is `bcftools query` over the benchmark VCF, with
a locus counted heterozygous when any of its reference positions carries a `0/1` record and
every-copy-non-reference when it carries `1/1` or `1/2` — the second because ng's classes count
non-reference **copies**, and a `1/2` site has no reference copy.
