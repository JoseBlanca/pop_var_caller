# The two STR genotype priors on GIAB HG002: not one genotype differs

*Measurement, 2026-08-18. Answers the first half of open question Q5 in
[`../ng/spec/calling_priors.md`](../ng/spec/calling_priors.md) §11 — does the marginalized
genotype prior genotype repeat tracts better than the plug-in one? Run on production's STR caller,
where both priors already exist behind an environment variable, so the answer could be had before
ng's STR path is built.*

---

## What was run

`benchmarks/ssr_hg002/` — the GIAB HG002 tandem-repeat benchmark v1.0.1 on GRCh38, 13,272 catalog
loci inside the Tier confident regions, truth from an **assembly** rather than from another caller.
HG002 is outbred and 72% of its truth genotypes are heterozygous, so the prior's inbreeding branch
is inert and the comparison isolates the frequency prior — the two faults that made the earlier
tomato-versus-HipSTR comparison unreadable.

The per-sample pileups (`.psp`) already existed from 2026-07-09 and were **not** rebuilt. Only
`ssr-call` was re-run, twice per coverage from the same file with the same binary, the arms
differing only by `PVC_SSR_MARGINALIZED_PRIOR=1`
([`driver.rs:287`](../../../src/ssr/cohort/driver.rs)). Coverages 50, 30, 20 and 15× from the
bundle's subsampled ladder; the native 288× was skipped because the reads swamp any prior there.

The plug-in arm reproduces the July run's record counts exactly at every coverage (50× 2873, 30×
2852, 20× 3572, 15× 5515), which is the check that the binary and catalog are the ones the bundle
was built with.

Scoring: `benchmarks/ssr_hg002/src/prior_genotype_accuracy.py`, derived from the bundle's own
`genotype_accuracy.py` so it inherits the period-aware truth matching and the anchor-extension that
the truth set requires (GIAB anchors an STR indel at the base *before* the tract).

## The result

Genotype accuracy is conditional on detection: of the loci a caller called length-variant, how many
carry HG002's exact allele lengths. `truth-variant` is 2,652 at every coverage.

| coverage | prior | detected | exact match | accuracy given detection |
|---|---|---|---|---|
| 50× | plug-in | 1,691 | 1,540 | 0.911 |
| 50× | marginalized | 1,679 | 1,535 | 0.914 |
| 30× | plug-in | 1,460 | 1,320 | 0.904 |
| 30× | marginalized | 1,435 | 1,307 | 0.911 |
| 20× | plug-in | 1,086 | 980 | 0.902 |
| 20× | marginalized | 1,060 | 965 | 0.910 |
| 15× | plug-in | 753 | 664 | 0.882 |
| 15× | marginalized | 732 | 652 | 0.891 |

**At every locus both priors emitted, the two genotypes are identical — 0 differences out of the 732 to 1,679 loci
both emitted, depending on coverage.** Checked twice, and by two routes: inside the scoring script, and by
joining the two VCFs on position and comparing the genotype field directly.

The entire difference is which loci are emitted at all. The marginalized prior withholds 13 to 32
loci per coverage — about **1 in 100** of those the plug-in emitted — and the marginalized arm adds
back at most one:

| coverage | withheld | of those, truth-variant | plug-in had them right | plug-in had them wrong | not truth-variant |
|---|---|---|---|---|---|
| 50× | 13 | 13 | 6 | 7 | 0 |
| 30× | 28 | 25 | 13 | 12 | 3 |
| 20× | 32 | 26 | 15 | 11 | 6 |
| 15× | 30 | 22 | 12 | 10 | 8 |

**The withheld loci are close to a coin flip.** At 30× it drops 13 correct genotypes and 12 wrong
ones; the accuracy-given-detection gain of 0.7 percentage points is that near-tie plus a handful of
loci with no truth variant at all. Nothing here separates the two priors.

## Why it is a zero, and what that means

**The single-sample failure the marginalized prior was built to fix does not exist on the STR
path**, and the reason is the seed, not the marginalization.

On the SNP path the plug-in prior's frequency estimate is regularised toward the **reference**
allele with a concentration of 10 against 0.01
([`posterior_engine.rs:107`](../../../src/var_calling/posterior_engine.rs)). With one sample that
pins the alternative frequency near 0.083 whatever the reads say, and the resulting prior favours
heterozygous over homozygous-variant by 22 to 1 — which is what called 214 true homozygous-variant
GIAB sites heterozygous at 5×.

The STR plug-in has no such pull, because its seed is **mode-centred**: `G₀` decays geometrically
away from the cohort's modal repeat count
([`allele_freq_prior.rs:25`](../../../src/ssr/cohort/allele_freq_prior.rs)), and with one sample
the cohort's mode is that sample's own mode. There is no reference allele being privileged, so
there is nothing for marginalization to undo.

**This confirms, from the other direction, the trap recorded in
[`calling_priors.md`](../ng/spec/calling_priors.md) §2.3:** the measured SNP gain came from setting
the reference concentration to 1 rather than 10, not from the act of integrating. Where the seed
was never reference-privileged, integrating changes nothing.

## What this does and does not settle

**Settles:** adopting the marginalized prior on ng's STR path costs nothing at the single-sample
end. That was worth knowing, because the only prior evidence — tomato, 51 samples, 30% fewer loci
emitted than HipSTR called polymorphic — looked alarming and had no way to tell a conservative
prior from a wrong one.

**Relocates that 30%:** since the single-sample end is neutral to four decimal places on genotypes,
the tomato effect is a **cohort** effect — the leave-one-out term doing something — and must be
judged in a cohort. It is not a property of the prior's form.

**Does not settle:** whether the marginalized prior helps in a cohort, which is the half of Q5 that
remains. At one sample the leave-one-out term is zero by construction
([`calling_priors.md`](../ng/spec/calling_priors.md) §6), so this run could not have tested it. The
cohort half needs a multi-sample panel with orthogonal truth, and the only multi-sample panel we
have is tomato, whose truth is another caller.

**A limitation of the instrument, stated plainly.** Production's STR caller detects 64% of
truth-variant loci at 50× and 28% at 15×, and its recall collapses below that (about 1 locus in 10
at 10×, essentially nothing at 5×). The loci it does detect are the well-covered ones, where the
read likelihood is sharp and any prior is a small term. So this benchmark cannot reach the
low-depth regime where a prior would matter most on a single sample — the SNP result was measured
at 5×, and there is no STR equivalent of that measurement to be had from this caller.

## Reproducing

```
./scripts/dev.sh bash benchmarks/ssr_hg002/src/run_prior_comparison.sh
cd benchmarks/ssr_hg002 && uv run src/prior_genotype_accuracy.py
```

VCFs at `benchmarks/ssr_hg002/results/ours_prior_comparison/{plugin,marginalized}/` (gitignored
with the rest of `results/`).
