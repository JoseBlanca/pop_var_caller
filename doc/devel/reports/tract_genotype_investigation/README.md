# Six investigations into ng's repeat-tract genotype errors

**Date:** 2026-09-02. These are the raw reports of six parallel investigations run to answer one
question the owner asked: *ng's repeat-tract genotypes are worse than they should be — how do we
improve them?* The conclusions drawn from them, and the runs that test those conclusions, are in
[`../ng_tract_genotype_improvement_2026-09-02.md`](../ng_tract_genotype_improvement_2026-09-02.md).

**Read that first.** These are the working; it is the answer.

| file | the question it was given |
|---|---|
| [`slippage.md`](slippage.md) | Measure the real slippage on HG002's reads, per stratum, and say how far HipSTR's shipped constants are from it. Produced [`fitted_slippage_hg002_30x.toml`](fitted_slippage_hg002_30x.toml), the rows the improvement report's winning run was given. |
| [`readmodel.md`](readmodel.md) | List every inherited constant in the repeat-tract read likelihood, measure what is measurable, and rank them by how much moving each could change a call. Found the outlier weight. |
| [`hetfhom.md`](hetfhom.md) | Characterise every case where ng calls a heterozygote and the truth says homozygous, and say what would remove them and at what cost. |
| [`prior.md`](prior.md) | Measure what the flat genotype prior costs, and what the true length spectrum looks like. |
| [`hipstr.md`](hipstr.md) | Score HipSTR — which fits its stutter model per locus — on exactly ng's ground with exactly ng's instrument, and say whether the gap is a candidate problem or a model problem. |
| [`unseen.md`](unseen.md) | Split the 268 true sequences no read carried between unobservable, an alignment loss, and a truth-set artefact. |

## What is and is not verified

**The numbers in these reports were not independently re-derived.** Three of their central claims
were tested by full runs of the shipped binary and did hold — the outlier weight's value, the
fitted slippage's value, and the direction of the genotype prior's effect (see the improvement
report §2 and §3). The rest stand on each investigation's own working, which is included beside
each report in the form the investigation left it.

**Two defects in existing tooling were found along the way** and are recorded in the reports rather
than fixed:

- `benchmarks/ssr_hg002/src/ng_tract_candidate_recall.py` keeps a truth record only where its REF
  span reaches the tract's start, so an insertion anchored at `start - 1` — where left alignment
  puts every repeat-length gain — is dropped. 35 of the 268 "no read carried it" cases have such a
  record, and 2,737 of the 19,613 clean tracts do. Widening the window takes the missing-sequence
  total from 434 to 564, so the window and the overlap rule need reworking together rather than a
  one-character fix (`unseen.md`).
- The same script's truth reconstruction misses records left-aligned before the tract, which
  `slippage.md` had to work around before it could fit anything; it reports 1,652 tracts called
  homozygous-reference that have a truth record within 30 bases.

## Scope

One sample (GIAB HG002), one chemistry, 30× and 50×, on the 50,000-interval tandem-repeat
benchmark. Nothing here says what any of it does on a cohort, or at three reads a position.
