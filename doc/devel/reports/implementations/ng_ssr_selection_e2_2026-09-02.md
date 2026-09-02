# E2 — what the shipped selector offers HG002, and which of §4.1 and §5 it reproduces

**Date:** 2026-09-02. **Plan:** [`candidate_alleles_ssr.md`](../../ng/impl_plan/candidate_alleles_ssr.md)
Milestone E step E2, and **Checkpoint E** (executed from
[`calling_loop_ssr.md`](../../ng/impl_plan/calling_loop_ssr.md) Milestone B).
**Design:** [`spec/candidate_alleles_ssr.md`](../../ng/spec/candidate_alleles_ssr.md) §4.1, §5, §13.
**Tools:** [`examples/ng_candidate_selection_probe.rs`](../../../../examples/ng_candidate_selection_probe.rs)
(its repeat-tract arm) and
[`benchmarks/ssr_hg002/src/ng_tract_candidate_recall.py`](../../../../benchmarks/ssr_hg002/src/ng_tract_candidate_recall.py).

---

## In one line

**Spec §4.1 reproduces; spec §5 does not.** Both true repeat *lengths* of an adjacent-length
heterozygote are offered at **97.8% (30×), 98.9% (50×), 98.9% (300×)** against the spec's
97–98%. Both *spellings* of one length are offered at **70.1%** at every depth, against the
spec's 86.1% and the 85.8% the shipped support share was expected to give. Two thirds of that
shortfall is upstream of selection.

## What was run

HG002 alone, the whole 50,000-region Tier BED, at 30×, 50× and 300×, through ng's own repeat
catalog (20,204 tracts inside that BED), ng's own merge, and the shipped `select_ssr`. Truth is
`benchmarks/ssr_hg002/truth/HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz`: both haplotype tract
sequences are reconstructed from the reference plus that VCF's phased edits inside the tract, the
window opened one base to the left so that a left-anchored indel — which is where left alignment
puts every repeat-length change — is applied. **20,006 tracts scored**; 198 dropped because a
truth record crosses the tract's end and the projection would be a guess.

## The numbers

All over the same 20,006 tracts. The last column is what spec §4.1 and §5 measured offline in
August, on the existing caller's Stage-1 pileup over its own 13,272-tract catalog.

| | 30× | 50× | 300× | 300× at 5 in 100 | spec, 300× |
|---|---|---|---|---|---|
| candidate sequences a tract | 1.416 | 1.416 | 1.409 | 1.448 | 1.26 |
| of which neither true nor the reference | 0.239 | 0.237 | 0.229 | 0.268 | — |
| **both true repeat lengths, one repeat apart** | **97.8%** | **98.9%** | **98.9%** | 99.0% | **97–98%** |
| both true repeat lengths, any distance | 92.4% | 94.1% | 94.1% | 94.3% | — |
| homozygous, both sequences offered | 99.6% | 99.6% | 99.7% | 99.7% | 99.8% |
| het, different length, both sequences | 90.3% | 91.8% | 91.9% | 92.1% | 86.0% |
| **het, same length, both spellings** | **70.1%** | **70.3%** | **70.1%** | 70.9% | **86.1%** |

**§4.1 and §5 measure different things and the difference is the whole story.** §4.1 scores
repeat *counts* — has nomination put both true lengths in play? §5 scores *sequences* — is each
true spelling in the candidate list? ng reproduces the first and falls 15 points short of the
second.

## The headline 1.615 was a denominator, and two thirds of the gap is that alone

A first run of the tract arm reported **1.615 candidate sequences a tract at 300×**. That is
taken over the 13,627 tracts the merge built a locus for — a set already enriched for tracts that
vary. The spec's 1.26 is over its catalog's whole tract set, and over ng's whole 20,204 the same
run gives **1.415**.

**And that figure is flat with depth where the enriched one is not**: 1.422 at 30×, 1.422 at 50×,
1.415 at 300×, against 1.971, 1.806 and 1.615 over the built subset. The swing in the second is
the merge's build rate, not the selector.

**Three things it is not.** The support share is not it and moves the wrong way: at the spec's
5 in 100 the figure *rises* to 1.455. The tract set's real content is not it either — HG002's own
variation averages 1.199 sequences a tract over ng's 20,204 tracts and 1.169 over the existing
catalog's 13,272, thirty thousandths apart. What is left is what the rule admits beyond the
truth: **0.229 sequences a tract**, against about 0.09 implied for the spec's figure on its own
tract set. Every one of ng's cleared the support rule on some sample's own reads.

## Where a true sequence is lost — and one third of it is selection's

At 300× with the shipped share, 387 true sequences are missing from a candidate list:

| where | misses |
|---|---|
| no read carried it, so the merge's table never held it | 233 |
| the merge's table held it, and the support bar refused it | 69 |
| it cleared the bar, and the per-sample top-`ploidy` cut dropped it | 57 |
| the merge refused the tract, so no locus was built | 28 |

**The first row is not selection's and it is the biggest.** Reading ten of the same-length-class
entries by hand: in every one the reads carry a sequence differing from the truth VCF's only in
*where* the change is placed inside the tract — the VCF writes a substitution, the reads show a
length change of one repeat. Some of those are two spellings of one haplotype rather than a lost
allele; the rest are a genuine difference in what ng's aligner assigns a read to.

**The spec names this ceiling for its own run and ours is lower.** §5 reports that 93.6% of
same-length truth sequences were shown by some read; here the merge's table holds the true
sequence at 76.4% of same-length tracts. That 17-point difference is the same-length shortfall,
and it sits in the read-to-tract-sequence layer, not in this module.

**The two rows that are selection's behave as designed.** The 69 the bar refused carried a median
9 reads in 100 of a sample's own against a 10-in-100 bar — thin, and refused at the bar's own
edge. The 57 the top-`ploidy` cut dropped are a diploid sample being asked for two lengths at a
tract where the truth has two; at 30× the same rows carry a median 26 reads in 100, so at low
depth the rung cut, not the bar, is what removes a true sequence.

## What the two tools now do

- **`examples/ng_candidate_selection_probe.rs`**, given a repeat catalog as its fourth argument,
  writes a per-tract dump (`NG_TRACT_DUMP`): every allele the merge's table held, with whether
  some sample's reads cleared the support bar and whether selection kept it. **Those two columns
  are what separates the four rows of the table above**; without the first, a sequence the rung
  cut dropped and one the bar refused read alike. `NG_TRACT_SHARE` moves the support share so the
  shipped value and the spec's swept one are comparable on one tract set.
- **`benchmarks/ssr_hg002/src/ng_tract_candidate_recall.py`** does the scoring. Handed a
  production `.cat` instead of a dump, it reports the truth floor and class sizes of *that* tract
  set, which is how two runs over different catalogs are made comparable at all.

## What this leaves open, and where it belongs

**The same-length ceiling is not this module's.** Three quarters of its losses are sequences no
read carried; the question is what ng's tract aligner assigns a read to when a substitution sits
inside a repeat, and it belongs with the STR generator and its aligner, not with step 6. The
selector's own two rows are 126 of 387 misses and both are the designed behaviour of a rule
whose constants are already measured.

**The spec's own figures are not restated as ng's.** They were taken by a different program over
a different tract set; §4.1's is reproduced here and §5's is not, and both facts are recorded
rather than either number being adopted.

## Validation

`cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` clean;
`cargo test --lib` 5,999 passed / 0 failed / 14 ignored, all in the container.
