# The hidden-duplication filter — D4: what it scores and flags, by kind of record

**Date:** 2026-09-07
**Plan step:** [hidden_paralog_filter.md](../../ng/impl_plan/hidden_paralog_filter.md) Milestone D, step D4
**Spec:** [hidden_paralog_filter.md](../../ng/spec/hidden_paralog_filter.md) §3.2, §8, §9
**Branch:** `ng-paralog-filter`

## The answer

**The coverage-only score at a repeat tract answers the prior and nothing else, over a band a copy
wide — and then turns on hard.** On HG002 it returns **one number for every tract between 0.20 and
1.27 times the one-copy scale**: 331 of the 580 print the identical ratio −0.6572 and 529 are within
0.001 of it. Above that band it is the steepest arm in the filter, and **the highest-scoring tract
in the run came within 0.68 σ₀ of the cut — about 17 more reads in a 364-read window.**

**So spec §8's tract-aware allele term is what any tract filtering would rest on, and the present
arm is not a safe zero either.** It cannot see most of the tract population at all, and where it
does see, it is close enough to the cut that a modestly deeper window would flag a tract on coverage
alone, with no allele evidence of any kind behind it.

**Production's filter agrees with ng wherever the two can be compared.** Over the 1,867 records both
callers wrote on the six-accession slice, production removed 17 and ng removed 27; **all 17 are
among ng's 27**, and production removed nothing ng kept. That is 17 of the 6,223 records production
removes over the whole of `regions.bed` — the comparison covers the two intervals ng ran, not
production's whole verdict.

## What each kind was scored and flagged

Both tables are from the tagging run, since a dropped record leaves no ratio in the file.
Percentiles are the value at sorted index `int(p·n)`. **`= max` marks a cell where the
ninety-ninth percentile is the maximum by construction** — `int(0.99 · n)` is the last index for
every n at or below 100, so on those rows the two columns are one measurement.

**D2 — six tomato accessions, one copy fitted between 3.02 and 29.24 reads a window, cut 6.2500:**

| kind | written | scored | flagged | lowest | median | p90 | p99 | highest |
|---|---|---|---|---|---|---|---|---|
| biallelic SNP | 2,217 | 2,217 | **36** | −104.93 | −10.91 | −5.04 | 11.12 | 34.76 |
| multiallelic | 27 | 27 | 0 | −92.12 | −14.17 | −3.26 | = max | −2.64 |
| insertion | 15 | 15 | 0 | −40.78 | −9.34 | −5.06 | = max | −4.49 |
| deletion | 22 | 22 | 0 | −42.32 | −9.16 | −6.53 | = max | −5.75 |
| equal-length substitution | 9 | 9 | 0 | −12.00 | −8.42 | 5.23 | = max | 5.23 |
| **repeat tract** | **21** | **21** | **0** | **−2.02** | **−2.00** | **−0.96** | **= max** | **0.47** |

**D3 — HG002, one copy fitted at 302.45 reads a window, cut 4.6500:**

| kind | written | scored | flagged | lowest | median | p90 | p99 | highest |
|---|---|---|---|---|---|---|---|---|
| biallelic SNP | 7,039 | 7,039 | **262** | −348.37 | −29.86 | −15.38 | 21.98 | 84.16 |
| multiallelic | 25 | 25 | 0 | −305.56 | −43.59 | −15.84 | = max | −6.44 |
| insertion | 253 | 253 | **4** | −182.23 | −29.07 | −10.81 | 13.32 | 16.84 |
| deletion | 257 | 257 | **15** | −253.09 | −21.34 | −0.25 | 13.71 | 47.92 |
| equal-length substitution | 88 | 88 | **3** | −308.69 | −33.97 | −5.42 | = max | 16.38 |
| **repeat tract** | **580** | **580** | **0** | **−0.66** | **−0.66** | **−0.66** | **−0.56** | **1.48** |

**Every record of every kind was scored on both runs; none was unscorable.** That is spec §3.2's
decision working: every record is scored, and the kinds are told apart by what evidence they hand
the scorer rather than by being skipped.

## The tract row, and what it says about spec §8

**The distinguishing feature is the spread.** A biallelic SNP's ratios run over **432 nats** on
HG002; a tract's run over **2.1**, and most of them sit at one value. Spec §3.2 states the reason: a
tract's spilled rows carry no read counts, so the allele term is zero under every genotype and every
carrier configuration, and the ratio rests on coverage alone.

**But "coverage says one copy, so the ratio is constant" is the wrong mechanism, and the right one
matters.** The ratio is that same constant for *every* copy number from about 0.2 up to about 1.1 —
a tract at half a copy scores exactly what a tract at one copy scores. What sets the width of that
plateau is **σ₀, not the coverage**: at HG002's σ₀ of 0.092 the carrier hypothesis sits at 1.5
copies, more than five standard deviations away, so its branch underflows and the ratio reduces to
the prior. Measured on the run's own window dump, the 469 tracts within 0.0001 of the plateau value
span raw window depths from **0.202 to 1.272 times the one-copy scale**.

**At tomato's larger σ₀ the branch does not underflow, and the arm does respond.** All 21 tomato
tracts have **distinct** ratios spanning 2.49 nats — the arm is measuring coverage there, it is
simply 5.78 nats below the cut. *(An earlier draft of this report claimed one pair coincided; none
does. It was a misreading of the analysis script's own output, which printed "1 of them share one
ratio" to mean the largest group had size one. The script now says so in words.)*

**How close the arm came, which is the number that decides §8.**

| | |
|---|---|
| copy number at which a zero-read record reaches HG002's cut of 4.65 | **1.372** — 4.04 σ₀ above one copy |
| the highest-scoring tract in the run | **1.310** — 3.37 σ₀ |
| the gap | **0.68 σ₀**, about 17 reads in a 364-read window |

The arm is flat and then abrupt: at σ₀ = 0.092 the ratio is −0.66 at one copy, 1.3 at 1.30, 3.5 at
1.35, 6.2 at 1.40, 11.9 at 1.5 and 56.7 at two copies. **The nat count makes the margin look safe
and the copy-number axis shows it is not.**

**What this does and does not settle.**

- **It settles that the arm cannot discriminate among tracts below about 1.1 copies**, which on
  HG002 is most of them. Whatever a tract's coverage is down there, it gets the prior.
- **It settles that the arm is not inert above that**, and that this run came 83% of the way to the
  cut on the axis the score actually reads. A tract flagged that way would be flagged on coverage
  alone, with the allele half switched off — which is the configuration spec §1 warns cannot tell a
  duplication from anything else that raises depth.
- **It does not settle whether tracts harbour collapsed duplications worth removing.** Nothing here
  measures that and neither benchmark has a tract truth set.
- **The owner's ruling of 2026-09-07 stands**: scoring a tract on coverage alone is a deferral, not
  a verdict on its evidence, taken because stutter makes the allele split unmodelled rather than
  uninformative. This is the size of what the deferral costs.

## Deletions, counted apart from tracts — and they behave nothing alike

The plan asks for this split because the two kinds' depth is biased in opposite directions. The
measurement justifies the split several times over:

| | D3 highest ratio | D3 flagged | D3 p90 |
|---|---|---|---|
| deletion | **47.92** | **15 of 257** | −0.25 |
| repeat tract | 1.48 | 0 of 580 | −0.66 |

**A deletion's ratios reach the far tail on this run and a tract's do not.** The observed contrast
is what justifies the split; **the causal reading needs care**, because a tract's coverage arm is
not capped — it reaches 56.7 at two copies, steeper than anything a deletion does. What separates
the two populations here is the copy numbers these tracts happened to have, plus the fact that a
deletion has a second arm and a tract does not. Lumping them — which a classifier testing
`len(ALT) <= len(REF)` would partly do anyway — would have averaged a live arm against a mostly
plateaued one and shown neither.

**Equal-length multi-base substitutions are counted as their own kind**, because they are neither
insertion nor deletion and a `len(ALT) <= len(REF)` test files them as deletions — an earlier draft
of D2's report did exactly that and reported 31 deletions where there are 22 deletions and 9 of
these. 3 of the 88 on HG002 were flagged.

**The name flatters them.** 8 of the 9 on tomato change exactly one base (`CA → TA`, `AC → AA`);
only 13,881,461 changes two. They are mostly SNPs written with a shared flanking base, so their
depth is not biased for the same reason a SNP's is not — which this step asserts rather than
measures.

## Where the flagged records are, by kind and by run

**On tomato at six accessions only biallelic SNPs were flagged** — 36 of them — and no other kind
came within a nat of the 6.2500 cut; the nearest miss is an equal-length substitution at 5.23,
which is 1.02 nats short.

**On HG002 four kinds were flagged**: 262 biallelic SNPs, 15 deletions, 4 insertions, 3 equal-length
substitutions. **The 19 indels among them are the ones D3's truth comparison found are real**: 17 of
19 are GIAB benchmark variants, because alignment bias puts a true heterozygous indel's alternative
fraction near a third, which is the model's `(T=3, m=1)` carrier value. That finding belongs to D3
and the owner has ruled it is studied later; what D4 adds is that the two indel kinds are
**19 of 284** removals here — 15 deletions and 4 insertions, the 3 equal-length substitutions being
neither — and **0 of 36 on tomato**, where indel depth is too low for the allele term to dominate.

## Production's filter over the same six accessions

Production has no tag mode — it drops what it flags and spec §9 records that as deliberate — so its
flagged set is recovered by calling the same six accessions twice, once with the filter and once
with `--no-paralog-filter`, and taking the difference.
[`tmp/d_milestone/d4_production.sh`](../../../../tmp/d_milestone/d4_production.sh) does it, over
production's own stored files from the benchmark, restricted afterwards to the two intervals ng ran
because production's `var-calling` takes no region list.

**It is not an oracle, and spec §7 says why**: production's depth is raw depth where ng's is
observation depth, its inbreeding coefficient is one cohort number where ng's is per sample, and it
scores biallelic SNPs only. The plan asks for the overlap, reported and not chased.

| | |
|---|---|
| records both callers wrote over the two intervals | 1,867 |
| production removed | 17 |
| ng removed | 27 |
| **removed by both** | **17** |
| removed by production and kept by ng | **0** |

**Production's flagged set is a strict subset of ng's.** Of the 19 records ng removed that
production did not, 9 are records production never called at all and 10 it called and kept.

**The two calibrations are not close, which makes the agreement worth more.** Production fitted a
duplication rate of **0.174464** and a cut of **2.9500**; ng fitted **0.020634** and **6.2500** on
the same six accessions. Two implementations disagreeing by a factor of eight on how common
duplications are, and by a factor of two on where to cut, still nominate the same seventeen records
— and fifteen of the seventeen lie in the two tight clusters D2 identified at 13,808,477–13,808,592
and 13,889,118–13,889,203.

## What this step adds

| file | what it is |
|---|---|
| [`scripts/ng_paralog_filter_by_kind.py`](../../../../scripts/ng_paralog_filter_by_kind.py) | the two tables above, from a tagging run: written, scored and flagged by kind, with each kind's ratio spread |

The classifier keys on the `STR` INFO flag first — which is what the filter itself keys on, since
`is_repeat_tract` decides which signals a record's samples carry — then on allele counts and
lengths, and it keeps equal-length substitutions apart from deletions.

No library code changed in this step.

## Deviations from the plan, recorded

1. **The production comparison is restricted to records both callers wrote**, 1,867 of them.
   Production wrote 1,879 records over the two intervals and ng wrote 2,311; comparing flagged sets
   over the union would count ng flagging records production never had the chance to.
2. **Production ran over the whole of `regions.bed` and was restricted afterwards**, because its
   `var-calling` has no region list. Its own calibration is therefore fitted over all of
   `regions.bed` rather than over the two intervals — one more reason it is a sanity comparison and
   not an oracle.
3. **The ratio distribution is given as five quantiles a kind, not as a full distribution.** The
   kinds with the fewest records — 9 equal-length substitutions on tomato — cannot support more.

## What D4 settles, and what it leaves

- **Spec §8's tract question has its number**: on HG002 the coverage-only arm returns one value for
  every tract between 0.20 and 1.27 times the one-copy scale — 331 of 580 printing it — and no tract
  was flagged on either run. **But it is not a safe zero**: the arm reaches the cut at 1.37 copies
  and this run's highest-scoring tract sat at 1.31, a gap of 17 reads in a 364-read window. Any
  tract filtering needs the tract-aware allele term, and leaving the coverage arm alone is itself a
  choice with a measured margin. **The two runs are not pooled**: they differ in species, cohort
  size, cut and depth, and they disagree about the plateau — tomato's 21 tract ratios are all
  distinct.
- **Deletions and tracts are different populations** and the plan was right to separate them: on
  HG002 a deletion's ratios reach 47.92 and 15 were flagged, against a tract's 1.48 and none.
- **Production and ng agree on every record production flags**, at cuts and fitted rates that are
  nowhere near each other.
- **Whether tracts need filtering at all is unmeasured**, and neither benchmark can answer it.
- **The indel behaviour on HG002 is D3's finding and is deliberately not acted on**; D4 records that
  it is confined to the high-depth run.
