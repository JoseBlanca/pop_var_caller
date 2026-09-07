# Fixes applied — ng_paralog_filter_d3

**Date:** 2026-09-07
**Review:** [ng_paralog_filter_d3_2026-09-07.md](ng_paralog_filter_d3_2026-09-07.md) (**1 Blocker**
/ 2 Major / 4 Minor, Request-changes; one sub-agent in an isolated worktree)
**Branch:** `ng-paralog-filter`

## What was done

**All seven applied to the report. No code changed, and the Blocker is deliberately not fixed in
code** — it is a measurement about the model, and the owner has ruled on what happens to it.

## The Blocker: the finding stands, and it is recorded rather than acted on

The report claimed the filter's calibration passed, on a biallelic-SNP restriction, and set the
all-kinds figure aside as method-uncertain. Normalising both sides against the same reference — the
fix for exactly that method concern — makes it worse:

| | records | in GIAB's HG002 v4.2.1 benchmark |
|---|---|---|
| substitutions removed | 265 | 2 — 0.8 in 100 |
| **indels removed** | **19** | **17 — 89 in 100** |
| everything removed | 284 | 19 — **6.7 in 100**, against a target of **1 in 100** |

**Re-derived by the orchestrator before being accepted**, with `bcftools norm -f <the analysis set>
-m -any` on both sides and `comm` over the normalised keys; the script is
[`tmp/d_milestone/d3_truth.sh`](../../../../tmp/d_milestone/d3_truth.sh). The removed indels'
alternative-read fractions were listed and sixteen of the twenty-two sit between 0.29 and 0.37,
against the half a real heterozygote gives — which is the alignment bias of an indel at 300× reads,
and is exactly the model's `(T=3, m=1)` carrier value.

**The owner's ruling, 2026-09-07:** *"don't worry right now, we'll study the giab problem and we'll
fix it if necessary. The main objective right now is to have a working filter so we can carry out
the real experiments and decide how to tweak it."*

**So neither spec §3.2 nor the code changes.** The report says what was measured, says it is not
being acted on and why, and gives the number to come back to. The step's headline was rewritten
from "the filter's calibration passes" to what the measurement actually shows: the target is met on
substitutions and missed sevenfold on the file, and all of the miss is indels.

## The two Majors

### The copy numbers were in the wrong units, for the third step running

The report divided window depth by `single_copy_scale` and called the quotient copies, concluding
the removed records were *under*-covered at 0.92 against the kept records' 0.99. The scorer divides
by `single_copy_scale × gc_multiplier(gc)`, and **the removed records sit at higher GC — median
0.473 against 0.406** (measured, on this run's own window dump), where Illumina depth falls. In the
model's own units they are at about **1.01 copies against 0.99**: marginally above, not 7% below.

**This is D1's defect and D2's defect a third time**, and the report now handles it differently
rather than restating a corrected number. It reports what it can measure — raw depth against the
fitted scale, and GC — says plainly that the two corrections push opposite ways, and states that
**the run prints the scale but not the GC curve, so a record's copy number cannot be recovered from
the output at all.** The load-bearing conclusion is unchanged and now rests on something weaker and
true: at about one copy there is no excess coverage, where a carrier is expected at 1.5.

**Making the copy number readable is D4's**, which needs it to explain ratios by record kind.

### The caveat about the truth set was false

The report warned the callset was "one step removed from the consortium's release". Its own header
reads
`view -R HG002_bench_azar_sorted_1000.bed … HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` —
**it is GIAB v4.2.1**, subset to the bench intervals and nothing else. Checked by the orchestrator
against the file. The false caveat is gone; the real one is in: 37 of the 265 removed-and-absent
substitutions sit in bench intervals with no truth call at all, where absence says least.

## The four Minors

| finding | what changed |
|---|---|
| the nats arithmetic was ~150 and ~20, "an order of magnitude" | **111.5 and 34.4, a margin of 3.2 to 1**, at the configuration the model actually scores best — `(T=4, m=1)`, not the `(3,1)` the draft assumed |
| "finding exactly the shape it was built to find" | withdrawn. Across the range the 284 occupy the grid values are never 0.042 apart, so every value lands near one; against a uniform null the clustering is 27% within 0.01 against 28% expected, and the 284's median fraction is a **quarter**, not the sixth the ten highest show. What is left is the separation from a half, which names no cause |
| the ten-highest table mixed `DP` and `AD` | the alternative fraction and the read counts are both from `AD`, which is what the scorer sums; the worst row printed 419 beside `232, 46` |
| "no repeat tract was flagged", asserted bare | now with the mechanism and the measurement: a tract is scored on coverage alone, **331 of 580 share one ratio exactly (−0.6572)**, the highest tract ratio is **1.4843** against a cut of 4.6500, and the deepest tract window is 1.42 times the fitted scale. So it is a fact about this run's coverage, and a tract cannot respond to the allele half that removed all 284 |

**One reviewer minor was not adopted.** It gave the kept heterozygotes' ninetieth-percentile
alternative fraction as 0.528 over 4,864 records; measured here it is 0.532 over 4,952. The gap is
which records count as heterozygous. Rather than take either number on authority, the report now
states its rule — a kept record whose two called alleles differ, alternative reads being `AD`
summed over the alternatives — and notes that a different selection moves the percentile by about
0.004.

## Validation

No library code changed, so the gate is `86f5fe86`'s, re-confirmed:

- `cargo test --all-features --lib --bins --tests` — lib **6,643 passed, 0 failed, 15 ignored**;
  one integration failure which is `main`'s.
- `cargo clippy --lib --bins --tests --all-features -- -D warnings` — 9 errors in three kinds.
- `cargo fmt --check` — 9 unique files, none this plan's.
- **The tomato oracle is unmoved** by the harness's new mode: the six-accession run in alignments
  mode still gives 2,311 records and sha256 `84ad19c22dd1…0590d`.
- **HG002's own baseline**, for later comparison: 8,242 records, sha256
  `0c991c8f43c30d3318049c974a0be3a8ed174463f3ae200aa3d89ea04876b7e6`.
