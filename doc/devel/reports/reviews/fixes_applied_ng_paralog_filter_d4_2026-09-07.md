# Fixes applied — ng_paralog_filter_d4

**Date:** 2026-09-07
**Review:** [ng_paralog_filter_d4_2026-09-07.md](ng_paralog_filter_d4_2026-09-07.md) (0 Blocker /
5 Major / 9 Minor, Request-changes; one sub-agent in an isolated worktree)
**Branch:** `ng-paralog-filter`

## What was done

**All five Majors and all nine Minors applied.** One defect was in the analysis script and is fixed
there; the rest were in the report's prose. **The corrected reading strengthens the step's
conclusion** rather than weakening it.

## The headline is now the opposite shape, and it argues §8 harder

The report's answer was "the coverage-only score at repeat tracts is inert" — 0 flagged of 601, most
of them at one value, and a tract would need coverage "far beyond anything in either run" to reach
the cut. Two of those three were wrong.

**What is actually true, measured:**

| | |
|---|---|
| copy number at which a zero-read record reaches HG002's cut of 4.65 | **1.372** — 4.04 σ₀ above one copy |
| the highest-scoring tract in the run | **1.310** — 3.37 σ₀ |
| the gap | **0.68 σ₀**, about **17 reads in a 364-read window** |

The arm is flat and then abrupt: −0.66 at one copy, 1.3 at 1.30, 3.5 at 1.35, 6.2 at 1.40, 56.7 at
two. **The run came 83% of the way to the cut on the axis the score actually reads**, and the nat
count hid it. So the coverage-only arm is *not* a safe zero: a modestly deeper tract window would be
flagged on coverage alone with the allele half switched off, which spec §1 says cannot tell a
duplication from anything else that raises depth. **That is a stronger case for §8's tract-aware
allele term than the one the report made.**

## The mechanism for the plateau was wrong, and the right one is the finding

The report said the ratio is constant "where a sample's coverage is at one copy". It is that same
constant for **every copy number from about 0.2 to about 1.1** — a tract at half a copy scores what a
tract at one copy scores. **σ₀ sets the plateau's width, not the coverage**: at 0.092 the carrier
hypothesis at 1.5 copies is more than five standard deviations away, its branch underflows, and the
ratio answers the prior.

**Confirmed directly from the run's own window dump** rather than by re-deriving the model: the 469
tracts within 0.0001 of the plateau value span raw window depths from **0.202 to 1.272 times the
one-copy scale**. One number, across a band a copy wide.

**And the tomato run disagrees with HG002, which is why it must not be pooled.** At σ₀ between 0.158
and 0.302 the carrier branch does not underflow, so the arm responds: **all 21 tomato tract ratios
are distinct**, spanning 2.49 nats. The report had claimed one pair coincided — see below.

## The script produced one of the errors, and is fixed

`ng_paralog_filter_by_kind.py` printed "1 of them share one ratio exactly (−0.9582)" for the tomato
run. That "1" is the size of the largest group — **no two tracts share a ratio at all** — and
−0.9582 is an arbitrary tie-break among 21 singletons. The report read it as a coincidence.

Three changes:

- when nothing repeats it now says **"no two tracts share a ratio: the arm responds to coverage over
  this cohort"** in words, instead of a count of one;
- when something does repeat it gives both the exact-print count and the count within 0.001, because
  the VCF prints four decimals and "exactly" cannot be established from the file — 331 print
  −0.6572 and **529** are within 0.001;
- the table marks **`= max`** wherever the ninety-ninth percentile is the maximum by construction,
  which is every row with 100 records or fewer, so two columns are not read as two measurements.

## The rest

| finding | what changed |
|---|---|
| the deepest tract window and the highest tract ratio quoted as one record | **two different records**: the deepest window is chr6:150,466,707 at 1.4156× the scale with a ratio of −0.6252, on the plateau; the ratio of 1.4843 is chr11:9,161,348 at 1.2033×. The GC multiplier decides — 0.295 against 0.457 — and this is the **fourth step running** to conflate raw depth over the scale with the copy number the model reads |
| "agrees wherever it speaks" | production removed **6,223** records over the whole of `regions.bed`; 17 of them — 0.27% — fall in the two intervals ng ran, and the agreement is over those |
| "over 440 nats" | **432.53** |
| "the indel arm is 22 of 284" | **19 of 284** — 15 deletions and 4 insertions; the 3 equal-length substitutions are neither, which is the report's own reason for separating them |
| "within a nat of the 6.2500 cut" | 5.23 is **1.02 nats** short |
| 601 tracts across two runs, as one denominator | **not pooled.** The runs differ in species, cohort size, cut and depth, and they disagree about the plateau — which is the reason not to add them |
| "the ninetieth percentile equal to the median equal to the lowest" | true at two decimals; at four, p90 is −0.6566 against −0.6572 |
| "equal-length multi-base substitutions … their depth is not biased at all" | **8 of the 9 on tomato change exactly one base**, so they are mostly SNPs written with a shared flanking base; the depth claim is now marked as asserted rather than measured |
| "a deletion's ratios reach the far tail **because** it carries the allele signal" | softened to what is measured. A tract's coverage arm is not capped — it reaches 56.7 at two copies — so the contrast is about the copy numbers these tracts had, plus the second arm a deletion has |

## Validation

No library code changed, so the gate is `4e061476`'s, re-confirmed after the script edit:

- `cargo test --all-features --lib --bins --tests` — lib **6,643 passed, 0 failed, 15 ignored**; one
  integration failure which is `main`'s.
- `cargo clippy --lib --bins --tests --all-features -- -D warnings` — 9 errors in three kinds.
- `cargo fmt --check` — 9 unique files, none this plan's.
- **Both by-kind tables re-run from the fixed script** and every cell in the report matches its
  output.
- **The standing oracles are unmoved**: tomato 2,311 records at sha256 `84ad19c2…0590d`, HG002 8,242
  at `0c991c8f…76b7e6`.
