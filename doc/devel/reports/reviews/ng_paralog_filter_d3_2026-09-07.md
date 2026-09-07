# Code Review: ng_paralog_filter_d3

**Date:** 2026-09-07
**Reviewer:** rust-code-review skill (orchestrator), one sub-agent in an isolated worktree
**Scope:** step D3 of the hidden-duplication filter plan — the filter on HG002 at three hundred
reads a position, and the harness's stored-file mode
**Status:** Request-changes

---

### 1. Scope

- **What was reviewed:** the D3 run and its implementation report, plus the psp mode added to
  [`ng_paralog_filter_runs.sh`](../../../../scripts/ng_paralog_filter_runs.sh). No library code
  changed in this step.
- **Out of scope:** D1's and D2's code, reviewed at their own steps; production; the copied
  statistics.
- **One sub-agent**, in its own detached worktree at `86f5fe86` with the two uncommitted harness
  scripts copied in and the release binary rebuilt there, re-running the whole measurement.

### 2. Verdict

**Request-changes.** One Blocker, two Major, four Minor.

**Every count, hash and fitted parameter reproduced bit for bit** — 8,242 records, the sha256, the
fit at 302.45 reads and σ₀ 0.092, π, the cut, the 2,768 saturating ratios and their sign, all seven
ratio percentiles, the kind counts, the flagged set's composition, the timing and the peak. **What
failed was the conclusion drawn from them**, and the Blocker inverts the step's headline.

### 3. Execution status

No library code changed, so the gate is `86f5fe86`'s: lib 6,643 tests pass, 0 failed, 15 ignored;
one integration failure which is `main`'s; clippy 9 errors in three kinds; `fmt` dirty on 9 unique
files. Re-confirmed by the reviewer after rebuilding.

### 4. Open questions and assumptions

1. **Spec §3.2's rule that an indel carries the allele signal is contradicted at this depth.**
   **Put to the owner and ruled on, 2026-09-07: not acted on now** — "we'll study the giab problem
   and we'll fix it if necessary. The main objective right now is to have a working filter so we can
   carry out the real experiments and decide how to tweak it." Recorded in the report and in
   `PROJECT_STATUS`; no spec or code change.

### 5. Findings

#### Blocker

**B1: the report's central conclusion is contradicted by its own truth set once the method is fixed**

**Confidence:** High. **Category:** the measurement.

The report restricted the comparison to biallelic SNPs, quoted **0.8 removed in 100 were real
against a target of 1 in 100**, called the calibration passed, and set the all-kinds figure aside as
method-uncertain: "an indel can be written more than one way, so an exact key match undercounts
agreement".

**Resolving that method concern makes the number worse, not better.** Normalising both sides against
the same reference (`bcftools norm -f <analysis set> -m -any`) and intersecting:

| | records | in the GIAB truth set |
|---|---|---|
| substitutions removed | 265 | 2 — 0.8 in 100 |
| **indels removed** | **19** | **17 — 89 in 100** |
| everything removed | 284 | 19 — **6.7 in 100** |

**Seventeen of the nineteen indels the filter removed are GIAB benchmark variants.** Their
alternative-read fractions cluster on 0.32 to 0.37 — the alignment bias of a true heterozygous
indel at 300×, where reads carrying the indel are harder to place — and the model's `(T=3, m=1)`
carrier configuration expects exactly one third. **The indel arm is not finding duplications; it is
finding true heterozygous indels.**

The knob targets the removed records, and over the removed records it is missed sevenfold.
Re-derived independently by the orchestrator before being accepted.

#### Major

**M1: "the removed records are under-covered" is measured in the wrong units, and the direction reverses**

**Confidence:** High.

The report divided each window's depth by `single_copy_scale` (302.45) and called the quotient
copies. The scorer divides by `single_copy_scale × gc_multiplier(gc)`
([`coverage_model.rs`](../../../../src/ng/paralog/coverage_model.rs)), and **the removed records sit
at higher GC** — median 0.473 against the kept records' 0.406, where Illumina depth falls. Inverting
the model over the ratios puts the removed records at about **1.01 copies against the kept records'
0.99**, marginally above rather than 7% below.

**The conclusion survives and the number does not.** At one copy there is still no excess coverage
where a carrier is expected at 1.5, so "the coverage half does no work" holds — as a 2% difference,
not as an inversion. This is the same defect D1 had and D2 had: dividing by a fitted parameter and
calling the result copies.

**M2: the caveat that the callset is not GIAB's own truth VCF is false**

**Confidence:** High.

The report warned it was "the benchmark's curated callset … one step removed from the consortium's
release". Its header records
`view -R HG002_bench_azar_sorted_1000.bed … HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`, and it
carries GIAB's own `platformnames`, `callsetnames` and `arbitrated` fields. **It is GIAB v4.2.1
subset to the bench intervals and nothing else.** The caveat errs safe and tells a reader to
discount the only external check in the plan for a reason that does not exist.

The caveat that should have been made instead: 37 of the 265 removed-and-absent substitutions sit
in bench intervals carrying no truth call at all, where absence says least — 14 in 100 of them.

#### Minor

- **Mi1 — the arithmetic offered for *why* is wrong in both terms.** For 72 alternative reads of 470
  the report claimed ~150 nats of allele advantage and ~20 of coverage penalty, "an order of
  magnitude". Measured against the model: the best configuration is `(T=4, m=1)` at **111.5** and
  **34.4**, a margin of **3.2 to 1**.
- **Mi2 — "finding exactly the shape it was built to find" is not earned.** Across the range the 284
  occupy, the grid values `1/8, 1/6, 1/4, 1/3, 1/2` are never more than 0.042 apart, so every
  possible value lands near one; against a uniform null the clustering is 27% within 0.01 against
  28% expected. And the 284's median fraction is 0.249 — a quarter — not the sixth the ten highest
  happen to show.
- **Mi3 — the ten-highest table mixed two depths in adjacent columns.** The depth column was FORMAT
  `DP` while the fraction came from `AD`, which the scorer uses; they disagree on four of ten rows,
  worst at chr9:131,452,079 where 419 is printed beside `232, 46`.
- **Mi4 — the tract ceiling was asserted without its mechanism.** "No repeat tract was flagged" is
  true and uninformative on its own: a tract is scored on coverage alone, so **331 of the 580 share
  one ratio exactly (−0.6572)**, the highest is 1.4843 against a cut of 4.6500, and none could have
  been flagged. The reviewer's own figure for the bar was in model-copy units the run does not
  print, so the report states the measured ratios instead.

**One of the reviewer's minors did not reproduce.** It reported the kept heterozygotes' ninetieth
percentile alternative fraction as 0.528 over 4,864 records; the orchestrator measures 0.532 over
4,952. The difference is which records count as heterozygous, not an error in either — so the report
now states its selection rule rather than adopting a figure from a different one.

### 6. What is right, and was checked

Twenty-three quantities reproduced exactly, including the sha256, every fitted parameter, the
saturating-ratio counts and their sign, the kind counts, the timing, the spill's size and the
13.4 GB peak — which the reviewer confirmed is the calling pass and not the filter, the filter-off
arm peaking at 13,428 MB against 13,428. The claim that this is the first real-data run of
`call-from-psps`' copy of the filter wiring is supported by the documentation record. The
one-record disagreement between the recorded cut and the flag reproduces exactly: 4.6098 against a
quoted 4.6500, 283 of 284 at or above it, and no kept record above it.
