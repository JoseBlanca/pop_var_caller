# The hidden-duplication filter — D1: one sample

**Date:** 2026-09-07
**Plan step:** [hidden_paralog_filter.md](../../ng/impl_plan/hidden_paralog_filter.md) Milestone D, step D1
**Spec:** [hidden_paralog_filter.md](../../ng/spec/hidden_paralog_filter.md) §4, §9
**Branch:** `ng-paralog-filter`

## The answer

**The filter runs at one sample, fits a coverage model for it, scores every record it wrote, and
removes one of them. Nothing fell back and no ratio came back non-finite.** Spec §9's second OPEN —
whether the copied precompute works at `N = 1` — is closed on real reads.

One tomato accession over the plan's two 100 kb intervals:

| | |
|---|---|
| records the run would have written without the filter | 217 |
| records scored | 217 |
| records that could not be scored | 0 |
| records removed at `--paralog-fdr 0.01` | 1 |
| fitted duplication rate `π` | 0.012594, **fitted from this run** — the `0.03` fallback did not fire |
| the cut | likelihood ratio 8.2500 |
| samples with a coverage model | 1 of 1 |
| **the fit itself** | **one copy at 5.22 reads a window, scatter σ₀ 0.300** |
| ratios past the histogram's ±100 ends | 0 |

The header says the calibration:

```text
##paralogFilter=target_fdr=0.0100;pi=0.012594;lr_cut=8.2500;em_converged=true;samples_with_coverage_model=1/1
```

## What was run, and on what

**One accession, not a cohort of one drawn from a cohort run.** `SRR7279481.p1.bench.cram` —
`bench.config.sh`'s own `SINGLE_CRAM`, sample `SRS3394712` — over the first two intervals of
`benchmarks/tomato1/regions.bed`, `SL4.0ch01:3,406,886-3,506,886` and
`SL4.0ch01:13,806,669-13,906,669`. Called with `--defaults` (no fitted parameters file) on four
threads, through `call-from-alignments`.

The cohort is called three ways every time, by
[`scripts/ng_paralog_filter_runs.sh`](../../../../scripts/ng_paralog_filter_runs.sh): the filter
off, the filter dropping what it flags, and the filter tagging it instead. **The tagging run is
where a flagged record can still be read** — a dropped record leaves no trace of its ratio in the
file — so everything counted per flagged record comes from that run.

**Depth.** The window each record's score read has a median mean depth of **8.49 reads** over the
217 records (mean 8.30, lowest 1.20, tenth percentile 3.29, ninetieth 12.83, highest 19.13).
**That is the depth at the 217 sites where a variant was written, not the accession's depth over
the 200 kb** — variant sites are where excess coverage collects, so it is a biased estimate and is
not comparable with `CLAUDE.md`'s three-reads-a-position figure for this cohort. The unbiased
number is the fit's own: **one copy is 5.22 reads a window**, which the run now prints.

**Percentiles throughout this report are the value at the sorted index `int(p·n)`** — a value some
record really has, not an interpolation between two. Recomputing the tables with `numpy.percentile`
gives different numbers at the fifth significant figure for p95, p99 and the depth p90.

**What the 217 records are:** 213 biallelic SNPs, 2 deletions, 2 repeat tracts. The sample is
called `1/1` at 202 of them and `0/1` at 15 — an inbred accession against the reference, which is
why homozygous-alternative records dominate every ranking below without that meaning anything
about the score.

## The record it removed, and the one it kept next to it

```text
SL4.0ch01  3471589  G→T  QUAL 378.4  FILTER hiddenParalog
           PARALOG_LR=21.5702  PARALOG_POST=1.000000
           GT 1/1  DP 15  AD 0,15
```

**Its window carries 19.01 reads against a fitted one-copy level of 5.22 — about 3.6 one-copy
levels before the GC correction, and close to the model's winsor cap of four copies.** That is the
duplication footprint, and at one sample it is the whole of the evidence: with no cohort there is
no *"the heterozygous samples are the over-covered ones"* coincidence to exploit, so the coverage
half does the work alone.

**The record ranked next was kept, and its window is marginally deeper still**: 19.13 reads, ratio
**8.1116** against the cut of 8.2500. The two are **13.4586 apart on the ratio axis** — a factor of
about e¹³·⁵ between the two likelihood ratios, not a near-tie. What is nearly a tie is the kept
record against the *cut*: it sits **0.1384** below it.

**Two things differ between them, and both push the same way.**

- **The read counts: `AD 0,15` against `AD 0,19`.** Both are fully homozygous-alternative, so under
  the real-variant story each read costs `ln(1 − ε) = −0.01005` at the model's error floor
  `ε = 0.01`, while the duplication story's best branch is `vaf = ½` and costs `ln 0.5 = −0.6931`.
  The model enumerates only the seven configurations `(T,m) = (3,1),(4,1),(4,2),(6,1),(6,3),(8,1),(8,4)`,
  so **`m/T` never exceeds ½** and no better branch exists. Four extra reads therefore cost the kept
  record `4 × 0.6831 = 2.73` nats — a fifth of the gap — before coverage is considered.
- **The GC their windows sit at: 0.453 against 0.367.** The fit's expected one-copy depth is the
  fitted level times a GC multiplier, so the same raw depth is a different copy number at the two
  places.

**How the remaining ~10.7 nats split is not readable from this run**, because the report prints the
fitted one-copy *level* and σ₀ but not the GC multiplier curve, so a record's exact copy number
cannot be recovered from the output. Reviewing this step, a sub-agent recovered it by inverting the
scorer over all 217 ratios: it puts the expected one-copy depth at **5.51 reads at GC 0.453 and
6.33 at GC 0.367**, giving copy numbers of **3.45 and 3.02**, and — the striking part — **had the
kept record sat at the flagged record's GC its ratio would have been about 19.4 and it would have
been flagged comfortably**. That reconstruction agreed with the run on the one number the run does
print: **σ₀ = 0.300, to three decimals.** It is a reconstruction and not a measurement, and it is
recorded here as such; what it says, if right, is that the calibration at one sample is soft **in
GC** rather than soft in general.

**The ten highest ratios**, with the window each score read. The multiple of one copy is against the
run's own fitted level of 5.22 reads and is **before** the GC correction, so it is an upper bound on
the copy number where the GC multiplier exceeds one:

| ratio | window depth | × the fitted one copy | GC | position | genotype |
|---|---|---|---|---|---|
| 21.570 | 19.01 | 3.64 | 0.453 | 3,471,589 G→T | 1/1, AD 0,15 |
| 8.112 | 19.13 | 3.66 | 0.367 | 13,874,388 C→T | 1/1, AD 0,19 |
| 3.790 | 16.85 | 3.23 | 0.433 | 3,470,328 G→A | 1/1, AD 0,16 |
| 2.590 | 13.95 | 2.67 | 0.453 | 3,467,331 A→G | 1/1, AD 0,13 |
| 1.127 | 13.97 | 2.68 | 0.423 | 3,468,385 T→C | 1/1, AD 0,8 |
| 1.107 | 9.65 | 1.85 | 0.292 | 3,429,960 **(repeat tract)** | 1/1, AD 0,9 |
| 1.087 | 10.59 | 2.03 | 0.423 | 13,881,841 A→C | 0/1, AD 6,4 |
| 0.789 | 15.34 | 2.94 | 0.418 | 3,468,706 A→T | 1/1, AD 0,13 |
| 0.767 | 16.63 | 3.19 | 0.367 | 13,874,104 G→A | 1/1, AD 0,18 |
| 0.362 | 12.36 | 2.37 | 0.477 | 3,467,229 C→T | 1/1, AD 0,11 |

The five lowest ratios sit at 1.00 to 1.92 times the fitted one-copy level (−12.704, −11.766,
−10.434, −10.403, −10.382).

**One of the ten is a repeat tract**, scored on coverage alone as spec §3.2 requires, and it lands
mid-table at 1.107. Two tracts in 217 records is too few to say anything about tracts; step D4 is
where that count is taken, on the six-accession run and on HG002.

## The whole ratio distribution

Over the 217 scored records:

| | ratio |
|---|---|
| lowest | −12.7036 |
| 25th percentile | −7.8277 |
| median | −6.2452 |
| 75th percentile | −4.5592 |
| 95th percentile | −0.0566 |
| 99th percentile | 3.7904 |
| highest | 21.5702 |

**Every one of the 217 is finite**, which is the `N = 1` question this step exists to answer: the
copied precompute's SFS grid degenerates to a single point at one sample (`[1/2N, 1 − 1/2N]` is
`[½, ½]`), and the run produced no `NaN` and nothing at either infinity. The run report's
`0 could not be [scored]` is the same fact counted by the code rather than read off the file.

**Nothing reached the histogram's ends.** The run printed no line about ratios past `[-100, 100]`
— the line exists and is printed only when the count is above zero, so its absence is a zero. The
concern C3's review recorded, that the ratio grows with the cohort while the histogram's range is
fixed, does not bite at one sample: the largest ratio here is 21.6, a fifth of the way to the edge.
That question is D2's and D3's.

**Depth does not order the ranking on its own.** Over the 217 records the correlation between the
window's mean depth and the ratio is **Pearson 0.456** — so a straight line through depth accounts
for **21%** of the variance in the ratio — and **Spearman 0.364**. Two things account for the rest,
and neither is the inbreeding coefficient: `ParalogScorePrecompute` is built once per pass from a
cohort-length slice of coefficients, so at one sample `F` is a single number applied identically to
all 217 records and **cannot move any record relative to any other**. What can:

- **the ratio is strongly convex in copy number, not linear**, so a linear correlation understates
  what is close to a deterministic relationship;
- **the read count rises with depth and pushes the ratio the other way**, by the 0.68 nats a read
  that the two records at the cut demonstrate.

## The three relations between the files, on real reads with something flagged

All three hold, checked by
[`scripts/ng_paralog_filter_relations.py`](../../../../scripts/ng_paralog_filter_relations.py):

- **the filter changed no record it did not flag** — the dropping run, with its two INFO keys taken
  back out, is the filter-off file with exactly the dropped record missing (256 lines compared);
- **dropping is tagging minus the tagged** — the 216 records the tagging run did not flag are, byte
  for byte, the 216 the dropping run wrote;
- **a flagged record is its off-run record with the verdict spliced in** — undoing the two INFO keys
  and the `FILTER` splice on the one flagged record gives back the off run's line exactly.

**This is the first time any of them has been checked with a record actually flagged.** C5 ran the
first two on the mode-equivalence fixture, where the cohort is too small for any sample to get a
coverage model, so both were assertions about the empty case (C5's report says so). The third did
not exist; this step's review found that the first two both work by *removing* the flagged records
from one side, so a flagged record's QUAL, INFO, FORMAT and sample columns were compared against
nothing at all.

## C5's missing oracle, closed on the way in

**Spec §10's first oracle holds: `--paralog-fdr 0` reproduces the pre-filter run byte for byte.**
Six accessions over the same two intervals give **2,311 records** and sha256
`84ad19c22dd14de583cd85805dcd2e5169e799d7a63691c979b7fa43d400590d` on everything but
`##commandline` — the plan's standing baseline, recorded on 2026-09-06 before the filter existed
and unchanged since.

This was the one of spec §10's four oracles C5 could not run, and the plan makes it the ground D
stands on, so it was run first.

**Matching it took one thing worth writing down: the output's own filename is inside the file.**
The header carries `##parametersFile=<basename>`, so a run written as `off.vcf` differs from the
baseline in that one line whatever the filter did — sha256 `9f6a0367…` on the first attempt, against
`84ad19c2…` when written under the baseline's own name. That one line is the only difference in the
whole file. `ng_paralog_filter_runs.sh` writes every run as `run.vcf` and moves it afterwards for
this reason.

## What this step changed in the caller

**The run report now says what each sample's coverage fit came to** — spec §3.1's outcome, which
plan step D1 asks for and which nothing printed before. One sample:

```text
  sample SRS3394712 fitted one copy at 5.22 reads a window, scatter 0.300
```

**This was added because the first draft of this report got the mechanism wrong without it.** With
no fitted level printed, the obvious stand-in for "one copy" is the median depth of the records the
run wrote — and that is 8.49 against the fit's 5.22, so the flagged record's window read as 2.24
copies when the model sees about 3.5. The number a reader carries away was wrong by a factor of
1.6, and nothing in the run's output could have corrected it.

- [`ParalogScoringContext::what_each_fit_came_to`](../../../../src/ng/run/paralog_filter/scoring_context.rs)
  returns one `WhatTheFitCameTo` per sample — the fitted one-copy depth and σ₀ — or `None` where the
  fit was refused. **A small owned summary rather than the models themselves**: a caller holding a
  `SingleCopyCoverageModel` could ask it for a copy number at a GC of its choosing and index the
  answer against the wrong sample, which is the shape of spec §6's trap 3, and C1's review
  deliberately made those accessors private. Two numbers per sample cannot be misindexed into a
  score.
- [`what_to_tell_the_operator`](../../../../src/ng/run/paralog_filter/finish.rs) names every fitted
  sample up to ten and **gives the spread past that** — the lowest, the median and the highest of
  both quantities. A run is up to several thousand samples (spec §4) and a report is read by a
  person; what it must not do is print the first ten and leave a reader believing that is the
  cohort.

**Two tests, and seven mutations run against them, all seven killed**: the ends of the spread
swapped, the median taken as the lowest, the cap raised past any cohort, the scatter printed where
the depth belongs, every line naming the first sample, a sample's depth and scatter swapped, and
the fit reporting σ₀ where the depth belongs. The fixtures give each sample a **different** fitted
depth for exactly this reason — a cohort fitted at one depth cannot tell a line that prints each
sample's own number from one that prints the first sample's number twelve times.

## What this step adds beside the caller

| file | what it is |
|---|---|
| [`scripts/ng_paralog_filter_runs.sh`](../../../../scripts/ng_paralog_filter_runs.sh) | calls one cohort three ways — off, dropping, tagging — and prints the records written, the filter's report lines, the header's calibration line, and each run's wall and peak resident |
| [`scripts/ng_paralog_filter_relations.py`](../../../../scripts/ng_paralog_filter_relations.py) | the three relations above, checked on real reads |

**Both were reviewed by being broken on purpose.** Sixteen deliberately corrupted inputs were fed
to the first version; nine were caught and seven were not. The four that mattered are fixed and the
repaired tool catches all ten cases it is now tested against:

- **neither relation could see the filter losing its own INFO keys or header lines** — both strip
  them before comparing, so a filter that stopped writing them passed clean. They are now asserted
  present on the two filtering runs and asserted **absent** on the off run, and a record carrying
  one of the two keys without the other is a failure.
- **a flagged record could be corrupted in any column undetected** — the third relation above is
  the fix.
- **records were matched by `(CHROM, POS, REF, ALT)`, which is not an identity.** Spec §3.2 says a
  repeat tract may share a position with the generic locus owning its anchor base, and two such
  records would make the check wrong in both directions — a false failure when one twin is flagged,
  a false pass when both are wrongly removed. It now walks the three runs positionally with a length
  check, which is exact under duplicates and shorter. **Not reachable on this data**: over 8,510
  records from six accessions there is no repeated `(CHROM, POS)` at all.
- **the harness's report section could go silently empty** — `sed` prints nothing and exits 0 when
  its marker is absent, so rewording the run report's first line would have taken π, the cut, the
  convergence flag and the drop count out of every report written from that output without a word.
  An empty report is now a failure.

Smaller ones applied: the harness names the binary it measured and when it was built; peak resident
converts units rather than printing a macOS figure a thousand times too large, and asserts the
`ru_maxrss` mark it reads is its own child's; a run that wrote no records says so rather than
printing three clean zeros; the `hiddenParalog` count matches whole `FILTER` ids rather than a
substring; and the off run's calibration line being present, or the filtering runs' being absent,
is a failure rather than a printed note.

## How to reproduce the depth numbers

The window each record's score read is **not in the VCF**. To get it, set
`NG_WINDOW_COVERAGE_FILE=<path>` on the run; the rows are
`sample \t contig \t position \t gc \t mean_depth`, and the last two are **`f32` bit patterns
written as decimal `u32`**
([`recorded_windows.rs`](../../../../src/ng/run/cohort_merge/recorded_windows.rs)). It must be the
**tagging** run's dump, since the dropping run is one record short.
[`tmp/d_milestone/analyse_d1.py`](../../../../tmp/d_milestone/analyse_d1.py) is the decoder every
table above came from.

## Deviations from the plan, recorded

1. **The reads were on this machine after all.** C5's report and the `PROJECT_STATUS` entry say
   `benchmarks/tomato1/crams/` is empty and that D1 to D4 wait on files this machine does not have.
   The directory is empty **in this worktree**, because `benchmarks/tomato1/crams/` is in
   `.gitignore` (line 41) and a git worktree carries no ignored files. The 63 CRAM slices and the
   repeat catalogue have been in the main checkout since June, and the container reaches them with
   `DEV_EXTRA_MOUNT=/Users/jose/devel/pop_var_caller/benchmarks`. Nothing was fetched and nothing
   was substituted.
2. **The caller gained a report line, which a measurement step would not normally do.** Plan step
   D1 asks for "the fit's outcome" and the run reported only how many fits were accepted; the cost
   of not printing it is measured above. Absorbed rather than escalated: it adds no behaviour, moves
   no VCF byte, and the standing oracle still gives 2,311 records at sha256 `84ad19c2…`.
3. **Wall and peak resident are read with `getrusage`, not GNU `time`.** The container image has no
   `/usr/bin/time`, and its `/bin/sh` is dash.
4. **Per-pass wall is not in this step.** D2's ask, and it needs the run to time itself. Whole-run
   figures: each of the three runs took **3.0 s**. Peak resident moved between **260.9 and 267.7 MB**
   over nine runs of the same three commands, with the dropping run highest every time — so the
   filter's cost at this size is inside the run-to-run noise, which is a statement about 217 records
   and one sample and not about the caller.

## What D1 leaves open

- **Whether the removed record is a real collapsed duplication is not knowable from one sample.**
  Its window is at about 3.5 times one copy's fitted depth, which is the footprint; at one sample
  there is no second signal to corroborate it. D2 is the first run where the coincidence the model
  is built on — the heterozygous samples being the over-covered ones — can exist at all.
- **The GC multiplier curve is still not printed**, so a record's exact copy number cannot be
  recovered from a run's output; only the level it is a multiple of. **For D4**, which counts
  flagged records by kind and wants their ratios explained, recording the relative copy number the
  scorer actually used — per record, beside the ratio — would remove the last reconstruction.
- **Two repeat tracts in 217 records** is not a measurement of anything about tracts. D4.
