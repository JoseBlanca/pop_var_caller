# What routing on ng's own floors recovered, measured

*2026-09-02. Step B4 of
[`run_ssr_observations.md`](../ng/impl_plan/run_ssr_observations.md), against the prediction
in [`ng_str_path_losses_2026-09-02.md`](ng_str_path_losses_2026-09-02.md) §5. Branch
`ng-ssr-observations`; the runs are ng's own, at `--defaults`, on the GIAB `per_sample`
benchmark.*

---

## 1. The headline

**Pooled over the three samples at 30×, SNP recall went from 0.935 to 0.974 and indel recall
from 0.818 to 0.946 — for one extra false positive in total.** At 50× it is 0.940 to 0.979 and
0.818 to 0.949. Nothing in the caller changed: the same reads, the same catalog file, the same
`--defaults` model. What changed is which stretches of reference the run sends to a locus
generator that does not exist yet.

The loss report predicted ≈0.97 and ≈0.94 at 30× and called it an upper bound. **The measured
figures are at that bound**, 2,008 true SNPs against the bound's 2,006 and 312 true indels
against its 311 — a difference of about two sites in two thousand, which says the bound's
arithmetic was right rather than that anything beat it.

**ng's indel recall now exceeds the production caller's on this benchmark.** 0.946 against
0.930 at 30× — the production figure quoted from the loss report's own run of
`src/var_calling/`, same samples, same regions, same day. Its SNP recall is still below:
0.974 against 0.987.

---

## 2. What was run

`benchmarks/giab/src/run_ng_per_sample.sh` at `30x` and `50x`, all three samples, with the
release binary built from this branch. Each sample is called alone, restricted to its own
100-interval confident BED, and scored against its own GIAB truth VCF.

The pre-change calls were kept: `results/per_sample/<cov>/ng_before_routing_change/` is a copy
of what was on disk before this branch touched anything, from the same script and the same
inputs. **That they predate the change is checkable rather than asserted**: their
`.parameters.toml` files carry no `[repeat_routing]` section, which every run this branch makes
now writes.

Scoring is `benchmarks/giab/src/score_ng_recall.sh`, written for this comparison and following
`accuracy_dashboard.py`'s method exactly: both truth and query restricted to the BED,
left-aligned and split with `bcftools norm -m -any`, filtered to one class, intersected on
POS + REF + ALT. **Allele concordance, not genotype concordance** — a 0/1 called 1/1 at a real
site is a true positive here. Truth is FILTER PASS; the query keeps its own FILTER.

---

## 3. The ground that moved

Routing does not depend on depth, so these are one set of numbers, identical at 30× and 50×.

| sample | bases asked for | repeat before | repeat after | ratio | satellite before → after | typed regions before → after |
|---|---:|---:|---:|---:|---:|---:|
| HG002 | 571,984 | 32,577 (5.7%) | 4,930 (0.9%) | 6.6× | 0 → 180 | 7,392 → 806 |
| HG003 | 400,643 | 22,562 (5.6%) | 3,096 (0.8%) | 7.3× | 0 → 172 | 5,145 → 597 |
| HG004 | 453,250 | 25,129 (5.5%) | 3,601 (0.8%) | 7.0× | 0 → 102 | 5,603 → 642 |

**The three ratios are the loss report's own, to the digit** — it predicted 6.6×, 7.3× and
7.0×, and HG002's two figures, 32,577 and 4,930, are the exact pair it named. That is the
routing change doing precisely what was computed from the catalog, on the run rather than on
paper.

**The satellite class grew from nothing**, which the report also predicted: the calling floors
cap a tract at 100 bp where the catalog file allows 500, so three stretches per sample that
were repeat tracts are now tandem arrays too long to type as callable. That is a **permanent**
refusal rather than an unbuilt path — 180, 172 and 102 bases, against the 27,647, 19,466 and
21,528 bases that moved to the SNP/indel caller.

---

## 4. Recall and precision, pooled over the three samples

Truth counts are the same in every row: 2,061 SNPs and 330 indels across the three confident
region sets.

| depth | class | | true positives | false negatives | recall | false positives | precision |
|---|---|---|---:|---:|---:|---:|---:|
| 30× | SNPs | before | 1,926 | 135 | 0.9345 | 16 | 0.9918 |
| | | **after** | **2,008** | **53** | **0.9743** | 17 | 0.9916 |
| 30× | indels | before | 270 | 60 | 0.8182 | 3 | 0.9890 |
| | | **after** | **312** | **18** | **0.9455** | 3 | 0.9905 |
| 50× | SNPs | before | 1,937 | 124 | 0.9398 | 16 | 0.9918 |
| | | **after** | **2,018** | **43** | **0.9791** | 17 | 0.9916 |
| 50× | indels | before | 270 | 60 | 0.8182 | 3 | 0.9890 |
| | | **after** | **313** | **17** | **0.9485** | 3 | 0.9905 |

**82 more true SNPs and one more false one, at 30×; 42 more true indels and no more false
ones.** Precision is flat to four decimals in three of the four cells and rises slightly in the
fourth. That is what says the recovered ground is being called and not merely emitted over: a
change that bought recall by loosening something would show here.

Per sample, at 30×:

| sample | SNP recall before → after | indel recall before → after |
|---|---|---|
| HG002 | 0.9347 → 0.9654 | 0.8417 → 0.9353 |
| HG003 | 0.9453 → 0.9906 | 0.8222 → 0.9889 |
| HG004 | 0.9234 → 0.9688 | 0.7822 → 0.9208 |

**HG003 is nearly perfect and HG004 is the weakest, and the spread is worth not over-reading**:
each sample has its own 100 regions and its own truth set, so these are three different
problems rather than three tries at one. The pooled figure is the one to quote.

---

## 5. What this does not say

- **Nothing here is a repeat tract being called.** The tract path is still unbuilt; what moved
  is ground that was *classified* as repeat and is ordinary sequence under ng's own floors.
  4,930 bases of HG002 remain unreachable, and that is what Milestone C and the calling-loop
  plan are for.
- **Depth barely matters.** SNP recall gains 0.005 from 30× to 50× and indel recall 0.003. The
  routing gain is 0.040 and 0.127. So the loss this closed was never an evidence problem, which
  is what the loss report said before the fix and what the fix confirms.
- **This is three human samples over about 1.4 Mb of confident regions at high depth.** It says
  nothing about a low-coverage panel, where the same re-routing sends the same ground to a
  generic caller that has fewer reads to work with. `design_principles.md` §0's hard case is
  untouched by this measurement.
- **Where the frontier belongs is still open.** These floors are ng's measured stutter onsets,
  and the period × length question — including whether period-1 tracts belong on the repeat
  path at all — is spec §8's deferral, answerable cheaply once both paths run.

---

## 6. Reproducing it

```
./scripts/dev.sh cargo build --release --bin pop_var_caller_exp
./scripts/dev.sh bash -c "NG_BIN=\$PWD/target-container/release/pop_var_caller_exp \
    ./benchmarks/giab/src/run_ng_per_sample.sh 30x"
./scripts/dev.sh bash -c "./benchmarks/giab/src/score_ng_recall.sh 30x ng"
./scripts/dev.sh bash -c "./benchmarks/giab/src/score_ng_recall.sh 30x ng_before_routing_change"
```

The container build is a Linux binary, so both the run and the scoring happen inside the
container; the host cannot execute it. A whole coverage tier — three samples, GRCh38, 100
regions each — takes about 40 seconds.
