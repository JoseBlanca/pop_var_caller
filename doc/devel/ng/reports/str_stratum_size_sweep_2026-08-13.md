# How few repeat tracts a stratum can hold: five thousand at three reads a site, a thousand at six

*Research report, 2026-08-13. Answers the question
[`parameter_prepass_joint_loci.md`](../spec/parameter_prepass_joint_loci.md) §6 question 1 reopened on
2026-08-13 and [`parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md) §6.2
left open: a later design reads one stratum at a time to bound memory, so it must cap how many tracts
a stratum keeps, and where the estimate starts to hurt had never been measured from below. One
program stands behind this: `examples/ng_joint_str_harness.rs`, `size` mode. Raw output in
`tmp/str_stratum_size_*.log`.*

---

## 1. What was asked and what came back

**The cap should be 5,000 tracts a stratum.** Below that the estimate degrades, and it degrades in a
particular order.

- **At tomato's three reads a site, the floor is 5,000 tracts.** There every one of the five fitted
  numbers is within 2.4% of the drawn truth on average and no more than 2.3% wide between draws.
  (*Three reads a site* means three reads from each of the 20 samples at each tract.)
- **At six reads a site the floor is 1,000 tracts**, where the widest of the five is 4.3% between
  draws. Doubling the depth moves the floor by a factor of five — but only because the numbers a read
  carries improve; the concentration improves not at all (§5).
- **Two of the five break first: how fast two-repeat slips fall off against one-repeat slips, and the
  concentration.** Going from 1,000 tracts to 250 at three reads a site, the fall-off's spread between
  draws goes from 5.5% to 13.2% and the concentration's from 2.8% to 7.9%, while the slippage level's
  only goes from 2.5% to 3.9%. Which of the two is worse depends on the depth: at three reads a site
  the fall-off is the widest of the five at every tract count, and at six reads the concentration is
  — because doubling the reads halves the fall-off's scatter and does nothing at all for the
  concentration's.
- **Nothing is badly biased, even at 50 tracts.** Of the sixty cells in the two tables below, 56 have
  a mean within 5% of the truth and the worst is 9.2%. **What a small stratum costs is not a wrong
  answer on average but an answer that moves from draw to draw** — at 100 tracts and three reads a
  site, individual draws put the fall-off anywhere from 0.139 to 0.338 against a truth of 0.250.

**Why 5,000 and not the floor itself.** The floor at six reads a site is 1,000 and at three reads it
is 5,000, so 5,000 is the floor at the depth the caller is aimed at rather than a margin above it. It
also costs nothing to take: a tract is about ten bytes a read group, so the largest section of a
records file is 50 kB a sample, and 50 MB across a thousand samples
([`parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md) §6.2). At that
cap tomato keeps 86,688 of its 462,701 STR loci and 8 of its 141 strata are capped at all
([`parameter_prepass_joint_loci.md`](../spec/parameter_prepass_joint_loci.md) §4.5).

---

## 2. What was measured

The program draws a stratum with a known truth, draws each tract's length frequencies, draws each
sample's genotype under a supplied inbreeding coefficient, draws reads, and fits. It was already
there; what is new is the `size` mode, which refits the same truth at falling tract counts with
several draws apiece.

**The five numbers are reported separately and never pooled.** They are fitted together but they are
not learned at the same rate, so a single accuracy figure would call a tract count adequate once the
fastest-learned of them had arrived:

| number | what it is | drawn truth |
|---|---|---:|
| slippage level | how often a read reports a tract length other than its allele's | 0.0800 |
| shorter-share | of the reads that slip, the share showing a **shorter** tract | 0.830 |
| fall-off | how fast two-repeat slips fall off against one-repeat slips | 0.250 |
| concentration | how monomorphic the stratum's tracts are — small means most tracts carry one length | 0.500 |
| commonest length | the share of chromosomes at the stratum's commonest length; one number standing for the whole length spectrum | 0.476 |

**The estimator is the Dirichlet-over-points description of a tract**: a tract's length frequencies
are drawn from the stratum's own spectrum, and the fit integrates over them on a fixed low-discrepancy
point set. That is the description
[`joint_str_estimator_2026-08-12.md`](joint_str_estimator_2026-08-12.md) recommended in place of the
spec's corners-and-edges support, and it is the only one of the candidates that fits the concentration
at all. Its cost does not grow with the number of length classes.

**The stratum:** three length classes, 20 samples, an inbreeding coefficient of 0.4 supplied rather
than fitted, and a concentration of 0.5 — a stratum where 30% of tracts carry a single length across
the panel and 18% carry three or more.

**Several draws per tract count, because at a few hundred tracts the scatter between draws is most of
the answer.** Twelve draws at 50, 100 and 250 tracts; eight at 1,000; five at 5,000; three at 20,000.
Each cell below gives the mean error across draws, then the standard deviation between draws, then
the range — all as a percentage of the truth.

### 2.1 Commands

Built and run in the development container, from a worktree on branch `ng-str-stratum-size`:

```sh
./scripts/dev.sh cargo build --release --example ng_joint_str_harness

# size <samples> <reads a site> <tract count, 0 for every count> <length classes> <concentration>
./scripts/dev.sh ./target-container/release/examples/ng_joint_str_harness size 20 3 0 3 0.5 \
  > tmp/str_stratum_size_depth3.log
./scripts/dev.sh ./target-container/release/examples/ng_joint_str_harness size 20 6 0 3 0.5 \
  > tmp/str_stratum_size_depth6.log

# One tract count on its own — how the two largest cells were re-run after a killed process
./scripts/dev.sh ./target-container/release/examples/ng_joint_str_harness size 20 3  5000 3 0.5
./scripts/dev.sh ./target-container/release/examples/ng_joint_str_harness size 20 3 20000 3 0.5
```

---

## 3. The sweep

### 3.1 Three reads a site — tomato's depth

| tracts | draws | slippage level | shorter-share | fall-off | concentration | commonest length |
|---:|---:|---:|---:|---:|---:|---:|
| 20,000 | 3 | −0.2% ± 0.4% | −0.0% ± 0.2% | −0.6% ± 1.5% | +0.5% ± 0.4% | −0.7% ± 0.4% |
| **5,000** | 5 | +0.6% ± 0.5% | −0.1% ± 0.5% | −2.4% ± 2.3% | +0.1% ± 1.2% | −0.3% ± 1.3% |
| 1,000 | 8 | +0.7% ± 2.5% | −0.5% ± 1.1% | −3.4% ± **5.5%** | −0.5% ± 2.8% | −0.4% ± 2.4% |
| 250 | 12 | +1.5% ± 3.9% | +0.2% ± 1.7% | −5.6% ± **13.2%** | +0.3% ± **7.9%** | −0.1% ± **5.8%** |
| 100 | 12 | +3.1% ± **7.1%** | +1.0% ± 3.3% | −0.2% ± **21.5%** | −4.3% ± **14.3%** | +0.0% ± **6.6%** |
| 50 | 12 | +2.8% ± **9.0%** | +1.4% ± **5.0%** | +4.1% ± **31.3%** | −4.3% ± **17.1%** | −4.0% ± **9.4%** |

### 3.2 Six reads a site

| tracts | draws | slippage level | shorter-share | fall-off | concentration | commonest length |
|---:|---:|---:|---:|---:|---:|---:|
| 20,000 | 3 | −0.0% ± 0.2% | +0.1% ± 0.1% | +0.6% ± 0.9% | −0.5% ± 0.3% | +0.2% ± 0.5% |
| 5,000 | 5 | −0.0% ± 0.4% | +0.1% ± 0.5% | −0.7% ± 0.9% | +1.0% ± 1.6% | +0.0% ± 0.8% |
| **1,000** | 8 | +0.4% ± 1.7% | +0.2% ± 0.6% | +1.6% ± 4.2% | +2.6% ± 4.3% | +1.7% ± 1.8% |
| 250 | 12 | −1.0% ± 2.1% | +0.4% ± 1.2% | +1.8% ± **8.1%** | +0.1% ± **8.4%** | +2.2% ± 4.3% |
| 100 | 12 | −2.0% ± 3.0% | +0.2% ± 1.8% | +4.8% ± **11.3%** | −5.3% ± **14.2%** | +3.9% ± **7.9%** |
| 50 | 12 | −1.1% ± 4.0% | +0.1% ± 2.8% | +9.2% ± **15.9%** | −3.2% ± **27.8%** | +7.1% ± **15.1%** |

**How to read these.** The bolded **row** in each table is the floor: the smallest count where every
one of the five is within a few percent of the truth both on average and between draws, taking "a few
percent" as **5% on the mean and 5% on the spread**. A bolded **spread** below that row is one over
the 5% bar. Reading each table downwards, the fall-off or the concentration is always the first to
cross it, the commonest-length share next; the slippage level does not cross until 100 tracts and the
shorter-share not until 50, and neither crosses anywhere at six reads a site.

---

## 4. Where the floor is, and what breaks first

**Read the two tables down the fall-off column.** At three reads a site it is 2.3% wide between draws
at 5,000 tracts, 5.5% at 1,000, 13.2% at 250, and 31.3% at 50. At 250 tracts the five spreads rank
fall-off 13.2%, concentration 7.9%, commonest length 5.8%, slippage level 3.9%, shorter-share 1.7%,
and **the fall-off is the widest of the five at every tract count in that table, the concentration
second at every count of 1,000 and below.**

**At six reads a site those two swap.** The fall-off halves with the extra reads and the concentration
does not, so from 5,000 tracts down it is the concentration that is widest: 27.8% against the
fall-off's 15.9% at 50 tracts, 14.2% against 11.3% at 100, 8.4% against 8.1% at 250.

**Either way it is the same two columns that go first, and the slippage level is never one of them.**
At 250 tracts, where both are already over the bar, the slippage level is still within 4% at three
reads a site and within 2.1% at six.

**Why those two.** They are the numbers with the least data behind them.

- The **slippage level** is carried by every read at every tract: 20 samples × 3 reads × 250 tracts is
  15,000 reads, and about 1,200 of them slip.
- The **fall-off** is carried only by the reads that slip by *two* repeats rather than one. At a
  fall-off of 0.25 those are about a fifth of the slipped reads, so roughly 225 reads at 250 tracts
  and 45 at 50 tracts. **That is the whole reason it goes first** — the count it rests on is the
  smallest, and it is smaller than the slippage level's by a factor of five.
- The **concentration** is not carried by reads at all. It is a statement about how unlike each other
  the tracts are, so a tract is one observation of it however deeply that tract is read.

**What a broken column looks like in the draws themselves.** At 100 tracts and three reads a site the
twelve draws put the fall-off at 0.139, 0.205, 0.216, 0.221, 0.240, 0.243, 0.246, 0.247, 0.284, 0.298,
0.314 and 0.338 against a truth of 0.250 — the two extremes are 44% low and 35% high. At 1,000 tracts
the same twelve-way scatter is a 0.224-to-0.267 spread, which is a number one can use.

---

## 5. Twice the depth buys back three of the five and not the other two

Comparing the two tables at the same tract count separates what more reads buy from what more tracts
buy. **They are not the same thing, and two of the five cannot be bought with reads at all.**

| at 100 tracts, spread between draws | three reads | six reads |
|---|---:|---:|
| slippage level | 7.1% | **3.0%** |
| shorter-share | 3.3% | **1.8%** |
| fall-off | 21.5% | **11.3%** |
| concentration | 14.3% | 14.2% |
| commonest length | 6.6% | 7.9% |

**The three read-driven numbers roughly halve when the reads double**, which is what doubling a count
does to its scatter. **The concentration does not move at all** — 14.3% against 14.2% at 100 tracts,
7.9% against 8.4% at 250. The commonest-length share does not halve either: 6.6% against 7.9% at 100
tracts and 5.8% against 4.3% at 250, which is scatter around no change. Both are counted in tracts —
a tract is one observation of how unlike its neighbours it is, however deeply it is read — so both can
only be bought with tracts.

That matters for the cap in a way the other columns do not: **a deeper run does not relieve the cap
for the concentration.** A cohort sequenced at 30× still needs the tracts, and the cap is the only
knob that supplies them.

---

## 6. What this cannot say

- **The truth is drawn, not real.** No real reads are in this, and which concentration a real tomato
  stratum has is still unmeasured. The repeat-tract records exist from real alignments
  ([`joint_records_on_real_alignments_2026-08-13.md`](joint_records_on_real_alignments_2026-08-13.md)),
  but the repeat-tract half of the estimator has never been run on them — only the ordinary-position
  half has ([`joint_fit_against_truth_2026-08-13.md`](joint_fit_against_truth_2026-08-13.md)). **Until
  that happens the cap rests on a drawn stratum**, and the concentration a real one carries is exactly
  the number that would move the floor.
- **One stratum shape.** Three length classes and a concentration of 0.5. A stratum whose tracts are
  more nearly monomorphic has less signal per tract for the concentration, and the floor there could
  be higher; a stratum with more length classes spreads the same reads over more classes, which
  should push the same way. Neither was run.
- **The floor is stated at 20 samples.** How the tract count and the sample count trade against each
  other was not swept — the panel is fixed at 20 throughout, which is the size the earlier comparisons
  used.
- **Three draws at the top of the range.** The 20,000-tract cells are three draws each, because each
  fit takes about three minutes there. They confirm that nothing is wrong at the top rather than
  measuring its scatter precisely; the counts that decide the cap have eight and twelve.
- **Nothing here says the cap is safe for the thin strata.** 68 of tomato's 141 strata hold fewer than
  a hundred tracts each ([`parameter_prepass_joint_loci.md`](../spec/parameter_prepass_joint_loci.md)
  §4.5), which is below every floor in this report. **A cap does not touch them** — they are already
  under it — and their answer remains borrowing from a neighbouring repeat count
  ([`parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) §4.3). What this report does add is
  a number for when borrowing is needed: **a stratum under 1,000 tracts cannot carry its own fall-off
  at three reads a site**, where the existing rule fires on a much thinner stratum than that.

---

## 7. Changes I would propose to the spec and architecture documents

None of these were made; they are what I would put in a follow-up.

1. **[`spec/parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md) §6.2** —
   the paragraph beginning *"What the cap should be is not settled"* ends *"Nobody has swept
   downwards, so where it starts to hurt is unknown, and a cap above a few thousand buys nothing
   anyone has measured."* Replace the last sentence with the measured floor: **5,000 tracts at three
   reads a site and 1,000 at six**, the two numbers that break first, and the observation that the
   concentration does not improve with depth so a deep cohort needs the same cap as a shallow one.

2. **[`spec/parameter_prepass_joint_loci.md`](../spec/parameter_prepass_joint_loci.md) §6 question 1**
   — close the part reopened on 2026-08-13. Its *leaning* was *"a cap of a few thousand, inside the
   measured range"*, and that turns out to be right with 5,000 under it. Record which of the five
   numbers sets it, because the reason is not the slippage level everyone would assume.

3. **[`spec/parameter_prepass_joint_loci.md`](../spec/parameter_prepass_joint_loci.md) §4.5** — the
   cap table lists 100, 500, 1,000, 5,000 and 20,000. Add what each now means for the estimate: 100
   and 500 are below the floor at both depths, 1,000 is the floor at six reads a site only, and
   5,000 is the floor at three. That turns a table of locus counts into a table a reader can choose
   from.

4. **[`spec/parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §4.1** — add one
   paragraph: the concentration is counted in tracts and not in reads, so it, and not the slippage
   numbers, is what sets how many tracts a stratum needs. §4.1 currently motivates the concentration
   entirely as the thing that makes *per tract or per stratum* one fitted number, and this gives it a
   second consequence that a run has to act on.

5. **[`spec/parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §11 question 10**
   — the lever *"a stratum is fitted on its own"* is listed as free of accuracy cost, which is true of
   the reading order but not of the cap that turns it into a bound. Note that the bound costs nothing
   at 5,000 tracts and costs the fall-off below 1,000.

6. **[`arch/parameter_prepass_joint_records.md`](../arch/parameter_prepass_joint_records.md) §2.2** —
   the per-stratum read handle is specified without saying what bounds one stratum's section, which is
   the memory guarantee the whole by-section shape exists for. Give the bound: the per-stratum cap,
   about ten bytes a tract a read group, so **50 kB a sample at a cap of 5,000**.
