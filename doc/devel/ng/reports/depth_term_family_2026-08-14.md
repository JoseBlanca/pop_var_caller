# Depth is not Poisson where the copy-number term would be used: its variance is 2.5 to 3.2 times its mean, and a Poisson term hands 1,000-to-1 odds of a duplication to 1 position in 74

*Research report, 2026-08-14. **Two recommendations. The term's `P` must be a negative binomial
with one dispersion number a sample, never a Poisson. The log-ratio is worth having in place of
the threshold — it ranks positions better at every depth measured — but only as a weight, never
as a yes-or-no test, because the point where it crosses zero flags three times as many positions
as the threshold did.** One program stands behind this,
`examples/ng_depth_term_family.rs`; eight tomato accessions from 2.8 to 30.2 reads a position,
raw output in `tmp/depthterm/`. Its threshold arm reproduces
[`locus_depth_vs_window_2026-08-13.md`](locus_depth_vs_window_2026-08-13.md) §4 to within 0.1 of
an enrichment at all eight accessions, which is what says the two sets of numbers can be
compared.*

---

## 1. The question, and the answer

Before calling anything, the caller estimates how heterozygous each sample is. A stretch of
genome a plant carries twice while the reference carries it once corrupts that estimate: both
copies' reads pile onto the one place the reference offers, and wherever the two copies differ
the position reads about half non-reference — which is what a heterozygote looks like. With
twenty-five samples or more the cohort recognises such a position for free. Below that, and at
one sample, the only evidence left is that the position carries about twice the reads the sample
normally has there.

That evidence is about to enter the estimator as a term in a likelihood,

```text
    ln P(d | 2m) − ln P(d | m)
```

where **`d` is the position's read count** and **`m` is the depth one copy is expected to give at
that position** — the sample's own depth at that position's GC content. Two things about it were
unsettled.

**Is `P` a Poisson?** No, and not by a small margin. At the depths where this term would be used
— 26.4 and 30.2 reads a position on the two deepest accessions — **the variance of depth is 2.5
and 3.2 times its mean**, measured inside a GC bin, over 7.5 million positions. A Poisson asserts
those two numbers are equal. The consequence is not a mis-tuned parameter: at 30.2 reads a
position a Poisson term gives **at least 1,000-to-1 odds of two copies to 136 positions in every
10,000 it scores**, one position in 74. A negative binomial fitted to the same data gives those
odds to 39 in 10,000. For scale, **5 positions in 10,000 read near half at all**.

**Does the log-ratio separate better than the threshold?** Its *ranking* does, at every accession
measured. Given the same number of positions to flag as the threshold arm took, the log-ratio
puts near-half positions **24.5 times above chance at 30.2 reads a position against the
threshold's 17.1**, and 14.8 against 12.4 at 26.4. **But its own zero is a bad place to cut.**
Scored as *positive means duplicated*, it flags 6.1 in every 100 positions where the threshold
flagged 1.8, and its enrichment falls to 9.0. At 2.8 reads a position it puts 58 positions in
100 above zero. So the term earns its place as a graded weight and not as a classification.

---

## 2. What is being counted, and what it is not

**Variance ÷ mean** over the read counts of single positions, inside one bin of GC content, on
one accession. It is 1 for a Poisson by definition, so the measurement needs no model and no
fit — it needs only that the positions compared are alike in what sets their expected depth.
That is why it is measured inside a GC bin: per-position depth on one accession runs from 23.7
reads at 16% GC to 33.6 at 30%, and a variance taken across that range would be measuring the GC
curve.

**The read counts are the exact ones the walk saw.** The joint records store depth as a five-bit
code standing for a range of counts — exact to eight reads, then eleven widening rungs — so a
dispersion measured on the stored code would be a property of the ladder. Nothing here reads the
ladder.

**Near half** means a position whose reads disagree with the reference in a fraction between 0.35
and 0.65, with at least two disagreeing reads. That is the artefact's signature, and the same
definition the earlier report used.

**Enrichment** is how many positions an arm flags *and* that read near half, over what the
flagging rate and the near-half rate together would give if the two were unrelated. 1.0 means the
arm knows nothing.

**There is no truth set.** Nobody has a validated list of the stretches these accessions carry
twice. Enrichment compares arms on the same positions; it is not a detection rate, and a flagged
position is not thereby a duplication.

**One difference in bookkeeping from the earlier report.** Mean depth here is averaged over the
positions the walk emitted a locus for; the earlier report divided by every generic base,
including those no read reached. The same accession is therefore 30.2 reads a position here and
28.7 there. Nothing else changed: the enrichments reproduce (§7).

---

## 3. Variance against mean, inside a GC bin

**SRR7279540, 30.2 reads a position** — GC bins two percentage points wide, bins under 5,000
positions dropped. The standard error is computed from each bin's own fourth moment, not assumed.

| GC% | positions | mean depth | variance | variance ÷ mean | with the unreached positions back at zero | dropping the near-half positions | deepest 1% of positions hold this share of the variance | a Poisson would put | a fitted negative binomial would put |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 16 | 7,057 | 23.67 | 78.8 | 3.33 ± 0.06 | 3.33 | 3.32 | 9.3% | 8.4% | 11.3% |
| 20 | 70,968 | 29.55 | 66.9 | 2.26 ± 0.02 | 2.63 | 2.26 | 11.5% | 8.2% | 9.8% |
| 24 | 209,708 | 31.29 | 65.1 | 2.08 ± 0.01 | 2.92 | 2.07 | 12.3% | 8.2% | 9.6% |
| 28 | 507,770 | 33.29 | 131.5 | 3.95 ± 0.03 | 5.23 | 3.93 | 40.7% | 8.2% | 11.0% |
| 30 | 659,149 | 33.61 | 129.3 | 3.85 ± 0.02 | 5.04 | 3.83 | 40.6% | 8.2% | 10.9% |
| 32 | 967,025 | 32.92 | 102.2 | 3.10 ± 0.01 | 4.51 | 3.09 | 33.3% | 8.2% | 10.4% |
| 34 | 1,084,406 | 32.43 | 85.7 | 2.64 ± 0.01 | 3.92 | 2.63 | 30.7% | 8.2% | 10.1% |
| 38 | 806,140 | 30.06 | 83.6 | 2.78 ± 0.01 | 4.03 | 2.77 | 24.0% | 8.3% | 10.3% |
| 42 | 423,497 | 25.28 | 108.4 | 4.29 ± 0.02 | 6.31 | 4.27 | 20.3% | 8.3% | 11.9% |
| 46 | 164,553 | 21.43 | 99.5 | 4.64 ± 0.02 | 7.64 | 4.64 | 11.9% | 8.4% | 12.7% |
| 50 | 115,941 | 20.19 | 80.9 | 4.01 ± 0.02 | 5.71 | 3.99 | 14.3% | 8.4% | 12.3% |
| 56 | 9,547 | 19.20 | 47.6 | 2.48 ± 0.04 | 5.06 | 2.48 | 5.4% | 8.4% | 10.8% |

Twenty-one bins survived the 5,000-position floor and twelve are shown; the nine omitted all fall
inside the same range. **Every bin is over 2.0, and the whole range is 2.08 to 4.64.** The
standard errors are 0.01 to 0.06, so no bin's ratio is in question.

**Per accession**, pooling the bins by how many positions each holds:

| accession | reads a position | variance ÷ mean, pooled | across its GC bins | with the unreached positions back at zero |
|---|---:|---:|---:|---:|
| SRR7279533 | 2.80 | **0.84** | 0.54 – 1.05 | 1.13 |
| SRR7279488 | 3.04 | **0.92** | 0.73 – 1.29 | 1.24 |
| SRR7279501 | 4.02 | 1.23 | 0.23 – 2.92 | 1.49 |
| SRR7279484 | 5.52 | 1.66 | 0.40 – 4.34 | 1.99 |
| SRR7279481 | 10.29 | 2.49 | 0.96 – 6.19 | 2.93 |
| SRR7279483 | 14.23 | 5.05 | 1.10 – 41.17 | 6.05 |
| SRR7279482 | 26.43 | **2.48** | 1.36 – 4.24 | 3.61 |
| SRR7279540 | 30.23 | **3.18** | 2.08 – 4.64 | 4.56 |

**The answer differs at the two ends of the depth range, and the end that decides the term is the
deep one, because that is the only place the term carries information.** At 2.8
and 3.0 reads a position the ratio reads *below* 1 — but that is an artefact of what the walk
emits, not a property of the data: **the walk emits no locus where no read landed**, and a Poisson
with its zeros removed is narrower than a Poisson. Putting those positions back at depth zero
raises the two to 1.13 and 1.24. Neither reading is exactly right there, because the positions no
read reached are a mixture: 10.6% of generic bases at 2.8 reads a position, of which about 6
points are the zeros a Poisson of that mean would produce and the remaining 4.5 are the same
structurally uncovered floor the deep accessions show. **The honest statement at the shallow end
is that depth is Poisson to within the measurement — the truth lies between 0.84 and 1.13.**

At the deep end the same ambiguity does not arise: a Poisson of mean 30 produces a zero once in
ten million million positions, so **the 5.1% of positions no read reached at 30.2× are unmappable
and not unlucky**, and the right reading there is the one over positions a read reached. That is
the 2.48 and 3.18 above.

**So overdispersion appears with depth.** That is not a paradox — the extra spread is
multiplicative (a position that is 20% easier to map collects 20% more reads at any depth), so at
3 reads a position it is buried under the counting noise and at 30 it is not. It also means the
term cannot be defended as *Poisson enough at low coverage*: it is Poisson enough exactly where it
carries no information.

---

## 4. Three things that are not causing it

**Not the width of the GC bin.** A bin two percentage points wide still holds a range of GC
contents, so some of what reads as spread at fixed GC could be the bin's own width. It is not:
halving and doubling the bin moves the pooled ratio by less than 0.03.

| GC bin width | SRR7279540 | SRR7279482 |
|---|---:|---:|
| 1 percentage point | 3.17 | 2.47 |
| 2 percentage points | 3.18 | 2.48 |
| 4 percentage points | 3.20 | 2.50 |

**Not the duplications themselves.** The positions the term exists to find are deep, so they
inflate the very variance being measured. Dropping every near-half position changes the ratio by
at most 0.024 in any bin of the deepest accession (3.950 → 3.927 at 28% GC). That is not a full
answer, because a duplicated stretch only reads near half where its two copies differ and its
other positions cannot be removed this way. The size of what remains can be bounded: if 1 position
in 100 were truly doubled — which is above the share the window arm flags at these accessions,
0.66% and 0.72% — it would add about 10 to a variance of 102, and the ratio would read 2.8 instead
of 3.10. **The overdispersion is not the signal.**

**Not the choice of one-copy level.** Reading `m` as each GC bin's median depth rather than its
mean — the median being far less moved by a heavy upper tail — changes the arms by under 5%
(17.86 against 17.05 at 30.2 reads a position). Nothing below turns on it.

---

## 5. Is it a wider spread or a heavy tail? Mostly a wider spread, with one library's exception

A position in a repeat-rich neighbourhood reads deep for a reason neither the GC curve nor the
copy-number term models. If that arrived as a few very deep positions rather than a uniformly
wider spread, it would change *which family* is right and not merely its parameter — a heavy tail
is a mixture, and no single negative binomial holds one.

The measure is what share of a bin's total squared deviation sits in its deepest 1% of positions,
printed beside what a Poisson and a fitted negative binomial put there (§3's last three columns).
**At the well-behaved accessions the tail is heavier than either family allows but does not carry
the finding**: at 26.4 reads a position the deepest 1% hold 18.6% of the variance at 32% GC
against a Poisson's 8.3% and a fitted negative binomial's 10.0%. Removing that entire 1% would
still leave variance ÷ mean at about **1.9**, nowhere near 1.

**One accession is different and it should not be averaged in.** SRR7279483, at 14.2 reads a
position, has four GC bins between 46% and 52% running **22 to 41**, with the deepest 1% of
positions holding 25% to 56% of the variance. Its other seventeen bins run 1.1 to 6.6. That is a
property of that library at high-GC positions, not of tomato: no other accession shows it. **A
per-sample dispersion number fitted over all bins would be set by those four bins** — the fit
returns 5.05 for that accession against 2.48 and 3.18 for the two deeper ones — which is an
argument for fitting the dispersion per GC bin rather than per sample, or for a robust fit that
those bins cannot dominate.

---

## 6. Which family, and what it costs in the fit

Two families fit an overdispersed count and they disagree about what happens as depth grows: a
**negative binomial**, whose variance is `m + m²/r` so that variance ÷ mean rises with depth, and
a **quasi-Poisson**, whose variance is `φ·m` so that the ratio is one number at every depth. The
GC bins inside one sample span a 1.75-fold range of mean depth on the deepest accession and a
2.6-fold range on the next, which is a lever to tell them apart. Neither family wins by much:

| accession | negative binomial, one size `r` | misses each bin's ratio by | quasi-Poisson, one `φ` | misses each bin's ratio by |
|---|---:|---:|---:|---:|
| SRR7279482, 26.4 reads a position | 19.2 | 0.42 | 2.48 | 0.36 |
| SRR7279540, 30.2 reads a position | 14.7 | 0.59 | 3.18 | 0.53 |

Both miss by 15% to 20% of the ratio they are fitting, and the quasi-Poisson misses slightly less.
The reason neither fits is visible in §3: **inside one sample the dispersion is not a function of
the bin's depth at all.** On the deepest accession the two deepest bins — 33.3 and 33.6 reads a
position, at 28% and 30% GC — have ratios of 3.95 and 3.85, while a bin at 31.3 reads has 2.08 and
the shallowest bin, 19.2 reads at 56% GC, has 2.48. The negative binomial predicts the ratio rises
with depth and it does not; the quasi-Poisson predicts it is flat and it is not.

**Recommend the negative binomial anyway, and the reason is decisive rather than a preference.**
A quasi-Poisson is not a distribution. It specifies a mean and a variance and stops; there is no
`P(d | m)` to take the logarithm of. It is usable for a variance-adjusted regression and it cannot
be a term in a likelihood. The negative binomial is the distribution with the same first two
moments, and matching them at the sample's own depth is what this term needs, since it is only
ever evaluated near `m` and `2m`.

### What the two families hand the fit

This is where the choice bites, and it is not about which positions get flagged. Both families'
log-ratios rise with `d` at a fixed `m`, so at one GC content they order positions identically.
What differs is **how loudly each argues**, and by three orders of magnitude.

At `m` = 30 reads a position with the fitted size `r` = 14.7:

| | Poisson | negative binomial |
|---|---|---|
| the log-ratio's slope | 0.69 nats a read | 0.18 nats a read |
| it crosses zero at | 1.44 `m` | 1.40 `m` |
| a genuinely doubled position, `d` = 2`m`, gets | 11.6 nats — **108,000 to 1** | 3.2 nats — **25 to 1** |
| 10-to-1 odds are first reached at | 1.55 `m` | 1.83 `m` |
| 100-to-1 odds are first reached at | 1.66 `m` | 2.26 `m` |

**The negative binomial's answer at a truly doubled position is 25 to 1, not certainty, and that
is the correct answer.** When depth varies three times as much as a Poisson allows, one position's
read count simply cannot separate one copy from two with confidence — no choice of family creates
information the data does not hold. The Poisson's 108,000 to 1 is that missing information,
invented.

Measured over every scored position rather than at one hypothetical, **how often each family hands
out strong evidence of a duplication**, in positions per 10,000 scored:

| accession | reads a position | near-half positions per 10,000 | Poisson ≥ 100:1 | negative binomial ≥ 100:1 | Poisson ≥ 1,000:1 | negative binomial ≥ 1,000:1 |
|---|---:|---:|---:|---:|---:|---:|
| SRR7279533 | 2.80 | 3.6 | 18.6 | 18.6 | 0.7 | 0.7 |
| SRR7279488 | 3.04 | 6.3 | 29.7 | 29.7 | 4.4 | 4.4 |
| SRR7279501 | 4.02 | 4.8 | 48.4 | 11.8 | 13.5 | 3.2 |
| SRR7279484 | 5.52 | 4.8 | 108.5 | 8.2 | 29.0 | 2.1 |
| SRR7279481 | 10.29 | 3.9 | 225.0 | 7.8 | 85.2 | 1.5 |
| SRR7279483 | 14.23 | 5.0 | 199.4 | 16.0 | 119.7 | 4.8 |
| SRR7279482 | 26.43 | 5.3 | 126.7 | 37.3 | 81.8 | 12.6 |
| SRR7279540 | 30.23 | 5.0 | 175.8 | 65.6 | 135.7 | 39.4 |

Read the 10.29 row: **a Poisson term would tell the fit, at 85 positions in every 10,000, that a
duplication is at least a thousand times likelier than not** — 1 position in 117, where 3.9 in
10,000 read near half at all. The negative binomial says that of 1.5. At the two shallowest
accessions the two columns are identical because the fitted size hit its ceiling: on the positions
the walk emitted there, depth is not overdispersed at all (§3), so the negative binomial *is* the
Poisson.

**On the doubled position's own spread**, a second modelling choice with a smaller consequence. If
the overdispersion is a local property the two copies share — mappability, a GC residual — the
doubled position keeps the same size `r` and its spread doubles with its mean. If the two copies'
read counts are drawn independently, it has size 2`r` and is relatively narrower. Both were
scored. At the two deep accessions they are almost the same object: they hand out 1,000-to-1
evidence to 12.6 against 12.6 positions per 10,000 at 26.4 reads a position, and 39.4 against 38.3
at 30.2. At 14.2 they are not — 4.8 against 0.0 — because that accession's fitted dispersion is set
by four aberrant bins (§5) and doubling the size undoes them. **Recommend the shared reading**,
because the causes of the overdispersion that can be named — mappability and residual GC bias —
act on the place and not on the copy, and because it is the one that stays conservative when the
dispersion estimate is poor.

---

## 7. The log-ratio against the threshold

Every arm scores the same positions: those at four reads or more whose GC bin holds enough
positions to have an expected depth. The first column is the earlier report's arm, rebuilt here.

**The reproduction.** Enrichment of the threshold arm, this program against
`locus_depth_vs_window_2026-08-13.md` §4's exact-read-count column:

| accession | 2.80 | 3.04 | 4.02 | 5.52 | 10.29 | 14.23 | 26.43 | 30.23 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| this program | 1.12 | 1.00 | 0.67 | 1.04 | 1.75 | 3.71 | 12.44 | 17.05 |
| the earlier report | 1.14 | 1.01 | 0.66 | 1.04 | 1.75 | 3.61 | 12.5 | 17.1 |

**The arms.** Enrichment, and in brackets the share of scored positions the arm flags:

| accession | reads a position | threshold, 1.6 ≤ `d`/`m` < 2.4 | threshold, no upper edge | log-ratio above zero | log-ratio, same number of positions as the threshold |
|---|---:|---|---|---|---|
| SRR7279533 | 2.80 | 1.12 (44.4%) | 1.53 (49.7%) | 1.38 (58.5%) | **1.61** (44.4%) |
| SRR7279488 | 3.04 | 1.00 (46.6%) | 1.50 (53.2%) | 1.44 (56.0%) | **1.63** (46.6%) |
| SRR7279501 | 4.02 | 0.67 (32.9%) | 1.08 (39.4%) | 1.07 (38.7%) | **1.19** (32.9%) |
| SRR7279484 | 5.52 | 1.04 (19.8%) | 1.55 (23.0%) | 1.38 (28.6%) | **1.70** (19.8%) |
| SRR7279481 | 10.29 | 1.75 (11.6%) | 2.38 (13.0%) | 1.75 (20.7%) | **2.62** (11.6%) |
| SRR7279483 | 14.23 | 3.71 (7.2%) | 5.04 (8.2%) | 3.27 (14.1%) | **5.59** (7.2%) |
| SRR7279482 | 26.43 | 12.44 (2.0%) | 13.76 (2.1%) | 4.86 (7.1%) | **14.82** (2.0%) |
| SRR7279540 | 30.23 | 17.05 (1.8%) | 20.33 (2.3%) | 9.00 (6.1%) | **24.54** (1.8%) |

The log-ratio columns are the negative binomial's. In the last column — the two arms given the
same number of positions — the Poisson's are within 3% of them at every accession, because at one
GC content the two families order positions identically (§6). In the *above zero* column they
part company, because they cross zero at different depths: the Poisson reads 5.0% and 10.71 at
30.2 reads a position where the negative binomial reads 6.1% and 9.00.

**Three findings, in the order that matters.**

1. **The log-ratio ranks better than the threshold, at every accession.** Given the same number of
   positions to flag it reaches 24.5 against 17.1 at 30.2 reads a position, 14.8 against 12.4 at
   26.4, and 2.6 against 1.8 at 10.3. Two things account for the gap. It has **no upper edge** —
   the threshold's band stops at 2.4 copies, and simply removing that edge already takes 17.05 to
   20.33. And it **ranks across GC bins by how much evidence a position carries**, not by a
   ratio: at a GC content where the sample is deep, twice the expected depth is stronger evidence
   than the same doubling where it is shallow, and only the log-ratio knows that. The first of the
   two is measured — removing the edge is a separate row in the table. The second is the mechanism
   that remains once the edge is accounted for, and it is an explanation rather than a second
   measurement.
2. **Its zero is not a place to cut.** Scored as *positive means duplicated* it flags 6.1 to 7.1
   positions in 100 at the deep accessions against the threshold's 1.8 to 2.0, and enrichment
   falls to 9.0 and 4.9 — below the threshold arm. The zero crossing sits at 1.40 `m` (§6), well
   inside the ordinary spread of a distribution whose variance is three times its mean. **What the
   term gives the fit is a weight, and the fit must use it as one.**
3. **At three reads a position it is the same flood the threshold was.** 58 positions in 100 above
   zero, at an enrichment of 1.4. This changes nothing about the earlier report's conclusion: the
   position's own read count is not a substitute for the coverage-by-window summary at the depths
   the tomato cohort runs at, and no choice of family or cut makes it one.

---

## 8. What this cannot say

- **No truth set.** Every number here is a comparison between arms on the same positions, never an
  accuracy. An arm at 24.5 is not finding 24.5 times more duplications.
- **The overdispersion and the thing being detected cannot be fully separated.** Some deep
  positions are duplications, and only those that also read near half can be removed. §4 bounds
  what remains at about 0.3 of the ratio; it does not eliminate it.
- **8 Mb of tomato chromosomes 1 to 12**, randomly placed, one species, one aligner, eight
  accessions. Unplaced contigs are excluded, and those are where an assembly's collapsed copies
  concentrate. **Nothing here has been measured on human data**, where mappability and GC bias
  behave differently and where a 30× whole-genome run is the ordinary case rather than the deepest
  sample available.
- **Single-base generic loci only.** Positions inside a locus wider than one base carry a depth but
  no alternative fraction, so they are outside both measurements. They are about 1 position in
  1,100 of the walk.
- **The dispersion is measured against a GC-only model of expected depth.** Any other systematic
  cause of depth variation — mappability above all — lands inside the measured spread. That is the
  right thing for this purpose, because a GC-only model is what the term will use, but it means the
  number is a property of the model and not only of the sequencing.
- **One dispersion number a sample fits the GC range to about 20%** (§6), and at one accession of
  eight it is set by four aberrant bins (§5). The measurement says the family; it does not say the
  fitting procedure is settled.

---

## 9. What I would change in the specifications — not made here

These are proposals; nothing under `doc/devel/ng/spec/` or `doc/devel/ng/arch/` was edited.

| document | proposed change |
|---|---|
| `spec/parameter_prepass_joint_fit.md` §2.2 | State the copy-number term's `P` as a **negative binomial** with a per-sample dispersion, and record why in one line: depth's variance is 2.5 to 3.2 times its mean at 26 to 30 reads a position (§3), and a Poisson form gives 1,000-to-1 odds of a duplication to 1 position in 74 where 5 in 10,000 read near half (§6). Say that the doubled position keeps the same dispersion, not twice it (§6). |
| `spec/parameter_prepass_joint_fit.md` §2.2 | Say explicitly that the log-ratio is a **weight and not a test**: its zero crossing at 1.40 `m` flags 6 to 7 positions in 100 (§7), so no branch of the fit may read *log-ratio > 0* as *this position is duplicated*. |
| `spec/parameter_prepass_joint_records.md` §4 | The per-sample coverage summary must carry **one dispersion number a sample** beside the depth-against-GC curve, or the term has no `r` to use. It is one `f32`. Note that it must be estimated robustly: at one accession of eight, four GC bins out of twenty-one would otherwise set it (§5). |
| `spec/parameter_prepass_joint_fit.md` §2.2 | Record the depth range the term is worth evaluating over. It carries nothing below about ten reads a position — 58 positions in 100 above zero at 2.8 (§7) — and the specification should say what the fit does there instead of leaving a term that returns confident noise. |
| `spec/parameter_prepass_generic.md` §4 | The expected depth `m` a position is scored against must come off a **per-position** GC curve, not a per-window one: single-position depth swings 1.4-fold across GC on one accession (23.7 to 33.6 reads), and a term that divides by the wrong `m` moves every position's log-ratio by the same error. |
