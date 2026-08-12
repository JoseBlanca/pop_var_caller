# The third class of site: it is there, coverage finds it, and it is thirty times smaller than the design assumed

*Research report, 2026-08-12. The measurement
[`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §2.2 and
[`parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md) §4 both stop and
wait for. **One program stands behind it**: `examples/ng_joint_duplication_probe.rs`, one walk over
one alignment, no cohort and no truth set. Raw output in `tmp/dup/`.*

---

## 1. The question, and the answer

The fit wants a third class of site: **a locus the sample carries more copies of than the reference
does**, collecting two copies' reads at one position, so every position where the copies differ
shows about half the reads disagreeing with the reference — indistinguishable from a heterozygote by
its read counts alone. Its discriminator was specified as **local relative coverage**, and the object
that would cost something is a per-sample coverage-by-window summary. Both documents gate that object
on one number: *on one tomato sample, what fraction of positions sit in windows near two copies and
show a near-half alternative fraction.*

**Three answers, and the third is the one that changes a decision.**

1. **The population is there and coverage separates it.** On tomato SRR7279482 at 25× depth,
   **1 position in 8,600** is both in a window carrying about twice the sample's normal coverage and
   reading between 35% and 65% alternative. Inside those windows the near-half rate is **1.26%**
   against **0.033%** in ordinary-coverage windows — a factor of 38, at matched read depth.
2. **The discriminator has nothing to do with depth and everything to do with reads per window.**
   At the same 500 bp window a sample at 3.6× shows a near-half enrichment of **1.3×** — no
   separation at all — where 25× shows **24×**. Widening the window to 5–10 kb restores it: the
   2.5× sample goes from 1.6× to **15×**. What the window needs is about **12,000 aligned bases**,
   which is depth times width.
3. **The class is about thirty times smaller than the specs say.** They put it at 1,700 to 8,400
   positions of two million. Measured across eight samples it is **150 to 590** — smaller than the
   count of near-half positions in ordinary-coverage windows, where the specs had it six to thirty
   times larger.

---

## 2. What was measured

`examples/ng_joint_duplication_probe.rs`. One pass over the reference gives every fixed window its
denominator — how many generic (non-repeat) positions it holds — and its GC fraction over exactly
those positions. One walk through the real locus generator then gives every generic position its
observation depth and, at single-base loci, its alternative-read count.

**Coverage tracks GC and the correction is fitted from the sample itself**: the median window depth
in each 2-point GC bin, falling back to the global median where a bin holds fewer than 50 windows. A
window's **relative coverage** is its mean depth over that curve at its own GC, rescaled so the
median window sits at 1.0. So 2.0 is a window carrying twice this sample's normal depth.

The correction is not cosmetic on tomato: median window depth runs from **16.2 reads a position at
20% GC to 29.0 at 36%**, a factor of 1.79 across the range — larger than the one-copy-to-two-copy
signal it would otherwise swamp.

**The data.** Eight of the 63 tomato bench alignments, spanning **2.5× to 28.7×** mean depth, over
`benchmarks/tomato1/regions.bed` — 80 spans of 100 kb picked at random from chromosomes 1 to 12
(`scripts/pick_regions.py`), of which 7,935,192 bases (99.2%) are typed generic by the repeat
catalog. A whole run takes 15 seconds.

**Depth here is observation depth**, the same quantity the generic path's accumulator sees: a read
that covered a position but anchored no border is not counted.

---

## 3. The population exists, and its shape is a duplication's

SRR7279482, 25.2× mean depth, 500 bp windows, positions at depth 8 or more.

| | windows | positions | near-half rate |
|---|---:|---:|---:|
| relative coverage below 0.6 | 1,835 | 901,227 | 0.106% |
| **0.6 to 1.4 — one copy** | 13,775 | 6,842,112 | **0.033%** |
| 1.4 to 1.6 | 213 | 105,866 | 0.161% |
| **1.6 to 2.4 — two copies** | **137** | **68,178** | **1.258%** |
| 2.4 to 3.5 | 15 | 7,474 | 1.289% |

**The gate number: 853 positions, 0.0116% of those scored, 1 in 8,594.** If window coverage and
alternative fraction were independent it would be 34 positions, so the joint cell holds **24.8 times**
what independence predicts.

**The shape is the one a duplication makes, not the one more depth makes.** The distribution of the
alternative fraction inside each band, as a percentage of that band's positions:

| band | 0.3–0.4 | 0.4–0.5 | 0.5–0.6 | 0.6–0.7 |
|---|---:|---:|---:|---:|
| one copy | 0.017% | 0.013% | 0.011% | 0.005% |
| **two copies** | **0.198%** | **0.555%** | **0.524%** | **0.083%** |

The two-copy band has a bump centred on a half that the one-copy band does not have: 44 times the
mass in 0.4–0.5, 50 times in 0.5–0.6.

**Two confounds, both ruled out by the same table.** The near-half rate inside two-copy windows is
1.2576% at depth 2 or more and 1.2579% at depth 16 or more — flat, so it is not an artefact of
alternative fractions being easier to measure where there are more reads. And **1.26% is what a
recent duplicate looks like**: it is the share of positions at which the two copies differ, which is
sequence divergence, not a rate anything in the caller sets.

**What the one-copy band's 0.033% is, is roughly heterozygosity.** 0.033% of 6.8 M positions is
2,279, which over two million positions would be 668 — the right order for a tomato accession, and
the reason it is the right comparison for §5.

**Across all eight samples**, at 5 kb windows and depth 4 or more:

| sample | mean depth | positions in two-copy windows | near-half rate, one copy | near-half rate, two copies | ratio |
|---|---:|---:|---:|---:|---:|
| SRR7279533 | 2.51× | 1.07% | 0.017% | 0.516% | 31× |
| SRR7279488 | 2.72× | 1.30% | 0.023% | 0.915% | 39× |
| SRR7279501 | 3.60× | 3.24% | 0.025% | 0.317% | 13× |
| SRR7279484 | 5.15× | 2.15% | 0.032% | 0.376% | 12× |
| SRR7279481 | 9.89× | 2.42% | 0.026% | 0.376% | 15× |
| SRR7279483 | 13.32× | 1.27% | 0.027% | 0.803% | 29× |
| SRR7279482 | 25.20× | 0.62% | 0.040% | 1.121% | 28× |
| SRR7279540 | 28.69× | 0.69% | 0.028% | 1.158% | 41× |

**Every sample separates, and the two-copy share sits near 1% of positions** — the same order as the
0.42% and 0.49% the histogram route's fit asks for on two of these samples
([`parameter_prepass_generic.md`](../spec/parameter_prepass_generic.md) §2.1), which is the first
independent confirmation that the two are the same phenomenon.

---

## 4. The window needs twelve thousand aligned bases, and that is a change to the spec

[`parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md) §4 fixes the
window at 500 bp. **On the cohort this caller is aimed at, 500 bp does not work.**

Enrichment of the joint cell over independence, 500 bp windows, positions at depth 4 or more:

| mean depth | 2.51× | 2.72× | 3.60× | 5.15× | 9.89× | 13.32× | 25.20× | 28.69× |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| enrichment | 1.6× | 1.6× | **1.3×** | 1.5× | 2.5× | 7.7× | 24.0× | 24.9× |

**At 3.6× a 500 bp window's mean depth is scatter.** Its two-copy band holds 9.7% of positions —
eleven times the deep samples' share — and those positions are no more likely to read near half than
any others.

**Widening the window recovers it**, and what matters is depth times width:

| window | SRR7279533, 2.51× | SRR7279501, 3.60× | SRR7279482, 25.20× |
|---:|---:|---:|---:|
| 500 bp | 1.6× | 1.3× | 24.0× |
| 1 kb | 3.6× | 2.3× | 31.1× |
| 2 kb | 8.5× | 5.2× | 24.2× |
| 5 kb | 14.0× | 6.6× | 21.3× |
| 10 kb | 15.1× | 5.9× | 25.6× |

2.51× at 5 kb is 12,550 aligned bases a window; 25.2× at 500 bp is 12,600, and the two return the
same enrichment. **Below about 12,000 the classification degrades, and the deep sample gains nothing
above it.**

**What to change, and it is not the stored width.** Store the summary at the fine grid — 500 bp,
1.6 MB a sample on tomato, as §4 already prices it — and let the **fit** sum adjacent windows up to
whatever width the sample's own depth requires. Storing at each sample's own width would break the
one property §4 leans on, that two samples' summaries are comparable by construction, and would put
a per-sample number into a value the identity check demands every sample agree on
([`parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md) §5). Summing is
exact and free; unsumming is not possible.

---

## 5. Whose duplication is it — the reference's or the individual's?

**The fit spec says the class's grain is the (locus, sample) pair**, on the argument that a
duplication carried by an individual is that individual's property, unlike a collapsed paralog which
mismaps in every sample. That argument decides whether the summary must be per sample at all, and it
had never been checked.

Eight samples on one 5 kb grid of 1,680 windows. **84 windows are called near two copies by at least
one sample**; of those, **40 by exactly one sample** and **11 by seven or eight**. So both components
are present, and the individual one is the larger.

The threshold-free version says the same thing. Taking one sample's two-copy windows and asking what
the other seven read there — 2.0 would mean the amplification belongs to the reference, 1.0 that it
belongs to the individual:

| windows called by | n | 481 | 482 | 483 | 484 | 488 | 501 | 533 | 540 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SRR7279482 | 10 | 2.01 | 2.05 | 1.99 | 2.10 | 1.93 | 2.01 | 1.94 | 1.36 |
| **SRR7279501** | 52 | **2.20** | **0.82** | 1.39 | **2.27** | 1.03 | 1.95 | 1.11 | **0.76** |

**SRR7279482's ten windows are everybody's** — every other sample reads between 1.36 and 2.10 there,
so those are the reference's own collapse. **SRR7279501's fifty-two are not**: SRR7279481 and
SRR7279484 read 2.20 and 2.27, while SRR7279482 and SRR7279540 read 0.82 and 0.76, which is one copy.
That is copy-number variation segregating in the panel, and a per-locus class would force one
accession's amplification onto samples that do not carry it.

Per-window near-half counts correlate between samples at **0.25 to 0.89**, with the highest values
inside the {481, 484, 501} group and the lowest between it and {482, 540} — the same split.

**So the fit spec's grain is right and the per-sample object is necessary.** What the cohort adds is
the ability to tell the two components apart, which is the bonus §2.2 already claims.

---

## 6. The size, and this is the finding that changes an argument

[`parameter_prepass_joint_loci.md`](../spec/parameter_prepass_joint_loci.md) §4.2 and
[`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §2.2 put this population
at **1,700 to 8,400 of two million positions**, and build on it the conclusion that *at most 3 in
every 100 positions carrying an alternative read are really heterozygous*. The measured figure is
**150 to 590 of two million**, across all eight samples.

**Where the old number came from, and why it was wrong by a factor of thirty.** It reads the
histogram fit's noisy-class weight — 0.42% of *sites* — as the count of duplicated positions showing
an alternative read. This measurement separates the two quantities that estimate conflates:

- **0.6% to 3.2% of positions sit in a window near two copies.** That matches the fitted weight, and
  it is what the fit's mixture is picking up.
- **Only 0.3% to 1.2% of those positions actually read near half**, because a duplication is silent
  wherever its two copies agree, which is 99% of their length.

The product is the artefact, and the product is small. **Against the same sample's near-half
positions in ordinary-coverage windows — 668 per two million on SRR7279482 — the artefact is about a
third, where the specs had it six to thirty times larger.**

**Two things follow, and they point in opposite directions.**

- **The class is still worth carrying.** It is 24 times concentrated where coverage says it should
  be, it is the correct sign, and the fit's only alternative for those positions is to call every
  sample heterozygous at a mid-frequency variant.
- **It is no longer the largest term in the heterozygosity budget.** The 6,000 error positions and
  2,500 noisy-class positions §4.2 counts are each an order of magnitude larger than it, so the
  paragraph concluding that duplication outweighs both is wrong and is corrected in place.

---

## 7. What this cannot say

- **There is no truth set**, so nothing here shows a two-copy window *is* a duplication. What it
  shows is that a population with a duplication's signature exists, at a size, and that relative
  coverage finds it.
- **8 Mb of chromosomes 1 to 12**, randomly placed. Unplaced contigs are excluded, and those are
  where an assembly's collapsed copies concentrate — so this is a floor rather than a genome-wide
  rate.
- **GC correction is measured to matter for the coverage curve and not measured to help the
  classification.** On SRR7279482 the enrichment is 24.8× corrected and **32.6× uncorrected**, with
  the correction adding 14,000 two-copy positions that carry no near-half signal. The curve's own
  1.79-fold range says the correction is doing something real; whether the fit wants it is open, and
  it is one flag in the accumulator either way.
- **The alternative-read fraction is only read at single-base loci.** Wider loci — an indel's
  reference span — contribute their depth to the window and nothing to the numerator; they are
  5,311 loci in 7.56 M.

---

## 8. What changed in the documents

| document | what changed |
|---|---|
| `spec/parameter_prepass_joint_fit.md` §2.2 | the open measurement is closed with the numbers; the 1,700–8,400 estimate replaced by the measured 150–590; the (locus, sample) grain now cites §5's evidence rather than only the argument |
| `spec/parameter_prepass_joint_records.md` §4 | the window stays 500 bp **stored** and gains the fit-side summing rule and the 12,000-aligned-base floor; the gating question is answered |
| `spec/parameter_prepass_joint_loci.md` §4.2 | the artefact term in the alternative-read budget recomputed; the "3 in 100" conclusion restated with the measured sizes |
