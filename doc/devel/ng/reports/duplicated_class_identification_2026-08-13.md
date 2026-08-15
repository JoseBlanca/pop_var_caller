# Ignoring a duplicated position costs a quarter of the inbreeding coefficient, and fifty samples find it without a coverage summary

*Research report, 2026-08-13. The measurement
[`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §11 question 11 asks for.
**One program stands behind it**: `examples/ng_duplicated_class_harness.rs`, which draws a cohort with
a known truth and fits the same drawn cohort three ways. Raw output in `tmp/dupclass/`.*

---

## 1. The two questions, and the answers

**Ignoring the class costs a quarter of the inbreeding coefficient, and the per-sample coverage
summary is not what recovers it — the cohort does, from about twenty-five samples up.**

1. **The two inflations do not cancel.** Both observed and expected heterozygosity rise, but not by
   the same amount, and not for the same reason. On a fifty-sample selfing panel at three reads a
   position, a fit with no third class puts `Hobs` **50.6%** above the drawn truth and `Hexp`
   **10.6%** above it, so `1 − Hobs/Hexp` reads **0.4471 where the truth is 0.5942**. At twenty-five
   reads it reads 0.4993 against 0.5855. Nothing else in the fit shows it: the two error rates and the
   noisy share come back within a percentage point either way.
2. **The cohort pattern alone recovers nearly all of it at fifty samples, and not at ten.** With the
   class identified only from the pattern across the cohort — no coverage reading anywhere — the
   coefficient reads **0.5807** against 0.5942 and `Hobs` is 3.0% high instead of 50.6%. At
   twenty-five samples `Hobs` is 6.8% high; **at ten it is 21.2% high** and the coefficient reads
   0.5498. The coverage summary holds `Hobs` to within 1.3% at every panel size.
3. **What the pattern misses is what matters least, and that is not a coincidence.** It finds nearly
   every duplication five or more of fifty samples carry — 59 of 64 — and 2 of the 37 carried by
   fewer. A duplication few samples carry inflates expected heterozygosity nearly as hard as it
   inflates the observed one, so it very largely cancels out of the ratio; a duplication the whole
   panel carries contributes
   **nothing** to expected heterozygosity and a heterozygote in every sample to the observed one. The
   pattern's blind spot and the harm are in different places.

**Recommendation: drop the requirement, keep the capability.** Make the third class identifiable from
the cohort pattern always — it costs nothing but the class, and it is where most of the benefit is —
and make the coverage-by-window summary **an input the fit uses when it has it** rather than an object
the design demands. The genome walk builds it for free in the run that goes straight from alignments
to a fit, so nothing is lost there. **In the two-phase run, do not spend a second full pass over every
sample's pileup rebuilding it by default**; spend it when the cohort is under about twenty-five
samples, when the panel is outbred, or when the run asks for it. The trade-off is explicit: on a
fifty-sample selfing cohort that pass buys about two percentage points on `Hobs` and one on the
inbreeding coefficient, against a fit with no third class at all, which is not survivable at any panel
size.

**Three things would change that recommendation, and two are measurable now.** The pattern's whole
evidence is that no sample is homozygous for the non-reference allele, and inbreeding is what makes
those common — **on an outbred panel the pattern finds only the duplications every sample carries**
(§4). It also becomes greedy where the class is large, calling 80 of 847 real variants duplicated on a
panel of overwhelmingly private duplications where the coverage fit calls 1 (§5). And **at one sample
the pattern has no power at all**, so a single-sample run keeps the third class only if it keeps the
coverage summary.

---

## 2. What was drawn, and what each fit was allowed to read

**The panel.** Fifty samples, 100,000 positions, drawn twice — once at **three reads a position**,
which is where the tomato archive sits, and once at **twenty-five**. Every position is one of three
kinds:

- **ordinary** — a population allele frequency from a monomorphic mass (0.990), a fixed-alternative
  mass (0.002) and a `Beta(0.3, 1.2)` over what segregates, then a genotype per sample under
  Hardy–Weinberg with inbreeding. This is the truth `examples/ng_joint_fit_harness.rs` already uses,
  unchanged;
- **duplicated** — a carrier frequency for the position, then a carrier indicator per sample. **A
  carrier gets twice the depth and about half its reads disagreeing with the reference**; a
  non-carrier is homozygous reference at ordinary depth. A carrier is *not* heterozygous, and the
  drawn truth does not count it as one;
- reads at every position carry sequencing error at one of two rates — 0.0019 at 99 positions in 100
  and 0.053 at the hundredth — the two-class noise model the fit already carries.

**How many duplicated positions.** The tomato measurement
([`duplicated_locus_probe_2026-08-12.md`](duplicated_locus_probe_2026-08-12.md) §6) counts, per
sample, **150 to 590 positions per two million** that are both in a window near two copies and reading
near half, against about **668 genuinely near-half positions per two million** in the same sample — so
the artefact is about **a third** of the real signal. **That ratio is what the drawn panel holds
fixed**, because the inflation it causes is a ratio: what a sample's observed heterozygosity gains is
the duplicated carrier rate divided by the heterozygous rate, and neither count alone decides it.

**Who carries a duplication.** Of 84 windows read near two copies by at least one of eight tomato
samples, 40 are read that way by exactly one sample and 11 by seven or eight (probe report §5). Those
counts are a *sample* of the carrier frequency rather than the frequency itself, so they are turned
into one: **9.0% of duplicated positions are carried by every sample** — the reference's own collapse,
which is what the seven-or-eight group is — and the rest have a carrier frequency drawn from
`Beta(1.19, 9.55)`, which reproduces the 40-against-33 split among eight samples to within a
thousandth. Fifty samples then see a duplication that eight samples would call private in one to five
of them, not in one.

**The coverage reading.** Every sample at every position also gets a **relative coverage** value — its
window's mean depth over its own depth-against-GC curve, rescaled so the median window sits at 1.0. It
is drawn around the sample's true copy number with the scatter a window of that many aligned bases
really has, calibrated to two points the probe report measured: at 12,600 aligned bases 86% of
single-copy windows land inside 0.6 to 1.4 and 0.86% reach 1.6 or above, and at 1,800 aligned bases
9.7% of positions land in the two-copy band. Those fix a standard deviation of
`√(copies² × 0.194² + copies × 313/aligned bases)`. **The default is 12,000 aligned bases a window**,
which is the floor the probe report found the classification needs and which the spec's summing rule
reaches at any depth.

**The three fits.** One drawn panel, three fits, everything else identical — same starting points,
same search, same eleven fitted numbers minus the three a fit without the third class does not have:

| fit | the third class | what its carrier state reads |
|---|---|---|
| `coverage` | yes | the sample's reads **and** its relative coverage |
| `pattern` | yes | the sample's reads only |
| `no-class` | no | — |

**The third class is an ordinary variant with one genotype removed.** At a position drawn from it,
each sample is either a carrier — half its reads disagreeing — or homozygous reference, with a carrier
frequency integrated over a fitted `Beta`, exactly as the ordinary class integrates over allele
frequency. What it has no room for is **a sample homozygous for the non-reference allele**. A real
variant at a frequency of a half leaves about a quarter of the panel there; a duplication leaves none,
and that difference is the entire evidence the `pattern` fit has.

**The coverage readings enter in one place only.** Under every class except the third, every sample is
single-copy, so its coverage reading contributes the same factor at every class and cancels out of the
position's likelihood entirely. What survives is the ratio *how much more likely is this reading if
the sample has two copies here rather than one*, multiplying the carrier branch of the third class. So
`coverage` and `pattern` differ by that one ratio and by nothing else.

**Neither fit reads the position's own depth.** A carrier really does get twice the reads, and both
fits condition on the depth they see rather than modelling it. That is what the spec requires at three
reads a position — per-base depth cannot tell a two-copy carrier at twelve reads from a single-copy
sample reading high (§2.2) — and at twenty-five reads it is conservative: it denies the `pattern` fit
a signal that in real data is swamped by GC content and mappability, neither of which this panel
draws.

**Nine panels were drawn across six settings, and each panel was fitted all three ways.** The first is
the headline and the rest say what the answer depends on:

| panel | what it changes | where it is read |
|---|---|---|
| **fifty samples, inbreeding 0.6** | the headline, at three reads and at twenty-five | §3, §4, §5 |
| fifty samples, **inbreeding 0** | an outbred panel instead of a selfing one | §4, §5 |
| fifty samples, **private duplications** | half of them carried by one sample of fifty | §4, §5 |
| **ten and twenty-five samples** | how many samples the cohort pattern needs | §5 |
| fifty samples, **1,500-aligned-base window** | the stored 500 bp grid at three reads, unsummed | §5 |

```text
ng_duplicated_class_harness <samples> <depth> <positions> <inbreeding> \
                            <aligned bases a window> <collapse share> <carrier a> <carrier b>
```

A depth of `0` runs three reads a position and twenty-five. Every run is deterministic from a fixed
seed: re-running the headline after the program grew an extra counter reproduced every reported number
exactly.

---

## 3. What ignoring the class costs

Fifty samples, 100,000 positions, an inbreeding coefficient of 0.6 — a selfing panel, which tomato is.
Every number is against the drawn truth of that same panel.

**At three reads a position.** The drawn panel holds 502 duplicated carrier positions per two million
against 1,304 heterozygous ones, and 101 duplicated positions carried by at least one sample: 13 by
one, 24 by two to four, 48 by five to twenty-four, 16 by twenty-five or more.

| fit | `Hobs` | against truth | `Hexp` | against truth | inbreeding coefficient | against truth |
|---|---:|---:|---:|---:|---:|---:|
| drawn truth | 6.518 × 10⁻⁴ | | 1.606 × 10⁻³ | | 0.5942 | |
| `coverage` | 6.581 × 10⁻⁴ | +1.0% | 1.588 × 10⁻³ | −1.2% | 0.5855 | −1.5% |
| `pattern` | 6.713 × 10⁻⁴ | +3.0% | 1.601 × 10⁻³ | −0.3% | 0.5807 | −2.3% |
| **`no-class`** | 9.819 × 10⁻⁴ | **+50.6%** | 1.776 × 10⁻³ | **+10.6%** | **0.4471** | **−24.8%** |

**At twenty-five reads a position.** A fresh draw: 396 duplicated carrier positions per two million
against 1,303 heterozygous, 90 duplicated positions carried by someone.

| fit | `Hobs` | against truth | `Hexp` | against truth | inbreeding coefficient | against truth |
|---|---:|---:|---:|---:|---:|---:|
| drawn truth | 6.516 × 10⁻⁴ | | 1.572 × 10⁻³ | | 0.5855 | |
| `coverage` | 6.520 × 10⁻⁴ | +0.1% | 1.534 × 10⁻³ | −2.4% | 0.5751 | −1.8% |
| `pattern` | 6.639 × 10⁻⁴ | +1.9% | 1.568 × 10⁻³ | −0.3% | 0.5766 | −1.5% |
| **`no-class`** | 8.497 × 10⁻⁴ | **+30.4%** | 1.697 × 10⁻³ | **+7.9%** | **0.4993** | **−14.7%** |

**Nothing else in the fit feels it.** The clean error rate, the noisy rate and the noisy share come out
within a percentage point of the truth in all three fits at both depths. The damage is confined to the
three numbers above, which is what makes it dangerous: a run with no third class returns error rates
that look right and an inbreeding coefficient that is a quarter too low.

---

## 4. The two inflations do not cancel, and inbreeding is why

**Both rise, and the observed one rises about five times as far.** At three reads, ignoring the class
puts `Hobs` **50.6%** above the truth and `Hexp` **10.6%** above it. The inbreeding coefficient is
`1 − Hobs/Hexp`, so what it sees is the gap: **0.4471 against a true 0.5942**, a quarter of its value
gone. At twenty-five reads the same two numbers are +30.4% and +7.9%, and the coefficient reads 0.4993
against 0.5855.

**The two inflations are nearly the same size in absolute terms, and that is exactly why the ratio
moves.** A duplicated position carried by a share `q` of the panel adds `q` to the heterozygote count
in `Hobs`, and `2q(1−q)` to the position's contribution to `Hexp`. Over the drawn carrier spread those
average 0.19 and 0.164 — within a fifth of each other. But **`Hobs` and `Hexp` are not the same size
to begin with**: at an inbreeding coefficient of 0.6 a sample is heterozygous at only 40% of the rate
random mating would give, so `Hobs` is 40% of `Hexp` and the same absolute addition is two and a half
times larger as a fraction of it. **The gap between the two inflations is therefore the panel's
inbreeding**, and it is widest exactly where this caller is aimed: an autogamous crop panel.

**And the damage is carried almost entirely by the duplications many samples share.** Compare `q` with
`2q(1−q)`: when few samples carry the duplication the second is nearly twice the first, so expected
heterozygosity gains about twice as much per position as observed heterozygosity does — which very
nearly offsets `Hobs` being the smaller of the two to begin with. When the whole panel carries it,
`2q(1−q)` is **zero**: the position makes every sample heterozygous and contributes nothing at all to
expected heterozygosity. **A duplication the reference itself collapsed is therefore the worst case
there is**, and 9 of every 100 duplicated positions on the tomato measurement are that.

The drawn panels bear it out. Redrawn so that duplications are overwhelmingly private — half of them
carried by exactly one sample of fifty, none by more than twenty-four — ignoring the class puts `Hobs`
43.0% high and `Hexp` 34.8% high, and the inbreeding coefficient reads **0.5725 against a true
0.5972**, down 4% instead of down a quarter. **So the class is worth carrying because of the
duplications the panel shares, not because of the ones one accession has.**

**On a panel with no inbreeding it does not cancel either — it goes negative, which is at least
loud.** The same panel drawn with every sample outbred, and the same one duplicated carrier position
for every three heterozygous ones:

| fit | `Hobs` | against truth | `Hexp` | against truth | inbreeding coefficient |
|---|---:|---:|---:|---:|---:|
| drawn truth | 1.569 × 10⁻³ | | 1.566 × 10⁻³ | | −0.002 |
| `coverage` | 1.560 × 10⁻³ | −0.6% | 1.536 × 10⁻³ | −1.9% | −0.016 |
| `pattern` | 1.802 × 10⁻³ | +14.8% | 1.811 × 10⁻³ | +15.7% | +0.005 |
| **`no-class`** | 2.179 × 10⁻³ | **+38.9%** | 1.995 × 10⁻³ | **+27.4%** | **−0.092** |

*(Holding the measured ratio fixed at zero inbreeding means holding the absolute rate high: this panel
carries 1,238 duplicated carrier positions per two million where tomato has 150 to 590, because an
outbred panel is two and a half times as heterozygous. The ratio is the quantity the measurement
pins; the absolute count is not tomato's.)*

The two inflations are much closer here — 39% against 27% — and the coefficient moves by **0.09 in
absolute terms** where the inbred panel moves by **0.15**. So it is smaller and it does not vanish.
**What it does instead is leave the interval:** −0.092 is a homozygote *deficit*, and
[`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §5.1 already refuses a
value outside `[0, 1]`. So on an outbred panel the missing class trips that refusal; **on a selfing
panel it does not — 0.447 is a perfectly legal number** — and that is the case this caller is aimed
at.

**Some of the excess is the model repeating itself.** Ignoring the class adds 38% to `Hobs` from the
carrier positions alone, and the measured excess is 51%. The remaining 13 points come from the
positions where a sample has no read: there the genotype posterior is the prior, whose heterozygosity
is `Hexp · (1 − F_hom_excess)` (spec §3.2), so an inflated frequency density raises `Hobs` a second
time at the fifth of positions that carry no evidence. **The class does not only add the positions it
occupies; it re-prices every position that has nothing to say.**

---

## 5. What the cohort pattern finds, and what it misses

**It finds nearly every duplication that five or more of fifty samples carry, and almost none carried
by fewer.**

At three reads a position, of the 101 duplicated positions somebody carries:

| carried by | drawn | `coverage` finds | `pattern` finds |
|---|---:|---:|---:|
| 1 sample | 13 | 10 | **0** |
| 2 to 4 | 24 | 23 | **2** |
| 5 to 24 | 48 | 48 | 43 |
| 25 or more | 16 | 16 | 16 |
| **all** | **101** | **97** | **61** |

At twenty-five reads, of 90: `coverage` finds 87 and `pattern` 52, with the same shape — 42 of 42 and
10 of 10 in the two crowded bands, **0 of 28 and 0 of 10** in the two sparse ones.

**Wrongly calling a real variant duplicated is not the failure mode.** `pattern` does it to **2 of 822**
genuinely segregating positions at three reads and 2 of 793 at twenty-five — one in four hundred.
`coverage` does it to none. Neither calls more than one monomorphic position duplicated.

**And what it misses costs little, because the positions it misses are carried by one sample each.**
The right count is not positions missed but **wrong genotypes left behind**: a position missed because
one sample of fifty carries it leaves one, where a position found because thirty carry it removes
thirty. At three reads the drawn panel holds **1,255 (position, sample) pairs where a sample really
carries a duplication**, and `pattern` leaves **113 of them** — 9 in every 100 — against `coverage`'s
**5**. That is why missing 40 of 101 positions costs `Hobs` three points where ignoring the class
entirely costs fifty-one: **`pattern` removes about 94% of the damage.**

**Depth barely matters to it.** Twenty-five reads a position buys `pattern` nothing over three: 52 of
90 against 61 of 101, and the same clean split at five carriers. The pattern is a statement about *how
many samples* show the position, not about how well any one of them is read — which is why the
low-depth archive is not the hard case here that it is for the coverage summary.

**The pattern works because the panel is inbred, and that is not a detail.** Its whole evidence is the
*absence* of samples homozygous for the non-reference allele, and inbreeding is what makes those
common: at a carrier frequency of 0.2, a real variant leaves 2 samples of 50 homozygous when nothing
is inbred and **7** when the inbreeding coefficient is 0.6. Drawn with no inbreeding, `pattern` finds
**42 of 219** duplicated positions instead of 61 of 101, and the split by carrier count is stark — 0
of 27 carried by one sample, 0 of 46 by two to four, **2 of 106 by five to twenty-four**, and 40 of 40
by twenty-five or more. **On an outbred panel it finds only what the whole panel carries.** `coverage`
on the same panel finds 210 of 219, including 21 of the 27 carried by one sample.

**How many samples it takes.** The same panel and the same three fits at three reads a position, with
the panel size the only thing changed. The column that matters is what is left of the artefact: the
(position, sample) pairs where a sample really carries a duplication and the fit does not call the
position duplicated.

| samples | `pattern` leaves | `pattern` `Hobs` | `coverage` `Hobs` | `no-class` `Hobs` | inbreeding coefficient: truth / `pattern` / `coverage` / `no-class` |
|---:|---:|---:|---:|---:|---|
| 10 | 93 of 238 — **39%** | +21.2% | +1.3% | +60.0% | 0.603 / 0.550 / 0.592 / 0.434 |
| 25 | 109 of 686 — **16%** | +6.8% | +0.0% | +64.2% | 0.606 / 0.599 / 0.613 / 0.444 |
| 50 | 113 of 1,255 — **9%** | +3.0% | +1.0% | +50.6% | 0.594 / 0.581 / 0.586 / 0.447 |

**Ten samples is not enough and twenty-five is.** The coverage fit is flat across the range — a
window's coverage says the same thing whoever else is in the run — while the pattern needs a cohort
before the absence of homozygotes means anything, which is the whole of the difference between the two
discriminators.

**A coverage window that collects too few aligned bases costs some of the advantage and not all of
it.** At three reads a position the stored 500 bp grid of
[`parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md) §4 holds 1,500
aligned bases, an eighth of what the probe report says the classification needs. Summed to that floor
the coverage fit finds 97 of 101 duplicated positions and leaves 5 of 1,255 carrier positions; left
unsummed at 1,500 it finds 85 and leaves 34 — still under a third of the pattern's 113, and 4 of the
13 private duplications instead of 10. **So the summing rule is worth having and it is not what
decides this question.**

**Its blind spot and the harm's location are the same place, and they cancel.** The positions
`pattern` misses are the ones one or two samples carry, and §4 is why those are the cheapest to miss:
a duplication few samples carry inflates expected heterozygosity nearly as hard as it inflates the
observed one, so it mostly cancels out of the inbreeding coefficient anyway. **What `pattern` finds is
exactly what does the damage.** That is not a designed property, and it is the reason its residual
error on the coefficient is 2.3% where ignoring the class costs 24.8%.

**Pushed the other way, it fails differently: it becomes greedy.** On the private-heavy panel — half
the duplications carried by one sample of fifty, which the eight-sample tomato counts rule out (that
truth would put 83 of every 100 doubled windows in exactly one of eight, where 48 is measured) — the
duplicated class is ten times heavier, and `pattern` starts calling real variants duplicated: **80 of
847** against `coverage`'s 1, plus 42 monomorphic positions against 14. `Hobs` then comes out 2.5%
*low* rather than high. So the failure mode reverses with the carrier spread, and **the coverage
summary's value in that regime is not finding more duplications — it is keeping the class from eating
the frequency spectrum.**

**One number neither fit recovers is the class's own weight.** The truth puts 1.07 duplicated
positions in every thousand; `coverage` fits 2.20 and `pattern` 1.86. Both are about twice the truth
while both sort the positions correctly, which says the weight is absorbing probability mass from
positions it only partly claims. **So the class's weight should not be emitted as a measurement of how
much duplication a sample carries** — the per-position posteriors are what carry that, and they are
what `Hobs` reads.

---

## 6. What this cannot say

- **There is still no truth set behind the tomato rate.** The 150 to 590 duplicated positions per two
  million, and the 668 real ones they are measured against, come from window coverage and
  alternative-read fraction on real alignments — not from a validated list of duplications
  ([`duplicated_locus_probe_2026-08-12.md`](duplicated_locus_probe_2026-08-12.md) §7). Everything here
  is as good as that ratio.
- **The carrier-frequency shape rests on 84 windows in eight samples.** Turning eight samples' counts
  into a frequency a fifty-sample panel would show is an extrapolation, and it is the input the
  cohort-pattern answer is most sensitive to: how well the pattern works is almost entirely a
  question of how many samples carry each duplication.
- **The coverage fit is handed a perfect coverage model.** It knows exactly the scatter the readings
  were drawn with, and each position's reading is independent of its neighbours'. Real windows share
  one reading across every position inside them, and the depth-against-GC curve is fitted from the
  sample rather than known. **So this arm is an upper bound on what a coverage summary buys**, which
  is the right bound for the decision: if the cohort pattern matched it, the summary would be
  droppable.
- **Both fits are handed the third class's true shape** — carriers independent given the carrier
  frequency, each reading exactly a half, the frequency itself from a `Beta`. That is generous to
  `pattern` in the same way the coverage model is generous to `coverage`, and it is why the fair
  reading of the whole table is the *comparison* between the two rather than either one's distance
  from the truth.
- **Depth is drawn as a plain Poisson**, with no GC content and no mappability. That is why neither
  fit is allowed to read a position's own read count: in this panel a doubled depth would be a much
  cleaner signal than it is in a real alignment, and giving it to either fit would measure the draw
  rather than the estimator.
- **100,000 positions, against the two million a real fit keeps.** All three fits see the *same*
  positions, so the differences between them are exact; what the locus budget adds is a common offset
  from the truth, shared by all three.
- **One drawn panel at each depth**, so a difference of one or two percent between two fits is not
  separable from the draw. The differences that carry the conclusions are ten times that.
- **A duplication carrier reads at exactly a half.** A sample carrying three copies would read a third
  or two thirds. The spread the probe report measured inside two-copy windows — about one position in
  five outside 0.4 to 0.6 — is what counting noise on 50 reads alone predicts, so a half plus counting
  noise fits the data; a panel with mixed copy numbers is untested.
- **The panel has no population structure.** Every ordinary position's genotypes are drawn from one
  frequency under one inbreeding coefficient. A panel that is landraces from several regions would put
  more positions at middling frequencies, which is where the third class competes for them — so the
  two-in-822 rate at which `pattern` calls a real variant duplicated is measured on the easy case.
- **Diploid throughout**, as the genotype prior is.

---

## 7. Changes this would propose to the specs

*Nothing under `spec/` or `arch/` is edited by this report. These are the changes I would make.*

**In [`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md):**

1. **§11 question 11 — close it.** The cohort pattern identifies the class from about twenty-five
   samples up and does not below; ignoring the class costs a quarter of the inbreeding coefficient on
   a fifty-sample selfing panel at three reads, and the two inflations do not cancel. The
   recommendation of §1 above is what I would write into the close.
2. **§2.2, the paragraph that files the cohort pattern as "a bonus rather than the mechanism"** —
   on a cohort of this size it is most of the mechanism. Replace with the measured split: every
   duplication five or more of fifty samples carry, none carried by fewer, 9 of every 100 carrier
   positions left behind, and one real variant in four hundred wrongly called duplicated.
3. **§2.2's "a real half-frequency variant produces that pattern with probability 0.5⁵⁰"** — true, and
   it is not the number that governs. The duplications that matter are not at a frequency of a half;
   they are carried by anything from one to fifty samples, and the pattern's power collapses below
   about five carriers. The sentence invites the reader to conclude the pattern is decisive, and what
   is decisive is the carrier count.
4. **§2.2's leaning that the pattern "does not suffice alone at three reads a site" is wrong in its
   reason.** Depth is nearly irrelevant to it — twenty-five reads a position buys the pattern nothing
   over three. What it needs is *samples*, not reads.
5. **§5.1's `[0, 1]` constraint on `F_hom_excess`** — add what it does and does not catch. On an
   outbred panel the missing class drives the coefficient to −0.09 and the constraint refuses it, so
   the failure is loud. **On a selfing panel it lands at 0.4471 against a true 0.5942, which is a
   perfectly legal value**, so the constraint is not a safeguard against this and nothing else in the
   fit is either.
6. **§3.2** — the reporting requirement (emit how many kept positions carried a read) gains a size.
   About a quarter of the excess `Hobs` from a missing third class comes from positions with no read
   at all, where the posterior is the prior and an inflated frequency density is charged a second
   time.
7. **§6.1, what the route does at one sample** — add that the third class goes with the per-locus
   site class *unless* the coverage summary is kept, since the pattern needs a cohort and the coverage
   reading does not.
8. **Anywhere the class's fitted weight is treated as a quantity** — it should not be emitted as a
   measurement of how much duplication a sample carries. Both fits recover it about twice too large
   while sorting the positions correctly.

**In [`parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md):**

9. **§4 — change the summary's status from required to used-when-available**, with the two-phase
   run's second pass over every sample's pileup conditional on cohort size rather than automatic,
   and the numbers of §1 beside it. §4's own decision that the summary is never stored already makes
   this a cost the fit chooses to pay or not; what changes is the default.
10. **§4's 12,000-aligned-base floor stays, and gains a size for what missing it costs**: left
    unsummed at 1,500 aligned bases the coverage fit leaves 34 carrier positions of 1,255 where the
    summed window leaves 5 — worse, and still four times better than no coverage at all.

**Neither document should change what it says about the class's grain.** The (position, sample) pair
is right, and this measurement is a second reason: the coverage discriminator is per sample and works
at any panel size, and the cohort pattern — which is per position — is exactly the one that fails on
the duplications a single accession carries.
