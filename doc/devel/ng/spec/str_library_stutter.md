# Each read group's own stutter: a multiplier fitted jointly with the repeat-class rates

**Status:** design, rewritten 2026-11-02 for the penalised fit; **approved by the owner the same
day**, with the answers of §11. It replaces
the design of 2026-10-10, which the owner approved but which was written for the climb, the fit
deleted on 2026-10-16. **What carries over unchanged from that design:** the problem (§1.1), the
owner's three rulings (§1.3), the model (§2), the pull toward one (§3.5) and how calling uses a
multiplier (§6). **What is new:** where the multipliers sit in the penalised fit (§3), what it costs
(§3.2, measured), the range of inputs (§4), the duplicated-tract rules (§5), the parameters file at
format version 8 (§7) and the renaming (§9). **The owner answered §11's three questions on
2026-11-02**: the name `ReadGroupStutter`; **one multiplier per read group and motif length**, not one
across every motif length, unless evidence shows the motif lengths behave alike; and kimura is free
for a couple of days for the measurements.

It amends the penalised fit ([`str_penalised_fit.md`](str_penalised_fit.md), whose §1.4 and §6 defer
this work here), the duplicated-tract rules ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md)
§4.6.1), calling's stutter lookup ([`read_likelihoods.md`](read_likelihoods.md) §4.4) and the
parameters file ([`parameters_file.md`](parameters_file.md)). Build order: §10 and
[`../../implementation_plans/str_library_stutter.md`](../../implementation_plans/str_library_stutter.md).

**Words.** This document says **read group** where the 2026-10-10 design said *library* (owner,
2026-11-02: "per read group (not per library)"). A read group is the grain the census already keeps
apart across a cohort: two samples claiming one read group are refused (`census.rs`). The code still
says `Library`; §9 renames it.

---

## 1. What this is

### 1.1 The problem

**`estimate-parameters` fits one stutter for every read group of the cohort.** It puts every read
group into one slippage group (`every_read_group_pooled`, `src/run/census_fit.rs:326`), the
repeat-tract fit refuses evidence declaring more than one (`fit_repeat_tracts_logging`,
`ssr_fit/penalised_fit.rs:1032`), and calling scores every read group's reads at its class's one slip
rate. *A class is every repeat tract with one motif length and one reference repeat count; the code
calls it a stratum.*

**Read groups do not stutter alike.**

- **Tomato, 51 accessions, one read group each** (2026-07): the share of a sample's reads at a length
  other than its most common one ran from 0.37% to 1.98% between samples, about fivefold; the spread
  was 66.6 times what a binomial allows. Within one motif length the slip level of the most and the
  least stuttering read group differed 4.9-fold in homopolymers and 5.6-fold in dinucleotides. It
  partly follows each read group's duplicate rate (Spearman ρ = +0.48), a sign of PCR.
- **The stuttering simulation**, where each sample's reads slip at the class rate times its own
  multiplier (log-normal, standard deviation 0.4 on the log scale; 0.33 to 3.45 over 50 samples):
  the pooled fit recovered the mean (median fitted over true 1.01 over 91 classes at 30 reads a
  sample) and no single sample. Calling with each sample's true stutter beat calling with the fit in
  all 16 runs of the first benchmark, by 0.1 to 6.1 points. *That was measured on the climb, before
  2026-10-16; §8 measures it again on the penalised fit.*
- **Possibly the same cause, 2026-11-02.** At tomato2 tracts where all 59 accessions read one length,
  the off-length reads vary from tract to tract 2.1 times (homopolymers), 2.7 times (dinucleotides)
  and 3.5 times (3 to 6 bases) as much as a binomial predicts; on HG002, one sample, 0.65 to 1.15.
  Read groups slipping at different rates would do this, because at 3 reads a sample each tract is
  read by a different mix of accessions. Whether they do is §8's last check.

### 1.2 What changes

**Each read group gets one stutter multiplier for each motif length, and its slip rate in a class is
the class's rate times its multiplier for that class's motif length.** The class rates and the multipliers are fitted together, so a class rate is
what the read groups share once each read group's own level is allowed for.

- **Only how often a read slips varies by read group.** Which way it slips (the shorter share) and
  how far (the fall-off) stay one per class. *An assumption with a known exception:* on tomato,
  dinucleotide read groups that slip more also slip with less of their mass at one repeat shorter
  (ρ = −0.69 between the two, after removing reads at lengths other accessions carry). §12 defers it.
- **One multiplier per motif length** *(owner, 2026-11-02)*: "I would lean to one per motif length
  unless we have evidence that they all behave the same." A read group whose PCR inflates
  homopolymer slips more than dinucleotide slips gets two different numbers. No measurement yet says
  whether read groups rank alike across motif lengths; §8.6 makes it, and if they do, one multiplier
  across motif lengths goes back to the owner as a simplification.
- **The substitution rate inside a tract stays pooled**, one a class (§12).

### 1.3 The owner's rulings (2026-10-10, standing)

1. **Not a shared stutter first and a multiplier after.** "I don't want to risk biasing towards the
   most common stutter observed in the libraries." The class rates are refitted with the multipliers
   in place, alternating until the multipliers settle (§3.1).
2. **Thin evidence is pulled toward a multiplier of one**, by a few pseudo-reads (§3.5).
3. **One multiplier per read group, not a full set of rates per read group.** A read group at 3 reads
   a position has too few slipped reads to fit about 100 classes of its own; one number is well
   determined. *Refined by the owner on 2026-11-02 to one a read group and motif length (§1.2): five or
   six numbers a read group, each fitted over all of that motif length's classes.*

### 1.4 Vocabulary

- **Class rate** `Lₛ`: the slip rate of class `s` for a read group whose multiplier is one.
- **Multiplier** `mᵣ`: read group `r`'s stutter relative to the class rates of one motif length. Each
  motif length has its own; where a formula is about one motif length, the motif length is left out
  of the symbol.
- **A read group's slip rate** in class `s`: `λᵣₛ = min(mᵣ · Lₛ, max(Lₛ, 0.95))`, `mᵣ` the read
  group's multiplier for the motif length of `s`. The 0.95 is
  `MULTIPLIED_LEVEL_CEILING` (`ssr_fit.rs:157`): a rate of one would say no read reports its allele,
  and the cap lifts no read group above it unless the class rate itself is above it.
- **Slippage group**: a set of read groups the run declares share one set of class rates. Today, and
  in this design, there is one; the multipliers live inside it (§12).
- **Row of reads**: one sample's reads at one tract, the count at each recorded length. The
  expectation step scores a row under every genotype.

---

## 2. The model

### 2.1 A read's slip

A read of read group `r` crossing a tract of class `s` reports a length other than its allele's with
probability `λᵣₛ`; given that it slipped, the direction and the distance follow the class's shorter
share and fall-off (`Slippage::read_probabilities`, `ssr_fit.rs:105`). The code already computes
`λᵣₛ` (`library_level`, `ssr_fit.rs:170`), and returns the class rate itself, to the bit, at a
multiplier of one.

### 2.2 What fixes the scale

`mᵣ · Lₛ` is unchanged when every multiplier of a motif length is doubled and every class rate of it
halved. **Each motif length's multipliers are scaled so that their geometric mean, each read group
weighted by its reads crossing that motif length's tracts, is one.** A class rate is then the rate of a read group of typical stutter, and a multiplier says how far
a read group is from typical, which is what a reader of the file needs and what the pull of §3.5 pulls
toward.

### 2.3 One read group alone

**Its multiplier is not fitted: it is one, written `defaulted`, with no error** (§4). With one read
group there is nothing to compare it with; its own rate is the class rate.

---

## 3. The fit

### 3.1 Where the multipliers sit

**The penalised fit runs one motif length at a time, in a loop of rounds** (`str_penalised_fit.md`
§2.2): each round weighs every sample's genotypes under the current numbers (the expectation step),
then draws every curve through the counts that gives (the maximisation step). **With a multiplier per
motif length, the multipliers become one more set of numbers that loop fits**, and each round
becomes:

1. **The expectation step, keyed by read group** (§3.2), at the round's class numbers and
   multipliers: today's counts, plus each read group's expected slipped and exactly-read reads in
   every class.
2. **The curves, with the multipliers held** (§3.3): every class number drawn as today, the slip
   rate's curve through each read group's counts at its own multiplier.
3. **The multipliers, with the curves held** (§3.4), from the same expected counts, then **scaled**
   (§2.2).

The first round starts at every multiplier one. **The loop stops as today**, when no number moves by
more than the tolerance (10⁻⁷) or the rounds run out (2,000), and the multipliers count among the
numbers that must stop moving. The duplicated-tract rules run as today, between fits of the motif
length (§5).

**Why this answers ruling 1.** The class rates and the multipliers are fitted in the same loop, from
the same expected counts, each round; at no point is a pooled stutter fitted to convergence and a
multiplier laid on top. Steps 2 and 3 each maximise the round's objective over one block of numbers
with the other held, which is expectation–conditional maximisation (ECM, Meng and Rubin 1993): every
round raises the objective, as plain EM does, and where the loop stops neither block can raise it.

**What a multiplier per motif length changes against the 2026-10-10 design.** No motif length waits
for another, so the outer rounds over whole fits that design needed are gone: each motif length is
fitted once, as now. **The cost is evidence**: a read group's multiplier for a motif length rests on
that motif length's slipped reads only. On tomato2 homopolymers hold 90% of the rows with reads and
pentanucleotides 0.04%, so a read group's pentanucleotide multiplier rests on almost nothing and the
pull of §3.5 holds it near one. That is the honest answer for it, and its error is withheld (§3.6).

### 3.2 The expectation step, keyed by read group

**Today a sample's read groups are pooled into one row before the row is looked up**
(`pool_reads_into` inside `DistinctRows::of`, `ssr_fit/expected_counts.rs:175`), and each distinct
row is scored once however many samples hold it. That is the first of the fit's two exact speed-ups
(`str_penalised_fit.md` §2.8). **Once read groups slip at different rates, two samples' rows are the
same only if their reads are the same and come from read groups at the same multiplier.**

**The row becomes a list**: for each read group of the sample that put a read on the tract, which
read group and its reads, in increasing read-group order (the evidence already holds exactly this,
`SampleTractReads::by_library`). A genotype's log-likelihood of the row is the sum over its read
groups of that read group's reads scored at its own rate, because reads are independent given the
genotype:

`ln P(row | g) = Σᵣ Σ_bucket readsᵣ(bucket) · ln P(a read of r reports bucket | g)`.

The log chance table (`ln_read_chances`) is built once a read group a round, from `λᵣₛ`; it is 190
genotypes by 17 recorded lengths (±8 repeats) at a ±9 allele window, 26 kB a read group, 56 MB at
2,169 read groups, one class at a time.

**What the expectation step returns, besides today's counts**: per read group, its expected
exactly-read reads `Nᵣₛ` and slipped reads `Sᵣₛ` in the class. Each read's share between "read
exactly" and "slipped" is taken under its own read group's rate. The direction and distance counts
stay pooled, because the shorter share and the fall-off are one a class. The pooled `exact` and
`slipped` stay too, as the sums over read groups, for the run log, the file's
`expected_slipped_reads` and the duplicated-tract rule.

**Keyed by read-group index, not by multiplier value.** Keying by the multiplier's value would pool a
sample's read groups in the first round, where every multiplier is one, but step 3 needs each read
group's own counts even when two share a multiplier, and pooling loses them. So **a run with one read
group keeps today's pooled table**, which keeps it today's to the bit, and a run with two or more
builds the keyed table in every round.

**What it costs, measured 2026-11-02** (`examples/ng_str_distinct_rows_by_read_group.rs`, the
evidence gathered as `estimate-parameters` gathers it, the depth rule applied):

| cohort | read groups | rows (tract × sample with reads) | distinct rows, pooled (today) | distinct rows, keyed by read group |
|---|---|---|---|---|
| tomato2, 59 accessions at about 3 reads | 59 | 773,836 | 9,199 (84 rows each) | 61,474 (12.6 rows each) |
| stuttering simulation, 50 samples at 3 reads | 50 | 127,804 | 9,723 (13 rows each) | 58,261 (2.2 rows each) |
| kimura's tomato cohort, 2,169 samples, median 11 reads | 2,655 | 119,072,017 | 4,795,352 (25 rows each) | 35,042,058 (3.4 rows each) |

*The kimura row (2026-11-02) was counted with the code just before main's reference checksum changed
(`b626ead2`'s parent), because its psps predate that change, and without the depth rule, which that
code lacks; on tomato2 the depth rule left out 240 tracts, so the counts would move little.* By motif
length the ratio runs from 6.6 (dinucleotides) to 11.3 (tetranucleotides); the largest class, a
homopolymer at 13 repeats, holds 2.69 million keyed rows against 332,000 pooled.

So **keying by read group scores 6.7 times as many rows a round on tomato2 and 6.0 times on the
simulation**, the same in every motif length to within 4.7 and 8.2 times. The expectation step is
most of a round's time (steps 2–3 of the penalised fit are 1–9%), so a round costs about six times
today's. Today's repeat-tract fit takes 38 s on tomato2 (Mac, the duplicated-tract rules included);
with every round keyed, about four minutes, if the loop needs as many rounds as today. *An estimate
from the row counts, measured in §8.7.*

**At thousands of read groups the cost is the same factor**: 7.3 on kimura's cohort, against 6.7 on
tomato2, because the pooled table there already repeats a row only 25 times, not 84. **What it
costs in memory**: a scored row holds one number a genotype, 190 at ±9, so the largest kimura class
holds about 4 GB of scored keyed rows against about 0.5 GB pooled, one class at a time. **What is
built**: the exact keyed table. **What is measured, on kimura** (free for a couple of days, owner
2026-11-02): the fit's time against the penalised fit's goal of one hour (`str_penalised_fit.md`
§1.3, goal 4), which at seven times today's rounds it may miss. *Owner, 2026-11-02:* measured as
one round on part of the genome, on kimura's existing psps, which are not regenerated; so the
timing runs on this branch's fit ported to the code those psps were walked with. **If it misses the goal**, the fallback is to round each multiplier to a grid of 1%
steps for the expectation step only, so that rows are keyed by the grid value and every read group
within 1% shares rows; a 1% step is well below a homopolymer multiplier's error at 3 reads (about 5%,
§3.6). That changes numbers by the rounding, so it is the owner's to approve, after the measurement,
not built now.

### 3.3 The class rates, with the multipliers held

**A class's slip rate is no longer a binomial of pooled counts**: each read group's reads slip at
`mᵣ · Lₛ`. The class's log-likelihood in its rate `p = Lₛ` is

`Σᵣ [ Sᵣₛ · ln(mᵣ p) + Nᵣₛ · ln(1 − mᵣ p) ]`,

plus the stutter prior of `str_penalised_fit.md` §2.10, unchanged: one read's worth at multiplier one,
`0.10` of it slipped. **The curve fitter already takes any per-class log-likelihood with its slope and
curvature** (`CurveData::evaluate`, `penalised_fit.rs:2378`), so this is one new kind of curve data,
not a new solver. On the logit scale `η` the fitter works on, with `p = logistic(η)`:

- slope: `p(1 − p) · Σᵣ [ Sᵣₛ / p − Nᵣₛ · mᵣ / (1 − mᵣ p) ]`;
- information: minus the second derivative, floored at a small positive number where the read groups
  with `mᵣ > 1` near the cap make it bend the wrong way, as the tabulated data already does
  (`penalised_fit.rs:2477`); the fitter's step halving keeps the climb honest there.

**At every multiplier one it is the binomial of today**: the slope becomes `Sₛ − p (Sₛ + Nₛ)`, the
existing one. The shorter share, the fall-off, F and the allele-length shares keep their curve data
unchanged.

### 3.4 The multipliers, with the classes held

With every class number held, read group `r`'s log-likelihood in its multiplier `m` is

`Σₛ [ Sᵣₛ · ln(m Lₛ) + Nᵣₛ · ln(1 − m Lₛ) ] + P · (ln m − m)`,

summed over every class of the motif length being fitted, the last term being the pull of §3.5. It is concave in
`m`, one number, solved by Newton's method inside `(0, 0.95 / maxₛ Lₛ)`. Each read group is solved on
its own, so thousands of read groups cost thousands of one-dimensional solves, about a millisecond in
all.

**The genotype a sample carries** enters through `Sᵣₛ` and `Nᵣₛ`: they are the expectation step's own,
each read shared between "slipped" and "read exactly" under the sample's genotype weights at the
round's numbers. Nothing else is needed; the 2026-10-10 design took the genotype from the climb's
likelihood, which is gone.

### 3.5 The pull toward one (ruling 2)

**`P` pseudo-reads that slipped exactly as often as a typical read group are added to every read
group**: the term `P · (ln m − m)`, whose maximum on its own is `m = 1`. A read group with many
slipped reads is barely moved by it; one with few is held near one. `P` is soft: 20, the per-locus
slip refit's own pull (`DEFAULT_LEVEL_PULL_BACK_SLIPPED_READS`), to be set from the simulation at 3
reads a sample (§8).

**What the pull does depends on the motif length**, because a read group's slipped reads do. On
tomato2 (59 accessions at about 3 reads, the expected slipped reads of the current fit's file summed
over each motif length's classes, divided by 59):

| motif length | slipped reads a read group | error of its multiplier, without the pull | how far 20 pseudo-reads pull it toward one |
|---|---|---|---|
| 1 | about 1,160 | about 3% | 2% of the way |
| 2 | about 160 | about 8% | 11% |
| 3 | about 9 | about 33% | 69% |
| 4, 5, 6 | about 0.3 to 1.2 | not placed | 94% or more |

So **at 3 reads a sample, a tomato read group's multipliers are measured for homopolymers and
dinucleotides, weakly for trinucleotides, and not at all for longer motifs**, where they stay near
one, which is today's pooled stutter. At 30 reads every row of the table has ten times the reads. A
deep cohort or a cohort of long reads moves the line down.

### 3.6 Standard errors

- **A multiplier's error** is the sandwich `H⁻¹ B H⁻¹` of §3.4's objective in `ln m`, with `B` summed
  over tracts (each tract's contribution to the read group's slope, squared), because samples at one
  tract are not independent (`str_penalised_fit.md` §2.6), and the pull's curvature in `H` and `B` as
  the stutter prior's is. Without the pull, and where the cap does not bind, it is about
  `1 / √(slipped reads)` on the log scale: about 8% at 160 slipped reads. **It is conditional on the
  class rates**: it leaves out their uncertainty, so it is narrower than the joint error. §8 measures
  how often it holds the simulation's truth.
- **A multiplier with fewer expected slipped reads than the pull's `P` gets no error** (`NotPlaced`,
  soft): its value leans on the pull at least as much as on its reads, and its error would measure the
  pull's width, not the read group. This follows the stutter prior's rule (`str_penalised_fit.md`
  §2.10.4) with the pull's weight in place of the prior's.
- **The class numbers' errors** are computed as today with every multiplier held at its fitted value;
  the slip rate's slope in the sandwich becomes §3.3's.

---

## 4. The range of inputs

| input | what the multipliers do |
|---|---|
| **one sample, one read group** (HG002 as benchmarked, to be checked) | not fitted: the multipliers' step does not run, the file writes one row a motif length, `value = 1`, `defaulted`, no error. Every number is today's to the bit. |
| **one sample, several read groups** | fitted. They are relative to each other only: the sample's genotypes are one set, so a read group's extra off-length reads at a tract its other read groups read exactly are stutter, not alleles. |
| **a cohort, one read group a sample, at about 3 reads** (tomato2, the simulation) | fitted; homopolymer and dinucleotide multipliers are measured, longer motifs' held near one by the pull (§3.5's table). The class shares, which every sample shares, are what keep a sample's multiplier apart from its own genotypes. |
| **a cohort at 30 reads and above** | fitted, errors narrow; the pull matters little except at the motif lengths that rarely slip. |
| **thousands of read groups** | one multiplier a motif length each: memory is six numbers and their errors a read group, plus the per-read-group counts, `Sᵣₛ` and `Nᵣₛ`, 16 bytes a read group a class, for one motif length's classes at a time. The cost is the keyed table (§3.2), measured on kimura. |

**A read group with no read crossing a motif length's tracts** gets no multiplier for it: `value = 1`,
`defaulted`, and the run log counts them. It weighs nothing in §2.2's mean.

---

## 5. The duplicated-tract rules

The depth rule reads no fitted number and does not change. **The stutter-tail rule changes in two
places** (`parameter_prepass_ssr.md` §4.6.1, `ssr_fit/tract_exclusion.rs`):

1. **The genotype each sample's reads favour** is taken from the row scored at each read group's own
   rate (the keyed table of §3.2), so a sample's reads from a stuttering read group are judged as that
   read group's.
2. **The rate the tract's off-length reads are tested against is the read-weighted mean rate**:
   `λ̄ = Σᵣ nᵣ λᵣₛ / N`, `nᵣ` the tract's reads from read group `r` and `N` all of them. So the
   expected off-length reads are `λ̄ N = Σᵣ nᵣ λᵣₛ`, exactly what the read groups' own rates predict, and
   the test is `K ≥ 2 λ̄ N` and `P(Binomial(N, λ̄) ≥ K) < 0.001`. **The binomial at the mean leans
   toward keeping a tract**: a sum of reads slipping at different rates is less spread than a
   binomial at their mean (Hoeffding 1956), so the true chance of so many off-length reads is lower
   than the test says, and fewer tracts are flagged than a per-read-group test would flag.

**When the rule runs.** As today, between fits of a motif length, until it flags none, against the
class rates and multipliers the last fit returned. Each refit after the rule starts its multipliers
where the last fit left them. A tract left out stays out, so the rule's rounds end.

**What it is expected to change.** On tomato2 the rule left out 324 homopolymer and 102 dinucleotide
tracts at the pooled rate (2026-11-02, with the genotype then taken from the class prior). A tract
read mostly by high-stutter accessions is the case a pooled rate flags wrongly; §8 counts how many
tracts the rule leaves out with and without the multipliers.

---

## 6. Calling

**A cell's frozen slip rate is its class's rate times its read group's multiplier for the class's
motif length**, `λᵣₛ` of §1.4,
set where calling gathers a tract's stutter: `TractScoringFits::gather_for_locus` and
`gather_for_bundle` (`calling/inference/repeat_tract_parameters.rs:305`, `:452`), which already loop
over every read group and look its slippage up per read group. The multiplied rate goes into the
stutter model and into the cell's frozen slippage, so the per-locus slip refit, where it runs, pulls
toward the read group's own rate.

- **`--defaults`, and a file of version 7 or older**: every multiplier one.
- **A read group with no row in the table**: the file is refused when it is read (§7), because the
  table covers every read group the file names. A read group the file does not name at all is already
  defaulted and counted (`slippage_defaulted_by_an_unknown_read_group`).

---

## 7. The parameters file, version 8

**A new table, `repeat_tracts.stutter_multiplier_by_read_group`**: one row a read group and motif
length, in the shape of `base_quality_calibration.by_read_group`:

```toml
stutter_multiplier_by_read_group = [
    { read_group = 0, period = 1, stutter_multiplier = { value = 1.42, standard_error = 0.04, warrant = "fitted_here" }, slipped_reads = 1160.2 },
    { read_group = 0, period = 4, stutter_multiplier = { value = 1.03, warrant = "fitted_here" }, slipped_reads = 0.3 },
]
```

- `warrant`: `fitted_here`; `defaulted` for a read group alone in its slippage group (§2.3) or with no
  read crossing the motif length's tracts (§4), with no error; a person may write `supplied`.
- `standard_error`: absent where §3.6 withholds it.
- `slipped_reads`: the read group's expected slipped reads over the motif length's classes, without
  the pull.
- **Validation**: a finite positive value; every read group of the file has exactly one row for every
  motif length the file holds a fitted class of.
- **Format version 8.** A version-7 build refuses a version-8 file with its existing message; this
  build reads a version-7 or older file with every multiplier one.
- **The file's note** says what a multiplier is, that a class rate is the rate of a read group of
  typical stutter (§2.2), that `mᵣ · Lₛ` is capped at 0.95, and that the shorter share and the
  fall-off are one a class.

---

## 8. How we know it works

1. **Drawn cohorts with known multipliers** (unit tests): class rates and multipliers drawn, reads
   drawn from them, the fit run end to end. At 3 and 30 reads, 2 to 50 read groups: every multiplier
   within two of its errors of the truth about 95 times in 100, and the class rates unbiased. One
   sample, two read groups, at 30 reads: both recovered relative to each other.
2. **One read group: today's fit to the bit** (unit test), and a cohort fitted with its multipliers
   held at one equal to today's to within rounding (the keyed table sums in another order).
3. **Each piece fails a test without it**: the keyed row (a two-read-group sample scored at two
   rates), §3.3's curve data (its slope against a finite difference, and the binomial at every
   multiplier one), §3.4's solve (the pull's limits: no reads gives one; many reads give the data's
   own ratio), the scaling, the tail rule's mean rate.
4. **The stuttering simulation** (`tmp/stutter_data/sim/` in the main checkout, `libraries.tsv`
   holding each sample's true multiplier): the fitted multipliers against the true ones sample by
   sample, at one sample 30×, 50 samples at 3× (inbred) and the PCR-free arm; and calling with the
   fit against calling with the true stutter (`ng_true_stutter`). **The bar**: the multipliers recover
   the truth within their errors, and the fit closes most of the gap to the true stutter. `P` is set
   from it. *The simulator gives each sample one multiplier for every motif length, so it tests the
   fit's recovery and not whether real read groups differ by motif length; that is item 6.*
5. **HG002 chr20, 30× and 5×**: first check it has one read group; then the file's numbers and the
   calls must not move.
6. **tomato2, and kimura's cohort**: the multipliers' spread against the 0.37–1.98% measured per
   accession; **whether the motif lengths behave alike** — each read group's homopolymer multiplier
   against its dinucleotide one, the two motif lengths with enough slipped reads at 3× (§3.5), as a
   correlation and as the share of read groups whose two multipliers differ by more than their errors
   allow. If they agree, one multiplier across motif lengths goes back to the owner (§1.2); the
   stutter-tail rule's
   tracts with and without; and **the binomial spread of §1.1 measured again with each read group's
   own rate**. If read groups' different stutter is its cause, the 2.1–3.5 falls toward one.
7. **Time and memory**: tomato2 and the simulation against today's fit; kimura against the one-hour
   goal (§3.2).
8. **The cross-platform checksums** (`src/cli/cross_platform_digests.rs`) re-recorded only after
   their move is explained; the four-accession oracle has four read groups and moves.

---

## 9. The renaming

**First build step, alone, with no number moving.** The code's names for the repeat-tract fit's
read groups say *library*. `ReadGroup` cannot be the new name of `Library`:
`read::input::read_groups::ReadGroup` already exists, and is a different thing (a read group as the
alignment header declares it). **Approved by the owner, 2026-11-02:**

| today | renamed to |
|---|---|
| `Library` (which read group, its slippage group, its multiplier) | `ReadGroupStutter` |
| `LibraryIndex` | `ReadGroupIndex` |
| `library_level` | `read_group_level` |
| `StratumEvidence::libraries` | `StratumEvidence::read_groups` |
| `StratumEvidence::check_libraries` | `check_read_groups` |
| `SampleTractReads::by_library` | `by_read_group` |
| `one_stutter_for_every_library` | `one_stutter_for_every_read_group` |

About 110 uses, in `ssr_fit.rs`, `expected_counts.rs`, `penalised_fit.rs`, `tract_exclusion.rs`,
`census_fit.rs`, the fitting bench and three examples. **Out of scope**: *library* in the SNP/indel
half (contamination, the error-rate fit, sequencing batches), which names a different grain and was not
part of the ruling. The parameters file's `library` key (the read group's declared library) is a fact
about the input and stays.

---

## 10. Build order

Each step its own commit, with a test that fails without it; the plan document's order, renamed.

1. **The renaming** (§9). Checksums unchanged.
2. **Per-read-group counts and the keyed table**, built and tested, used by nothing: every number
   unchanged; a keyed table at every multiplier one gives the pooled counts to within rounding, and
   each read group's counts sum to the pooled ones.
3. **The estimators**: §3.3's curve data and §3.4's solve with its pull, scaling and error, tested
   against drawn truth with the classes held at their true values.
4. **The multipliers in each motif length's loop** (§3.1) and the tail rule's mean rate (§5), wired
   into `estimate-parameters`; the refusal of more than one slippage group stays. **Numbers move.**

Before step 2, the counting tool of §3.2 runs on kimura's cohort, so the keyed table's size at
thousands of read groups is known before the code that holds it is written.
5. **The parameters file, version 8** (§7).
6. **Calling reads the multipliers** (§6). **Calls move.**

Then §8's measurements, and a report.

---

## 11. Decisions (owner, 2026-11-02)

1. **The name: `ReadGroupStutter`** for today's `Library`, and the rest as §9's table.
2. **One multiplier per read group and motif length**: "I would lean to one per motif length unless
   we have evidence that they all behave the same." The lead had recommended one across motif
   lengths. §3.1 is built this way; §8.6 gathers the evidence that could reverse it.
3. **Kimura is free for a couple of days** for the measurements: the keyed table's size (§3.2, before
   build step 2) and the fit's time (§8.7). The exact keyed table is built; the 1%-grid fallback only
   if kimura misses its hour, and then back to the owner.

---

## 12. Deferred

- **Extrapolating the multipliers the reads cannot place** *(owner, 2026-11-02)*: "for the ones that
  we can't calculate we could take the estimated ones and do some kind of extrapolation. But to decide
  how to do it we need estimates of real libraries to check the behaviour." At 3 reads a sample that
  is the 3- to 6-base multipliers (§3.5); they could be predicted from a read group's homopolymer and
  dinucleotide multipliers instead of pulled to one. How, is decided once §8.6 has fitted multipliers
  of real read groups on tomato2 and kimura's cohort, and on any deep cohort where longer motifs are
  placed too. Until then the pull holds them near one.

- **The shorter share and fall-off per read group**, or as a function of its multiplier: tomato's
  dinucleotides show high-stutter read groups slipping with less mass at one repeat (§1.2). Only if
  §8 shows calls lost to it.
- **The substitution rate per read group**: pooled, one a class.
- **Slippage groups other than one**, declared by the run (`slippage_group_by_read_group` already
  exists in the file); the multipliers would live inside each group, scaled within it.
