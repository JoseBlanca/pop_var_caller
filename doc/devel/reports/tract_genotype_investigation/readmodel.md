# Every other constant in the repeat-tract read likelihood — what each does, and which is worth moving

*Measured 2026-09-02 on GIAB HG002, the 20,204 benchmark tracts of
`tmp/tract_qual/ground/tier.bed`, at 30×, one sample, one read group. Everything below is a
fact about that corner: one high-coverage human, short reads, novoalign, 2014 HiSeq chemistry.
Nothing here has been checked at 3× or on a cohort.*

Two lines are already being worked and are excluded here: **fitting the slippage numbers** (the
four direction shares and the two one-step shares, today `StutterModel::hipstr_shipped`) and
**fitting the genotype prior's length spectrum** (today flat). What follows is everything else.

---

## 1. The inventory — every constant a repeat tract is scored under that nobody fitted

A `--defaults` run scores a tract under the numbers in the third column. Every one of them is a
constant compiled into the binary.

| # | constant | where | value | what it decides |
|---|---|---|---|---|
| 1 | `DEFAULT_OUTLIER_WEIGHT` | `src/ng/calling/likelihood/ssr.rs:83` | 0.01 | the share `λ` of a tract's reads the model expects came from somewhere no candidate explains. Taken from the existing production caller (`src/ssr/cohort/em.rs`) and declared inherited |
| 2 | the shape of the junk distribution `U` | spec §4.5, built in `ssr_emission::fill_reachable_lengths` | uniform over the reachable tract lengths | what such a read is assumed to show. Chosen by the owner on 2026-08-19 to be a property of the locus rather than the cohort; the *uniform* shape was never measured |
| 3 | `DEFAULT_SSR_SUBSTITUTION_RATE` | `src/ng/calling/inference/repeat_tract_parameters.rs:130` | 0.001, *defined as* the SNP path's `DEFAULT_ERROR_RATE` | the per-base rate at which a read's tract letters differ from the allele's, wherever the pre-pass fitted nothing — which on a `--defaults` run is every cell |
| 4 | the wrong-base divisor | `alignment::emission::FlatEmission`, spec §3.5 | 3 | a mismatching base costs `ε/3`, not `ε`. Physical (three bases to go wrong into), not a free parameter |
| 5 | `PART_REPEAT_SHARE_OF_WHOLE` | `src/ng/calling/likelihood/stutter_rates.rs:48` | 0.05 | on the **fitted** route, how much stutter mass goes to length changes that are *not* a whole number of motif copies, as a fraction of the whole-repeat mass. Production's `OUT_FRAME_REL`; its own comment calls it a Step-4 placeholder that Step 5 never replaced |
| 6 | the part-repeat shares inside `hipstr_shipped` | `src/ng/alignment/stutter.rs:313`–`315` | 0.01 and 0.01, against whole-repeat shares of 0.05 | the same quantity on the **defaults** route, where it works out at **0.20**, four times the constant in row 5 |
| 7 | `part_repeat_one_step_share` tied to the whole-repeat one | `stutter_rates.rs`, `stutter_rates_for` | one value for both | how fast part-repeat changes fall away with size. HipSTR keeps the two independent (0.95 against 0.80 in its own EM start); ng ties them, declared a placeholder |
| 8 | `MAX_WHOLE_REPEAT_SLIP` | `src/ng/alignment/stutter.rs:65` | 10 repeats | past this a whole-repeat change scores exactly zero and the read falls to the junk term. Production's number, whose own comment calls it "a provisional choice" |
| 9 | `MAX_PART_REPEAT_SLIP` | `src/ng/alignment/stutter.rs:92` | 10 re-indexed steps | the same for part-repeat changes. Split out from row 8 so the two are independently settable; both still 10 |
| 10 | `GEOM_MIN` / `GEOM_MAX` | `src/ng/alignment/stutter.rs:104`, `:106` | 0.01 / 0.99 | clamps on a one-step share, and `GEOM_MIN` doubles as the floor under the derived same-length share |
| 11 | `MIN_SHARE_FROM_THIS_INDIVIDUAL` | `src/ng/calling/likelihood/ssr.rs:182` | 1 × 10⁻¹² | floors `1 − λ − c` positive so `ln` cannot see a negative |
| 12 | equal weights over slip placements | `ssr_emission::letters_over` | `1 / placements` | in an interrupted tract a slip could have landed in any run; the model averages the realisations with equal weight rather than weighting by run length |
| 13 | the unreachable mass is reported, never applied | `SsrScoringContext::unreachable_mass`, and no consumer | — | the stutter distribution loses mass it cannot place (2 parts in 100 at period 1, where the part-repeat branch is unreachable by arithmetic). The row carries the number and never divides by it |

**Row 6 is an inconsistency rather than a placeholder, and it is worth stating on its own.** The
two routes into the stutter distribution disagree about the same quantity by a factor of four: a
tract scored from a fit puts part-repeat changes at 0.05 of the whole-repeat mass, and the same
tract scored from `hipstr_shipped` — which is what every `--defaults` run gets — puts them at
0.20. Neither number was measured. §2.2 below measures it.

**One constant that is *not* in play at a tract, named so nobody looks for it.** The base-quality
multiplier (`DEFAULT_ERROR_PROBABILITY_MULTIPLIER`, 1) never reaches this path: spec §4.3 rules
that the tract substitution term takes a fitted per-*base* rate and never the read's own per-*read*
quality, and the code follows.

---

## 2. What the reads say

### 2.1 The per-base error rate inside a tract — measured

**Method.** 13,302 of the 20,204 benchmark tracts have no GIAB truth record within 5 bases of
them, so HG002 is homozygous reference there and every read should match the reference exactly.
Every aligned base of every read over those tracts was compared with the reference
(`select_homref_tracts.py`, `tract_error_rate.py`; 611,789 reads, no MAPQ filter).

| | bases compared | mismatched | rate | mean reported base quality |
|---|---|---|---|---|
| inside the tract, all periods | 5,828,111 | 11,189 | 0.00192 — **1 in 521** | 35.3 |
| inside, period 1 | 2,992,639 | 3,959 | 0.00132 — 1 in 756 | 35.3 |
| inside, period 2 and above | 2,835,472 | 7,230 | 0.00255 — **1 in 392** | 35.2 |
| the 50 bases either side, all periods | 38,374,419 | 269,860 | 0.00703 — 1 in 142 | 34.0 |

**The flank number is not the control it looks like, and reporting it as one would mislead.**
These tracts sit inside GIAB's tandem-repeat benchmark regions, so the 50 bases either side of one
tract are usually inside another repeat. The rate is flat across those 50 bases — 0.0079 at 41–50
bases out against 0.0109 at 1–5 bases out for homopolymers (`flank_profile.py`), so it is not an
edge effect of the tract — but it is a measurement of *repeat-adjacent* sequence, not of ordinary
genome. **The comparison that settles the question is the tract-internal rate against the
constant, and against what the instrument claims.**

**The instrument's claim is the honest control, and it is out by a factor of six.** A mean
reported base quality of 35.2 inside a tract is an error probability of 3.0 × 10⁻⁴. The reads
mismatch at 1.92 × 10⁻³. **So inside a tract the sequencer's own quality understates the observed
error rate 6.4-fold.**

**On the doc comment this was meant to settle.** `DEFAULT_SSR_SUBSTITUTION_RATE`'s documentation
says base quality inside tracts is systematically worse than outside them, and that 0.001 is
therefore "very likely optimistic at a tract". **The premise is not true on this data** — reported
quality inside a tract is 35.3 against 34.0 in the surrounding 50 bases, i.e. very slightly
*better*. **The conclusion holds anyway, and by a smaller factor than the wording suggests**:
0.001 against a measured 0.00132 at homopolymers (1.3-fold) and 0.00264 at period 2 and above
(2.6-fold). The sentence should be re-grounded on the measurement rather than on the premise.

**The rate is not one rate — a tail of reads carries a third of the mismatches**
(`per_read_inside.py`, reads that span a tract end to end and carry no indel inside it):

| | period 1 | period 2 and above |
|---|---|---|
| reads spanning a homozygous-reference tract | 225,744 | 126,713 |
| carrying an indel inside it | 6,822 — 1 in 33 | 886 — 1 in 143 |
| of the indel-free rest: 0 mismatches inside | 216,114 — 0.9872 | 121,905 — 0.9688 |
| 1 mismatch | 2,534 — 1 in 86 | 2,707 — 1 in 46 |
| 2 mismatches | 220 — 1 in 995 | 629 — 1 in 200 |
| **3 or more mismatches** | **54 — 1 in 4,054** | **586 — 1 in 215** |
| pooled rate | 0.00126 | 0.00264 |
| pooled rate with the 3-or-more reads set aside | 0.00126 | **0.00168** |

At period 2 and above, **1 read in 215 carries 3 or more mismatches inside the tract, and those
reads carry 2,301 of the 6,266 mismatches — 37 in 100 of them.** That is the shape the model's
third term exists for: a bulk at 1 in 594 plus a small population of reads that are simply not
from this tract. A single flat rate has to average the two and lands at 1 in 378.

### 2.2 The part-repeat share — measured

At a homozygous-reference tract the read came from the reference allele, so the change between
what the read shows and the reference **is** the stutter distribution, measured
(`length_changes.py`; reads spanning the tract with at least 3 bases of flank either side).

| | period 1 — 212,749 reads | period 2 and above — 120,420 reads |
|---|---|---|
| showed the reference length | 205,609 — 0.9664 | 119,527 — 0.9926 |
| a whole repeat longer | 2,613 — 0.01228 | 332 — 0.00276 |
| a whole repeat shorter | 4,527 — 0.02128 | 501 — 0.00416 |
| part of a repeat longer | 0 by arithmetic | 23 — 0.00019 |
| part of a repeat shorter | 0 by arithmetic | 37 — 0.00031 |
| **part-repeat mass ÷ whole-repeat mass** | **not defined** | **0.072** (60 reads against 833) |
| what a defaults run assumes | 0.20 | 0.20 |
| what a fitted run assumes | 0.05 | 0.05 |

**At period 1 the part-repeat branch cannot be reached at all**, because every length change is a
whole number of one-base repeats. The two shares still cost 2 parts in 100 of the distribution's
mass, which `unreachable_mass` reports and nothing renormalises (row 13).

**At period 2 and above the measured ratio is 0.072**, which sits between the two constants the
caller uses — 0.05 on the fitted route and 0.20 on the defaults route. **So the defaults route is
2.8 times too generous to part-repeat changes and the fitted route is 1.4 times too stingy**, and
§3 measures that neither matters for a call.

### 2.3 Reads nothing explains — a bound on the outlier weight

Two things put a read beyond every candidate: letters no substitution rate makes plausible, and a
length past a slip cutoff. Counted at the same homozygous-reference tracts:

| | period 1 | period 2 and above |
|---|---|---|
| 3 or more mismatches inside the tract | 54 of 218,922 | 586 of 125,827 |
| a length change past the cutoffs (10 repeats / 10 steps) | 39 of 212,749 | 16 of 120,420 |
| **together, as a share of reads** | **0.00043 — 1 in 2,300** | **0.00479 — 1 in 209** |

**Read literally, the shipped 0.01 is about right at period 2 and above and about 20 times too
high at homopolymers.** §3 finds that read literally is not how the number earns its keep.

---

## 3. What moving each one does to a genotype call

**The measurement.** The repeat-tract row was re-implemented over the evidence ng itself produced
— `tmp/attrib/tier_30x_candidates.tsv`, which carries every distinct tract sequence HG002's reads
showed at each benchmark tract and how many reads showed it — and scored exactly as spec §2.1,
§4.2, §4.3 and §4.5 specify, against GIAB's truth lengths (`constant_sweep.py`). **8,686 tracts
have both a truth genotype and a candidate set: 5,538 at period 1, 3,148 at period 2 and above.**

**This is not ng and its absolute accuracy is not ng's.** The candidate set is the dump's, reads
that ran off their own end inside the tract are absent, and the prior is flat. It calls 0.8629 of
the 8,686 right — 0.8472 at period 1, 0.8904 at period 2 and above, against ng's own 0.886 and
0.903. What it is for is the **difference** between one setting and another over identical
evidence.

**The harness moves when a number that matters moves, which is what makes the flat rows below
mean something.** Replacing the shipped slip level of 0.10 with the 0.007 measured at period-2+
homozygous-reference tracts costs 137 calls in 8,686 (accuracy 0.8629 → 0.8471). So a row that
does not move is a number that does not matter, not a harness that cannot see.

| what was moved | calls right of 8,686 | accuracy | against the shipped setting |
|---|---|---|---|
| **shipped: λ = 0.01, ε = 0.001** | 7,495 | 0.86288 | — |
| outlier weight λ = 0.001 | 7,491 | 0.86242 | −4 calls |
| λ = 0.05 | 7,506 | 0.86415 | +11 calls |
| λ = 0.10 | 7,515 | 0.86519 | +20 calls |
| **λ = 0.20** | **7,523** | **0.86611** | **+28 calls** |
| λ = 0.40 | 7,508 | 0.86438 | +13 calls |
| substitution rate ε = 0.0005 | 7,495 | 0.86288 | 0 |
| ε = 0.00264 (measured, period 2+) | 7,497 | 0.86311 | +2 calls |
| ε = 0.005 | 7,498 | 0.86323 | +3 calls |
| ε = 0.05 | 7,501 | 0.86357 | +6 calls |
| part-repeat share 0.20 → 0.072 (measured) | 7,495 | 0.86288 | **0** |
| part-repeat share 0.20 → 0.05 (fitted route) | 7,495 | 0.86288 | **0** |
| part-repeat share 0.20 → 0.00 | 7,495 | 0.86288 | **0** |
| whole-repeat slip cutoff 10 → 3 | 7,496 | 0.86300 | +1 call |
| whole-repeat slip cutoff 10 → 20 | 7,495 | 0.86288 | 0 |
| one-step share 0.95 → 0.5 (a slippage number, for scale) | 7,495 | 0.86288 | 0 |
| λ = 0.20 **and** ε = 0.005 together | 7,521 | 0.86588 | +26 calls |

**The headroom, so the sizes above have something to be a share of.** Of the 1,191 tracts called
wrong at the shipped settings, **778 have both truth lengths among the candidates the row was
scored over** — 9 tracts in 100 of the 8,686. Those are the calls a better likelihood could
reach; the remaining 413 are lost before the likelihood sees them, in candidate generation. **So
every constant in this document put together reaches 28 of a possible 778.**

**What raising the outlier weight actually does, and it is not what the name says.** `λ · U` is a
floor under every emission, and `U` is one over the number of tract lengths the cutoffs reach —
median 22 at period 1 and 35 at period 2 and above. So:

| | the floor under every emission | the shipped stutter mass at 1 / 2 / 3 repeats shorter |
|---|---|---|
| λ = 0.01, period 1 | 4.6 × 10⁻⁴ | 4.75 × 10⁻² / 2.38 × 10⁻³ / 1.19 × 10⁻⁴ |
| λ = 0.20, period 1 | 9.1 × 10⁻³ | the same |

**At the shipped 0.01 the floor sits between a 2-repeat slip and a 3-repeat one, so a read three
or more repeats from a candidate already carries no information about the genotype. At 0.20 the
floor sits between a 1-repeat slip and a 2-repeat one.** Raising λ is therefore a statement that
*a read more than one repeat away tells you almost nothing*, and on this data that is a better
statement than the one the stutter model makes. It is the same device as freebayes'
read-dependence factor of 0.9 and GATK's Phred-45 cap on how much one read may discriminate
between alleles (spec §3.1) — a bound on any single read's pull — and ng has no other.

**It is not a patch for the shipped slippage numbers.** Re-run with the slippage numbers measured
in §2.2 instead of HipSTR's, the optimum is still at λ = 0.20: 7,516 right against 7,467 at
λ = 0.01, +49 calls. So the two lines do not collide — whoever fits slippage will not make this
finding go away.

**Where the 28 calls come from** (λ 0.01 → 0.20): the call changes at **96 of the 8,686 tracts,
52 wrong-to-right against 24 right-to-wrong**, the remaining 20 wrong either way. By depth: 43
gained and 20 lost between 10 and 29 reads, 7 gained and 4 lost below 10 reads, 2 gained and 0
lost at 30 reads or more. **So it works on tracts with middling depth and it is a 2-to-1 bet, not
a free win** — a fifth of the tracts it touches it makes worse.

---

## 4. The ranking

**1. `DEFAULT_OUTLIER_WEIGHT` (0.01) — the one to move.** It is the only constant in this document
that changes calls at all: +28 of 8,686 at λ = 0.20 (0.8629 → 0.8661), and the gain survives
replacing the shipped slippage numbers with measured ones. The mechanism is legible — the floor
`λ · U` bounds how far one read may pull a genotype, and at 0.01 it bounds it two repeats out
instead of one. Two cautions, both from the same measurement: read literally as *the share of
reads nothing explains* the number should be **0.0005 at homopolymers and 0.005 at period 2 and
above** (§2.3), the opposite direction; and 24 of the 96 tracts whose call moves get worse. **So
what the evidence supports is not "0.01 is wrong, use 0.2" but "this number is doing a job nobody
named it for, it is worth about 0.3 points of accuracy, and it should be swept and set
deliberately — per period, since the two periods want opposite corrections."** It is already
settable from the parameters file (`stated_constants.repeat_tract_outlier_weight`), so the sweep
needs no code change.

**2. The shape of the junk distribution `U` — unmeasured, and it is half of row 1.** The floor is
`λ · U`, and every conclusion about λ above is really about the product. `U` is uniform over 22 to
35 lengths by a decision that was never measured, so the floor moves by a factor of 1.6 between a
homopolymer and a period-4 tract for reasons that have nothing to do with either. Nothing here
measures it; measuring it means asking what lengths the reads that no candidate explains actually
show, which the candidate dump can answer and this run did not.

**3. `DEFAULT_SSR_SUBSTITUTION_RATE` (0.001) — measurably wrong, barely consequential.** The
measured rate inside a tract is 0.00132 at period 1 and 0.00264 at period 2 and above, so the
constant is optimistic by 1.3-fold and 2.6-fold; the bulk rate once the 1-in-215 junk reads are set
aside is 0.00168 at period 2 and above. Correcting it moves 2 to 3 calls in 8,686. **Its real cost
is in confidence rather than in calls**: the evidence one mismatching base gives for one allele
over another is `ln(3(1 − ε)/ε)`, which is 34.8 Phred at 0.001 and 30.6 Phred at 0.00264 — the
model is **4.2 Phred over-confident per mismatching base**, which lands on QUAL and GQ, not on the
genotype. Worth correcting because it is cheap and because the doc comment that flags it is
wrong about *why* (§2.1), not because it will move the accuracy number.

**4. The part-repeat shares — an inconsistency worth closing, a number not worth fitting.** The
two routes disagree by a factor of four about the same quantity (0.05 fitted against 0.20
defaulted, §1 row 6) and the measured value is 0.072. Setting it to the measured value, to the
fitted-route constant, or to **zero** changes **0 calls in 8,686**. At period 1 the branch cannot
be reached at all. **Recommendation: make the two routes agree, and do not commission the
per-period part-repeat estimator the specs keep deferring** — this measurement says it would buy
nothing at a genotype. What it might still buy is calibration, which was not measured here.

**5. `MAX_WHOLE_REPEAT_SLIP` / `MAX_PART_REPEAT_SLIP` (10 and 10) — settled, leave them.** They
discard 39 reads in 212,749 at period 1 and 16 in 120,420 at period 2 and above (§2.2), and moving
the whole-repeat cutoff anywhere from 3 to 20 moves at most 1 call in 8,686. The outlier floor
truncates the distribution long before the cutoff does, which is why: at λ = 0.01 the floor bites
at 3 repeats and the cutoff at 10. **These two constants are not worth a measurement; they are
worth a sentence in the spec saying that the outlier floor, not the cutoff, is what ends the
distribution.**

**6. `GEOM_MIN` / `GEOM_MAX` (0.01 / 0.99) and `MIN_SHARE_FROM_THIS_INDIVIDUAL` (10⁻¹²) — inert.**
Neither clamp binds at any value in play: the one-step shares are 0.95 and `1 − λ − c` is 0.99.
They exist for parameter rows that do not occur. No measurement is owed.

**7. The unreachable mass reported and never applied (row 13) — cannot change a genotype, can
change a QUAL.** At period 1 every candidate loses the same 2 parts in 100, so it cancels within a
locus; between candidates at one locus the difference is about 6 parts in a million at the shipped
rates. It does not cancel *between loci*, which is where the data likelihood is compared — so it
belongs to whoever owns emission and QUAL, not to genotyping.

**8. Equal weights over slip placements (row 12) — untouched here.** It only bites in interrupted
tracts and this harness compares lengths rather than sequences, so it could not see it. Recorded
as unmeasured rather than as unimportant.

---

## 5. What this does not say

- **One sample, 30×, one read group, human, short reads.** Nothing here was measured at 3×, at
  300×, or on a cohort. The outlier weight's gain is concentrated between 10 and 29 reads a
  tract, which is the middle of this benchmark's depth range and the whole of tomato's — so if it
  transfers anywhere it transfers there, but that is a prediction and not a measurement.
- **The sweep in §3 is a re-implementation, not ng.** It shares the formula and the constants and
  not the candidate generation, the censored-read term or the prior. Its differences are
  trustworthy; its absolute accuracy is 2.4 points below ng's and should not be quoted.
- **Calibration was not measured at all.** Every number in §3 is a hard call against truth.
  Several constants that move no call — the substitution rate above all — move confidence, and
  spec §4.1 records that calibration is what separated the two candidate models in the first
  place.

## 6. The scripts

All under `tmp/agent_readmodel/`, all runnable from the host with `uv run --no-project python`:

| script | what it does |
|---|---|
| `select_homref_tracts.py` | picks the 13,302 tracts with no truth record within 5 bases, and writes the region files the other three read |
| `tract_error_rate.py` | the per-base mismatch rate and mean base quality, inside the tract against the 50 bases either side, by period (§2.1) |
| `flank_profile.py` | the same mismatch rate by distance from the tract, which is what shows the flank rate is flat rather than an edge effect (§2.1) |
| `per_read_inside.py` | the per-read mismatch histogram inside a tract — the bulk-plus-tail shape (§2.1, §2.3) |
| `length_changes.py` | the five shares of the stutter distribution measured against truth, and the mass past each cutoff (§2.2) |
| `constant_sweep.py` | the genotype sweep of §3 |

The three that read the BAM take a headerless SAM stream on stdin:

```
cd tmp/agent_readmodel
samtools view -M -L homref_regions.bed \
  /Users/jose/devel/pop_var_caller/benchmarks/ssr_hg002/bam/30x/HG002_TR_v1.0.1_Tier_30x.bam \
  | uv run --no-project python tract_error_rate.py 0
```

Each script's output as run for this report is saved beside it as `*_output.txt`, so a claim
above can be checked without re-reading the BAM.
