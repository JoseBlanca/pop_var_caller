# Making ng's repeat-tract genotypes better: what worked, and what did not

**Date:** 2026-09-02. **Asked for by the owner** after
[the QUAL experiment](ng_tract_qual_experiment_2026-09-02.md) §5 reported the genotype accuracy.
**No change to the caller** — everything here is parameter values a run can already be given, and
measurements of what they are worth.

---

## 1. The answer

**Two parameter changes, no code, and they compose: +0.6 points of genotype accuracy at 30× and
+0.4 at 50×.** Fitting the slippage numbers from HG002's own reads instead of running on HipSTR's
shipped constants, together with raising the repeat-tract outlier weight from 0.01 to 0.10:

| | 30× homopolymer | 30× period 2+ | 50× homopolymer | 50× period 2+ |
|---|---:|---:|---:|---:|
| shipped (`--defaults`, what ships today) | 0.8856 | 0.9033 | 0.8950 | 0.9132 |
| fitted slippage + outlier weight 0.10 | **0.8907** | **0.9095** | **0.8997** | **0.9160** |

Scored the way STR callers conventionally are, on the two **repeat lengths** rather than the two
sequences — the number every caller here already writes as `REPCN`:

| | 30× homopolymer | 30× period 2+ | 50× homopolymer | 50× period 2+ |
|---|---:|---:|---:|---:|
| shipped | 0.9147 | 0.9127 | 0.9294 | 0.9237 |
| fitted slippage + outlier weight 0.10 | **0.9185** | **0.9182** | **0.9325** | **0.9249** |

**And the headline claim I made before this work is wrong.** I recommended fitting the slippage as
*the* lever on the strength of a simulator result. On real reads a *flat* change to any of the
three stutter numbers is worth nothing: the shipped values already sit at the optimum. What the
fit is worth comes entirely from its per-stratum shape, and it is a third of the gain. The larger
half comes from a constant nobody was looking at.

---

## 2. What was tried, and what each was worth

Every row is a full run of `call-from-alignments` over GIAB's HG002 tandem-repeat benchmark
(50,000 Tier intervals, 20,204 typed repeat tracts) at 30×, scored by
[`tract_qual_experiment.py`](../../../benchmarks/lib/tract_qual_experiment.py). Accuracy is over
the comparable tracts; 3 to 4 in 100 cannot be laid out and are excluded from both sides.

| setting | homopolymer | period 2+ | het called for a hom truth | hom called for a het truth |
|---|---:|---:|---:|---:|
| slip share 0.02 | 0.8805 | 0.9034 | 120 | 13 |
| slip share 0.05 | 0.8842 | 0.9048 | 98 | 20 |
| **slip share 0.10, the shipped value, supplied as rows** | **0.8851** | **0.9037** | **88** | **27** |
| slip share 0.15 | 0.8829 | 0.8976 | 82 | 42 |
| slip share 0.20 | 0.8781 | 0.8941 | 80 | 58 |
| slip share 0.30 | 0.8681 | 0.8850 | 71 | 96 |
| slip share 0.40 | 0.8384 | 0.8527 | 61 | 195 |
| shorter share 0.65 (shipped 0.50) | 0.8849 | 0.9035 | 87 | 30 |
| shorter share 0.80 | 0.8841 | 0.9031 | 89 | 31 |
| fall-off 0.15 (shipped 0.05) | 0.8857 | 0.9037 | 87 | 27 |
| fall-off 0.30 | 0.8854 | 0.9029 | 91 | 25 |
| slip rising with tract length, best of three shapes | 0.8859 | 0.9048 | 83 | 28 |
| **outlier weight 0.10 (shipped 0.01)** | **0.8892** | **0.9059** | **73** | **31** |
| outlier weight 0.20 | 0.8887 | 0.9051 | 70 | 37 |
| fitted per-stratum slippage | 0.8868 | 0.9053 | 90 | 21 |
| **fitted slippage + outlier weight 0.10** | **0.8907** | **0.9095** | **78** | **23** |

Three readings.

**The slippage share is a dial that trades one error for the other and does not move the total.**
Over a twenty-fold range it changes accuracy by half a point at most, while the two error classes
swing by a factor of nine and three in opposite directions — 120 spurious heterozygotes at 0.02
against 61 at 0.40, and 13 collapsed heterozygotes against 195. The shipped 0.10 sits at or beside
the peak. **The direction split and the fall-off do nothing at all**: every setting tried is within
0.1 points of the shipped one, which is a stronger statement than "they do not help".

**A per-stratum fit is worth about a sixth of a point on its own** (0.8868 against 0.8851), and
what it buys is a *different* error: it cuts collapsed heterozygotes from 27 to 21 and leaves
spurious ones roughly where they were, where the outlier weight does the opposite. That is why
the two compose rather than overlapping.

**The outlier weight is the larger half, and nobody was looking at it.** It is the share of reads
at a tract that came from somewhere the stutter model cannot explain — inherited from the existing
caller at 0.01 and never measured. `λ · U` is a floor under every read's emission, so the number
is really a bound on how far one read may pull a genotype, which is the job freebayes does with a
read-dependence factor and GATK with a Phred-45 cap, and which ng has nothing else doing. Raising
it from 0.01 to 0.10 removes 15 spurious heterozygotes of 88 and costs 4 collapsed ones.

---

## 3. What was ruled out, with the numbers that rule it out

**The genotype prior.** Its length spectrum is flat at `--defaults` and the truth is nothing like
flat — 79 chromosomes in 100 sit at the reference length for homopolymers against the 11 in 100 a
flat shape over nine offsets asserts, and it is strongly stratified (0.97 at a 6–8-base
homopolymer, 0.51 at 21+). But **fitting it reaches 10 of 648 errors and puts 77 correct calls at
risk**, because the dominant error is a truth homozygote called heterozygous — 130 of the 241
errors the prior could touch — and a prior peaked on the reference length makes exactly that
worse. Worth about +0.003 accuracy against seven correct calls risked for each one reached.

**A stricter candidate bar.** Every rung loses more true alleles than it removes spurious
heterozygotes: 2 reads and 15% removes 7 homopolymer errors and destroys 10 correct calls; 3 reads
and 10% removes 11 and destroys 50. The ladder saturates at 37 of 137 removed for 617 true alleles
lost.

**A GQ floor.** GQ 30 withdraws 43 of 86 wrong homopolymer calls and 413 of 3,113 correct ones —
ten good for one bad, and the bad become no-calls rather than right calls.

**An allele-balance collapse rule.** Best net anywhere on the curve is +2 tracts in 3,515; the
majority allele is the truth only 42 times in 86.

**Every other constant in the read model.** The part-repeat shares disagree by a factor of four
between the two ways of expressing the shipped model and setting them to anything, zero included,
changes no calls. The slip-size cutoffs change at most 1 call in 8,686 anywhere between 3 and 20 —
the outlier floor ends the distribution long before the cutoff does. The per-base substitution
rate inside a tract is measurably wrong (the reads give 1 mismatch in 756 at homopolymers against
the model's 1 in 3,300) and is worth 2 to 3 calls in 8,686.

---

## 4. The outside bar: HipSTR, which fits its stutter model per locus

On the 2,044 period-2-or-more tracts both callers reach:

| | comparable | right | accuracy |
|---|---:|---:|---:|
| ng 30× | 1,974 | 1,832 | 0.928 |
| HipSTR 30× | 1,874 | 1,733 | 0.925 |
| ng 50× | 1,973 | 1,852 | 0.939 |
| HipSTR 50× | 1,890 | 1,792 | 0.948 |

**A caller that fits its stutter model per locus is 0.3 points behind ng at 30× and 0.9 ahead at
50×** — 18 tracts of 1,973. And HipSTR's own fitted slip level has a median of 0.04 against ng's
fixed 0.05, so **the median locus does not need fitting**; only its top decile (0.07 up to 0.33)
differs. That is a third independent route to the same conclusion as §2's sweep.

**What this comparison does not control**, and it is a lot: different tract catalogs, different
candidate rules, different priors, `--min-reads 5`, and HipSTR emitting homozygous-reference
records where ng does not. HipSTR's region file holds **no period-1 loci at all**, so ng's
homopolymer numbers have no comparator here.

---

## 5. Where the accuracy actually is, and it is not in the model

At 30×, of 648 comparable genotype errors, **407 are a truth sequence ng never offered as a
candidate**. Every parameter in §2 and §3 together reaches perhaps 40 of the other 241. So the
model is not where the remaining accuracy lives.

Where the missing sequences go, over the same ground at 30× — 434 missing true sequences:

| | |
|---:|---|
| **268** | no read carried it, so the merge's allele table never held it |
| 61 | it cleared the support bar; the per-sample top-`ploidy` cut dropped it |
| 59 | the merge refused the tract, so no locus was built |
| 46 | the merge's table held it and the support bar refused it |

And the 268 split further:

| | |
|---:|---|
| 121 | absent from the reads even at 300× with a median 115 spanning reads — not recoverable |
| **67** | **an alignment loss: reads carry the sequence and ng's table does not** (46 unambiguous at 30×) |
| 66 | every base the truth needs is in the reads, but no read spells the tract that way |
| 14 | the allele is longer than a 150 bp read can span with 20 bp of flank |

**The 46 unambiguous alignment losses are the finding nobody suspected**, and three are verifiable
by hand: at `chr3:33,877,690` an 11 bp poly-A where 10 of 23 reads carry `CAAAAAAAAAA` and ng's
table holds only the bare A-run — the leading `C` is dropped and the tract is reported one repeat
short; the same mechanism at `chr3:37,126,860` (10 of 11 reads at 30×, 120 of 138 at 300×) and at
`chr11:37,147,255`, where the only 17-base allele ng holds is a one-read sequencing-error spelling
rather than the true one carried by 12 reads of 14.

**And a number that may make three quarters of this moot: 174 of the 268 are the reference
length** — the truth record inside the tract is a substitution, not a repeat-count change — and
195 of 268 have the missing allele's length already in ng's table. Scored on repeat length, which
is what every caller here emits and how the field scores STRs, most of this bucket is not a
bucket. §1's second table is that measure.

---

## 6. What to do

0. **Done, 2026-09-03 (owner's decision): the shipped `DEFAULT_OUTLIER_WEIGHT` is 0.05**, taken
   at the conservative end of the plateau rather than the 0.10 that scored best, because 0.05
   moves least far from what the reads themselves suggest while taking about three quarters of the
   gain. Measured end to end through `--defaults` rather than through a parameters file, at 30× on
   this benchmark: homopolymer genotype accuracy **0.8856 → 0.8881**, period 2+ **0.9033 →
   0.9059**; scored on repeat length, 0.9147 → 0.9166 and 0.9127 → 0.9150. **Its warrant stays
   `Defaulted`** — it is a stated constant chosen by a sweep, not an estimate of the share it is
   named for, and `likelihood/ssr.rs` says so at the constant. **What it still owes**: a sweep per
   motif period, and a check at three reads a position on the tomato panel, neither of which this
   benchmark can give.

1. **Adopt the two parameter changes** — but note that the fitted slippage came from HG002's own
   reads, so it is a fit for this sample and this chemistry, not a shipped default. What it argues
   for is the **fit-mode command** (spec §3.4, deferred): the machinery to produce a parameters
   file per run already reads back correctly, and this measures what it is worth. The outlier
   weight is different — 0.10 beat 0.01 on both period classes at both depths, so it is a
   candidate for the shipped constant, and it should be swept per period before it is set.
2. **Look at the tract aligner, not the model.** 46 tracts of 20,204 lose an allele the reads
   plainly carry, with a legible mechanism — a base at the tract's edge dropped, an interruption
   inside it discarded. That is 2 in 1,000 tracts, and it is the only class in §5 that is both
   large and clearly a defect rather than a limit of the data.
3. **Settle whether a tract genotype is scored on its sequence or its repeat length**, because it
   changes the headline by 3 points and it changes which of §5's buckets is worth anything.
4. **Do not build allele discovery for this.** It is aimed at the 61 of 434, against 268 that are
   not there to be found and 46 that are lost upstream of selection.

---

## 7. How this was measured, and what it does not cover

Every number in §2 is a full run of the shipped binary at that parameter setting, scored by one
instrument. **Two controls, and they say different things.** A run given the `--defaults` run's own
parameters file unedited is **byte-identical** to it, which is what says the harness changes
nothing by itself. A run given the *shipped numbers written out as slippage rows* — §2's `slip
share 0.10` row — differs from `--defaults` by 2 genotypes in 3,648 (0.8851 against 0.8856 at
homopolymers, 0.9037 against 0.9033 at period 2+): a supplied row rebuilds the part-repeat shares
as a twentieth of the whole-repeat mass where the shipped model states them as 0.01 each. So §2's
rows are comparable with each other, and the ±0.05-point offset against §1's baseline is that
difference and not a result.

§3's and §5's numbers come from six independent investigations run in parallel, whose raw reports
are kept in [`tract_genotype_investigation/`](tract_genotype_investigation/). **Their figures are
not independently re-derived here** except where a run in §2 corroborates them: the outlier
weight's value, the fitted slippage's value, and the direction of the prior's effect all were.

**One sample, one chemistry, two depths.** HG002 at 30× and 50×, and nothing here says what any of
it does on a 63-accession cohort at three reads a position, where the outlier weight's floor and
the candidate bar both bind differently.
