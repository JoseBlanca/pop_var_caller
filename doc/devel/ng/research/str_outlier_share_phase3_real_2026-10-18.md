# Outlier share, phase 3: the truth-free fit on HG002's real reads

**Date:** 2026-10-18. **Plan:** [`str_outlier_share_plan.md`](str_outlier_share_plan.md) §4,
phase 3, run on real data before any calling test. **Sample:** HG002, GIAB tandem-repeat
benchmark regions, 5x / 10x / 30x, every tract's complete reads from the psp (phase 1's dumps).
**Code and output:** `experiments/str_outlier_share/real_fit.py`,
`experiments/str_outlier_share/phase3_real_fit_{5x,10x,30x}.tsv`.

## What was run

Phase 2's prototype fit, unchanged, on each stratum — one repeat-unit length and one reference
repeat count — holding at least 100 tracts (42 to 44 strata a depth, up to 3,000 tracts each). It
fits the slip rate, its direction and fall-off, the stray share and the allele-length shares,
summing over each tract's unknown genotype. It sees only reads a whole number of repeats from the
reference, as today's stutter fit does; reads a non-whole number off and reads with a flank deletion
are counted beside it. No truth is used.

## What was found

**Pooled over a period class, the share lands where phase 1's truth-based count put it, at every
depth.** Read-weighted over the class's strata:

| | 5x | 10x | 30x | phase 1, against truth (30x) |
|---|---|---|---|---|
| period 1 | 0.0031 | 0.0039 | 0.0027 | 0.0021 two or more repeats off, at clean tracts |
| period 2 | 0.0045 | 0.0044 | 0.0032 | 0.0006 |
| period 3–6 | 0.0002 | 0.0004 | 0.0002 | 0.0016 |

The phase 1 column counts reads two or more repeats off at tracts it classed as clean; the fit's
share is spread uniformly over 25 lengths, so the two are related but not the same quantity. The
fitted value sits above the count at periods 1 and 2 and below it at 3–6. The shortfall at 3–6 is
not explained here: reads three or more repeats off are 0.0013 of reads at clean tracts there, and
the fit sees them.

**Per stratum, the share is not one number.** At homopolymers it rises steadily with tract length
(30x): 0.0005 at 8 repeats, 0.0019 at 10–13, 0.004–0.008 at 16–20, 0.012–0.015 at 21–26. The
fall-off fitted alongside it is small, 0.04–0.25, so stutter's chance of slipping three or more
repeats is small, and the share is probably taking the long tail of stutter at long tracts rather
than reads from elsewhere. At dinucleotides it swings between neighbouring strata — 0.0000, 0.0202,
0.0017, 0.0107, 0.0000 from 10 to 14 repeats — so with 150–250 tracts a stratum the share and the
stutter tail are not separated.

**The reads the fit does not see are larger than the share it fits.** As shares of all reads at the
class's tracts:

| | non-whole number of repeats off | flank deletion |
|---|---|---|
| period 1 | — | 0.010–0.014 |
| period 2 | 0.0018–0.0022 | 0.0038–0.0062 |
| period 3–6 | 0.0039–0.0050 | 0.0018–0.0041 |

Flank-deletion reads include real alleles: since the held edge fix, a deletion allele that touches
the tract's edge is spelled with one. Phase 1 counted 0.0024 of reads at clean homozygous
homopolymer tracts carrying one; the 1–1.4% here is all tracts. At calling, a read whose flank
deletion differs from its allele's gets no stutter probability and falls entirely to the outlier
term, so a share fitted for calling has to account for these reads too.

## What it means

- **A fitted share near 2–4 reads in 1,000 at periods 1 and 2 is what the reads say**, at each of
  the three depths, and it agrees in size with the truth-based count.
- **Fitting it as one number per period class would be wrong at long homopolymers**, where it
  climbs to 15 in 1,000 together with a steep fitted stutter tail. Either the share is fitted per
  stratum and smoothed along repeat count, as the slip rate already is, or the stutter shape needs a
  heavier tail and the share stays flat. Which of the two is right is the next question.
- **The fit has to see the reads calling sends to the outlier term**: non-whole lengths and
  mismatched flank deletions.

## What this does not show

- Whether calling with these values beats `--defaults` at 0.07 (phase 4).
- Whether the rise with tract length is stutter's tail or real stray reads; the truth-based counts
  of phase 1 were not split by repeat count.
- One sample only; the cohort case is simulated (phase 2), not measured.

## Addendum: a smoothed per-stratum share, a heavier stutter tail, and calling with them

Same data; scripts `smooth_fit.py`, `tail_fit.py`, `far_placement.py`, outputs `phase3_*` in
`experiments/str_outlier_share/`.

**Where the far reads sit, against truth.** At clean truth-homozygous homopolymers (30x), reads two
repeats off are about a tenth as common as one repeat off, but three repeats off are half as common
as two (at 18–22 repeats: 20.9 against 35.3 per 10,000), which one geometric fall-off cannot give.
Reads six or more repeats off, which no stutter reaches, rise with tract length: 2.7, 6.7, 23.5 and
40.7 per 10,000 at 8–12, 13–17, 18–22 and 23 or more repeats.

**A share smoothed along repeat count beats both alternatives on held-out tracts.** Fitted on even
tracts and scored on odd ones, a second-difference penalty on the share's logit beat independent
strata and one share per period in all four cells (10x and 30x, periods 1 and 2). The fit here also
counts the reads calling sends to the outlier term — non-whole lengths, and flank-deletion spellings
under 20% of their tract's reads — at the outlier term's probability. With them, homopolymer shares
run from 0.003 at 8 repeats to 0.05–0.10 at 22–26.

**Calling with the fitted per-stratum shares ties `--defaults`.** HG002, the fitted model from the
30x psps, the held edge fix in every arm, correct tract genotypes over common tracts:

| | 5x | 10x | 30x |
|---|---|---|---|
| `--defaults`, 0.07 | 2,873 | 4,709 | 7,225 |
| fitted, per-stratum shares | 2,873 | 4,709 | 7,210 |
| fitted, one share 0.003 | 2,876 | 4,710 | 7,201 |
| fitted, one share 0.03 | 2,874 | 4,715 | 7,228 |

At 30x the fitted model's losses against `--defaults` are 30 homopolymer heterozygotes against 9
gained; of the 49 heterozygotes where the reads favour the heterozygote and the prior overrules
them, 22 have two non-reference alleles, the genotype class the fitted prior under-predicts about
three-fold (`str_outlier_share_phase1_2026-10-18.md`).

**A heavier stutter tail does not earn its place.** A mixture of two geometric fall-offs against
one, with the smoothed share, held out on odd tracts (negative favours the mixture):

| | 10x homopolymers | 30x homopolymers | 10x dinucleotides | 30x dinucleotides |
|---|---|---|---|---|
| with outlier-only reads | +1.2 | −2.9 | −22.7 | +27.7 |
| without them | −0.8 | +9.4 | — | — |

With the outlier-only reads the mixture's slow share fits to about zero everywhere: those reads can
only be explained by the stray share, which they hold high enough to absorb the tail. Without them
the mixture does take a slow tail (20–50% of slips at about half the homopolymer strata) and lowers
the stray share by 10–30%, but the share still climbs about seven-fold from short tracts to long
ones (30x: 1.5–2.3 in 1,000 at 9–15 repeats, 10–15 at 21–26), and the held-out score does not
improve. The rise is reads from outside stutter's reach, as the truth count above shows.

**What follows.** The stray share is a real, length-dependent quantity worth fitting per stratum
and smoothed; a heavier stutter tail is not needed. The largest input to it at homopolymers is the
minority flank-deletion spellings (about 1 in 100 reads), and whether those are reads from
elsewhere or a spelling artefact of the mapper at the tract's edge is not yet known.

## Addendum 2: the minority flank-deletion reads are not reads from elsewhere

Script `flank_reads.py`. At homopolymer tracts (30x, psps generated with the held edge fix), reads
spelled with a flank deletion whose spelling carries under 20% of its tract's reads — what calling
sends to the outlier term — number 1,967. For each, its total length change (the tract's change
minus the flank bases it deletes) was compared with the tract's true allele lengths:

| truth at the tract | total length is a true allele's | total length is no true allele's |
|---|---|---|
| homozygous reference | — | 684 (32 in 10,000 of those tracts' reads) |
| heterozygous | 558 | 170 |
| homozygous deletion | 461 | 81 |
| homozygous insertion | — | 13 |

In every one the deleted flank base differs from the repeat's base (907 on the left, 1,060 on the
right). 10x gives the same picture at a third of the counts.

- **About half are a true allele spelled at the tract's edge**: the deletion is real, and the mapper
  placed it on the flank base — the mechanism of the spelling defect behind the held edge fix.
- **The other half are one-base length errors beside the tract**, mostly at homozygous-reference
  tracts, where the tract reads at full length and the neighbouring base is lost. They are errors of
  stutter's size placed outside the tract, not reads from another locus.

**So they should not raise the stray share.** Both kinds are a question of how a read is compared
with an allele: scored by its total length, each would be a match or a one-step slip. That belongs
to the edge fix's spec amendment, not to the outlier share; the share fitted without them is the
one to carry forward (30x homopolymers: 1.5–2.3 in 1,000 at 9–15 repeats, 10–15 in 1,000 at 21–26).

## Correction (2026-10-18, same day): no read here is "outlier-only"

The two addenda above, and `smooth_fit.py`, treat reads a non-whole number of repeats off and
minority flank-deletion spellings as reads calling sends wholly to the outlier term. **That is
wrong for both.** Calling compares a read and an allele with different flank deletions over the
wider of the two extents (`str_edge_deletion_alleles.md` §5.6; `fill_emissions` in
`calling/likelihood/ssr.rs`): same length after respelling is priced at the substitution rate, a
different length as an indel that is not whole repeats. Non-whole lengths are priced by the stutter
model's part-repeat regime (`read_likelihoods.md` §4.2). So:

- **the share to carry forward is the one fitted without those reads** (addendum 1's "without
  them" row; 30x homopolymers 1.5–2.3 in 1,000 at 9–15 repeats, 10–15 in 1,000 at 21–26);
- **the calling comparison in addendum 1 used per-stratum shares fitted with them**, which run to
  5–10% at long homopolymers, and has to be repeated with shares fitted without them;
- addendum 2's split of the flank-deletion reads stands — half a true allele spelled at the edge,
  half a one-base error beside the tract — and under §5.6 both are already scored against the
  alleles rather than sent to the outlier term. What decides how they score is whether the read
  records its flank deletion at all, which is the held edge fix (a deletion that *ends at* the
  tract's edge, not only one that crosses it).
