# What the flat tract prior costs on HG002, and what the true length spectrum looks like

ng calls a repeat tract's genotype from a stutter likelihood times a prior. With
`--defaults` — which is every benchmark run — the prior's shape is the tract
ladder's bottom rung: **flat over the candidate lengths, at concentration 1.0**
(`HG002_30x.raw.parameters.toml`: `length_spectrum_by_stratum = []`,
`length_spectrum_by_period = []`,
`fallback_length_spectrum_concentration = { value = 1.0, warrant = "defaulted" }`).
Nobody had measured what that costs. This is the measurement.

Two halves. **A** asks what HG002's tract lengths actually do — the empirical
length spectrum a fit would produce. **B** asks whether ng's genotype errors are
the kind a non-flat prior would fix.

## How it was measured

Both halves use `benchmarks/lib/tract_qual_experiment.py`'s own comparison —
`prepared_vcf`, `TractGround`, `tract_reference_bases`, `records_by_tract`,
`haplotypes_over_tract` — imported rather than rewritten, because that
comparison is subtle and three earlier versions of it were wrong. Both sides are
left-aligned by `bcftools norm` against GRCh38 and cut to the confident regions
(`tmp/tract_qual/ground/tier_sorted.bed`) before anything is compared; the tract
ground is `tmp/tract_qual/ground/tier.bed`, 20,204 tracts with their motif
period.

**A tract's genotype is compared as the tract's own two sequences**, rebuilt on
the reference window, exactly as the scorer does it. **An allele's offset is
that sequence's length difference from the reference tract, divided by the motif
period** — so `-1` means one repeat unit shorter than the reference carries.
Measuring offsets from the rebuilt sequences rather than from ng's `REPCN` keeps
the two sides on one measure; `REPCN` exists only on ng's side and the truth set
has no equivalent.

Scripts, all in this directory and all run with `uv run --no-project python`:

| script | what it does |
|---|---|
| `spectrum_and_errors.py` | the two measurements; writes `spectrum_<depth>.json` and `errors_<depth>.json` |
| `report_tables.py` | half A's tables, and half B's called-against-true offset grid |
| `report_extra.py` | candidate-set width; the fitted prior scored per stratum; the tracts the genotype denominator cannot see |
| `report_control.py` | the counterweight — what the fitted prior would do to the calls ng gets right |
| `report_levers.py` | the prior's other two knobs, the concentration and the inbreeding coefficient |

Captured output: `tables_30x.txt`, `tables_50x.txt`, `extra_*.txt`,
`control_*.txt`, `levers_*.txt`.

---

## A. The true length spectrum, and it is nothing like flat

**Denominator: 18,046 tracts** — every tract of the ground lying wholly inside
the confident regions, so **36,092 chromosomes**. 1,956 tracts were dropped for
straddling a confident-region edge, 188 for a truth genotype the sequence
rebuild refuses, 14 for a truth no-call. Where the truth set writes no record at
a tract, both chromosomes are counted at offset 0 — that is what a hom-ref tract
means, and leaving those out is what would make the spectrum look flat.

### Chromosomes by offset in whole repeats

| class | −4 | −3 | −2 | −1 | **0** | +1 | +2 | +3 | +4 | \|off\|>4 | partial | total |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| homopolymer | 57 | 111 | 358 | 1,504 | **16,922** | 1,562 | 360 | 144 | 68 | 288 | 0 | 21,374 |
| period 2+ | 132 | 182 | 421 | 876 | **10,931** | 837 | 396 | 194 | 120 | 536 | 93 | 14,718 |
| period 2 | 105 | 136 | 297 | 554 | 7,478 | 534 | 261 | 146 | 90 | 475 | 52 | 10,128 |
| period 3 | 8 | 13 | 27 | 82 | 1,514 | 73 | 37 | 18 | 11 | 25 | 6 | 1,814 |
| period 4 | 19 | 30 | 91 | 194 | 1,482 | 200 | 86 | 27 | 17 | 33 | 27 | 2,206 |
| period 5 | 0 | 2 | 4 | 34 | 350 | 20 | 11 | 3 | 2 | 2 | 4 | 432 |
| period 6 | 0 | 1 | 2 | 12 | 107 | 10 | 1 | 0 | 0 | 1 | 4 | 138 |

`partial` is a chromosome whose length change is not a whole number of repeat
units — 93 chromosomes in 36,092, about 3 in 1,000, so nothing here turns on
them.

**How far from flat.** Flat over the nine offsets −4..+4 puts 11 chromosomes in
100 at each. The truth puts **79 in 100 at offset 0 for homopolymers** (16,922
of 21,374) and **74 in 100 for period 2 and above** (10,931 of 14,718) — seven
times what flat says. Against its own nearest neighbour the reference length is
**11 times more likely for homopolymers** (16,922 against 1,504 at offset −1)
and **12.5 times for period 2+** (10,931 against 876).

**But "flat" in ng is not flat over nine offsets.** The prior is flat over
*this locus's candidate lengths*, and at the 6,058 tracts scored in half B the
candidate set holds:

| candidate lengths offered | tracts |
|---:|---:|
| 1 | 458 |
| 2 | 4,460 |
| 3 | 1,138 |
| 4 | 2 |

So flat usually means **50/50 between two lengths**. Renormalise the homopolymer
spectrum onto a two-candidate set `{0, −1}` and it says **92/8**. That is the
real gap the bottom rung leaves: at a two-candidate tract the fitted shape would
move the odds by about **11:1, roughly 10 phred**, in favour of the reference
length.

### The spectrum is strongly stratified by reference repeat count

| period class | reference repeats | p(0) | p(−1) | p(+1) | chromosomes |
|---|---|---:|---:|---:|---:|
| homopolymer | 6–8 | **0.971** | 0.009 | 0.013 | 3,526 |
| homopolymer | 9–12 | 0.870 | 0.050 | 0.053 | 8,730 |
| homopolymer | 13–20 | 0.682 | 0.107 | 0.114 | 7,300 |
| homopolymer | 21+ | **0.510** | 0.140 | 0.125 | 1,818 |
| period 2+ | ≤5 | 0.888 | 0.056 | 0.028 | 392 |
| period 2+ | 6–8 | **0.918** | 0.025 | 0.027 | 7,766 |
| period 2+ | 9–12 | 0.676 | 0.092 | 0.081 | 2,992 |
| period 2+ | 13–20 | 0.454 | 0.112 | 0.108 | 2,556 |
| period 2+ | 21+ | **0.267** | 0.098 | 0.097 | 1,012 |

A short tract is nearly monomorphic — 97 chromosomes in 100 at the reference
length for a 6–8-base homopolymer — and a long one is barely peaked at all: at a
period-2+ tract of 21 repeats or more, only 27 chromosomes in 100 sit at the
reference length. **One pooled spectrum per period would be wrong at both ends.**
This is §4.4's "what it gives up: the repeat-count trend within a period",
measured: it is the difference between 0.97 and 0.27.

### Restricted to tracts the truth calls a variant, the peak nearly vanishes

| class | −4 | −3 | −2 | −1 | 0 | +1 | +2 | +3 | +4 | total |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| homopolymer | 57 | 111 | 358 | 1,504 | 2,222 | 1,562 | 360 | 144 | 68 | 6,674 |
| period 2+ | 132 | 182 | 421 | 876 | 1,549 | 837 | 396 | 194 | 120 | 5,243 |

At a tract where the truth carries any record, only **1 chromosome in 3** is at
the reference length (2,222 of 6,674 for homopolymers). **This is the population
half B lives in**, because the genotype-accuracy denominator is exactly "tracts
the truth calls". The marginal spectrum's 79-in-100 peak is real, and most of it
is contributed by tracts where there is no genotype decision to get wrong.

---

## B. What it costs: the errors are not the kind a length prior fixes

Scored on `HG002_30x.raw.vcf` against the same ground: **6,058 tracts both sides
call, 5,410 genotypes right, 648 wrong** — 0.893 overall, which is the standing
0.886 / 0.903 split. (50x: 6,112 both call, 5,517 right, **595 wrong**.)

**Of the 648 errors, 407 are not the genotyper's to fix**: a sequence the truth
carries was never among the alleles ng's records offered, so no prior over that
set could have been right. That is candidate selection's. **The prior's whole
reachable pool is the remaining 241.**

### The errors do not move away from the reference length

For each error, the distance of the called genotype from the reference length
(the two called offsets' absolute values, summed) against the same for the truth:

| | n | call further from reference | call closer | same distance |
|---|---:|---:|---:|---:|
| errors, homopolymer | 402 | 114 | 176 | 112 |
| errors, period 2+ | 246 | 119 | 99 | 28 |
| **errors, all** | **648** | **233** | **275** | **140** |
| errors, truth's alleles were offered | 241 | 90 | 97 | 54 |
| correct calls (control) | 5,410 | 0 | 0 | 5,410 |

**They are symmetric.** 233 of the wrong calls sit further from the reference
length than the truth does and 275 sit closer — if anything ng errs slightly
*towards* the reference, not away from it. A prior that favours the reference
length would fix a few of the first group and break an equal number of the
second.

The called-against-true offset grid for the 241 reachable errors says the same
thing (rows are the truth's offset, columns ng's, counted per chromosome):

| true \ called | −4 | −3 | −2 | −1 | 0 | +1 | +2 | +3 | +4 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| −3 | 0 | 9 | 3 | 0 | 2 | 0 | 0 | 0 | 0 |
| −2 | 0 | 2 | 18 | 6 | 2 | 1 | 0 | 1 | 0 |
| −1 | 1 | 1 | 10 | 56 | **20** | 1 | 0 | 0 | 0 |
| 0 | 1 | 2 | 4 | **17** | 67 | **29** | 3 | 1 | 1 |
| +1 | 0 | 0 | 0 | 0 | **28** | 82 | 6 | 1 | 0 |
| +2 | 0 | 0 | 0 | 0 | 2 | 10 | 24 | 0 | 0 |
| +3 | 0 | 0 | 1 | 0 | 0 | 2 | 4 | 10 | 1 |

Read the four bold cells. **46 chromosomes whose truth is the reference length
were called one repeat off** (17 + 29) — a reference-favouring prior would fix
those. **48 chromosomes that are really one repeat off were called at the
reference length** (20 + 28) — the same prior makes those worse. The two are the
same size.

### Scoring the fitted prior directly

Take half A's per-stratum spectrum, seed it as ng would (`α_j = concentration ·
w_j`, `w` renormalised over the locus's own candidates), and compute for each
error the phred by which it lifts the truth's genotype over the one ng called,
relative to what flat gave. Compare that swing against the GQ ng reported — the
margin to its own second-best genotype, so a swing smaller than GQ cannot flip
the call. Run backwards on the correct calls: a correct call is at risk when the
same prior lifts some rival genotype (one repeat away on either copy, within the
locus's candidates) by more than its GQ.

| concentration | errors the fitted prior reaches | errors it pushes further from the truth | correct calls put at risk |
|---:|---:|---:|---:|
| 1.0 | **10** | 127 | **77** |
| 2.0 | 11 | 127 | 90 |
| 5.0 | 13 | 137 | 99 |

At 50x: 9 reached, 135 pushed further, 35 correct calls at risk.

**The fitted spectrum is net negative here — it risks about seven calls for
every one it could reach.** Both numbers are upper bounds on their side (GQ is
the margin to the best rival, not specifically to the truth or to the rival
tested), which does not change the ratio.

### Why: the dominant error runs against the reference-length peak

Of the 241 reachable errors, split by the hom/het shape of truth against call:

| | 30x | 50x |
|---|---:|---:|
| truth homozygous, ng called heterozygous | **130** | 123 |
| truth heterozygous, ng called homozygous | 46 | 23 |
| both homozygous, wrong length | 31 | 36 |
| both heterozygous, wrong lengths | 34 | 48 |

**ng calls a heterozygote where the truth is homozygous 2.8 times as often as
the reverse (130 against 46).** And of those 130, **55 are a truth homozygous at
a non-reference length that ng called heterozygous *including* the reference
length**, against 26 where the truth was homozygous at the reference length.

A prior peaked on the reference length makes that dominant case worse, not
better. Worked, at a homopolymer of 13–20 repeats with two candidates `{−1, 0}`,
concentration 1, truth homozygous at −1 and ng calling het(0, −1):

```
flat    α = {−1: 0.500, 0: 0.500}   log P(hom −1) −0.288   log P(het) −0.693   hom favoured by +1.8 phred
fitted  α = {−1: 0.136, 0: 0.864}   log P(hom −1) −1.868   log P(het) −1.448   hom favoured by −1.8 phred
```

Flat already prefers the homozygote by 1.8 phred; the fitted spectrum reverses
that and prefers the heterozygote by 1.8 phred. **The flat prior is accidentally
better aimed at ng's commonest tract error than the true spectrum is**, because
what that error needs is more weight on homozygotes, and the fitted spectrum
buys its sharpness at the reference length by making the off-reference
homozygote cheap.

### The prior's other two knobs, on the same instrument

The hom/het balance is set by the concentration and by the inbreeding
coefficient (`inbreeding_coefficient = { value = 0.0, warrant = "defaulted" }`),
not by the shape. Both were swept, flat shape held:

| concentration (F=0) | errors reached | pushed worse | unchanged | swing too small | correct calls at risk |
|---:|---:|---:|---:|---:|---:|
| 0.05 | 18 | 46 | 65 | 112 | 105 |
| 0.1 | 12 | 46 | 65 | 118 | 69 |
| 0.2 | 9 | 46 | 65 | 121 | 35 |
| 0.5 | 0 | 46 | 65 | 130 | 0 |
| **1.0 (as shipped)** | 0 | 0 | 241 | 0 | 0 |
| 2.0 | 0 | 130 | 65 | 46 | 0 |

| inbreeding F (concentration 1) | errors reached | pushed worse | unchanged | swing too small | correct calls at risk |
|---:|---:|---:|---:|---:|---:|
| 0.0 (as shipped) | 0 | 0 | 241 | 0 | 0 |
| 0.25 | 0 | 46 | 65 | 130 | 0 |
| 0.5 | 1 | 46 | 65 | 129 | 6 |
| 0.75 | 10 | 46 | 65 | 120 | 46 |
| 0.9 | 18 | 46 | 65 | 112 | 105 |

Both point the same way and both hit the same wall: every setting that reaches
18 errors puts 105 correct calls at risk, **about six risked for one reached**.
The reason is in the "swing too small" column — 112 to 130 of the 241 errors
move in the right direction and not far enough. **ng's tract errors are not
close decisions the prior loses by a hair; they are decisions the likelihood
gets wrong by more than any prior of this shape can pay for.**

### One thing this denominator cannot see

At **144 tracts** inside the confident regions the truth carries no record at
all and ng called a variant anyway (100 homopolymer, 44 period 2+; 149 at 50x).
Those are false positives, not genotype errors, so they are outside the 648
entirely. A reference-favouring prior would act on them — but 144 is small
beside 648, and this measurement does not say how many of them are close enough
to flip.

---

## Verdict

**Fitting the length spectrum is not the lever, and it should not be the next
thing built for genotype accuracy.**

Three findings, in the order that decides it:

1. **The true spectrum is far from flat, and the ladder's bottom rung really is
   losing information.** 79 chromosomes in 100 sit at the reference length for
   homopolymers, against the 11 in 100 flat-over-nine-offsets asserts, and the
   peak runs from 0.97 at a short tract to 0.27 at a long one. Half A is
   unambiguous.

2. **None of that information helps at the tracts where genotypes go wrong.**
   The genotype denominator is tracts the truth calls variant, and there only 1
   chromosome in 3 sits at the reference length. ng's errors are symmetric about
   the reference length — 233 further out, 275 closer in, of 648.

3. **Dropping the fitted spectrum in would lose more than it gains.** 10 of the
   648 errors reachable at 30x, against 77 correct calls put at risk; the
   dominant error class (truth homozygous, called heterozygous, 130 of the 241)
   is one the sharper prior actively worsens.

**How many of the 648 could a prior plausibly reach? About 10, and no more than
about 18** — that is 1 error in 65 at best, worth roughly +0.003 on a genotype
accuracy of 0.893, and only if the 77 to 105 correct calls it endangers all
happened to survive. The honest ceiling on the whole prior, shape and knobs
together, is around **1 error in 65 reached against 6 to 7 risked.**

**Where the accuracy actually is:**

- **407 of the 648 errors are candidate selection's** — the truth's sequence was
  never offered, so 63 errors in 100 are decided before the prior is consulted.
  That is the biggest single block and it is 4 times the whole prior's reach.
- **Of the 241 that are the genotyper's, 130 are a homozygote called
  heterozygous.** That asymmetry — 2.8 to 1 against its reverse — is the stutter
  likelihood's shape, not the prior's: it says stutter reads from a homozygous
  tract are being read as evidence for a second allele. 112 to 130 of these move
  the right way under every prior tried and none of them move far enough.

**What the spectrum measurement is still worth.** Half A is a real result and
should be kept: it says the tract ladder's second rung — one pooled spectrum per
motif period — would be wrong at both ends of the repeat-count range (0.97
against 0.27), which is §4.4's stated open question 2 answered with numbers. If
the spectrum is ever fitted, it must be per stratum and not per period. It is
just not the thing to fit next.

Everything here is one sample, one benchmark, GIAB HG002 at 30x and 50x on
GRCh38 tandem repeats. A cohort changes half A's arithmetic — the spectrum
becomes a population's rather than one diploid's — and nothing measured here
says what the flat rung costs at 63 accessions at 3×.
