# The SNP/indel path charges an indel error at the substitution rate

**Date:** 2026-10-21. Branch `generic-indel-error`, from main bbcc9f46. Prompted by a colleague's
GIAB chr22 report on pvc bbcc9f46 (`tmp/inform_to_pop_var/README.md` in the main checkout,
sections 2, 3, 8 and 9).

**Status (2026-10-21, later the same day): the owner ruled for §7's amendment.** It is in
`read_likelihoods.md` §3.5 and built (commit f800df67): an indel error is charged without the
calibration scale. The prototype switch this report describes, `PVC_EXP_INDEL`, is gone; the
built rule gives calls byte-identical to its `nomult` arm on all 15 GIAB configurations and on
tomato, so every `nomult` figure below is the shipped behaviour's. The `table`, `max` and `floor`
arms existed only in commit 18dd5712.

**One difference from §7's proposed text: partial reads keep the scale** (spec §5.3). Deciding that
a fragment is one indel from an allele would need an alignment the partial rule deliberately avoids,
and the prototype — which never touched partials — already restored every lost indel. The built
rule also states that it holds at any scale, below one as well as above.

**One follow-up from the code review, not done here.** At a locus holding only the reference and
one indel, with no partial reads, the scale now enters no genotype's likelihood, yet the call's
reported warrant (`weakest_warrant_at_the_locus`, `summarise_condition.rs`) still folds in the
calibration's provenance, so it can name a weaker warrant than the call rests on. It is a label, not
a number; whether to narrow the fold is a small decision for later.

The colleague measured two defects on the SNP/indel path (the *generic path*: every locus that is
not a catalogued repeat tract):

- **Defect A — fitting the parameters loses indels.** Calling with parameters fitted on the
  sample's own data, instead of `--defaults`, cuts false SNPs at every depth but misses more
  indels, on generic ground only.
- **Defect B — `0/1` at homozygous indels.** At ten 300× sites where 80 to 88 reads in 100 that
  cross the indel carry it, pvc writes `0/1` and DeepVariant `1/1`.

The working hypothesis was that both come from one gap: the generic path has no indel error rate
of its own. **That is right for Defect A and wrong for Defect B.**

## Summary

**Defect A: the fit makes indel errors several times too cheap, and one number does it.** The fit
sets a per-read-group multiplier on every read's error probability — 6.2 to 6.4 on HG002 at 10×
and 30×, 14 at 300× — so that the average charge matches the rate at which reads disagree with the
truth, a rate set almost wholly by base substitutions. The multiplier is applied to indel errors
too, which are much rarer. On generic ground a read shows an indel the sample does not carry 0.4
to 13 times in 10,000 (§2); with the multiplier the model charges it about 30 to 45 times in
10,000 at 30× and 60 to 95 at 300×. So a real indel's reads become cheap to explain away as errors and the indel is missed.
**Setting that one multiplier back to 1, and changing nothing else in the fitted file, gives back
every missed indel and every SNP gain at once** (§3).

**The prototype fix keeps both gains.** Charging an indel the genotype cannot explain at the
read's own error, without the multiplier, restores generic-path indel recall to the defaults' at
every depth tested and leaves every fitted SNP number unchanged (§4). Over three GIAB samples it
takes generic indels missed from 61 to 50 at 5×, 26 to 18 at 10×, 10 to 7 at 30× and leaves 300×
at 4; SNP false calls stay at the fitted 17, 13, 5 and 2 against the defaults' 24, 34, 8 and 11.
It moves nothing in a `--defaults` run, and needs no new parameter or fit.

**Defect B: the reads that make `0/1` win are not indel errors, so an indel error rate does not
fix them.** Under `--defaults` a reference-looking read at an indel locus is charged about 5 in
10,000 — within a factor of two of the measured rate on short-repeat ground (§2). To make `1/1` win
at the colleague's sites each such read would have to be charged at least 1 in 45, 17 times
the highest rate measured on generic ground. The reads come from two other places (§5):

- **reads that end close to the indel, which the aligner writes without it.** At truth `1/1`
  indels of 7 bases or more, reads running only 5 to 9 bases past the indel show the reference 34
  times in 100; reads running 25 bases or more past it, 0.6 times in 100;
- **reads from a second copy.** Four of the colleague's five generic-path sites are in
  low-mappability or segmentally duplicated sequence, and three of her ten positions recur in two
  samples, each with 12 to 18 reference reads per 100.

The first has a direct fix — require a read to extend further past an indel before it counts as
reference evidence, a sibling of 7a1e2417 — which this branch measures but does not build.

**One decision was needed** (§7): whether to amend the read likelihood so that the calibration
multiplier applies to substitution errors only. The owner ruled yes.

---

## 1. What an indel costs today

**A read either shows an allele the genotype carries, and is charged only for which copy it came
from, or it does not, and is charged for being wrong** (`read_likelihoods.md` §3.3). The charge for
being wrong is, per read,

```text
ε_read × scale(read group) / m
```

- **`ε_read`** is the worse of two numbers the walk mints for each read: the **lowest base quality
  over the locus window**, and the mapping quality (`minted_ln_read_error`,
  `src/locus_generation/pileup/open_record.rs:3607`). An insertion or deletion inside the window
  contributes the lowest quality over itself plus one base each side, freebayes' `l + 2` window
  (`indel_bq_proxy_insertion` / `_deletion`, `decompose.rs:181`). **So an indel is charged by base
  quality, which describes substitutions; nothing in the number is about indels.**
- **`scale`** is the calibration multiplier (§3.2 there). Under `--defaults` it is 1 for every read
  group (`parameters_file/defaults.rs:441`). Fitted, it is the pre-pass's fitted error rate divided
  by the geometric mean of the minted errors, one number per read group, applied to every wrong
  read whatever kind of error makes it wrong (`likelihood/mod.rs:503`, `charged_error`).
- **`m`** is 3 when the read differs from every allele the genotype carries by one substitution,
  and 1 otherwise — so 1 for every indel (§3.5 there; `fill_error_spreads`, `generic.rs:102`).

**Who sets the fitted rate.** The pre-pass fits one error rate per read group, plus a per-sample
"noisy site" class, from `(depth, alternative-read count)` histograms over generic sites, where a
read over a reference copy shows another base with probability `ε` (`parameter_prepass_generic.md`
§2, §2.1). The rate the read likelihood calibrates against is the share-weighted average of the
clean and noisy rates (`read_likelihoods.md` §3.7). That spec says outright that one flat rate
"understates the alternative-read rate wherever the non-uniform contributors pile up and overstates
it everywhere else" (`parameter_prepass_generic.md` §2) — this report measures how that reaches
indels.

**So, for a read the genotype explains only by an indel error:**

| | charged per read | HG002 5 Mb set |
|---|---|---|
| `--defaults` | lowest base quality over the window, floored by mapping quality | about 5 in 10,000 |
| fitted | the same × the multiplier | × 6.2 at 10×, 6.4 at 30×, 14.1 at 300× |

**What the prior adds.** SNPs and indels share one diversity parameter (`calling_priors.md` §4.2,
Q1 open). It is not the lever here: the fitted expected heterozygosity on HG002 at 30× is 0.00097
against the default's 0.001, and §3 shows the multiplier alone accounts for the whole change.

## 2. The real indel error rate on generic ground

**Measured on invariant ground at 300×.** Every reference stretch inside the samples' GIAB
confident regions with no truth variant within 15 bases was classed by the longest exact tandem
repeat covering it — homopolymer length, or the repeat unit's length and copy count, or no repeat.
For each stretch I counted reads (mapping quality 20 or more, as pvc admits) that cross it with at
least 5 aligned bases each side, and among them the reads carrying an insertion or deletion inside
it. The sample carries the reference there, so each such read is an error. Beside it, the error
pvc would charge such a read if it were the reference read at an indel locus: the geometric mean of
the lowest base quality over the stretch plus one base, floored by mapping quality — the
uncalibrated `ε_read` of §1. Tool: `benchmarks/giab/src/indel_error_by_repeat.py`.

Reads in 10,000 crossing the stretch that show an indel there, HG002's 1,000-region set (5.1 Mb),
with HG003 and HG004 (100 regions each) beside it:

| context | stretches | HG002 | HG003 | HG004 | charged under `--defaults` | charged fitted, 30× |
|---|---:|---:|---:|---:|---:|---:|
| no repeat | 1,738,648 | 0.4 | 0.5 | 0.4 | 4.4 | 28 |
| homopolymer 3 | 186,502 | 1.1 | 1.2 | 1.1 | 5.0 | 32 |
| homopolymer 4 | 63,220 | 1.6 | 1.6 | 1.7 | 5.1 | 33 |
| homopolymer 5 | 21,530 | 3.2 | 3.4 | 2.9 | 5.3 | 34 |
| homopolymer 6 | 6,266 | 5.2 | 4.3 | 6.2 | 5.3 | 34 |
| homopolymer 7 | 2,401 | 8.3 | 8.8 | 9.0 | 5.2 | 33 |
| dinucleotide, 2 copies | 109,543 | 0.9 | 1.0 | 0.7 | 5.5 | 35 |
| dinucleotide, 4 copies | 1,204 | 3.2 | 1.6 | 1.9 | 6.2 | 40 |
| dinucleotide, 5 copies | 272 | 12.8 | 10.6 | 12.8 | 6.8 | 44 |
| *homopolymer 8 (STR path)* | 816 | 22 | 13 | 26 | 5.4 | — |
| *homopolymer 12 (STR path)* | 364 | 292 | 293 | 306 | 37 | — |

The charged columns are HG002's; "fitted, 30×" multiplies by that run's 6.44. The full table, every
period to 6 and copy count to 8, is the tool's output.

**Read the table in two halves.** Below the STR path's routing floors — homopolymers up to 7,
dinucleotides up to 5 copies, which is the generic path's repeat ground — **the default charge is
within a factor of two of the truth**: four to ten times too high where there is no repeat, about
right at 5- and 6-base homopolymers, 1.6 to 1.9 times too low at 7-base homopolymers and 5-copy
dinucleotides. **The fitted charge is 3 to 70 times too high everywhere on that ground.** Above
the floors the default charge falls far behind (an 8-base homopolymer slips 22 times in 10,000,
charged 5) — but that ground is the STR path's, which has its own slippage model.

**The three samples agree** to within the counting noise of the rarer contexts, as expected of one
sequencing batch.

**The STR path's fit already estimates this curve, without a truth set.** Its fitted share of reads
that slip, from this one HG002 sample, against the truth-masked measurement (per 10,000 reads):

| homopolymer length | 8 | 9 | 10 | 11 | 12 |
|---|---:|---:|---:|---:|---:|
| STR fit, 10× | 19 | 50 | 121 | 203 | 274 |
| STR fit, 300× | 23 | 56 | 137 | 194 | 243 |
| measured on invariant ground | 22 | 50 | 130 | 194 | 292 |

So the estimator for an indel-error curve exists and works at one sample and at 10×; it stops at
the routing floor because the census only gives it catalogued tracts.

**One figure in the STR pre-pass spec is not an indel error rate.** `parameter_prepass_ssr.md`
quotes ordinary indel error at "roughly half a read in 100 at every repeat count" at dinucleotide
tracts. That table was measured on the tomato archive with no truth set, counting every read whose
tract length differs from the reference by a non-whole number of copies — real variants included.
On truth-masked GIAB ground the same kind of movement at dinucleotides is 0.7 in 10,000, about 70
times lower. The spec's use of it — to say where whole-copy slippage outgrows other movement — is
unaffected; it should not be quoted as an error rate.

## 3. Defect A reproduced, and the one number behind it

**Data.** HG002's 1,000 GIAB confident regions (5.1 Mb, 7,975 truth variants of which 1,120
indels), from the full-depth CRAM in `benchmarks/human_genome_bottle/` and downsampled (seed 42) to
10× and 30×; the GIAB v4.2.1 truth, which those regions sit wholly inside. And the per-sample panel
(HG002, HG003, HG004, 100 regions each) at 5×, 10×, 30× and 300×. Each sample was walked once,
fitted on its own psp, and called with `--defaults` and with `--parameters`; the hidden-duplication
filter was off (`--hidden-dup-fdr 0`) so that a concurrent change to it cannot move these numbers.
Scored with RTG vcfeval 3.13 inside the confident regions, twice: genotype match and allele match
(`--squash-ploidy`).

**Error words used below.** *Missed*: a truth variant vcfeval cannot match even ignoring the
genotype. *False*: a call matching no truth allele. *Wrong genotype*: matched by allele, not by
genotype. *Generic ground*: the variant sits in no repeat, a homopolymer below 8 bases, or a tandem
repeat below the STR path's copy floors; everything else is *STR ground*. These are classes of the
reference sequence at the variant, not of which path called it.

**The colleague's pattern reproduces**, smaller because the regions are smaller. On the 5 Mb set at
10×, generic-ground indels missed go from 48 with defaults to 67 fitted (no repeat 17 → 22,
homopolymers of 4–7 bases 9 → 14, short tandem repeats 22 → 31); STR-ground indels barely move
(78 → 81). At 30×, 15 → 20. On the three-sample panel, 50 → 61 at 5×, 18 → 26 at 10×, 7 → 10 at
30×, 4 → 4 at 300×.

**The lever is the multiplier, and nothing else.** The fitted file with only
`error_probability_multiplier` set to 1.0 — prior, inbreeding, STR slippage and everything else
still fitted:

| 5 Mb HG002 | `--defaults` | fitted | fitted, multiplier 1 |
|---|---:|---:|---:|
| 10×: generic indels missed | 48 | 67 | **48** |
| 10×: SNP false / missed | 92 / 396 | 21 / 484 | **92 / 396** |
| 30×: generic indels missed | 15 | 20 | **15** |
| 30×: SNP false / missed | 38 / 58 | 17 / 71 | **38 / 58** |

Every SNP count and the generic indels missed return exactly to the defaults' (one generic
wrong genotype at 10× stays: 4 with defaults, 5 in both fitted arms). **The changes that remain are on
STR ground** — wrong genotypes 4 → 8 at 30× — and come from the STR side of the fit; nothing in
this report touches them.

**It is not a one-sample artefact.** Fitting the GIAB trio's three samples together over their
shared 100 regions gives multipliers of 6.53, 7.13 and 8.99 at 30× against 6.51, 7.21 and 8.91 from
each sample alone; at 300×, 13.4, 15.0 and 19.2 against 13.8, 15.5 and 19.6. The multiplier is a
property of the read group, and a cohort fit measures the same thing. (The trio's indel error
counts are too small to compare call by call: 1 to 5 per arm.)

**Two things noticed and not pursued.** The multiplier roughly doubles from 30× to 300× (6.4 to
14.1 on HG002's 5 Mb set; 6.5 to 13.8 on the trio's HG002) though it is meant to be a property of
the library, which suggests the fitted rate absorbs something that grows with depth. And for SNPs
the multiplier is itself a trade at low depth: on the 5 Mb set at 10× it removes 71 false SNPs and
adds 88 missed ones; at 30× and 300× it is a net gain (104 → 95 and 74 → 55 SNP errors). The
colleague's chr22 numbers show the same shape (10×: 1,872 fewer false, 1,342 more missed).

## 4. Fixing Defect A: the options, the prototype, and the ends of the range

**The prototype** (commit 18dd5712, `PVC_EXP_INDEL`; unset, every output is byte-identical to main —
checked on HG002 30×, defaults and fitted). For a genotype that can explain an observed allele only
by an insertion or deletion — every allele it carries differs from the observed one by one pure
indel — the wrong read is charged by one of:

- **`nomult`** — the read's minted error without the multiplier;
- **`table`** — a fixed rate by the repeat the indel sits in, from §2's HG002 measurement (the
  repeat found in the carried allele with 80 bases of reference either side);
- **`max`** — the larger of those two; **`floor`** — the larger of the calibrated charge and the
  table.

`max` gave the same calls as `nomult` and `floor` the same as main's fitted run (the table never
reaches the multiplied charge on generic ground), so the tables show the two that differ.

| | generic indels missed / false / wrong genotype | SNPs missed / false / wrong genotype |
|---|---|---|
| **5 Mb HG002, 10×** | | |
| `--defaults` | 48 / 0 / 4 | 396 / 92 / 49 |
| fitted | 67 / 0 / 5 | 484 / 21 / 52 |
| fitted, `nomult` | **48 / 0 / 5** | 483 / 21 / 51 |
| fitted, `table` | **40 / 3 / 4** | 483 / 21 / 52 |
| **5 Mb HG002, 30×** | | |
| `--defaults` | 15 / 4 / 1 | 58 / 38 / 8 |
| fitted | 20 / 3 / 1 | 71 / 17 / 7 |
| fitted, `nomult` | **16 / 4 / 1** | 71 / 17 / 7 |
| fitted, `table` | **14 / 5 / 1** | 71 / 17 / 7 |
| **5 Mb HG002, 300×** | | |
| `--defaults` | 7 / 2 / 0 | 43 / 27 / 4 |
| fitted | 9 / 2 / 0 | 43 / 8 / 4 |
| fitted, `nomult` | **7 / 2 / 0** | 43 / 8 / 4 |
| fitted, `table` | **7 / 2 / 1** | 43 / 8 / 4 |
| **3-sample panel, 5× / 10× / 30× / 300× — generic indels missed** | | |
| `--defaults` | 50 / 18 / 7 / 4 | |
| fitted | 61 / 26 / 10 / 4 | |
| fitted, `nomult` | **50 / 18 / 7 / 4** | SNPs identical to fitted at every depth |
| fitted, `table` | **49 / 16 / 7 / 4** | SNPs identical to fitted at every depth |

**Neither costs SNPs**: the SNP columns of both are the fitted run's to within one call, because
neither touches a substitution. **`table` on `--defaults`** (not shown) recovers the same few
indels — 48 → 40 at 10× on the 5 Mb set — for the same 3 extra false generic indels.

**Tomato at about 3×** — the tomato2 oracle's four accessions over its first 20 regions, fitted
jointly (multipliers 2.7, 9.7, 9.8, 4.9). There is no truth set, so this measures only what moves.
`--defaults`: `nomult` moves nothing; `table` adds 7 generic indel records and changes 22 sample
genotypes, out of 7,652 records. Fitted: `nomult` adds 36 generic indel records, removes 10 and
changes 43 genotypes, out of 6,888; `table` adds 67, removes 23, changes 70. Nothing but generic
indels moves in any arm.

**The three options at the ends of the range:**

| | one sample → thousands | 3× → hundreds | estimated without a truth set by |
|---|---|---|---|
| **`nomult`** — no multiplier on indel errors | identical: nothing is estimated | identical arithmetic at every depth; at 3× the error charge decides calls (`read_likelihoods.md` §3.8), and this one is the instrument's | nothing to estimate |
| **fixed table** | identical | identical | nothing — the numbers are one HG002 Illumina library's |
| **fitted curve per read group** — the STR slippage curve carried below the routing floor, plus a flat term where there is no repeat | per read group, so one sample has one; pooled over every generic site the census holds, so thousands of samples cost one curve each | pooled over thousands of sites per read group, so 3× is thin per site and rich in total; §2 shows the STR fit recovering hp8–12 from one sample at 10× | the census and the joint fit, summing over the unknown genotype as the STR fit already does — no catalogue of known variants |

- **`nomult` is the right first step.** It removes the defect, costs nothing measurable, adds no
  parameter and needs no fit. Its weakness is §2's: a base quality is not an indel rate, so it
  overcharges indels outside repeats about tenfold and undercharges at 7-base homopolymers and
  5-copy dinucleotides by up to twofold. Both are the defaults' behaviour today.
- **A fixed table is rejected as a default.** It is one library's numbers. A PCR-amplified library
  slips more than this PCR-free one, and tomato's libraries are not human Illumina at all; a table
  fitted here would be confidently wrong there, which is the corner `CLAUDE.md` says a design must
  not take.
- **A fitted curve is what the table approximates, and it is a separate spec.** What `table` buys
  over `nomult` is the size of that spec's prize: 8 fewer missed generic indels at 10× on 5 Mb for 3
  more false ones, 2 fewer for 1 more at 30×, none at 300×. That is real but small, and it needs the
  census to count indel reads at generic sites by repeat context — a change to what the pre-pass
  stores. Recommended to defer until the clearance fix of §5 is in, because the two act on the same
  reads.

**What the prototype does not cover.** Reads that saw only part of the locus are still charged the
calibrated error when no allele fits them (`score_partials`); an amendment should treat them the
same way. And a genotype carrying alleles that differ from the observed one in different ways (one
by an indel, one by substitutions) keeps main's charge.

## 5. Defect B: where the reference reads come from

**It does not reproduce locally.** At 300× on the 5 Mb set, pvc writes no generic-ground indel
with the wrong genotype, with defaults or fitted; the colleague's group 2 is 10 sites among about
19,000 truth indels over three samples, so about 0.3 would be expected here. So this section works
from the reads.

**What the colleague's sites show** (`het_for_hom_reads.tsv`, group 2): of the five generic-path
sites, four are low-mappability or segmental duplication; the fifth, HG003 chr22:47088098 `T>TAAA`,
slides over 3 bases, where §2 measures about 1 indel error in 10,000 reads, yet 50 of its 318
crossing reads are reference. Three positions appear twice, in HG002 and HG003 (chr22:23344056,
42616310, 49449094) — a reference share that recurs at one position in two people is a property of
the sequence or the mapping, not of chance slippage.

**What `1/1` would need.** At 51 reference reads among 288, `0/1` scores 288 × ln ½ = −199.6 for
the reads. For `1/1` to beat it, the 51 must be charged at least about 1 in 45 each (reads only; the
prior is negligible at this depth). That is 17 times the highest generic-ground rate in §2
(dinucleotides at 5 copies, 1 in 780) and 44 times what defaults charge. A rate that large
would also make a heterozygote's reads cheap to dismiss under `0/0` at low depth.

**Where reference reads at true homozygous indels actually come from — measured.** At the 423
truth `1/1` indels in HG002's 5 Mb set with at least 20 crossing reads, most show almost no
reference reads (190 none; 134 under 2 in 100); 18 show 10 or more in 100. Split by how far each
read runs past the indel's slide before it ends (`benchmarks/giab/src/hom_alt_ref_reads_by_clearance.py`):

| indel length | 1–2 bases past | 3–4 | 5–9 | 10–14 | 15–24 | 25 or more |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 15.2 | 6.7 | 2.8 | 1.5 | 1.7 | 1.4 |
| 2–3 | 8.0 | 8.7 | 5.3 | 3.0 | 1.4 | 0.9 |
| 4–6 | 17.2 | 19.1 | 8.8 | 3.3 | 0.5 | 0.3 |
| 7 or more | 66.7 | 30.3 | 34.2 | 18.6 | 11.1 | 0.6 |

Reference reads in 100. **A read that ends close to an indel is often written by the aligner
without it** — matching the inserted bases as mismatches or clipping them costs less than opening a
gap near the end — and that read then shows the reference although it came from the alternative
allele. For short indels the effect is gone by about 10 bases past the slide; for indels of 7 bases
or more it lasts to about 25. The colleague's own rule for an "informative" read is 5 aligned bases
each side, so her group-2 reference counts include exactly these reads.

**pvc already escapes most of these, and the reason points at the fix.** At the 16 high-reference
sites, every indel locus pvc built (14) is called `1/1`, by wide margins (for example chr2:138668090, `1/1` −56
against `0/1` −235), because its loci there are wide: a read ending near the indel does not cover
the whole locus and counts as a *partial* read, which is compatible with every allele. Where a
locus is no wider than the slide — the case 7a1e2417 widened to the slide's last base — a read
needs only one base past the slide to count as reference. The same artefact also makes false SNPs:
four bases before a truth 8-base insertion at chr11:84218398, 32 of 197 reads carry the insertion
written as mismatches, and pvc writes a `0/1` G>A at chr11:84218394.

**The fix this suggests is evidence shaping, not an error rate**: a read counts as complete at an
indel locus only if it runs at least about 10 bases past the slide for indels up to 6 bases, and
about 25 for longer ones; otherwise it is a partial read. That is 7a1e2417's rule with a margin
measured from real alignments instead of one base. It is not realignment: every read keeps its
placement. It costs information, not bias: of the reads crossing a truth `1/1` indel in the 5 Mb set, 1 in
10 runs less than 10 bases past a 1-base indel and would become partial, 1 in 14 at 4–6 bases. Not built here: records open at their left-most anchor and
cannot reach left, so a clean implementation needs a decision about how the walk marks a read as
partial without widening the record — its own spec.

**The second-copy sites are the hidden-duplication filter's ground**, which a parallel session is
working on, and are not addressed here.

## 6. What would confirm this on chr22

The colleague's own files settle both claims at her scale, without a new walk:

- **Defect A.** Her fitted parameters files, one per sample and depth: read
  `error_probability_multiplier` (expected about 6–9 at 10–30×, 13–20 at 300×). Then, over her
  fitted-run psps, call twice more with the same file: once with the multiplier edited to 1.0, once
  with this branch's binary (f800df67 or later). Expected: generic-path indel errors back to
  the defaults' level (her section 8: homopolymer 67 → 36, other tandem repeat 96 → 55, extreme GC
  21 → 6, rest 62 → 25 at 30×), SNP errors at the fitted level.
- **Defect B.** Her 300× CRAMs and psps at the five generic group-2 sites
  (`het_for_hom_reads.tsv`): run `hom_alt_ref_reads_by_clearance.py` on them, and `--explain-loci`
  to see the locus spans. Expected: the reference reads pvc counts end within about 10–25 bases of
  the slide at the mappable site, and at the four low-mappability sites the reference share holds at
  any clearance (a second copy, not alignment).

## 7. Recommendation and the decision

**Decision: amend the read likelihood so that the calibration multiplier applies only to errors
that are substitutions** — option `nomult`, step one of §4. I recommend it.

- **What changes.** Fitted runs only: generic-path indels the genotype cannot explain are charged
  the read's own minted error, as `--defaults` charges them today. Measured gain: every generic indel
  the fit lost comes back (§4), at every depth from 5× to 300×. Measured cost: none on GIAB in SNPs
  or in indel false calls; on tomato at 3×, 46 generic indel records and 43 genotypes move among
  6,888 records, with no truth to say which way.
- **What does not change.** `--defaults` output is byte-identical. Nothing new is fitted or stored;
  the parameters file is unchanged.
- **If not adopted**, fitted runs keep losing generic indels: on the colleague's chr22, 2,165 →
  2,630 indel errors at 10× and 610 → 753 at 30×.

**Two follow-ups, each its own spec, in this order:** (1) the clearance rule of §5 for Defect B;
(2) a fitted indel-error curve per read group (§4's third option), if the clearance rule leaves the
prize at the size measured here.

### Proposed amendment — `read_likelihoods.md` §3.5, a new closing paragraph

> **AMENDED (proposed 2026-10-21, `ng_generic_indel_error_2026-10-21.md`): the calibration scale
> applies to substitution errors only.** §3.2's scale is the ratio of a fitted rate to the mean
> minted error, and the fitted rate is learned from reads that disagree with the reference at
> generic sites — almost all of them by a substitution. Indel errors are rarer: on HG002 at 300× a
> read shows an indel the sample does not carry 0.4 to 13 times in 10,000 on generic ground, while a
> fitted scale of 6 to 14 charges it 30 to 95 times in 10,000, so a fitted run explains real indels
> away and misses them. **So when every allele a genotype carries differs from the observed
> allele by one insertion or deletion, the wrong read is charged its minted error unscaled —
> `ε̄` without `scale` — and `m = 1` as before.** The same holds for a partial read that no allele
> fits, when every carried allele is an indel away from what it showed. A genotype whose carried
> alleles differ from the observed one in mixed ways keeps the scaled charge. `--defaults` runs, whose
> scale is 1, are unchanged.
>
> This is a stopgap with a known shape, not an indel model: a minted error is a base quality, which
> overstates indel error outside repeats about tenfold and understates it at 7-base homopolymers
> and 5-copy dinucleotides up to twofold. The model that replaces it is a fitted indel-error curve
> per read group, of which the STR path's slippage curve is the upper half (§10).

And in §10 (deferred), a new row: *a fitted indel-error curve on the generic path, by repeat unit
and copy count, per read group — the STR slippage curve extended below the routing floor with a
flat term for no repeat; home: `parameter_prepass_generic.md` and the census; prize measured in
`ng_generic_indel_error_2026-10-21.md` §4.*

## How to reproduce

All scratch lives under the worktree's `tmp/` (not committed); the scripts there are the ones used.

- The rule as built: this branch from f800df67. The four prototype arms: commit 18dd5712 with
  `PVC_EXP_INDEL=nomult` (or `table`, `max`, `floor`) on `call-from-psps`; unset, output equals
  main.
- Error-rate curve: `uv run --with pysam --with numpy python -I
  benchmarks/giab/src/indel_error_by_repeat.py SAMPLE ALIGNMENT BED TRUTH_VCF REFERENCE OUT.tsv`
  (also writes `OUT.tsv.hom_alt.tsv`, the per-site table §5 reads).
- Clearance table: `benchmarks/giab/src/hom_alt_ref_reads_by_clearance.py OUT.tsv.hom_alt.tsv
  ALIGNMENT REFERENCE ENDS.tsv`.
- The 5 Mb evaluation set: `benchmarks/human_genome_bottle/HG002_bench_azar_merged_1000.bed` with
  the GIAB v4.2.1 HG002 benchmark VCF (downloaded from GIAB's FTP into `tmp/hg5mb/`); its
  downsamples are `samtools view -s 42.098764` (30×) and `42.032921` (10×) of
  `benchmarks/human_genome_bottle/crams/HG002_reads_selected_1000_rg.cram`.
