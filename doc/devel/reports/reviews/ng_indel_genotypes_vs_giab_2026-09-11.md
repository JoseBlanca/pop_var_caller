# Why ng calls homozygous insertions heterozygous

**Date:** 2026-09-11
**Branch:** `ng-indel-gt`; the calls examined were produced on `ng-accuracy`
(worktree `/Users/jose/devel/pop_var_caller-ng-accuracy`) on 2026-09-11.
**Data:** `benchmarks/giab/per_sample` — HG002, HG003 and HG004, each on its own
random 100 confident intervals of GRCh38 (450–570 kb a sample) against its own GIAB
v4.2.1 truth VCF, at 5×, 10×, 15× and 30×. ng run as
`generate-psps` → `estimate-parameters` → `call-from-psps` with the fitted
parameters, every caller gated at QUAL ≥ 30.
**Scoring:** `benchmarks/lib/score_against_truth.py` and
`benchmarks/lib/genotype_disagreements.py`.

---

## The answer in one paragraph

**The genotyping is doing what its input says; the input is wrong.** At the 17
truth-homozygous indels ng genotyped as heterozygous at 30×, ng credits 83 reads to
the reference allele. Realigning every read at those sites against both haplotypes
directly — no gap scoring, just which of the two sequences the read's bases fit —
**one** read fits the reference better, 24 fit the two equally, and the rest fit the
alternative. At 16 of the 17 sites **not a single read supports the reference**. The
credit is manufactured by the rule that decides what a read showed, and the reason is
one line of code: **an insertion's locus covers one reference base, the anchor, so a
read is asked only whether it agrees with the reference at that one base.** A read
that agrees there and contradicts the reference twenty bases downstream is recorded as
an exact reference observation.

**AD cannot mislead while the likelihoods are right, because they are the same
number.** `SampleReadCounts::allele_reads` — the `AD` column — is filled by summing
`row.support.num_reads` over the merge's `supported` rows
([`run/records.rs:232`](../../../../src/ng/run/records.rs#L232)), and the read
likelihood is built by walking the same rows through the same allele remapping
([`calling/evidence_shaping.rs:273`](../../../../src/ng/calling/evidence_shaping.rs#L273),
`GenericObservation::fill_from_supported_alleles`). A read reaches the genotype model
as *an allele identifier and a count*, never as a sequence, so there is no second
opinion for the likelihood to hold.

---

## 1. The truth set is right, and the scoring is right

Both possibilities the investigation was asked to rule out first were checked against
the GIAB VCF directly, and both are ruled out.

At **HG002 chr1:243535155** — a 19-base insertion — the truth record is `1/1` with
`GQ 314`, and GIAB's own read counts are `ADALL=1,261`: across every platform, one
read on the reference against 261 on the insertion. At **HG003 chr15:96140584** — a
24-base insertion — the truth record is `1/1` with `GQ 404` and `ADALL=84,346`.

The second site is the more useful of the two, because GIAB's own short-read count
there is not clean either: 84 reads in 430 sit on the reference, one in five. So
reference-looking reads at a long homozygous insertion are a real property of
short-read data, not an ng artefact. **What is an ng artefact is the size.** ng's own
count at that site is `AD=13,13` — half its reads — where GIAB's pipeline, on the same
kind of data, gets one in five, and freebayes on exactly these reads gets zero.

---

## 2. What separates read assignment from genotyping

Take the two records the investigation started from. ng's output at 30×:

```
HG002  chr1:243535155  T -> TATTTAAA , TATTTTAAAATATATTTAAA   QUAL 390.6  0/2:99:37:12,5,20
HG003  chr15:96140584  C -> CTGATTGGTCCATTTTACAGATGGT          QUAL 109.6  0/1:99:26:13,13
```

Now look at what the aligner actually put there. At **chr15:96140584**, the reads
overlapping the anchor split into two groups by their CIGAR:

```
  96140442  148M          NM:i:1      \
  96140445  148M          NM:i:1       |  aligned as reference
  96140461  144M4S        NM:i:1       |  (12 reads, NM 1 to 12)
  96140483  148M          NM:i:5       |
  96140491  148M          NM:i:5      /
  96140509  76M24I48M     NM:i:25     \
  96140521  64M24I60M     NM:i:24      |  carrying the 24-base insertion
  96140545  40M24I84M     NM:i:24      |  (14 reads)
  96140575  27M24I97M     NM:i:24     /
```

The reference around the anchor is an imperfect 24-base tandem repeat
(`CTGATTGGTCCATTTTACAGAGTG` `CTGATTGGTCCGTTTTACAGAGTG` `CTGATTGGTGCGTTTACAAACC`), so
BWA can place a read carrying the extra copy as plain `148M` with four or five
mismatches rather than opening a 24-base gap. **Realigning each read against the two
haplotypes says which is right**: the `148M` reads at 96140461–96140491 sit 4–5
mismatches from the reference and **1** from the insertion. They carry the insertion.

Four reads — those ending within nine bases of the anchor — fit both equally, and they
are the only reads at the site that genuinely cannot tell. **Zero fit the reference
better.**

ng credits 13 of the 26 to the reference.

### The same comparison over all 12 insertions

Realignment is ungapped against each full haplotype window, so the indel is inside the
sequence being matched rather than in a gap whose cost has to be chosen
([`benchmarks/giab/src/read_haplotype_support.py`](../../../../benchmarks/giab/src/read_haplotype_support.py),
added with this report):

| sample | site | inserted | ng credits REF | ng credits ALT | fits REF | fits ALT | fits both |
|---|---|---:|---:|---:|---:|---:|---:|
| HG002 | chr1:169767410 | 1 | 3 | 20 | 0 | 22 | 2 |
| HG002 | chr1:179298130 | 5 | 3 | 23 | 0 | 27 | 0 |
| HG002 | chr1:243535155 | 19 | 12 | 20 | 0 | 32 | 5 |
| HG002 | chr2:29156875 | 4 | 4 | 22 | 0 | 27 | 0 |
| HG003 | chr8:4310529 | 4 | 6 | 27 | 0 | 33 | 1 |
| HG003 | chr15:96140584 | 24 | 13 | 13 | 0 | 22 | 4 |
| HG003 | chr11:123865478 | 1 | 2 | 11 | 0 | 13 | 1 |
| HG003 | chr4:141191914 | 10 | 8 | 19 | 0 | 22 | 6 |
| HG004 | chr15:65223392 | 17 | 6 | 16 | 0 | 21 | 2 |
| HG004 | chr9:11024723 | 2 | 5 | 16 | 1 | 23 | 2 |
| HG004 | chr13:95122827 | 1 | 4 | 18 | 0 | 22 | 0 |
| HG004 | chr7:50121974 | 2 | 5 | 37 | 0 | 42 | 1 |
| **total** | **12 sites** | | **71** | **242** | **1** | **306** | **24** |

**Of the 71 reads ng puts on the reference allele, 1 fits the reference better than
the insertion, 24 fit the two equally, and 46 fit the insertion better.**

---

## 3. The line of code

A read's observation at a locus is **the bases it showed across the locus's reference
positions**, and a read is credited to an allele by byte equality against that allele's
sequence — `SequenceObservation::matches_reference` is a slice comparison
([`locus_generation/mod.rs:416`](../../../../src/ng/locus_generation/mod.rs#L416)).

How many reference positions a locus covers is decided by the event that opened the
record:

```rust
// src/ng/locus_generation/pileup/decompose.rs:55
pub fn footprint_span(&self) -> u32 {
    match self {
        ReadEvent::Match { .. } => 1,
        ReadEvent::Insertion { .. } => 1,
        ReadEvent::Deletion { deleted_len, .. } => *deleted_len + 1,
    }
}
```

**An insertion's footprint is one reference base however long the inserted sequence**,
because the anchor is the only reference base the event touches. That is a faithful
description of the *event*. It is the wrong window for the *question* — *did this read
come from a haplotype carrying the insertion?* — because the answer lives in the bases
after the anchor, and the comparison never reaches them.

The specification saw the hole and left it open.
[`read_likelihoods.md` §5.4.1](../../ng/spec/read_likelihoods.md) says an insertion has
no partial reads at all — *"An insertion's reference span is its anchor base, however
long the inserted sequence — width stays 1, so every read covering the anchor is
complete"* — and parks the rest in a parenthesis: *"(A read that ran out inside a long
inserted sequence is a different problem — a truncated allele sequence, not a partial
witness — and it belongs to read preparation, not here.)"* Read preparation never
took it, and the measurement above says the parenthesis is understated twice over: the
reads that ran out are the *smaller* half (24 of 71), and the larger half is reads that
reached well past the insertion and disagree with the reference there.

---

## 4. How big it is, and how it grows

### With insertion length

Reads credited to the reference as a share of the site's reads, at the
**truth-homozygous** indels ng found, three samples pooled, QUAL ≥ 30, 30×. The
freebayes column is the same sites, the same reads, its own `AD`:

| inserted bases | sites | ng | freebayes |
|---:|---:|---:|---:|
| 1 | 39 | 0.039 | 0.008 |
| 2–3 | 14 | 0.052 | 0.003 |
| 4–6 | 4 | 0.133 | 0.000 |
| 7–12 | 2 | 0.154 | 0.000 |
| 13 and over | 3 | 0.388 | 0.210 |

**The prediction holds and the comparison is the point.** ng's reference credit at a
homozygous insertion rises from 4 reads in 100 at one base to 39 in 100 above thirteen.
freebayes, which compares a read against a haplotype window rather than against a
single anchor base, credits essentially nothing below 13 bases — and at 13 and over,
where a 150-base read stops being able to settle it, it climbs to 21 in 100 and
freebayes gets 2 of those 3 sites wrong too.

**That one quantity is the whole indel-genotype gap between the two callers.** Rerun
here from `genotype_disagreements.py` over all three samples at 30×: ng gets the
genotype wrong at **21 of the 297** truth indels it found, freebayes at **9 of 301**,
the production caller at **64 of 307**. *(The brief this investigation started from
gave ng's denominator as 276; 297 is what the documented command returns, and the
other two callers' figures reproduce exactly. It moves ng's error rate from 7.6% to
7.1% and changes nothing else.)*

### With depth: flat

Same measure, by coverage (the number in brackets is how many sites the cell pools):

| inserted bases | 5× | 10× | 15× | 30× |
|---:|---:|---:|---:|---:|
| 1 | 0.043 (26) | 0.039 (35) | 0.037 (26) | 0.040 (41) |
| 2–3 | 0.075 (10) | 0.058 (15) | 0.022 (10) | 0.060 (15) |
| 4–6 | 0.250 (3) | 0.148 (3) | 0.158 (3) | 0.133 (4) |
| 7–12 | — | 0.222 (2) | 0.154 (2) | 0.154 (2) |
| 13 and over | 0.250 (1) | 0.353 (2) | 0.357 (2) | 0.388 (3) |

**Depth does not touch it**, which is what a systematic misclassification looks like
and what sampling noise does not: the miscounted reads scale with coverage exactly as
the correctly counted ones do. The absolute number of wrong genotypes is flat too —
12, 17, 11, 17 across 5× to 30× — while the number of truth indels ng finds grows from
103 to 297.

### How few wrong reads it takes

Over all four depths, 363 truth-homozygous indels ng found, ordered by the share of
reads it credited to the reference:

| reference-credited share | sites | not called 1/1 |
|---|---:|---:|
| below 0.02 | 244 | 3 |
| 0.02–0.05 | 17 | 0 |
| 0.05–0.08 | 23 | 0 |
| 0.08–0.10 | 13 | 0 |
| 0.10–0.12 | 15 | 6 |
| 0.12–0.15 | 12 | 9 |
| 0.15–0.20 | 10 | 10 |
| 0.20 and over | 29 | 29 |

**The genotype turns at about one read in eight.** Below one in ten a homozygote
survives (297 sites, 3 lost); above one in six it never does (39 sites, 39 lost).

That threshold is not itself a defect — it is what any genotype likelihood does, and
freebayes' is in the same place. A read the genotype cannot produce is charged its own
error probability, tens of Phred; a read it can produce buys only `log 2`, 3.0 Phred,
for the homozygote over the heterozygote. **So the turn comes at about one read on the
reference for every seven on the insertion**, which is what the table shows and what
the model predicts: this sample's fitted calibration multiplies each read's reported
error by 6.4, putting a typical Q30 base at an error of 6 in 1,000, and
`ln(0.5 / 0.006) / ln 2` is 6.3. With a correct read assignment that trade is the
right one; the defect is the assignment, not the trade.

---

## 5. Deletions: a second, smaller route, through partial reads

At a deletion the footprint is `deleted_len + 1`, so the comparison does cover the
deleted bases and the length trend disappears: the reference-credited share is 2 in 100
at one base, 3 in 100 at two or three, and **zero** at four and above, at every depth.

Five of the 17 wrong genotypes are nonetheless deletions, and they split three ways.
**Three are the same defect on a shorter window** — `chr1:1388565` (`AD=6,22`),
`chr1:68898298` (`AD=4,21`) and `chr19:39771337` (`AD=2,12`), all short deletions
inside a repeat, where three or two reference positions are still not enough ground to
tell the aligner's placement from the truth. **One is the repeat-tract path** and
belongs to §6. **And one comes a different way entirely.**

`HG002 chr1:106818720` is a 7-base deletion called `0/1` with **`AD=0,15`** — no
reference reads at all. Its 19 reads are 15 carrying the deletion and 4 that end
between one and seven bases into the locus. A read that stops inside a locus becomes a
*partial* observation, scored by whether each allele could have produced the bases it
showed
([`calling/likelihood/generic.rs:388`](../../../../src/ng/calling/likelihood/generic.rs#L388),
`allele_is_compatible_with_partial`). For a read anchored at the left border and not
reaching the right one, the test is `allele.starts_with(bases)`.

**At a deletion that test can only ever come out one way.** The alternative allele is
the anchor base alone — one byte — so a partial holding three or seven bases is never a
prefix of it, while it is always a prefix of the reference allele. Every multi-base
partial at a deletion is therefore reference-only evidence, and is charged against the
homozygous-deletion genotype at its full read quality.

**Those four reads are uninformative, not reference reads.** The deleted run is
`ACTTAAA` and the reference immediately after it is `ACTTAAT` — a tandem duplication —
so a read that stops within seven bases of the anchor shows the same bases on either
haplotype. The realignment says so: 0 fit the reference better, 3 fit both equally, 16
fit the deletion. Four uninformative reads outvote fifteen informative ones, and `GQ`
is 47.

This is the case [`read_likelihoods.md` §5.4.1](../../ng/spec/read_likelihoods.md)
predicted in words — *"A deletion has [partial reads], and they carry the reference
allele preferentially"* — without noticing that the compatibility rule makes it
unconditional rather than merely likely.

---

## 6. Repeat tracts: the error runs the other way, and is smaller

Splitting the 21 disagreements by the ground ng's own region typing assigns them
(`examples/ng_typed_region_dump`, catalog floors — the floors a calling run uses):

| ground | class | direction | sites |
|---|---|---|---:|
| ordinary sequence | insertion | truth 1/1 → called 0/1 | 12 |
| ordinary sequence | deletion | truth 1/1 → called 0/1 | 3 |
| repeat tract | deletion | truth 1/1 → called 0/1 | 2 |
| repeat tract | deletion | truth 0/1 → called 1/1 | 3 |
| ordinary sequence | insertion | truth 0/1 → called 1/1 | 1 |

**15 of the 17 homozygote-as-heterozygote calls are on ground ng types as ordinary
sequence and sends down its generic SNP/indel path.** That is a difference from the
production caller, whose version of this failure the report of 2026-07-04 put at 86% in
homopolymers and short tandem repeats. ng's region typing does not call most of these
sites repeats — but the local sequence at nearly all of them is repetitive enough for
the aligner to place the alternative haplotype as reference, which is the property that
matters here and is not the same property.

**The reverse direction is small, is at tracts, and is the same rule mirrored.** Of the
4 truth-heterozygous indels ng called homozygous, 3 are at tracts where ng credits
*fewer* reads to the reference than support it — `HG003 chr10:5596734` has `AD=0,8`
where 5 reads of 18 fit the reference, `HG004 chr8:24760294` has `AD=0,28` where 6 of
31 do. The reference-supporting reads became partial observations and stopped counting.
The fourth, `HG004 chr10:4848078`, is a compound site: truth carries a 1-base deletion
*and* a 4-base insertion, ng emitted only the deletion, and calling it homozygous
follows. That one is a missing candidate allele, not this defect.

---

## 7. What a fix looks like

The fix cannot live in calling. A psp stores folded observations, not reads: the bases
a read showed beyond the locus footprint were never written down, so `call-from-psps`
has nothing to re-examine. **The comparison has to be widened where the observation is
built**, in the pileup walk.

Two shapes, both measured on the 30× BAMs by rebuilding the walker's own classification
over a wider window
([`benchmarks/giab/src/indel_comparison_window_sim.py`](../../../../benchmarks/giab/src/indel_comparison_window_sim.py),
added with this report).

**Read the predictions as the direction and the ordering, not as the caller's own
numbers.** Run at a window of one base it reproduces today's rule, and the control
says how close that is: it gives 0.099 / 0.073 / 0.170 / 0.186 / 0.360 across the five
length buckets where ng measures 0.040 / 0.052 / 0.133 / 0.154 / 0.388. The four
longer buckets land within a quarter of ng's figure; the one-base bucket is 2.5 times
ng's, because the simulation applies none of ng's read filters — no mapping-quality
floor, no base masking, no depth cap — and those remove most of its extra reads
exactly where the credit is one or two reads a site.

**Shape A — widen an insertion's record footprint to the inserted length**, so
`footprint_span` returns `inserted_len + 1` the way a deletion returns
`deleted_len + 1`. Reference credit at truth-homozygous insertions falls from
0.099/0.073/0.170/0.186/0.360 to **0.087/0.028/0.075/0.051/0.116** — **the length
trend goes away**, which is the property to look for, because the trend is the
signature of the defect. Across the 12 insertion sites that got a wrong genotype the
credit falls from 71 reads to 24, and 5 of the 12 lose it entirely. At truth
**heterozygous** insertions the reference share moves from 0.49–0.63 to 0.42–0.49 —
still either side of a half, and far above the one-read-in-eight turning point, so no
true heterozygote is put at risk. The cost is visible in the VCF: `REF` at an insertion
stops being the anchor base alone, which
[`run/callers.rs:5939`](../../../../src/ng/run/callers.rs#L5939) currently pins as a
test and three documents rest on, and every stored psp has to be regenerated.

**Shape B — a fixed 30-base confirmation window.** Reference credit goes to **zero** at
11 of the 12 sites (the twelfth is a 19-base insertion in a `TTTTAAAATATA` repeat where
the reads genuinely cannot settle it). But it turns a quarter to a third of the reads
into partial observations, and at truth heterozygous insertions the reference share
falls to 0.24–0.40. At 5× that is close enough to one-in-eight to start losing real
heterozygotes, and 5× is inside the range this caller has to work over.

**Recommendation: Shape A.** It removes two thirds of the false reference credit and
flattens the length trend, its window is a property of the locus rather than a
constant somebody has to tune, it reuses the widening machinery deletions already
have, and it leaves true heterozygotes at the allele balance they should have. Shape B
buys the last third by spending read depth, and read depth is the wrong currency for a
caller that has to work at three reads a position.

**What Shape A does not fix**, and should be said before it is chosen: it leaves 24 of
the 71 false reference reads in place, so the 13-and-over bucket still sits near 12 in
100 and some long insertions will still flip. Reaching those needs a read-to-allele
alignment over flanks — the machinery the repeat-tract path already has in
`src/ng/alignment` — rather than a byte comparison over any window. That is a larger
piece of work and it is separable: Shape A does not close the door on it.

**The deletion route needs its own fix and it is small.** `allele_is_compatible_with_partial`
should compare a partial's bases against the allele's sequence *continued into the
reference beyond the locus*, not against the allele's own sequence alone. Today a
7-base partial can never be a prefix of a 1-base deletion allele, so it is scored as
reference evidence with no possibility of any other verdict. One of the 17 wrong
genotypes at 30× is that, and the same site is wrong at 15× for the same reason; the
change touches one function and needs no psp regenerated, which makes it the cheapest
item here and worth taking whatever is decided about Shape A.

---

## 8. What this did not measure

- **20×, 50× and 300× were not on disk** when the investigation ran, and were still
  being generated when it was written. The depth table above covers 5× to 30×; the
  pattern over those four depths is flat, and 300× is the cell most likely to differ,
  because at 300× a long insertion's flanking reads are numerous enough for a
  haplotype-aware caller to separate and for ng's single-base comparison still not to.
- **Only three samples, only human, only 450–570 kb each.** The bucket at 13 bases and
  over holds 3 sites. The direction of the length trend is solid across every depth;
  its slope above 13 bases rests on very little.
- **Whether widening the footprint fragments the candidate allele set** was not
  measured. With a wider window every sequencing error downstream of the anchor makes
  its own byte string, so the alternative allele's reads scatter across several
  candidates and selection has to pool or drop them. The simulation counts a read as
  "not reference" without asking which candidate it lands on, so it predicts the
  reference credit correctly and says nothing about the allele table. **That is the
  first thing to measure if Shape A is taken.**
