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

## The answer

**The genotype arithmetic is right. What it is given to work with is wrong.**

There are 17 places where GIAB says the sample carries an indel on both copies of its
chromosomes and ng called it on one copy only. At those 17 places ng decided that 83 of
the reads came from a chromosome without the indel.

**They did not.** I took every read at those sites and asked which of the two possible
chromosomes it fits: the reference, or the reference with the indel written into it.
The read's bases are slid along each of the two sequences and the mismatches counted,
so the indel sits inside the sequence being matched and no gap cost enters the answer.
Of the 83 reads ng put on the reference, **one** fits the reference better than the
indel, 24 fit the two equally well, and the rest fit the indel. **At 16 of the 17
sites no read at all fits the reference better.**

**The reason is one line of code.** ng decides which allele a read supports by
comparing the bases the read showed across the stretch of reference the variant
occupies. For an insertion that stretch is a single base — the one the inserted
sequence hangs off — however long the insertion is. So the only question ever put to
the read is *do you agree with the reference at this one base?* A read that agrees
there and contradicts the reference twenty bases further on is filed as a reference
read.

**One thing this rules out.** It would be possible for `AD` to be a misleading summary
while the genotype model itself saw something better — but not here, because they are
the same number. The `AD` column is filled by summing `row.support.num_reads` over the
merge's `supported` rows
([`run/records.rs:232`](../../../../src/ng/run/records.rs#L232)); the read likelihood
is built by walking those same rows through the same allele renumbering
([`calling/evidence_shaping.rs:273`](../../../../src/ng/calling/evidence_shaping.rs#L273),
`GenericObservation::fill_from_supported_alleles`). A read arrives at the genotype
model as an allele number and a count, never as a sequence. There is no second opinion
available to it.

---

## 1. The truth set is right, and the scoring is right

Before blaming the caller, the cheaper explanations have to go: the truth set could be
wrong about the genotype, or the scoring could be matching the wrong pair of records.
Both were checked against the GIAB VCF directly, and both are out.

**HG002 chr1:243535155**, a 19-base insertion. GIAB calls it `1/1` with `GQ 314`. Its
`ADALL` field — the reads GIAB's own pipeline saw, pooled over every sequencing
platform it used, PacBio included — is `1,261`: one read on the reference against 261
on the insertion.

**HG003 chr15:96140584**, a 24-base insertion. GIAB calls it `1/1` with `GQ 404`, and
`ADALL` is `84,346`.

The second site is the more useful of the two, because GIAB's own counts there are not
clean either: 84 reads in 430 sit on the reference, one in five. **So reference-looking
reads at a long homozygous insertion are a real feature of short-read data, not
something ng invents.** What ng gets wrong is the size of the effect. Its count at that
site is `AD=13,13` — half its reads on the reference — where GIAB's pipeline gets one
in five, and freebayes, on exactly the reads ng was given, gets none.

---

## 2. One site in full

Take the two records the investigation started from. ng's output at 30×, where the
last field is `genotype : genotype quality : depth : reads per allele`:

```
HG002  chr1:243535155  T -> TATTTAAA , TATTTTAAAATATATTTAAA   QUAL 390.6  0/2:99:37:12,5,20
HG003  chr15:96140584  C -> CTGATTGGTCCATTTTACAGATGGT          QUAL 109.6  0/1:99:26:13,13
```

Now look at what BWA actually put at **chr15:96140584**. The 26 reads covering the
anchor fall into two groups, and the CIGAR says which: a read whose CIGAR contains
`24I` was placed as carrying the insertion, one that is all `M` was placed as
reference. `NM` is how many bases of that read disagree with the reference where it
was put.

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

**BWA's placement here is wrong, and the sequence says why.** The reference around the
anchor is an imperfect 24-base tandem repeat —
`CTGATTGGTCCATTTTACAGAGTG` `CTGATTGGTCCGTTTTACAGAGTG` `CTGATTGGTGCGTTTACAAACC` — and
the insertion is one more copy of that unit. Adding a copy to a repeat leaves the local
sequence nearly unchanged, so BWA can lay a read carrying the extra copy flat against
the reference at the cost of four or five mismatches, which is cheaper than opening a
24-base gap. That is what it did for the twelve reads in the first group.

**Comparing each read against the two chromosomes says which group is right.** The
`148M` reads starting at 96140461 to 96140491 sit four or five mismatches from the
reference and **one** from the insertion. They carry the insertion. Only four reads at
the site — the ones that stop within nine bases of the anchor, before the sequences
diverge — fit the two equally, and **none fits the reference better.**

ng puts 13 of the 26 on the reference.

### The same comparison at all 12 insertions

The comparison slides the read along each of the two chromosome sequences and counts
mismatches; no gaps, because the indel is already inside the sequence being matched
([`benchmarks/giab/src/read_haplotype_support.py`](../../../../benchmarks/giab/src/read_haplotype_support.py),
added with this report). "ng puts on REF" and "ng puts on ALT" are the two numbers of
its `AD` field.

| sample | site | bases inserted | ng puts on REF | ng puts on ALT | fits REF | fits ALT | fits both |
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

**How ng decides which allele a read supports.** When the walk passes a read over a
locus it keeps the bases that read showed across the locus's reference positions, and
nothing else. It then compares that string of bases, byte for byte, against each
allele's sequence; the allele it equals is the one the read is counted for. The
comparison against the reference allele is a plain slice equality,
`SequenceObservation::matches_reference`
([`locus_generation/mod.rs:416`](../../../../src/ng/locus_generation/mod.rs#L416)).

**How wide that stretch of reference is** is decided when the record is opened, by the
kind of event that opened it:

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

**An insertion gets one reference base, however long it is.** That is a true statement
about the *event*: an insertion adds sequence between two reference bases without
covering any of them, so the anchor really is the only reference base it touches.

But it is the wrong window for the *question ng is asking*, which is **did this read
come from a chromosome carrying the insertion?** A read that carries the insertion
shows the same base at the anchor as a read that does not. What separates them is what
comes after, and the comparison never reaches it.

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

### It grows with the length of the insertion

The table below takes every place where GIAB says the sample carries an insertion on
**both** copies of its chromosomes. At such a place no read should look like the
reference. The number in each cell is the fraction of reads at the site that the
caller nevertheless put on the reference allele — the first number of `AD` in its VCF.

Sites are grouped by how many bases the insertion adds. HG002, HG003 and HG004 are
pooled, at 30× coverage, counting only records that passed the QUAL ≥ 30 gate.
freebayes is run on the same sites and the same reads, and its column is read off its
own `AD`.

| bases the insertion adds | sites | ng | freebayes |
|---:|---:|---:|---:|
| 1 | 39 | 0.039 | 0.008 |
| 2–3 | 14 | 0.052 | 0.003 |
| 4–6 | 4 | 0.133 | 0.000 |
| 7–12 | 2 | 0.154 | 0.000 |
| 13 and over | 3 | 0.388 | 0.210 |

**Down the ng column**: at a one-base insertion, 4 reads in 100 are put on the
reference. At an insertion of 13 bases or more, 39 in 100 are. That is the prediction
the investigation set out to test, and it holds.

**Across to freebayes**: fewer than 1 read in 100 goes on the reference until the
insertion passes 13 bases. freebayes compares a read against a window of haplotype
rather than against one anchor base, so a read that carries the insertion disagrees
with the reference somewhere inside that window and is not counted as a reference read.
Past 13 bases a 150-base read stops being able to settle the question for either
caller, and freebayes climbs to 21 in 100 and gets 2 of those 3 sites wrong itself.

**This one quantity is the whole indel-genotype gap between the two callers.** Rerun
here from `genotype_disagreements.py` over all three samples at 30×: ng gets the
genotype wrong at **21 of the 297** truth indels it found, freebayes at **9 of 301**,
the production caller at **64 of 307**. *(The brief this investigation started from
gave ng's denominator as 276; 297 is what the documented command returns, and the
other two callers' figures reproduce exactly. It moves ng's error rate from 7.6% to
7.1% and changes nothing else.)*

### It does not grow with depth

The same measure, split by coverage instead of pooled. The number in brackets is how
many sites that cell averages over.

| inserted bases | 5× | 10× | 15× | 30× |
|---:|---:|---:|---:|---:|
| 1 | 0.043 (26) | 0.039 (35) | 0.037 (26) | 0.040 (41) |
| 2–3 | 0.075 (10) | 0.058 (15) | 0.022 (10) | 0.060 (15) |
| 4–6 | 0.250 (3) | 0.148 (3) | 0.158 (3) | 0.133 (4) |
| 7–12 | — | 0.222 (2) | 0.154 (2) | 0.154 (2) |
| 13 and over | 0.250 (1) | 0.353 (2) | 0.357 (2) | 0.388 (3) |

**Every row is flat.** More coverage does not dilute the wrongly-placed reads, because
it brings more of them in the same proportion. That is what a systematic
misclassification looks like, and it is what sampling noise does not: noise would
shrink as the square root of the depth. The number of wrong genotypes is flat too —
12, 17, 11, 17 from 5× to 30× — while the number of truth indels ng finds over the
same ground grows from 103 to 297.

### How few wrongly-placed reads it takes to change the genotype

Every truth-homozygous indel ng found, at all four depths — 363 of them — sorted by
what fraction of the site's reads ng put on the reference allele. The right-hand
column counts how many of those sites ng then failed to call homozygous.

| fraction of reads put on the reference | sites | not called 1/1 |
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

## 5. Deletions: the window grows with the variant, so the same defect stays small

**A deletion is compared over the bases it deletes, and that is the whole difference
from an insertion.**

An insertion's window is one reference base whatever its length, so the longer the
insertion, the more of the read's evidence falls outside the window and the more reads
are wrongly put on the reference. A deletion's window is the deleted bases plus one: a
seven-base deletion is compared over eight bases of reference. A read that really
carries the deletion has eight bases' worth of chance to disagree with the reference,
and it takes it. **The window grows with the deletion where it does not grow with the
insertion**, so the trend runs the other way and dies out:

| bases the deletion removes | ng puts this fraction of reads on the reference |
|---:|---:|
| 1 | 0.020 |
| 2–3 | 0.028 |
| 4 and over | 0.000 |

Those figures hold at every depth from 5× to 30×.

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
the 71 wrongly-placed reads in place, so insertions of 13 bases and over still put
about 12 reads in 100 on the reference, and some of them will still be genotyped
heterozygous. Those are the reads a 150-base read length cannot settle — at
`chr1:243535155` the insertion is a repeat unit copied into a `TTTTAAAATATA` tract, and
a read that stops 21 bases past the anchor matches the reference exactly on either
chromosome. **Nothing that compares sequence can recover those**, because there is no
sequence difference to find; they need a longer read.

**Neither shape realigns anything, and neither should.** The generic SNP/indel path
takes the mapper's placement as given — the only thing it changes is the leftmost
spelling of an indel the mapper already found
([`read/left_align.rs`](../../../../src/ng/read/left_align.rs): *"Bases and qualities
are copied through untouched; only the CIGAR changes"*) — and realignment stays on the
repeat-tract path, where every algorithm in `src/ng/alignment` lives (owner, 2026-09-11).
Both shapes keep that: they change only **how much of what the read already showed gets
compared**, at the position the mapper put it. The realignment in this report is a
measuring instrument — it establishes which haplotype each read fits independently of
any caller — and it belongs in the benchmark tree, not in the walk.

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
