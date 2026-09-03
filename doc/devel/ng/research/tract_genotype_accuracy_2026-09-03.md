# ng at repeat tracts: how we measure it, what moves it, and what to try next

**Date:** 2026-09-03. **Written to be read by somebody starting fresh**, so it repeats things the
last conversation took for granted. It covers three questions the owner asked for: how the numbers
are produced (§2, and **§3 is the trap** — read it before trusting any genotype figure), what has
been tried and what it was worth (§4, §5), and what to do next (§6).

**Branch** `ng-ssr-calling-loop`, worktree `/Users/jose/devel/pop_var_caller-ssr-calling-loop`,
nine commits ahead of `main`, working tree clean, 6,026 library tests green.

---

## 1. What ng is doing at a repeat tract, in one paragraph

ng is the from-scratch SNP/indel/STR caller in `src/ng/`. A **repeat tract** is a microsatellite —
a stretch like `ATATATATAT` — and ng routes those to their own path: candidate alleles chosen by
`select_ssr`, evidence shaped by `shape_ssr_locus`, and genotypes from the same calling loop the
SNP path uses but under a **stutter model** instead of a base-error model. Stutter (also called
slippage) is a read reporting a different number of repeat copies than the allele it was copied
from; it is the dominant error at a tract and everything here turns on how it is priced.

Every benchmark run of ng today is `--defaults`, which means **nothing is fitted**: the stutter
model is HipSTR's shipped constants, one pair of numbers for every class of tract, and the
genotype prior's length spectrum is flat. The run's own report says so — *"0 of 7 groups the file
says were fitted"*.

---

## 2. How the numbers are produced

### 2.1 The three grounds

| ground | what it is | why it exists |
|---|---|---|
| **`benchmarks/ssr_hg002`** — GIAB's HG002 tandem-repeat benchmark, 50,000 Tier intervals over 6.09 Mb, 36,497 assembly-based truth records, BAMs at 5×–300× | ng types **20,204 repeat tracts** in it and writes **6,351 tract records at 30×** | **This is the ground that can answer anything.** One sample only, so it says nothing about cohorts |
| **`benchmarks/giab/per_sample`** — the GIAB trio's 100 random confident intervals a sample, at 5×–300× | ng writes **149 tract records pooled over three samples at 30×** | Where every standing ng number was measured, so a figure here is comparable with C4's. **Far too small for anything about tracts** — its calibration table is a column of zeroes |
| **`examples/ng_tract_simulator`** — tracts whose genotypes we chose, sequenced under a slippage we set | exact truth, settable slippage, and the only place the fitted-versus-defaulted split has two sides | **It cannot measure candidate-selection losses**: its "truth allele never offered" count is 0 at every tract, because its alleles are drawn within three repeats of the reference and its reads carry nothing else |

**A trap in the third one.** The simulator says a wrong stutter model costs a lot (period-2+
genotype accuracy 0.932 at a true slippage of 0.25 against an assumed 0.10, rising to 0.990 when
the true model is supplied). **On real reads that effect is an order of magnitude smaller** (§4).
The simulator over-states the stutter lever because it contains no other failure — no mapping
error, no interruption, no allele the aligner mis-reads. A recommendation made from it alone was
wrong once already.

### 2.2 The instrument

[`benchmarks/lib/tract_qual_experiment.py`](../../../../benchmarks/lib/tract_qual_experiment.py)
takes a truth VCF, a caller's VCF, a confident-region BED and a **tract-ground BED**, and writes
three tables. Read its module docstring — it records why each rule is what it is.

- **`calibration.tsv`** — emitted records binned by QUAL against the share that really are at a
  variant tract. **Unit: the tract.** A record counts as truly variant when any truth record falls
  inside the tract it sits at.
- **`sweep.tsv`** — precision and recall as a QUAL threshold sweeps. **Unit: the allele**, on
  `benchmarks/giab/src/score_ng_recall.sh`'s rule (contig, position, REF and ALT equal after both
  sides are left-aligned and split), so its numbers are comparable with the standing ones. Carries
  `fp_with_no_called_copy` — false alleles listed in ALT that no sample was given a copy of, which
  is a different failure from a sample genotyped wrong.
- **`genotype.tsv`** — **the one that matters most, and §3 is about how to read it.** Where truth
  and caller both call a tract, are the caller's two alleles the truth's two. Reported twice: as
  the two tract **sequences** and as the two **repeat lengths**.

Everything is split by motif period — homopolymer against period 2 and above — because slippage
rises steeply as the period falls.

### 2.3 The tract ground, and which floors define it

The BED is built from `examples/ng_typed_region_dump <reference> <regions.bed> calling`, keeping
the `ssr_locus` rows as `chrom start end period`. **Pass `calling`, not `catalog`.** `catalog` is
the floors the catalog file was *stored* at; `calling` is what `call-from-alignments` actually
routes on, and the two differ by a lot. The dump's own doc comment still claims `catalog` is what
the run uses, which stopped being true when the routing floors changed.

Built ground lives at `tmp/tract_qual/ground/` — `tier.bed` (20,204 tracts) and `tier_sorted.bed`
(the confident regions). `benchmarks/lib/run_tract_qual_experiment.sh tandem_repeat_tier` rebuilds
both.

### 2.4 Running things

The benchmark data is gitignored and exists **only in the main checkout** at
`/Users/jose/devel/pop_var_caller/benchmarks`. From the worktree:

```
DEV_EXTRA_MOUNT=/Users/jose/devel/pop_var_caller/benchmarks ./scripts/dev.sh <command>
```

mounts it read-only. `dev.sh` forwards only `CARGO_*`/`RUST*` environment, so anything else must
be set **inside** the container command: `./scripts/dev.sh bash -lc 'FOO=bar ...'`. Every ng run
needs `NG_REFERENCE_CHECK=skip`. The container has `python3`, `bcftools` and `samtools` and no
`uv`; the host is the other way round.

- ng over the tandem-repeat ground: `benchmarks/ssr_hg002/src/run_ng_coverages.sh 30x 50x`
  (51 s and 61 s; writes `results/ng/HG002_<cov>.raw.vcf` and its `.parameters.toml`).
- One parameter setting, run and scored: `benchmarks/ssr_hg002/src/sweep_tract_parameters.sh`.

---

## 3. ⚠ The genotype-coding trap, and it has cost five wrong measurements

**This is the most important section here.** Repeat-tract genotype accuracy has now been reported
wrong five times, and every wrong version looked plausible. Four erred against the caller and the
fifth — the diagnosis in §3.4b — erred for it. If you change the comparison, read §3.4c first: it
is the standing check that catches this class of mistake.

### 3.1 What GIAB's truth does that ng does not

**A heterozygote carrying two different non-reference alleles is written as TWO records**, at the
same position, phased, one per haplotype:

```
chr1 4416244 . C CGTGT     ... GT=0|1
chr1 4416244 . C CGTGTGT   ... GT=1|0
```

ng writes the same genotype as **one multi-allelic record**:

```
chr1 4416244 . C CGTGT,CGTGTGT ... GT=1/2
```

These agree perfectly. **1,412 of 6,303 tracts on this benchmark are that shape.** Any comparison
that keeps one record a side, or that treats two records at a tract as two independent edits to
compose, gets them all wrong.

### 3.2 Four ways to get it wrong, each measured

1. **Compare allele strings as written, after `bcftools norm`.** Two records describing one event
   over different spans — `AGT → A` against `AGTGTGT → A,AGTGT` — share no string, because
   normalisation trims each record against *its own* ALT column, so a record with two ALTs trims
   less than one with a single ALT. **324 of 6,303 tracts read as genotype errors that were only a
   difference of spelling.** Reported accuracy 0.771 and 0.628, where the same callsets score
   0.877 and 0.867 under the comparison as it now stands.
2. **Pad both records to their union span.** Fixes that pair, and still fails wherever the two
   sides put their records at different places in the tract, because the span between them is in
   neither record's REF.
3. **Keep one record a tract a side.** Throws away 1,711 tracts of 6,303, 1,412 of them §3.1's
   shape — the class most worth measuring.
4. **Reconstruct haplotypes but refuse overlapping records.** Two truth records at the same
   position overlap as *edits* and do not as *haplotypes*: each copy takes a non-reference allele
   from exactly one of them. Refusing on spans alone scored 3,354 of 3,648 homopolymer tracts
   incomparable.

**And one more, in the plumbing:** applying the confident-region masks **before** left-alignment.
A truth record already written on an anchor base one base outside its interval is dropped, while
the query record describing the same event has not moved there yet, is kept, and then moves onto
the anchor with nothing to match. `chr1:69,233,430` is such a site — the truth carries
`TATAATAATA → T`, the Tier interval starts at 69,233,431, and ng's identical call scored a false
positive at QUAL 922. **Left-align first, mask after**, so both sides are treated alike. This is a
deliberate departure from `score_ng_recall.sh`, which masks first.

### 3.3 The rule that works

For each tract: take **all** of a side's records, read the tract's reference bases (one batched
`samtools faidx` for the whole ground), and for each haplotype copy lay on it only the records
whose allele for that copy is **non-reference**. Compare the resulting sequences as a multiset. A
side whose phase is stated (GIAB writes `0|1`) yields one pair; a side whose phase is open (ng
writes `0/1`) yields every assignment, and the two agree if any pair is shared.

**Records are collected from ten bases either side of the tract and the two rebuilt copies are
compared over the tract's own bases** — §3.4b for why those are two different numbers, and why
setting both to one base was wrong twice over.

**A record is reduced to the change it makes before any of this**, common suffix then common
prefix trimmed off, and then left-aligned against the reference. A record's written span is not
its change: `bcftools norm` trims a record carrying two alternate alleles less than one carrying a
single allele, so two files describing one event arrive anchored differently. That is what made
324 of 6,303 tracts read as genotype errors when allele strings were compared as written.

**7 tracts of 6,694 cannot be laid out** — a copy needing two overlapping changes at once — and
are counted apart and excluded from the denominator, because a tract the instrument cannot compare
is neither right nor wrong. Under the old rule it was 245 of 6,303, and most of those were tracts
the rule refused rather than tracts that are genuinely ambiguous.

**Seven hand-checkable shapes pin the implementation**, in the scorer itself:
`benchmarks/lib/tract_qual_experiment.py --self-test`, which needs no files and no tools and runs
in under a second. Two phased truth records against one multi-allelic call agree; a genuinely
different call does not, and lands in the counter naming what would have to change; one event
written over two spans agrees; a caller SNP one base past the tract's end leaves the tract right
and the neighbourhood wrong; a truth record five bases out is collected; left-alignment moves a
change written inside a repeat back to the repeat's start; and an event the two sides anchor on
opposite sides of the tract's edge disagrees on the tract and agrees on the neighbourhood. Each
was checked against a mutation of the code it covers.

### 3.4 Sequence or repeat length — settled, and it turned out to be a half-point question

Scored on the two tract **sequences**, base for base, ng is at 0.877/0.867 at 30×. Scored on the
two **repeat lengths** — which is how the STR field scores, and what ng, HipSTR and the existing
caller all emit as `REPCN` — it is at **0.881/0.873**.

**The 3-point gap this section used to report was the §3.4b defect, not a real difference between
spelling and counting.** With the comparison running one base past the tract's end, a substitution
in that one base broke the letter-for-letter score and left the length score alone; the gap was
2.9 points at homopolymers. Confined to the tract's own bases it is **0.4 points**.

**Settled 2026-09-03: the headline is the sequence, letter for letter, with the length beside
it.** Half a point is what the stricter question costs, and the stricter question is the one that
sees an interruption or a substitution *inside* a tract — the class §6.2's aligner cases sit in.
Both are in `genotype.tsv` (`genotype_accuracy` and `length_accuracy`), and a third column,
`neighbourhood_accuracy`, scores the tract with ten bases of flank.

### 3.4b ⚠ A fifth way it was wrong — two settings, wrong in opposite directions

**Found 2026-09-03 and fixed the same day.** Every number in §4 and §5 below has been re-scored
from the same callsets. **The diagnosis first written here was itself wrong, and the correction
is in §3.4c**, which is the fifth time this quantity has been reported wrong.

The comparison rebuilds each side's two chromosome copies from the reference and its records.
Two settings govern that, and until 2026-09-03 both were "the tract plus one base":

1. **How far from the tract records are collected.** One base was too few.
2. **How much DNA the two rebuilt copies are compared over.** One base was too many.

**Comparing too much charged a tract for a variant outside it.** At `chr1:9,955,404-9,955,422`, a
19-base poly-A, the truth writes `C → CAAA` at 9,955,403 and ng writes the identical record — plus
a SNP at 9,955,423, one base past the tract's end. ng's account of the tract is exactly the
truth's, and the tract was scored wrong. **46 of the 648 errors were that shape.**

**Collecting too little rebuilt the truth's own haplotype from a fragment of its own claim.** At
`chr1:150,329,038-150,329,047`, a 10-base poly-T, the truth describes the stretch with four
records — at 150,329,033, 35, 36 and 37 — and the one-base reach collected only the last. ng
writes it as two records that reconstruct the identical DNA, `AAACTTTTTTTTTTTTTTTTTGAGAT`, so ng
is exactly right and was scored wrong. **173 of the 648 wrong tracts have a truth record within
ten bases that the collection left out, against 208 of the 5,410 right ones** — at period 2 and
above, 38 tracts in 100 among the errors and 4 in 100 among the correct calls.

**The rule now.** Collect every record within ten bases of the tract, whichever path wrote it;
compare the two rebuilt copies over the tract's own bases. Two alternatives were measured and
rejected. Refusing any tract where a change reaches outside it raises the headline by 4.7 points
and **corrects not one verdict** — it refuses 503 tracts of which 57 in 100 were wrong against a
base rate of 11 in 100. Comparing over a real flank instead of the tract turns 62 tracts from
right to wrong for every 13 it turns from wrong to right at ten bases, and the ratio worsens as
the flank widens (33 against 6 at five bases, 119 against 26 at twenty) — because it charges the
tract for every variant the caller writes *near* it, which is `chr1:9,955,404` again further out.
That comparison is still reported, as `neighbourhood_accuracy`, because it is the only one that
can see a boundary-straddling variant ng genuinely misses.

### 3.4c ⚠ The diagnosis in §3.4b was wrong twice, and this time it erred *for* the caller

**§3.4b as first written named two tracts as instrument artefacts and inferred a size from them.
Both descriptions were wrong, and so was the size.**

- `chr1:14,722,151-14,722,162` was called "identical but for the base outside the tract". **ng is
  genuinely wrong there.** The truth's variant copy deletes the tract's last `A` and the `G` after
  it; ng's deletes two `A`s inside the tract and keeps the `G`. Both copies lose two bases, so the
  lengths agree and the sequences do not. No window rule makes them agree — all six were tried.
- `chr1:150,329,038-150,329,047` was called a variant the pad wrongly admitted. It is an
  instrument error, but by the other mechanism: three of the truth's four records were never
  collected (above).

**And the count.** §3.4b said §4's 407 "never offered" was overstated by up to 111 — the tracts
where the merge's table held the truth's spelling and selection kept it. That count reproduces
exactly (648 errors, 407 never-offered, 179 with the right length and the wrong spelling, 111 of
those kept as candidates). **The inference does not: at most 10 of the 111 turn right under any of
the six rules tried, and none under the closest to the old one.** They are real disagreements. A
kept candidate the record does not use is not a kept candidate the record cannot express.

**The standing rule this leaves.** Five measurements of tract genotype accuracy have now been
wrong; four erred against the caller and the fifth for it. Every rule that compares *more* tracts
gives a *lower* number, because the tracts a rule refuses are the ones ng gets wrong. **Treat any
proposed correction that raises the headline as suspect until it is shown to correct verdicts and
not merely to drop tracts** — the cross-tabulation in §3.4b is the check, and it is cheap.

### 3.5 Two defects in existing tooling, found and not fixed

- `benchmarks/ssr_hg002/src/ng_tract_candidate_recall.py` keeps a truth record only where its REF
  span reaches the tract's start, so an insertion anchored at `start − 1` — where left alignment
  puts every repeat-length gain — is dropped. 35 of the 268 cases have such a record, and 2,737 of
  the 19,613 clean tracts do. Widening the window takes the missing-sequence total from 434 to
  564, so the window and the overlap rule need reworking together.
- The same script's truth reconstruction misses records left-aligned before the tract; the
  slippage investigation had to work around it and reports 1,652 tracts called homozygous-reference
  that have a truth record within 30 bases.

---

## 4. Where ng stands, and what the errors are

**All re-scored 2026-09-03** with the corrected comparison (§3.4b), from the same callsets. The
run labelled "shipped" is `--defaults` **as it stood when the callset was written**, with the
repeat-tract outlier weight at 0.01; the shipped default has since moved to 0.05, whose row is in
§5.2.

**Genotype accuracy, comparable tracts, tandem-repeat benchmark:**

| | 30× homopolymer | 30× period 2+ | 50× homopolymer | 50× period 2+ |
|---|---:|---:|---:|---:|
| shipped at the time (outlier weight 0.01) | 0.8771 | 0.8665 | 0.8916 | 0.8767 |
| **fitted slippage + outlier weight 0.10** | **0.8823** | **0.8725** | **0.8954** | **0.8788** |
| shipped, scored on repeat length | 0.8812 | 0.8728 | 0.8954 | 0.8830 |
| both changes, scored on repeat length | **0.8852** | **0.8778** | **0.8977** | **0.8844** |
| shipped, scored over the tract and ten bases of flank | 0.8432 | 0.8597 | 0.8479 | 0.8655 |

Tracts scored at 30×: 3,864 homopolymer and 2,823 period 2+, with 7 the instrument cannot compare.
The old rule scored 3,515 and 2,543 and refused 245.

**The 852 errors at 30×, partitioned by what would have to change to fix them:**

| | homopolymer (475) | period 2+ (377) |
|---|---:|---:|
| a truth sequence was **never offered** as a candidate | **245** | **219** |
| called heterozygous, truth homozygous | 141 | 101 |
| called homozygous, truth heterozygous | 37 | 19 |
| wrong some other way, over a set that held the right alleles | 52 | 38 |

**Five errors in ten are decided before the model is consulted.** That is 464 of 852, against the
old rule's 407 of 648 — the class grew in count and shrank in share, because the corrected rule
scores 629 more tracts.

**⚠ The follow-through below has not been re-derived and should not be quoted.** It was produced
by `benchmarks/ssr_hg002/src/ng_tract_candidate_recall.py`, which §3.5 records two defects in, and
its 434 was keyed to the old comparison's 407. §3.4c is the specific warning: the one inference
that *was* re-derived — that 111 of these are an instrument artefact — turned out to be wrong.

<details><summary>Superseded: where the 434 missing sequences went under the old comparison</summary>

| | |
|---:|---|
| **268** | no read carried it, so the merge's allele table never held it |
| 61 | it cleared the support bar; the per-sample top-`ploidy` cut dropped it — **the only class a discovery round is aimed at** |
| 59 | the merge refused the tract, so no locus was built |
| 46 | the merge's table held it and the support bar refused it |

and the 268 split:

| | |
|---:|---|
| 121 | absent from the reads even at 300× with a median 115 spanning reads — not recoverable |
| **67** | **an alignment loss: reads carry the sequence and ng's table does not** (46 unambiguous at 30×) |
| 66 | every base the truth needs is in the reads, but no read spells the tract that way |
| 14 | the allele is longer than a 150 bp read can span with 20 bp of flank |

</details>

---

## 5. The levers: what is settable, what was tried, what it was worth

### 5.1 What can be changed today without touching code

A run's own `.parameters.toml` is written beside every output and can be edited and fed back with
`--parameters`. Two things in it matter here:

- **`repeat_tracts.slippage_by_stratum_and_group`** — per `(period, reference_repeats,
  slippage_group)`: `share_of_reads_that_slip`, `shorter_share`, `fall_off`. Empty at `--defaults`,
  which means HipSTR's shipped 0.10 / 0.50 / 0.05 everywhere.
  **Rows must be appended at the END of the file** with the empty
  `slippage_by_stratum_and_group = []` line deleted — an array-of-tables closes the table it sits
  in. Write them with `benchmarks/ssr_hg002/src/tract_slippage_rows.py`.
  **A row is needed for every repeat count a *candidate* can reach**, not just every reference
  length, because the lookup is by the candidate's count.
- **`stated_constants.repeat_tract_outlier_weight`** — 0.01, inherited from the existing caller
  and never measured. A value somebody typed must carry `warrant = "supplied"`.

`benchmarks/ssr_hg002/src/sweep_tract_parameters.sh` does all of this and scores the result. **Its
control is byte-identical to the `--defaults` run**; run it after any change to the script.

### 5.2 The sweep, on real reads, 30×

**Re-scored 2026-09-03** with the corrected comparison, from the same 22 callsets. The two error
columns are homopolymer only, as they were before.

| setting | homopolymer | period 2+ | het called for a hom truth | hom called for a het truth |
|---|---:|---:|---:|---:|
| the run's own defaults (outlier weight 0.01) | 0.8771 | 0.8665 | 141 | 37 |
| slip share 0.02 | 0.8718 | 0.8669 | 184 | 23 |
| slip share 0.05 | 0.8757 | 0.8682 | 156 | 30 |
| slip share 0.10 (shipped, supplied as rows) | 0.8766 | 0.8668 | 143 | 37 |
| slip share 0.15 | 0.8743 | 0.8612 | 133 | 52 |
| slip share 0.20 | 0.8694 | 0.8571 | 129 | 68 |
| slip share 0.30 | 0.8594 | 0.8477 | 110 | 108 |
| slip share 0.40 | 0.8301 | 0.8184 | 96 | 211 |
| shorter share 0.65 (shipped 0.50) | 0.8766 | 0.8670 | 139 | 40 |
| shorter share 0.80 | 0.8759 | 0.8663 | 142 | 41 |
| fall-off 0.15 (shipped 0.05) | 0.8773 | 0.8668 | 143 | 37 |
| fall-off 0.30 | 0.8771 | 0.8661 | 147 | 35 |
| slip rising with tract length, best of three shapes | 0.8775 | 0.8682 | 139 | 38 |
| outlier weight 0.05 (the shipped default since 2026-09-03) | 0.8796 | 0.8692 | 129 | 40 |
| **outlier weight 0.10** | **0.8808** | **0.8692** | **124** | **41** |
| outlier weight 0.20 | 0.8806 | 0.8685 | 118 | 48 |
| outlier weight 0.30 | 0.8802 | 0.8677 | 108 | 56 |
| fitted per-stratum slippage | 0.8789 | 0.8687 | 147 | 31 |
| **fitted slippage + outlier weight 0.10** | **0.8823** | **0.8725** | **131** | **33** |

**What worked. Every one of these is smaller than the earlier scoring made it, and the ordering
is unchanged.**

- **The outlier weight, 0.01 → 0.10** — the largest single move, +0.37 points at homopolymers and
  +0.27 at period 2 and above. It is the share of reads at a tract the stutter model cannot
  explain; `λ · U` is a floor under every read's emission, so the number is really **a bound on how
  far one read may pull a genotype** — the job freebayes does with a read-dependence factor and
  GATK with a Phred-45 cap, and which ng has nothing else doing. It removes 17 spurious
  heterozygotes of 141 at homopolymers and costs 4 collapsed ones.
- **A per-stratum slippage fit**, +0.18 alone, and it fixes a *different* error (collapsed
  heterozygotes 37 → 31), which is why the two compose to +0.52 rather than to +0.37.

**What did not work, and this contradicts the recommendation that started the work:**

- **Any flat change to any of the three stutter numbers.** Over a twenty-fold range of the slip
  share, homopolymer accuracy moves half a point while the two error classes swing two-fold and
  nine-fold in **opposite** directions. It is a dial that trades one error for the other, and the
  best flat setting (0.05, at 0.8757/0.8682) is within 0.1 points of the shipped 0.10 at
  homopolymers and 0.14 above it at period 2 and above. `shorter_share` and `fall_off` do nothing
  at all — every setting is within 0.1 points of the shipped one.
- **A slippage that rises with tract length**, hand-set: +0.04 at best over the run's own defaults.

### 5.3 What the six investigations found

Raw reports in [`../../reports/tract_genotype_investigation/`](../../reports/tract_genotype_investigation/),
with a README saying which claims a run corroborated.

**⚠ All six were run under the old comparison.** Where a verdict counts *tracts scored wrong*, its
denominator has moved: 852 errors at 30× rather than 648, and 242 spurious heterozygotes rather
than 137. Where a verdict counts *reads* or *candidate table entries* it is untouched. The column
below says which, and the HipSTR row has been re-scored because it is the only outside bar and the
correction changed its direction.

| investigation | verdict |
|---|---|
| **slippage** — fit it from HG002's reads | Real slip share is **0.0039 at an 8-repeat homopolymer rising to 0.088 at 30**, 0.027 tract-weighted: the shipped 0.10 is **too high by 3.7-fold**, and is the value at the very top of the measured range. `shorter_share` is 0.58–0.73 against the shipped 0.50 and `fall_off` 0.16–0.45 against 0.05 — both too **low**. Produced the rows the winning run used. |
| **readmodel** — every inherited constant | The outlier weight is the only one that changes calls. The part-repeat shares can be set to anything, zero included, and change nothing. The slip-size cutoffs change at most 1 call in 8,686. The per-base substitution rate inside a tract is measurably wrong (reads give 1 mismatch in 756 at homopolymers against the model's 1 in 3,300) and worth 2–3 calls in 8,686. **Its doc comment's premise is false**: base quality inside a tract is 35.3 against 34.0 outside — better, not worse. |
| **hetfhom** — the spurious heterozygotes | Counted under the old comparison, so **the shares hold and the counts do not**: 110 of *its* 137 are a length change, 78 of those exactly one repeat, concentrated at long tracts — 7 in 100 at homopolymers of 25+ repeats against 1 in 100 below 13. **The class does not shrink with depth, and that survives re-scoring**: 242 at 30× against 236 at 50×, a fall of 2 in 100, while total errors fall from 852 to 774, 9 in 100. An error more reads do not buy down is a model error. Every proposed filter costs more than it buys (below), on the old counts. |
| **prior** — the flat length spectrum | The truth is nothing like flat — **79 chromosomes in 100 sit at the reference length** against a flat shape's 11 in 100, and it is strongly stratified (0.97 at a 6–8-base homopolymer, 0.51 at 21+). **But fitting it reaches 10 of the old comparison's 648 errors and risks 77 correct calls**, because the dominant error is a homozygote called heterozygous and a reference-peaked prior makes that worse. **Do not fit it.** Half A answers `population_diversity.md` §4.4's open question: a per-period pooled spectrum is wrong at both ends; it would have to be per stratum. |
| **hipstr** — the outside bar | **Re-scored 2026-09-03**, and the correction moved it: on the period-2+ tracts both reach, ng **0.8998** against HipSTR **0.8806** at 30× and ng 0.9096 against HipSTR 0.9102 at 50×. Under the old comparison ng was level at 30× (0.928 against 0.925) and 0.9 points behind at 50× (0.939 against 0.948) — **ng was the arm the old rule penalised, because ng also writes SNP-path records beside a tract and HipSTR writes none, so only ng was charged for its own neighbours.** HipSTR's own median fitted slip level is **0.04 against ng's fixed 0.05** — the median locus does not need fitting; only its top decile does. Still uncontrolled: different catalogs, candidate rules, priors, `--min-reads 5`. **HipSTR's region file holds no period-1 loci**, so ng's homopolymer numbers have no comparator. HipSTR writes `.` in QUAL, so no calibration or sweep number for it is producible. |
| **unseen** — the 268 no read carried | The split §4 now carries as superseded, three alignment losses verifiable by hand, and the observation that 174 of 268 are at the reference length. **Its three hand-verified alignment losses stand** — they were read off the reads and the candidate table, not off the genotype comparison — and they are §6.2's starting point. |

**Ruled out with numbers** (from `hetfhom`):

- **A stricter candidate bar** — 2 reads and 15% removes 7 homopolymer errors and destroys 10
  correct calls; 3 reads and 10% removes 11 and destroys 50; the ladder saturates at 37 of 137
  removed for **617 true alleles lost**.
- **A GQ floor** — GQ 30 withdraws 43 of 86 wrong homopolymer calls and 413 of 3,113 correct ones,
  and the bad become no-calls rather than right calls.
- **An allele-balance collapse rule** — best net anywhere on the curve is +2 tracts in 3,515.
- **Lowering the support bar to feed discovery** — 10% to 5% supplies the truth to 2 more of the
  242 never-offered tracts and hands 137 more tracts a candidate the truth does not carry.

---

## 6. What I would do next, and why

Ordered by expected value, with what I am confident about and what I am guessing at marked.

### 6.1 ✅ Done — the comparison is fixed and the scoring convention settled

The two open questions this section used to hold are answered. The comparison now collects records
from ten bases out and compares over the tract's own bases (§3.4b), and the headline is the
sequence letter for letter with the repeat length beside it (§3.4). Every number in §4 and §5 has
been re-derived from the same callsets, and seven hand-checkable shapes pin the scorer
(`--self-test`).

**What is still owed here**: the follow-through in §4 — where the missing sequences go — has not
been re-derived, and `benchmarks/ssr_hg002/src/ng_tract_candidate_recall.py`, which produced it,
carries the two defects §3.5 records. Fixing that script and re-running it is the next measurement,
because §6.2 and §6.5 are both sized from its output.

### 6.2 The tract aligner: three hand-verified losses, and a count that needs re-deriving

**Confident this is a real defect.** Three cases with coordinates:

- `chr3:33,877,690` — an 11 bp poly-A, heterozygous. 10 of 23 reads at 30× carry `CAAAAAAAAAA`;
  ng's table holds `AAAAAAAAAAA` and a 10-base run. **The leading `C` is dropped and the tract is
  reported one repeat short.**
- `chr3:37,126,860` — homozygous non-reference. 10 of 11 reads at 30×, and 120 of 138 at 300×,
  spell `AAAAAAAAAAAAAGAAA`; ng's table holds a bare 13-base A-run. **The `GAAA` tail inside the
  tract span is discarded.**
- `chr11:37,147,255` — same mechanism, and the only 17-base allele ng holds is the one-read
  sequencing-error spelling rather than the true one carried by 12 reads of 14.

**The three cases stand** — they were read off the reads and the candidate table, not off the
genotype comparison, so the 2026-09-03 correction does not touch them. **The 46 does not**: it came
from the follow-through §4 now marks superseded. The mechanism is legible: something at the tract's
edge or an interruption inside it is being discarded by the tract realigner
(`src/ng/alignment/ssr_best_path_*`). It is the only large class in §4 that is a defect rather than
a limit of the data. **Start by reproducing one of the three by hand** — pull the reads with
`samtools view`, run the tract locus generator over that interval, and see where the sequence is
lost.

### 6.3 Set the outlier weight deliberately, per period

**Confident in the direction, not the value.** The shipped default moved from 0.01 to 0.05 on
2026-09-03; re-scored, 0.05 is worth +0.25 points at homopolymers and +0.27 at period 2 and above
over 0.01, and 0.10 is worth +0.37 and +0.27. Two things say the remaining choice needs care: read
literally as *the share of reads nothing explains*, the measured value is **1 in 2,300 at
homopolymers and 1 in 209 at period 2+** — the opposite ordering to what the sweep wants — and
raising it trades one error for another, 141 spurious heterozygotes at homopolymers falling to 108
at weight 0.30 while collapsed heterozygotes rise from 37 to 56. **So sweep it per period, and
check it at 3 reads a position on the tomato cohort before moving the shipped default again**,
because a floor under every emission behaves very differently when there are three reads rather
than thirty. That is the project's stated range commitment and nothing here has tested it.

### 6.4 The fit-mode command

**Confident it is worth building; less confident it is next.** §5.2 says a per-stratum fit is worth
+0.18 points, and it is a fit for *this sample and this chemistry*, so it cannot ship as a default —
it argues for the deferred fit-mode command (spec `calling_loop_ssr.md` §3.4) that produces a
parameters file per run. The machinery to read one back works. **But +0.18 points is small**, and
HipSTR — which does exactly this per locus — is **1.9 points behind ng at 30× and level at 50×**
on the tracts both reach (§5.3), so do not expect more from fitting than the sweep measured.

### 6.5 Milestone E (allele discovery): I would leave it parked

**Confident.** The plan's Milestone E builds a discovery round that admits tract lengths hiding
under stutter. Its decision half is built, tested and committed (`src/ng/calling/inference/discovery.rs`,
12 tests, five mutations each predicted and each caught). **But it is aimed at 61 of 434 missing
sequences** — a count from the superseded follow-through in §4, so re-derive it before acting on
it — and the class it can only enlarge, a heterozygote called for a homozygote, is now **242** at
30x rather than 137. Two further things came out of building it and should be recorded before anyone
resumes:

- **ng's retrace needs no posteriors.** The tract locus generator already realigns every read, so
  an observation's bases *are* the sequence HipSTR's alignment retrace would imply. The eligible
  set is a function of the observations and the candidate table alone, and a test asserts that a
  second round over the same evidence admits nothing. **So discovery is a pre-pass, not a round
  wrapped around the loop** — no second convergence, and no append-only emission store is needed.
  That contradicts the premise of spec §4.1's decision to "look against the converged posteriors",
  and the spec should be amended rather than the code bent to it.
- **The remaining wiring belongs inside `select_ssr`, not the loop.** A discovered allele's reads
  have to move out of the record's "no written allele explains these" column and into its `AD`,
  and that bookkeeping — the merge-index-to-candidate map and the per-sample leftover — is built
  inside selection.

### 6.6 Things I suspect and have not tested

**These are guesses. Treat them as hypotheses with a stated test, not as findings.**

- **The one-repeat spurious heterozygote at long tracts may be a per-tract problem, not a
  per-stratum one.** 78 of 137 are exactly one repeat away and they concentrate at long tracts,
  which is where a *locus's own* slippage differs most from its stratum's. The spec already has a
  per-locus slippage re-fit designed and unbuilt (`calling_em_loop.md` §5.1, three pull-back
  settings). **Test:** score those 137 tracts alone under a re-fit and see whether the class moves;
  if the per-stratum fit already got most of it, drop the idea.
- **The outlier term's uniform shape may matter more than its weight.** `λ · U` spreads the junk
  mass evenly over every reachable length; real junk is probably concentrated near the called
  allele. **Test:** replace `U` with a shape falling off with distance and re-run the same sweep.
  If the weight's optimum moves toward the measured 1-in-209, the shape was the problem.
- **`bcftools norm` may be doing damage nobody has checked.** Every number in §4 passes both sides
  through it. **Test:** score once with a comparison that never normalises (§3.3's reconstruction
  does not need it) and see whether anything moves.
- **The 66 "ambiguous" cases in §4 — every base the truth needs is in the reads, but no read spells
  the tract that way — smell like the same aligner defect as the 46.** If they are, §6.2's prize is
  112 tracts rather than 46.

### 6.7 What nothing here has touched

**One sample, one chemistry, two depths.** Every number is HG002 at 30× and 50×. Nothing says what
any of it does on the 63-accession tomato cohort at three reads a position, where the outlier
weight's floor, the candidate bar and the prior's concentration all bind differently — and that is
half of the range the caller is committed to (`CLAUDE.md`, `design_principles.md` §0).

---

## 7. State of the plan and the tree

**Plan:** [`doc/devel/ng/impl_plan/calling_loop_ssr.md`](../impl_plan/calling_loop_ssr.md).
A1–D2 ✅, E1–E3 ☐. **Checkpoint D is a hard pause the owner has not yet cleared** — the QUAL report
exists and `calling_quality_ssr.md`, which takes the decision from it, is unwritten.

**Commits on `ng-ssr-calling-loop`:**

| | |
|---|---|
| `314cb3da` | D1 — the QUAL experiment's instrument |
| `8eef9887` | D2 — the QUAL experiment and Checkpoint D |
| `aafc7d54` | the genotypes at a tract (**its numbers were wrong; superseded**) |
| `e20abef0` | E1 — what a discovery round would admit (decision half only) |
| `2efa47da` | the genotype measurement corrected — 0.886, not 0.771 (**also superseded**) |
| `65df5c24` | two parameter values lift the tract genotypes by 0.6 points (**numbers superseded**) |
| `073b678c` | the repeat-tract outlier weight ships at 0.05, not 0.01 |
| `f5bd2d41` | freebayes does not cap a read's evidence, and the outlier weight's doc said it did |
| `5939608a` | a fifth way the tract genotype comparison is wrong (**its diagnosis was wrong; see §3.4c**) |

**The QUAL half of the instrument was not touched.** The calibration and sweep tables re-derive
byte-identically from the same callsets — 68 and 52 rows — so every QUAL number below and in the
Milestone D report stands as measured.

**Reports:**
[`ng_tract_qual_experiment_2026-09-02.md`](../../reports/ng_tract_qual_experiment_2026-09-02.md) —
the QUAL answer, and §5 the genotypes;
[`ng_tract_genotype_improvement_2026-09-02.md`](../../reports/ng_tract_genotype_improvement_2026-09-02.md) —
the sweep and what moved;
[`tract_genotype_investigation/`](../../reports/tract_genotype_investigation/) — the six raw
investigations.

**The QUAL answer, for completeness**, since it is the milestone that is formally open: ng's tract
QUAL is as well calibrated as its SNP QUAL on the same benchmark (above QUAL 200, wrong 5 times in
2,882 against ordinary sequence's 19 in 10,272) and is **not** a usable gate — at a homopolymer
tract precision peaks at 0.850 at QUAL 50 and falls back to 0.831 by QUAL 200 having shed more
than half the recall.
