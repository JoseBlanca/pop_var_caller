# ng — a deletion across an STR segment's edge, called as one allele of one record

*Design, 2026-10-10; §6 amended 2026-10-29 so that pvc writes no `*` allele. The direction is the
owner's decision of 2026-10-10; §4–§8 are the implementing session's answers to the questions the
draft left open, awaiting the owner's confirmation. Reads on:
[`locus_generation_ssr.md`](locus_generation_ssr.md) §3 (how a read is measured across a segment,
including the edge-deletion rule this builds on), [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md)
§5 and §7, [`read_likelihoods.md`](read_likelihoods.md) §4, and the measurements in
[`../../reports/implementations/str_accuracy_2026-10-10.md`](../../reports/implementations/str_accuracy_2026-10-10.md).*

---

## 1. The problem

**An STR segment and the sequence beside it are called by two different paths.** The STR path
calls the repeat from each read's spelling of it; the generic path calls the flanking sequence
position by position and stops at the segment's edge. A deletion that runs from the repeat into a
flank is therefore cut in two at the edge.

**Since commit c0b55125 both halves are usually called**: a read whose own alignment carries a
deletion across the edge counts as crossing the segment, so the repeat's half becomes an STR
candidate, and the generic path calls the flank's half. Written as two records, they spell the
sample's sequence. What is left:

- **Accuracy.** On the stuttering benchmark at one sample, main gets 16–17 of the 27 pairs carrying
  such a deletion right at 15× and above, where GATK gets 26–27. At 50 samples (PCR, 15×), 347 of
  1,085 carrier pairs are still not right; 289 of them are at alleles whose sequence no combination
  of the written records can spell.
- **Phase.** A sample heterozygous at both records may carry the two halves on one chromosome (one
  long deletion) or on two (two short ones). Two records do not say which.
- **Representation.** Where a truth set writes the deletion as one record, an allele-exact
  comparison counts the two records as a false positive and a miss although the sequence is right:
  on HG002's tandem-repeat benchmark about 70 such records a depth after c0b55125.

GATK gets these right by calling one region holding the flank and the repeat, scoring each read
against whole haplotypes. The aim here is the same result on ng's two-path design.

## 2. The decision

**A read that crosses a segment and carries a deletion into a flank records, beside its spelling of
the repeat, how many flank bases the deletion removed on each side. The STR path's alleles become
the repeat's spelling plus that flank deletion, so the whole deletion is one allele, phased by the
reads. When the calls are gathered — before the filters and the VCF — a record whose alleles include
a flank deletion has its reference widened over those flank bases, and the generic records that
overlap it are merged into it, one record whose alleles are haplotypes (§6).**

Every other read is measured exactly as today. Where no sample carries a flank deletion, nothing
changes: the generic records stay in the flanks and the STR record stays on the repeat, with the
loci as they are defined now.

## 3. Why not the alternatives

**A fixed flank margin in every segment's spelling.** Proposed first (100 bases each side), and
rejected by measurement: the STR path calls a segment only from reads that cross its whole window,
and with 150-base reads every base of margin costs crossing reads. On the simulated sample at 15×
(3,000 segments):

| flank added each side | median crossing reads a segment | segments with fewer than 3 |
|---|---|---|
| 0 (today) | 10 | 28 |
| 10 | 8 | 98 |
| 20 | 6 | 342 |
| 50 | 0 | 2,660 |
| 100 | 0 | 3,000 |

Too few crossing reads is already the largest loss left at one sample, so any fixed margin trades
the main weakness for this one. Carrying flank bases only on the reads that have a flank deletion
costs no read.

**Combining the generic and STR loci after both are built.** Rejected because the evidence cannot
be joined there: the STR path stores per-read spellings of the repeat, the generic path stores
per-position counts, and nothing at calling time says which reads carry both halves. A combined
locus could only multiply the two independent answers, which is what two records already are, and
could not recover phase.

## 4. The flank deletion, and how two reads of one sequence get one key

**A flank deletion is two numbers: how many flank bases are missing just left of the segment, and
how many just right.** It belongs beside a spelling of the repeat, and together they say what the
read's sequence is over a stretch wider than the segment: with `L` bases missing on the left and
`R` on the right, the read's sequence over the reference from `L` bases before the segment to `R`
bases after it is the spelling, nothing more. `(0, 0)` is every observation today.

**The same sequence can be cut into a spelling and a flank deletion in more than one way**, because
the mapper chooses where to put a deletion inside a run of equal bases. If the flank base beside
the segment is an `A` and the repeat starts with an `A`, a deletion of that flank base and of an
`A` of the repeat is the same sequence. Two reads of one chromosome can then come with different
pairs, and they must not become two alleles.

**So every pair is put in its smallest form before it is stored.** While the left count is above
zero and the spelling begins with the first of the deleted flank bases, that base is given back to
the flank: the spelling loses its first base and the left count drops by one — the same sequence,
with the deletion placed one base further right. The right side is the mirror image: while the
spelling ends with the last deleted flank base. The result does not depend on where the mapper put
the deletion, and it has a property the design leans on throughout: **a deletion that can be placed
wholly inside the segment comes out as `(0, 0)`** — a plain repeat allele, measured and called
exactly as today. Stutter, which adds or removes whole units inside the repeat, is always such a
deletion. Only a sequence whose flank really differs keeps a flank deletion.

**The smallest form is the allele's identity; it is not where the generic path puts the
deletion.** The generic path left-aligns every read before its pileup (`read/left_align.rs`, a
port of GATK's `leftAlignIndels`), which slides a deletion as far *left* as it goes — into the
left flank wherever the flank's last bases repeat what the deletion removed. `AG` deleted from
`…CG|AGAGAG…` is a plain one-unit contraction in its smallest form, and the generic path, sliding
it left, writes the flank's `G` as deleted. So where the calls are gathered (§5.8, §6), each
allele's flank extent is computed again from its sequence, as the generic path would place it,
and the record is widened over that.

*Changed during implementation (2026-10-10).* The first build slid every read's deletion left at
measurement, so that the counted flank bases were the generic path's. That made every stutter
contraction at a segment like `…CG|AGAG…` an allele with a flank deletion, no longer explained by
the plain alleles' stutter. It was replaced before any measurement used it.

**Tried and not adopted: counting a flank base the read shows substituted, right against the
repeat, as part of the extent** (2026-10-11). The mapper can write one chromosome two ways — a
deletion taking the flank's first base is also, base for base, a deletion inside the repeat with
that flank base substituted — and the second spelling reads as a plain repeat allele, a second key
for the same chromosome. On one sample at 15× (PCR, segment `seg02581`) that split six reads four
to two and the genotype called the two spellings as the sample's two chromosomes. Taking the
substituted flank bases into the spelling gave both spellings one key, and on all 16 runs it
changed the edge-deletion pairs right by between −2 and +3 a run while adding up to 7 false calls
a run at 50 samples (PCR, 5×: 1,448 against 1,441). The split is real but rare; the measurement
does not pay for the extra observations every flank sequencing error then makes.

On the stuttering benchmark (50 samples, PCR, 15×) the truth has 97 alleles of the kind
`deletion_across_edge`. 74 are a deletion and nothing else, each taking 1 to 8 flank bases, all
from one side; the other 23 also change the repeat's length.

## 5. What changes, path by path

### 5.1 Read measurement

**A read spelled from its own alignment records the flank deletion of every deletion in that
alignment that crosses a segment edge**, in its smallest form (§4). That covers both kinds of read
spelled this way: the reads the rule of c0b55125 makes complete, and the complete reads the
delimiter anchors whose own alignment is used for the spelling (`locus_generation_ssr.md` §3).
Every other read — partial reads, reads spelled by the delimiter, reads with no crossing deletion —
records `(0, 0)`.

**A deletion that touches the edge counts as crossing it** (amended 2026-10-18): one lying in the
flank and ending at the segment's first base, or starting just after its last, records its flank
bases as one that crosses the edge does, and makes the read complete by the same rule. The case is
a misread base, not a variant the mapper split. Where a run of one base is followed by a different
base, a read that misreads that base as the run's shows a run one base too long, and minimap2
writes it as the run at full length and the neighbouring base deleted — a deletion that ends at the
edge without crossing it. Before the amendment such a read was spelled over the segment's positions
alone and counted as the reference length: at HG002's homozygous one-base-deletion homopolymers at
30×, 565 of the 926 reads at the reference length carried such a deletion
(`../../reports/ng_str_edge_touching_2026-10-18.md`). Recorded, the read keeps
its flank deletion and is compared with each allele over the wider extent (§5.6): against the
deletion allele it differs by the misread base alone, and is priced at the substitution rate.

*This is not the variant §4 tried and did not adopt.* That counted a flank base the read shows
**substituted** as part of the extent, which turned every flank sequencing error into an extra
observation; this counts only a base the read's alignment **deletes**, which a substitution error
does not produce.

**A flank deletion longer than the flank stored with the segment is not recorded**: the read is
measured as today. The record can only be widened over flank bases the locus holds.

**Two reads with the same spelling and different flank deletions are different observations**:
the per-sample tally keys on spelling, flank deletion, how much of the segment the read saw, and
read group.

### 5.2 The psp file

**The flank deletions travel as one new field at the end of each record**: a list of
`(observation, left, right)` for the observations whose flank deletion is not `(0, 0)`, empty
otherwise — one byte a record in the usual case. The psp format goes from 1.1 to 1.2.

A field at the end of the record is the one change the format already allows without breaking
readers (`psp_file_format.md` §4.3, the precedent of 1.1): a reader of 1.1 skips it and calls as
today, and a reader of 1.2 given a 1.1 file reads every flank deletion as `(0, 0)`, which is
exactly what a 1.1 file measured. A new field per observation would have cost a byte on every
generic observation too, and a major version.

### 5.3 The parameter census

**Observations with a flank deletion are left out of the stutter census.** They are a different
sequence, not a slip of the polymerase, and counting them as a repeat-length change would charge
them to stutter. (The fit is not run in this work; this keeps it right for when it is.)

### 5.4 Joining the samples' alleles

**The cohort's list of alleles at a segment keys on the pair** (spelling, flank deletion), where
today it keys on the spelling alone. The reference is `(reference tract, (0, 0))` as before.

### 5.5 Candidates

- **Selection keys on the pair throughout.** A rung of the repeat-count ladder is still the
  spelling's repeat count; two alleles on one rung that differ only in their flank deletion are
  two sequences on that rung, as two spellings are today.
- **A flank deletion earns a candidate under the same bar as any sequence**: reads from at least
  two crossing reads of one sample, and at least a tenth of that sample's crossing reads
  (`candidate_alleles_ssr.md` §5). The smallest form (§4) is what makes the bar safe to reuse:
  stutter, which adds or removes whole units inside the repeat, can always be placed inside the
  segment and so never produces a flank deletion. A flank deletion seen by one read is not
  enough, as no sequence seen by one read is (the one-read candidate was measured and rejected on
  2026-10-10: at 5×, one sample, 216 newly right against 252 newly false).
- **A flank-deletion allele that is one misread base from another allele, and a small minority
  beside it, is not a candidate** (2026-10-11). Written out over the wider of the two extents,
  it differs from the reference or an admitted plain allele by one substituted base, and in every
  sample that shows it, it has under a third of the two alleles' reads together. That is the
  shape of a systematic read error beside a run: on HG002 at 30×, 4 reads against 13 showed the
  `G` after 25 `T`s misread, earned a candidate, and made a false heterozygote. Its reads still
  count for the allele they are one base from, at the substitution rate. Without the minority
  test the rule also removed the commonest real edge deletion on the stuttering benchmark — one
  flank base, which written out is one base from a plain contraction — and lost 36 to 77 carrier
  pairs a run at 50 samples. The third is a judgement between a misreading's share (under a
  quarter there) and a heterozygous carrier's (about a half), not a fitted number.
- **The periodicity check never counts a read with a flank deletion as off the grid** — it is the
  read that check's edge-deletion exception was written for, now known exactly instead of guessed
  from the spelling (`candidate_alleles_ssr.md` §7). *(The check itself was removed on 2026-10-22,
  owner, so nothing is counted off the grid any more.)*

### 5.6 Read likelihood

**A read and an allele over the same flank extents are scored as today, on their spellings.**
Where the extents differ, both are respelled over the wider of the two, each adding the reference
flank bases it leaves out:

- **if they come out the same length, they differ only by substitutions**, and today's emission
  prices them at the substitution rate — a flank base read wrong, or one allele showing a flank
  substitution the other has not;
- **if they do not, one has flank bases deleted that the other has not.** That is not stutter,
  which changes whole units inside the repeat; if it is anything the read model knows it is a
  sequencing indel across the edge. The repeat spellings are scored as ever, times the model's own
  price for an indel that is not whole units (`read_likelihoods.md` §4.2, the part-repeat
  regime), sized by the flank bases. A read of an allele that lost a flank base to such an error
  then still counts for that allele, at that price, instead of for nothing: on HG002 at 30×, two
  reads of an 11-`C` allele missing the base before the run had stopped counting for it, and the
  allele lost to the reference. *(The first build scored this case zero.)*

*Changed during implementation (2026-10-11).* The first build scored every difference in extent
zero. A deletion taking a flank base can respell to the same length as a plain allele — `GTT` ×11
then `G`, against `GTT` ×11 then the flank's `A` — and such a pair differs by one base, which the
substitution rate prices, not by a deletion. The measured build is the one that prices it.

**Reads that end inside the segment: compared with the allele and the flank it leaves.**
*(AMENDED 2026-10-21, owner. This paragraph said such reads are "scored as today", compared with
the allele's spelling alone, and named the change below as a refinement left out of the first build.
Measured since: [`../../reports/ng_generic_indel_error_2026-10-21.md`](../../reports/ng_generic_indel_error_2026-10-21.md)
§5.2.)*

A read that ends inside the segment — a *partial read* — held one flank and ran towards the other
until it stopped. It records no flank deletion of its own: it never reached the far flank. **When the
allele it is compared with deletes flank bases on the side the read runs towards, the read is
compared with the allele's spelling followed by the flank that allele leaves on that side.** Every
other partial read is compared with the allele's spelling as before.

**Why: without the flank, a homozygous edge deletion is called heterozygous.** An allele that deletes
flank bases has a spelling shorter than the reference tract, so a partial read that stops near where
the reference tract ends has, under that allele, already run past the allele's tract and into the
flank beyond the deletion. Compared with the spelling alone, those extra bases can only be explained
as a longer tract — slippage — and are priced as such. Measured on HG002 chromosome 20 at 60×:

- at chr20:53879889 (truth `GAGACAGAC>G` at 53879901, `1/1`; tract `GAGA…GA`, 16 bases; flank
  `CAGACAGAC…`) the allele is `GAGAGAGAGAGAG` with 5 right-flank bases deleted. Seven partial reads
  held left and spelled `GAGAGAGAGAGAGAGA`. Under the allele that is its 13 bases followed by `AGA`,
  the start of the flank it leaves — the allele and the reference read the same there, so the reads
  cannot tell them apart. Compared with the spelling alone they scored 0.003 each under the allele
  against about 1 under the reference; reads a base or two longer scored under one in a million;
- those dozen partial reads outweighed 36 complete reads, each 0.87 likely under the allele and
  3 × 10⁻¹¹ under the reference, and the site was written `0/1`. chr20:34924682 (`TTTG>T`) and
  58447656 (`GTGTGTGCA>G`) fail the same way.

**The rule, precisely.** For a partial observation that held the left flank and a candidate `a` whose
flank deletion takes `R > 0` bases on the right, let `F` be the reference right flank without its
first `R` bases. Then

```text
Lr(partial of ℓ bases | a)  =  Σ   P_stutter(Δ | a) · subst( the ℓ bases | first ℓ bases of ( a stretched by Δ, then F ) )
                               Δ : (L + Δ) + |F| ≥ ℓ
```

summed over every length change `Δ` the stutter distribution reaches, where `L` is the length of
`a`'s spelling. Two things differ from [`read_likelihoods.md`](read_likelihoods.md) §5.2's form:

- **the condition is on the stretched spelling plus `F`, not on the stretched spelling alone.** A
  stretching shorter than the read is now admissible — the flank supplies the rest — and it is the
  stretching that explains these reads (`Δ = 0`, the allele unchanged). A stretching whose spelling
  plus `F` is still shorter than the read contributes nothing; with `F` the tract's stored flank (15
  bases in the chromosome 20 psps) and a read needing only 5 aligned flank bases to be complete
  (`MIN_ALIGNED_FLANK_BASES`, `src/locus_generation/ssr.rs:1371`), that does not arise from a read
  the generator classes as partial;
- **there is no factorised shortcut.** §5.2 takes the letters out of the sum on a pure tract because
  every stretching agrees on its first `ℓ` bases; with `F` appended they do not — a stretching by one
  unit puts a repeat base where `Δ = 0` puts a flank base — so the sum is always taken in full.

A slip is placed as §5.2 places it: on a pure tract at the end away from the held flank, on an
interrupted tract in each run with equal weight, a part-repeat change resized at that end. `F` is
appended after the stretched spelling. A read that held the **right** flank is the mirror image:
the candidate's flank deletion on the **left** (`L_del > 0`), `F` the reference left flank without
its last `L_del` bases, and the read, the spelling, the motif and `F` all written back to front, as
§5.2 already does for a right-held read. A whole-tract witness with no recorded side is scored from
both sides, each half under its own side's rule.

**This is what the STR-bundle path already does.** A bundle compares a partial read with the allele
followed by `flank_beyond_the_hull` — the reference flank on the side the read ran towards, less the
allele's own deletion (`src/calling/likelihood/ssr_emission.rs:2165`). A segment's partial reads
now see the same flank under the same definition.

**Unchanged:**

- complete reads (the paragraphs above);
- a candidate with no flank deletion on the side the read runs towards — the reference, every plain
  repeat allele, and a flank-deletion allele whose deletion is on the side the read held;
- STR bundles.

**Not adopted: the flank for every candidate.** Comparing every partial read with every allele
followed by its flank — the reference included, which a read that runs a base or two past the
reference tract also currently scores as slippage — was measured beside the rule above. It gains
more: on chromosome 20 at 30×, 71 indels with the wrong genotype against 83 for the rule above and
84 before; on HG002's tandem-repeat benchmark at 10×, 6,500 indel errors against 6,510 and 6,518.
But it adds 7 SNP errors on that benchmark at 30× (13,640 against 13,633), all on repeat ground, and
changes 35 STR genotypes on the tomato oracle's four accessions at about 3× with no truth to score
them. It changes how every partial read is scored at every tract, which is
[`read_likelihoods.md`](read_likelihoods.md) §5.2's rule rather than this design's, and needs its
own study (§9, deferred).

**One known cost.** HG003 chr17:28427006, truth `CA>C` (`0/1`) inside an `A` run, fitted at 30×: an
allele deleting the flank `C` before the run now beats the one-`A` contraction. Partial reads that
never reach the `C` used to count against the flank-deletion allele and now count for both, and the
complete reads do not separate the two. It is the only call the rule moved on the per-sample GIAB
panel (three samples at 5×, 10×, 30× and 300×).

*Measured on the build (2026-10-21, branch `str-partial-flank`).* That cost was the prototype's, and
came from its placing every slip at the tract's end. The tract there,
`AAAAAAAAAAAACAAAAAAACAAAAAAAAAAAAA`, has three runs of `A`; with each slip placed in each run, as
above, the one-`A` contraction wins again and the panel keeps every genotype of main. The cost the
build does carry is on HG002's tandem-repeat benchmark at 30×: at chr11:26874157 the truth is a
compound heterozygote, one copy deleting `AAAG` from a tract followed by `AAAG` × 7, the other
carrying `G>A` in the flank; the deletion is now called `1/1` and the generic `G>A` record, which it
covers, is dropped (§6). At chr11:44023063 the same mechanism removes a false SNP beside a true
homozygous flank deletion. Report: `../../reports/ng_generic_indel_error_2026-10-21.md` §5.2.

**At the ends of the range.** The rule is per read and estimates nothing, so it behaves the same at
one sample and at a thousand. It matters in proportion to how many partial reads a tract has, which
grows with depth: at 3× a tract has few, and on the tomato oracle (four accessions at about 3×) the
rule moved no call; at 60× a dozen partial reads were enough to flip a site.

**How we know it works.**

- Output is byte-identical to the build before this amendment wherever no candidate deletes flank
  bases: on the tomato oracle's slice, on HG002's 5 Mb set at 10×, 30× and 300×, and on the GIAB
  panel except chr17:28427006. *(Measured on the build: no genotype moves on any of the three,
  chr17:28427006 included; records whose candidates delete flank bases move in QUAL, allele
  frequency or GQ — 0 to 8 a call set there, 11 and 8 on the tomato oracle.)*
- chr20:34924682, 53879901 and 58447656 are written `1/1` at 60×.
- On HG002's tandem-repeat benchmark (`benchmarks/ssr_hg002/`, `--defaults`), indel errors (missed
  + false + wrong genotype) fall from 8,605, 6,518 and 5,003 to 8,597, 6,510 and 4,999 at 5×, 10× and
  30×, and no SNP count moves.
- On chromosome 20, indels with the wrong genotype fall from 37 to 34 at 60× and 84 to 83 at 30×
  (`--defaults`), from 35 to 32 and 67 to 65 fitted; no SNP count moves.
- Unit tests: a left-held partial read reading into the flank an allele leaves scores within a few
  per cent of the reference at a site where the two read the same (the 53879889 fixture), and well
  above the slippage-only score; the right-held mirror with a left flank deletion; a candidate
  without a flank deletion on the read's far side scores bit-identically to before; an interrupted
  tract places its slips per run before `F` is appended; a stretching whose spelling plus `F` is
  shorter than the read contributes zero.

### 5.7 The prior

**An allele with a flank deletion sits on the rung of its spelling's repeat count**, which is
where the spelling of the same reads sits on main today. Its stutter stratum and its place in the
contamination seed are the spelling's too. Nothing in the prior changes. Two alleles on one rung
each take the rung's full weight, as two spellings on one rung already do (`calling_priors.md`
§5.2, Q3).

### 5.8 Writing the STR record: widened over the flank

When the record is assembled, after the alleles no sample carries are dropped, each written allele
is given its **reach into each flank as the generic path would place it**:

- **on the right**, its flank deletion's right count — the smallest form already places a deletion
  as far left as it goes on that side;
- **on the left**, how many of the left flank's bases, counted back from the segment, the generic
  path would show as changed: the allele's sequence with its flank (the reference's left flank
  without the deleted bases, then the spelling) is compared with the reference from the flank's
  first base, and every flank base after the point where they first differ is reached. For a
  plain allele that is zero unless its deletion slides into the flank, as `AG` from `…CG|AGAG…`
  does (§4).

If any written allele reaches into a flank, the record's reference is widened by the largest reach
on each side, and every allele is respelled over the widened reference: the reference's flank
bases, less the allele's own flank deletion, around its spelling. The repeat count written per
sample (`REPCN`) is still that of the spelling, not of the widened sequence.

An allele that comes out empty is padded by [`vcf_output.md`](vcf_output.md) §5's anchor rule,
over the widened record: with the reference base **left** of the widened span. While the reach
leaves some of the left flank, that base is the flank's own; where it spends the whole flank, the
base lies just before the flank the locus holds and is read from the reference. Only where the
widened span starts at the contig's first base is it the base to the **right**. Taking the base to
the right whenever the flank was spent wrote a right-padded record in the middle of a contig,
which the writer refuses; that stopped the 2,174-sample tomato run on 2026-10-24.

## 6. Gathering: one record over the STR record and the generic records it overlaps

*(AMENDED 2026-10-29: pvc never writes a `*` allele. Until this amendment the pass rewrote each
overlapping generic record, writing `*` for every copy the STR allele deletes; that design, and why
it was replaced, is at the end of this section, §6.9. The user guide explains what `*` does to
population-genetics tools: `doc/user/overlapping_variants.md`.)*

**This runs on the finished records, in genome order, before the filters and the VCF writer** —
before the paralog filter's spill, which in psp mode already encodes each record — as a pass over
the emitted records, never a coupling between loci being called (`run_streaming.md` §4.3). It is
`src/run/edge_deletion_gathering.rs`, wrapped around both drivers, so the records are prepared for
the writer on the calling thread as the pass releases them.

**Records are held** while a later record could still be widened over them: anything ending more
than 256 bases before the newest record's start is released, in the writer's order (a generic
record before a tract where two start together, `vcf_output.md` §5). Generic records after a
widened STR record are merged as they arrive.

**Every generic record that overlaps a widened STR record is merged into it: one record over the
union of their spans, whose alleles are the haplotypes the samples carry over that span.** This is
the representation `vcf_output.md` §3 already uses for generic loci chained by overlap, and it
restores that section's rule that two records never claim the same base. The merged record is
still the STR record: it keeps `STR`, `RU`, `PERIOD` and `REPCN`, and the generic record is not
written.

The rest of this section says how each part of the merged record is built. Two facts about the
inputs shape all of it:

- **The STR record's alleles are whole events, phased by the reads that crossed the segment.** An
  allele that deletes flank bases already spells the flank over the bases it reaches.
- **Nothing links a generic call to an STR call read by read.** The STR path keeps each read's
  spelling of the segment and the generic path keeps counts per position; neither keeps which
  read it was, in the psp file or in memory (§3 rejected a combined locus for the same reason). So
  where a sample is heterozygous at both records, the run holds no read that says which STR allele
  carries which flank allele. §6.2 says what the merge does there.

How often the case arises, measured on main's build at 5ffb2444 with the benchmark suite's
outputs (`benchmark_results/pvc/5ffb2444/` in the benchmarks worktree):

| run | generic records written with `*` | samples needing a pairing the reads cannot give (§6.2) |
|---|---|---|
| stuttering benchmark, 1 sample, each of 8 runs (PCR and PCR-free, 5× to 60×), true stutter | 0 to 2 | 0 |
| stuttering benchmark, 50 samples, each of 8 runs, true stutter | 9 to 36 | 30 to 40 a run |
| tomato, 8 accessions (about 3×), `--defaults` | 6 | 0 |
| tomato, 62 accessions (about 3×), `--defaults` | 97 | 15, at 3 records |

The run does not print how many generic records the pass dropped (records whose genotypes named
only the reference and `*`); the VCF cannot show them. They are the cases §6.7 leaves unchanged.

### 6.1 The span and the alleles

**The merged record spans the union** of the widened STR record and every generic record that
overlaps it. Generic records are disjoint among themselves, so no further record can overlap the
union. Several generic records are merged one at a time, the record built so far playing the STR
record's part.

**Each allele is a haplotype: an STR allele, with a generic allele written over the flank bases
the generic record covers.** Which generic allele goes with which STR allele is decided sample by
sample (§6.2), and the record's alleles are the haplotypes the written genotypes name:

- an STR allele that **reaches** the generic record — its deletion takes, or its sequence changes,
  a base the generic record covers, as today's `reach_per_allele` says — is written as it is: its
  own sequence already says what the sample carries there;
- an STR allele that does not reach it is the reference over those bases, and is written with the
  generic allele it was paired with in their place. Paired with the generic reference, it is the
  STR allele unchanged.

**Order:** the reference; then the STR record's alleles, in their order, each unchanged; then each
STR allele paired with a generic alternative, ordered by STR allele and then by generic allele. An
allele no genotype names is dropped, its reads handled as §6.4 says, as assembly does.

**Bases every allele shares at the ends of the union are trimmed, as far back as the STR record's
own span.** A generic indel record carries a base before its event that both its alleles keep —
`AA`→`A` at the base before a homopolymer — so the union starts there, and every merged allele
then begins with that base. Trimmed, the record starts and ends where the STR record did and says
the same; a base trimmed from the front is the one an allele left empty is padded with. On the
stuttering benchmark (50 samples, PCR, 15×) the untrimmed merge moved 37 STR records one or two
bases left for nothing else.

**Not cut by the STR path's cap of 32 alleles.** The cap bounds the cost of calling — a locus of
`A` alleles is scored against `A(A+1)/2` genotypes per sample — and the merged record is not called
again. At most the record has its STR alleles times the generic record's alleles; measured, the
largest merged record would carry 12 alleles where its STR record carries 7 (50 samples, PCR, 5×),
and 9 where it carries 5 (tomato, 62 accessions).

### 6.2 Which STR allele each flank allele goes with

For each sample with both calls, its copies are paired in three steps:

1. **Each STR copy that reaches the generic record takes the generic copy whose allele its own
   sequence spells there** — the flank half of its deletion, or the substitution the mapper wrote
   instead (§4) — and where no copy's allele is spelled that way, a reference copy, then any. This
   is the rule the pass has used since 2026-10-11, and its reason stands: on one sample at 15× the
   mapper wrote a one-base flank deletion as a substitution for four of six reads. **Where the
   generic call did not agree, the STR record wins** (owner, 2026-10-10): the generic copy is
   overridden, and the run counts the copies taken in agreement and the copies overridden.
2. **The STR copies left pair with the generic copies left.** Where either side is homozygous, or
   one copy is left, there is one pairing and the reads were never needed.
3. **Where both are heterozygous, the pairing comes from the cohort.** Every sample whose pairing
   the first two steps settled is a count for each (STR allele, generic allele) pair it carries.
   With `h(s, g)` that count plus one half, the pairing `(s₁ with g₁, s₂ with g₂)` is chosen over
   the other when `h(s₁, g₁) · h(s₂, g₂)` is larger, and its probability, `p`, is that product
   over the sum of both. At one sample, or where no sample settled a pairing, both are one
   quarter and `p` is one half; the pairing is then the first by allele order, and the genotype's
   quality says it is a guess (§6.3).

**Measured against the simulator's haplotypes** (stuttering benchmark, 50 samples, all 8 runs,
true stutter): of the 238 samples where one of the two pairings matches the truth, the cohort rule
chose that one 238 times. In a further 31 neither pairing matches, because one of the two
genotypes is wrong; 23 of those are at 5×. **The simulator flatters the rule**: each flank SNP it
places rides on one STR allele only, so every carrier pairs it the same way. In a real population
an STR mutates often enough that the same SNP sits on several lengths, and `p` will be lower; it
is not measured on real data, which has no truth for it.

**A sample no-called at one of the two records:**

- **no-called at the generic record, called at the STR record**: it keeps its STR genotype, each
  STR allele written as it is — with the reference over the flank where it does not reach. That is
  no new claim: the widened STR record already spelled the reference there for that sample, and a
  consumer of the two records today reads the same. Its `GQ` is the STR record's. Measured: 5 and
  14 samples on the two 50-sample runs at 5×, none at 15× and above; 70 on tomato (62 accessions).
- **no-called at the STR record**: `./.` in the merged record, with its `DP` and `AD` written
  (§6.4). A haplotype needs the STR half, and no reference stands in for a repeat length nobody
  called. Where its generic call named an alternative, that call is lost: 2 samples on the
  stuttering benchmark (50 samples, PCR-free, 5×), none on the other runs, and 111 on tomato
  (62 accessions) — of 6,014 sample genotypes at those 97 records (97 × 62).

### 6.3 `GQ`

`GQ` is the probability that the written genotype is wrong. The merged genotype is wrong when the
STR genotype is, or the generic copies kept are, or the pairing is.

- **Where the generic genotype adds nothing** — every copy it keeps is the reference and step 1
  found it in agreement, or it is a no-call — `GQ` is the STR record's. That covers every sample of a generic record
  the replaced pass dropped.
- **Where it adds an alternative, or was overridden in step 1**, `GQ` is the lower of the two
  records' `GQ`s, as the replaced pass already wrote for an overridden copy. The lower and not the
  combined error of the two because the two calls come from largely the same reads, so their errors
  are not independent.
- **Where the pairing came from the cohort** (§6.2, step 3), with `e` the error of the lower of
  the two `GQ`s: `GQ = −10 log₁₀(1 − (1 − e) · p)`, capped at 99. At `p` one half that is at most
  3: a coin's chance of being wrong, written as one.

### 6.4 `AD`, `DP` and the mapping qualities

**From the STR record; the generic record's reads are not written.** They are largely the same
reads, and the generic record counts them per base, not per haplotype, so adding them would count
a read twice and split none.

- **`AD`**: an STR allele's reads go to the one allele built on it, where the record has only one;
  otherwise to the one allele built on it that the sample's genotype names; otherwise — none named,
  or two, as at a sample homozygous at the STR and heterozygous in the flank — to the reads no
  written allele explains (`DP − ΣAD`). The reads cannot be split there: no read was compared with
  a haplotype.
- **`DP`**: the STR record's, unchanged.
- **The mapping-quality pools** (`MQREF`, `MQALT`) follow `AD`: an allele's pool is the reads `AD`
  gives it, at its STR allele's mean mapping quality.

### 6.5 `REPCN`

**Each allele's `REPCN` is its STR allele's**, counted from the allele less its flank bases, as
since §5.8. Two alleles with the same `REPCN` that differ in the flank are two sequences on one
rung, as two spellings of one length already are (`vcf_output.md` §3 keeps them for that reason). A
consumer that reads lengths sees the same `REPCN` as before the merge, whichever pairing §6.2 chose.

### 6.6 `QUAL`, `FILTER`, `AF` and `AC`

- **`QUAL`** is the larger of the two records', with that record's `ABPEN` and `SPPEN`. The merged
  record holds a variant when either did, so the chance that it holds none is at most the smaller
  of the two; the larger `QUAL` says no more than that.
- **`FILTER`** is the STR record's verdict, unless it is `PASS` and the generic record's is not.
- **Both count a generic record only where it adds an allele of its own** — where some written
  haplotype carries one of its alternatives that the STR allele does not spell. A generic record
  whose every called alternative is the flank half of an STR allele measured the same event as
  the STR record and adds nothing to the site.
- **`AF`**: an STR allele's fitted copies are shared among the alleles built on it in proportion to
  the copies the written genotypes give each. **`AC` and `AN`** are counted from the written
  genotypes, as on every record.
- **The filters after the gathering judge the merged record as one record** (owner, 2026-10-29).
  The hidden-duplication filter scores it like any multi-allelic record, so its verdict can differ
  from the one either record would have had alone. On tomato (62 accessions, about 3×, fitted
  parameters) that happens at 5 sites in 169,115 records: two STR records dropped alone are
  written merged (one scored 0.9998 alone, 0 merged), and three generic alleles that were dropped
  alone are written inside a passing STR record. *Rejected:* dropping a generic record before the merge where the filter
  would drop it alone — it would make the gathering depend on a filter that runs after it.

### 6.7 The records the replaced pass dropped

**They are merged like any other, and the result is the STR record unchanged.** Every sample's
generic genotype named only the reference and the flank half, so §6.2 pairs each STR allele with
the generic reference or takes the flank half into the STR allele that spells it, and §6.3 keeps
the STR record's `GQ`. The one difference from the replaced pass is a sample whose generic copy
the STR record overrode: the lower `GQ` the replaced pass wrote went with the dropped record, and
now it is the merged record's (§6.3).

### 6.8 At the ends of the range

- **One sample.** The cohort has no pairing to offer, so a sample heterozygous at both records
  gets `p` one half and `GQ` at most 3, with `REPCN` intact. On the stuttering benchmark's one-sample
  runs no sample reaches this case.
- **Thousands of samples.** The pairing counts are one pass over a record's samples, and the
  merged record's alleles are bounded by the two records' alleles multiplied, not by the samples.
  The held records are a flank's worth, as before.
- **Low depth.** No-calls are commoner at 3× than at 15×, and §6.2 turns a no-call at the STR
  record into a no-call of the merged record: on tomato at about 3× that loses the 111 flank
  genotypes above.

### 6.9 The pass this replaced (2026-10-10 to 2026-10-29)

For every generic record overlapping a widened STR record, sample by sample, each copy the STR
allele reaches became `*` — "deleted by an overlapping record", as GATK writes it — on the generic
copy chosen as §6.2's step 1 chooses; the STR record won where the generic call disagreed, with the
lower `GQ`; and a generic record left naming only the reference and `*` was dropped. A generic
record that kept another sample's SNP was written beside the STR record, overlapping it. Replaced
because population-genetics tools mostly do not understand `*`: read as missing, it makes
missingness follow the deletion's carriers; recoded as the reference, it invents reference copies
nobody carries; counted as an allele, it counts one deletion once per site it covers
(`doc/user/overlapping_variants.md`).

**Also rejected, 2026-10-29:**

- **Dropping the generic record's alternatives where they meet a widened STR record.** No `*`,
  and simpler, but it deletes real flank SNPs: on the stuttering benchmark the 30 to 40 samples a
  50-sample run that §6.2 pairs from the cohort are all carriers of one.
- **Not widening an STR record where a flank variant would fall inside the widening.** That
  returns those segments to the split deletion of §1, the problem this design exists to fix.

## 7. Limits known in advance

- **A deletion that eats the whole flank anchor stays invisible**: the read rule needs five aligned
  bases beyond the repeat on each side (`MIN_ALIGNED_FLANK_BASES`).
- **Insertions across the edge** are not covered by the c0b55125 rule and not by this design;
  whether they occur in the benchmark, and how often, is unmeasured.
- **A generic allele that carries the flank half and another variant** — a SNP chained into the
  same generic locus — is not the flank half alone, so the STR copy that reaches it does not
  recognise it (§6.2, step 1) and overrides a copy: the STR record wins, and the other variant is
  lost from that sample's haplotype. The run counts the copies overridden. *(Until 2026-10-29 the
  sample could instead end with the deletion spelled twice.)*
- **A widened record covers flank bases where other samples may carry a SNP.** Since 2026-10-29
  that SNP is an allele of the merged record (§6), on the STR allele the sample's pairing gives it.
  *(Until then a sample whose STR allele has no flank deletion was spelled with the reference over
  those bases by the STR record while a generic record gave it the SNP, and a comparison applying
  whole alleles saw two edits on one base. 17 of the 97 benchmark segments with an edge deletion
  also have a SNP beside the segment; the measurement reports them separately, §8.)*
- **Cohort size.** At one sample the evidence for a flank deletion is that sample's two or more
  crossing reads — the same bar as any allele, so a heterozygous edge deletion at 5× will often go
  uncalled, as other alleles do there. At a thousand samples the reconciliation touches only
  segments with a written flank-deletion allele, and costs one pass over the samples of each
  generic record within a flank's reach; the held records are a flank's worth (at most the locus's
  stored flank, tens of bases) of generic records, so memory does not grow with the genome.
- **`--explain-loci` rows** describe the records before the merge.
- **Two widened STR records over one stretch of flank are not merged with each other.** Where the
  right widening of one segment and the left widening of the next would share bases, a generic
  record overlapping both merges into the first only, and the two STR records are written
  overlapping. Not seen so far: no two records share a base in any of the 16 stuttering-benchmark
  VCFs or tomato's (62 accessions), against 97 pairs on main's tomato VCF, all of them `*`
  records (2026-10-29). Segments that close are usually one STR bundle.

## 8. Building and measuring it

**Built as the only behaviour, with no switch** (owner, 2026-10-10): the branch bets on this design.
If it wins on both benchmarks it is merged and is the only way forward; if it loses, it is
abandoned. Each step is measured against main's binary at 95e9b916 and against the step before it.

**Decisions of 2026-10-10 on the questions this section's predecessors left open:** where the two
paths disagree, the STR record wins (§6, item 3); and the stuttering benchmark's scorer trims from
each allele the bases it shares with the reference before checking whether two records' alleles
touch one base, so that a widened record over a flank SNP is not scored as a collision (§7). The
scorer change is shown on its own — GATK and main rescored — before any number of this work uses
it.

Three steps, each measured against the one before:

1. **Measurement and storage** (§4, §5.1–§5.3): the reads record their flank deletion; nothing is
   called differently yet. Checked by counting, on the benchmark psps, the observations with a
   flank deletion at the segments whose truth has one, and with the VCF unchanged.
2. **Calling and the widened record** (§5.4–§5.8). Measured as below; the generic halves are still
   written, so carriers can show the deletion twice — this step's number is the STR record's own.
3. **Reconciliation** (§6). Measured as below.

The measurement, as for the 2026-10-10 changes:

- **The stuttering benchmark, all 16 runs**, true stutter, with `scripts/ng_str_accuracy/`
  (`copy_run.sh`, `call_and_score.sh`, `compare.py`): right genotypes and false calls overall, and
  separately the pairs carrying a deletion across the edge (`edge_deletions.py`) and the pairs
  with a SNP beside the segment (truth change `snp_beside`).
- **HG002's tandem-repeat benchmark** at 30× and 50× (`hg002_two_arms.sh`): right genotypes by
  period class, the per-tract verdict flips, and the allele-exact false positives, which should
  fall where the two-record spelling becomes one record.
- The edge-deletion pairs at one sample should approach GATK's 26–27 of 27 at 15× and above.

## 9. Deferred, with a recommended home

- **Comparing every partial read with every allele followed by its flank** — the reference and every
  plain allele, not only alleles that delete flank bases (§5.6, "Not adopted"). Measured: more
  indels right (chromosome 20 at 30×, 71 wrong genotypes against 83), 7 more SNP errors on HG002's
  tandem-repeat benchmark at 30×, 35 tomato genotypes moved with no truth. **Home:**
  [`read_likelihoods.md`](read_likelihoods.md) §5.2, whose rule it changes for every tract; it wants
  the stuttering benchmark's 16 runs and the per-period breakdown of the SNP cost before a decision.
  The prototype is commit 707d725f on branch `generic-indel-error` (`PVC_EXP_PARTIAL_FLANK=all`).
