# ng — a deletion across an STR segment's edge, called as one allele of one record

*Design, 2026-10-10. The direction is the owner's decision of 2026-10-10; §4–§8 are the implementing
session's answers to the questions the draft left open, awaiting the owner's confirmation. Reads on:
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
overlap it are reconciled with it (§6).**

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
  from the spelling (`candidate_alleles_ssr.md` §7).

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

**Reads that end inside the segment are scored as today.** They see one flank only; the delimiter
anchored them there, so they record no extent and are compared with the allele's spelling as now.
Scoring them against the allele's extent on the side they anchored is a possible refinement, left
out so that its effect is not mixed into this one.

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

## 6. Gathering: reconciling the generic records

**This runs on the finished records, in genome order, before the filters and the VCF writer** —
before the paralog filter's spill, which in psp mode already encodes each record — as a pass over
the emitted records, never a coupling between loci being called (`run_streaming.md` §4.3). It is
`src/run/edge_deletion_gathering.rs`, wrapped around both drivers, so the records are prepared for
the writer on the calling thread as the pass releases them.

**Records are held** while a later record could still be widened over them: anything ending more
than 256 bases before the newest record's start is released, in the writer's order (a generic
record before a tract where two start together, `vcf_output.md` §5). Generic records after a
widened STR record are reconciled as they arrive.

For every generic record overlapping a widened STR record, sample by sample:

1. **Each copy of the sample whose STR allele reaches a base of the generic record — changes it
   where the generic path looks — becomes `*` there**, which is how VCF writes "deleted by an
   overlapping record" and how GATK writes it. The generic copy it takes is the one whose allele
   the STR allele's own sequence spells over that record: the flank half of the deletion, or the
   substitution the mapper wrote instead (§4), whichever form the generic path called. Where no
   copy's allele is spelled that way, a reference copy, then any.
2. **Where the generic call did not agree, the STR record wins** (owner, 2026-10-10): the copy is
   rewritten anyway, and the sample's `GQ` there becomes the lower of the two records' `GQ`s. The
   run counts the copies taken in agreement and the copies rewritten.
3. **A generic record whose genotypes then name nothing but the reference and `*` is dropped.**
   Otherwise it stays — another sample's SNP in the flank, say — and an allele no genotype names
   any more is dropped, its reads moved to the reads no written allele explains, as assembly does.
4. A sample with no STR call, or no generic call, is left as the generic record has it.

*Changed during implementation (2026-10-11).* The first build recognised the flank half only as a
deletion and turned that allele into `*` for every sample. On one sample at 15× the mapper had
written a one-base flank deletion as a substitution for four of six reads; the generic path called
the substitution, the rule gave the `*` to the sample's reference copy, and the substitution stayed
on the wrong chromosome.

## 7. Limits known in advance

- **A deletion that eats the whole flank anchor stays invisible**: the read rule needs five aligned
  bases beyond the repeat on each side (`MIN_ALIGNED_FLANK_BASES`).
- **Insertions across the edge** are not covered by the c0b55125 rule and not by this design;
  whether they occur in the benchmark, and how often, is unmeasured.
- **A generic allele that carries the flank half and another variant** — a SNP chained into the
  same generic locus — is not the flank half alone and is not turned into `*`. That sample can end
  with the deletion spelled twice; the run counts these.
- **A widened record covers flank bases where other samples may carry a SNP.** A sample whose STR
  allele has no flank deletion is spelled with the reference over those bases by the STR record
  while a generic record gives it the SNP. The VCF reading is that the STR allele makes no change
  there, but a comparison that applies whole alleles sees two edits on one base. 17 of the 97
  benchmark segments with an edge deletion also have a SNP beside the segment; the measurement
  reports them separately (§8).
- **Cohort size.** At one sample the evidence for a flank deletion is that sample's two or more
  crossing reads — the same bar as any allele, so a heterozygous edge deletion at 5× will often go
  uncalled, as other alleles do there. At a thousand samples the reconciliation touches only
  segments with a written flank-deletion allele, and costs one pass over the samples of each
  generic record within a flank's reach; the held records are a flank's worth (at most the locus's
  stored flank, tens of bases) of generic records, so memory does not grow with the genome.
- **`--explain-loci` rows** describe the record before reconciliation.

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
