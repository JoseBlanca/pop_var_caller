# ng — a deletion across an STR segment's edge, called as one allele of one record

*Design draft, 2026-10-10. The direction is the owner's decision of 2026-10-10; the details in §4
and §6 are for the implementing session to settle and bring back. Reads on:
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
overlap it are reconciled with it (§5).**

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

## 4. What changes, path by path — to be settled by the implementation

- **Read measurement** (`locus_generation/ssr.rs`, `classify_read`). The edge-deletion rule of
  c0b55125 already identifies the reads and spells their repeat part from their own alignment. It
  must also record the flank bases removed on each side. **Open:** how the observation carries it —
  an extended spelling over the deleted flank positions, or a separate field — and what that does
  to the psp format and its version.
- **Candidates** (`allele_candidates/ssr.rs`). Two observations with the same repeat spelling and
  different flank deletions are different alleles. Selection, the support bar and the
  edge-deletion exemption of the periodicity check (§7 of `candidate_alleles_ssr.md`) must key on
  the pair. **Open:** whether a flank deletion seen only by a sample's few reads can earn a
  candidate under the same bar.
- **Read likelihood** (`calling/likelihood/ssr*.rs`). A read's emission against an allele whose
  flank deletion differs from the read's own is not stutter: the stutter term covers whole-unit
  changes inside the repeat. **Open:** what a mismatch in flank deletion costs — an indel error
  rate, or zero, and whether partial reads (which cannot see the far flank) are scored as today.
- **The prior** keys alleles by repeat count today (`calling_priors.md` §5). **Open:** where an
  allele with a flank deletion sits on that ladder.
- **Gathering the calls** (§5).

## 5. Gathering: widening the record and reconciling the generic calls

Where an STR record has an allele with a flank deletion, the record's reference is widened to cover
the deleted flank bases, and every allele of the record is respelled over the widened reference.
The generic records overlapping the widened span are then reconciled:

- **The generic allele that is the flank half of the same deletion is removed**, for every sample,
  since the STR record now carries the whole event.
- **Other generic records there stay** — another sample's SNP in the flank, say — with the samples
  carrying the widened deletion given the spanning-deletion allele `*` at those positions, which is
  how VCF writes overlapping variants and how GATK writes them. Removing these records outright
  would delete other samples' calls.
- **Where the two paths disagree** — the STR record says a sample carries a flank deletion and the
  generic record over that flank says it does not — **open**: which wins, and whether the sample's
  genotype there is written missing.

This runs where the calls are gathered into the call stream, before the filters (including the
hidden-duplication filter, which scores STR records from their read split) and the VCF writer.

## 6. Limits known in advance

- **A deletion that eats the whole flank anchor stays invisible**: the read rule needs five aligned
  bases beyond the repeat on each side (`MIN_ALIGNED_FLANK_BASES`).
- **Insertions across the edge** are not covered by the c0b55125 rule and not by this design;
  whether they occur in the benchmark, and how often, is unmeasured.
- **Cohort size.** At one sample the read evidence for a flank deletion is a few reads; at a
  thousand samples the gathered record's `*` handling is per sample over every overlapping generic
  record. Both ends need their cost and behaviour stated by the implementation.

## 7. How we know it works

The same protocol as the 2026-10-10 changes, each run against the build before it:

- **The stuttering benchmark, all 16 runs**, true stutter, with `scripts/ng_str_accuracy/`
  (`copy_run.sh`, `call_and_score.sh`, `compare.py`): right genotypes and false calls overall, and
  separately the pairs carrying a deletion across the edge (`edge_deletions.py`) and the pairs with
  a SNP beside the segment (truth change `snp_beside`), which the reconciliation of §5 touches.
- **HG002's tandem-repeat benchmark** at 30× and 50× (`hg002_two_arms.sh`): genotype accuracy by
  period class, the per-tract verdict flips, and the allele-exact false positives, which should
  fall where the two-record spelling becomes one record.
- The edge-deletion pairs at one sample should approach GATK's 26–27 of 27 at 15× and above.
