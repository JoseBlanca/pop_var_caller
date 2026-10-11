# A deletion across an STR segment's edge, called as one allele of one record

*Implementation report, 2026-10-11, branch `str-edge-deletion` (from main at 95e9b916). Design:
[`str_edge_deletion_alleles.md`](../../ng/spec/str_edge_deletion_alleles.md). Every ng run here
used the stuttering benchmark's true stutter, one group per library; no fit was run (the fit pools
every library into one group, a defect being fixed elsewhere). Scratch, not committed: `tmp/` in the
branch's worktree.*

---

## 1. What changed

Before this work, a deletion running from an STR segment into a flank was written as two records:
the STR record carried the repeat's half and the generic path's record the flank's half. They lost
phase, and an allele-exact comparison counted them as a miss and a false positive.

Now:

- **A read records the flank bases its deletion took**, left and right, beside its spelling of the
  repeat, in the smallest form the sequence allows (a deletion that fits inside the repeat records
  none). The psp format is 1.2: one field at the end of every record. 1.1 files read as having none,
  and appending to one is refused.
- **The merge, candidate selection and the read model key on spelling plus that extent.** A read and
  an allele over different extents are compared over the wider one: differing by substitutions,
  they are priced at the substitution rate; differing by flank bases deleted, at the read model's
  price for an indel that is not whole repeat units.
- **An allele one misread base from another is not a candidate** when, in every sample that shows it,
  it has under a third of the two alleles' reads. That is the shape of a systematic read error beside
  a run (on HG002, 4 reads against 13 misread the `G` after 25 `T`s, the same way).
- **The STR record is written widened** over the flank bases its alleles reach as the generic path,
  which left-aligns every read, places them. `REPCN` still counts the repeat.
- **Where the calls are gathered**, before the filters and the VCF, each sample's copy whose STR
  allele deletes bases a generic record covers is written `*` there — "deleted by an overlapping
  record" — on the generic copy the STR allele's own sequence spells. Where the two paths disagree,
  the STR record wins (owner, 2026-10-10) and the sample's `GQ` there becomes the lower of the two. A
  generic record left with nothing of its own is dropped. Records near a segment are held and
  released in the writer's order.

Tried and not adopted, each measured: sliding every deletion left at measurement (it turned every
stutter contraction at a segment like `…CG|AGAG…` into a flank-deletion allele); recording
substituted flank bases beside the repeat (−2 to +3 edge-deletion pairs a run, up to 7 more false
calls a run at 50 samples); folding every allele one base from another (it removed the commonest
real edge deletion and lost 36 to 77 carrier pairs a run at 50 samples); leaving homopolymers out
(within 2 tracts on HG002, up to 16 edge-deletion pairs a run fewer on the benchmark than the
version adopted).

## 2. Two scorer fixes

- **The stuttering benchmark's scorer** trims from each allele the bases it shares with its
  record's reference before checking whether two edits touch one base (owner, 2026-10-10). Before
  it, main's correct two-record answers at a segment's right edge were scored wrong, because the
  generic record's padding base is the segment's last base: under the fix main gains 2 to 8 pairs a
  run at one sample and 132 to 296 at 50 samples, GATK at most 27. So main's start, at one sample
  (PCR, 15×), was 25 of 27 edge-deletion pairs right, not the 16–17 quoted before this work; the gap
  was real at PCR-free (16).
- **The HG002 scorer** reads `*` as no change on that copy, as VCF defines it; it read it as bases.
  Main writes no `*`, so main's numbers do not move.

## 3. The stuttering benchmark: four callers

Pairs called right (one sample at one STR segment) / reference pairs given a false call. Each run
has 1,125 variant and 1,613 reference pairs at one sample, and 56,621 and 80,279 at 50 samples; the
262 hexanucleotide segments of motif AAAAAC/AAAAAG/AAAAAT are left out, as in the 2026-10-10 report.
freebayes ran with defaults (reference, the run's targets, `--ploidy 2`, every sample's reads
together) and its records were kept at QUAL 20 or more (owner, 2026-10-11); GATK's calls are the
first pass's.

| run | main | this branch | freebayes, QUAL ≥ 20 | GATK |
|---|---|---|---|---|
| PCR, 1 sample, 5× | 508 / 14 | 508 / 14 | 321 / 1 | 529 / 9 |
| PCR, 1 sample, 15× | 1,030 / 17 | 1,031 / 17 | 801 / 6 | 1,025 / 18 |
| PCR, 1 sample, 30× | 1,113 / 3 | 1,116 / 3 | 1,018 / 4 | 1,095 / 21 |
| PCR, 1 sample, 60× | 1,121 / 1 | 1,124 / 1 | 1,100 / 5 | 1,107 / 23 |
| PCR-free, 1 sample, 5× | 548 / 2 | 548 / 2 | 349 / 1 | 541 / 2 |
| PCR-free, 1 sample, 15× | 1,034 / 1 | 1,041 / 1 | 848 / 1 | 1,044 / 0 |
| PCR-free, 1 sample, 30× | 1,112 / 0 | 1,119 / 0 | 1,049 / 0 | 1,112 / 1 |
| PCR-free, 1 sample, 60× | 1,119 / 0 | 1,123 / 0 | 1,108 / 0 | 1,120 / 0 |
| PCR, 50 samples, 5× | 40,041 / 1,449 | 40,074 / 1,441 | 31,893 / 6,844 | 25,558 / 2,004 |
| PCR, 50 samples, 15× | 53,476 / 338 | 53,500 / 336 | 49,948 / 10,649 | 49,890 / 7,309 |
| PCR, 50 samples, 30× | 55,961 / 76 | 55,978 / 76 | 53,456 / 9,864 | 53,290 / 9,445 |
| PCR, 50 samples, 60× | 56,486 / 14 | 56,496 / 14 | 54,575 / 9,585 | 53,952 / 9,961 |
| PCR-free, 50 samples, 5× | 41,764 / 853 | 41,790 / 851 | 34,159 / 964 | 27,532 / 172 |
| PCR-free, 50 samples, 15× | 54,592 / 81 | 54,612 / 81 | 52,459 / 1,260 | 52,627 / 546 |
| PCR-free, 50 samples, 30× | 56,366 / 11 | 56,381 / 11 | 55,445 / 639 | 55,809 / 404 |
| PCR-free, 50 samples, 60× | 56,573 / 0 | 56,581 / 0 | 56,266 / 338 | 56,307 / 229 |

**The pairs carrying a deletion across the segment's edge** (27 a run at one sample, 1,085 at 50):

| run | main | this branch | freebayes, QUAL ≥ 20 | GATK |
|---|---|---|---|---|
| PCR, 1 sample, 15× / 30× / 60× | 25 / 23 / 23 | 26 / 26 / 26 | 22 / 27 / 27 | 26 / 26 / 27 |
| PCR-free, 1 sample, 15× / 30× / 60× | 16 / 19 / 22 | 23 / 26 / 26 | 21 / 24 / 26 | 26 / 26 / 27 |
| PCR, 50 samples, 15× / 30× / 60× | 998 / 1,024 / 1,028 | 1,019 / 1,040 / 1,038 | 1,026 / 1,067 / 1,080 | 1,015 / 1,062 / 1,069 |
| PCR-free, 50 samples, 15× / 30× / 60× | 996 / 1,029 / 1,042 | 1,014 / 1,041 / 1,050 | 991 / 1,060 / 1,070 | 1,042 / 1,081 / 1,085 |

The pairs with a SNP beside the segment, which the gathering of the calls touches, are never fewer
than main's in any run (3,556 a run at 50 samples; at most 3 more right).

What the grid says:

- **Against main**, this branch is never worse: +0 to +7 right pairs a run at one sample, +8 to +33
  at 50 samples, false calls equal or up to 8 fewer.
- **At 50 samples ng is far ahead of both haplotype callers** on everything but edge deletions: with
  PCR at 15× and above, freebayes and GATK each make 7,000 to 10,600 false calls a run against
  ng's 14 to 336, and at 5× ng gets 7,600 to 14,500 more pairs right a run.
- **At one sample** ng and GATK are within 21 pairs a run; freebayes at QUAL 20 makes the fewest
  false calls (0 to 6) and gets 187 to 230 fewer pairs right at 5× and 15×.
- **Edge deletions at 50 samples, 30× and above, are still where freebayes and GATK lead**: by 19 to
  42 pairs a run of 1,085. Both score each read against whole haplotypes that span flank and repeat;
  ng's remaining losses are concentrated at deletions taking one or two flank bases.

## 4. HG002's tandem-repeat benchmark, 30× and 50×

Main against this branch, one binary each, scored by `benchmarks/lib/tract_qual_experiment.py` with
`*` read as no change.

| | main, 30× | branch, 30× | main, 50× | branch, 50× |
|---|---|---|---|---|
| right genotypes, homopolymers (3,995 tracts) | 3,503 | 3,503 | 3,595 | 3,595 |
| right genotypes, period 2 and above (2,980 tracts) | 2,557 | 2,566 | 2,603 | 2,613 |
| allele-exact false positives, homopolymers | 594 | 452 | 587 | 441 |
| allele-exact false positives, period 2+ | 574 | 495 | 584 | 509 |

At homopolymers the unchanged count is tracts trading places: 8 gained and 8 lost at 30×, 6 and 6
at 50×. Allele-exact true positives rise by 56 to 81 per depth and class.

## 5. What is left

- **Edge deletions at 50 samples, 30× and above**, 19 to 42 pairs a run behind freebayes and GATK.
- **The share in the misread-allele rule (a third)** is a judgement, not a fitted number.
- **HG002 homopolymers** trade tracts both ways; the 8 lost at 30× were not inspected one by one.
- **Insertions across a segment's edge** are not covered (spec §7).
