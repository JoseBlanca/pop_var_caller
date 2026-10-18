# A deletion that touches a repeat tract's edge counts as a flank deletion

**Date:** 2026-10-18. **Spec:** `doc/devel/ng/spec/str_edge_deletion_alleles.md` §5.1, amended.
**Code:** `deletion_takes_flank_bases_at_a_tract_edge` and `flank_bases_deleted` in
`src/locus_generation/ssr.rs`.

## The defect

**Some reads of a deletion allele at a homopolymer were counted as the reference length.** Where a
run of one base is followed by a different base (`…CTTTTTTTTTG…`), a read that misreads the
neighbouring base as the run's base shows a run one base too long. minimap2 writes such a read as
the run at full length and the neighbouring base deleted: a one-base deletion that **ends at the
tract's first base, or starts just after its last, without crossing it**.

The STR path recorded a flank deletion only for a deletion that crossed the tract's edge. A touching
one was ignored, and the read was spelled over the tract's positions alone, at the reference length.
At a tract homozygous for a one-base deletion, those reads became evidence for a reference allele
the sample does not carry.

**Measured on HG002 at 30x** (the GIAB tandem-repeat benchmark, minimap2 alignments, per-read dump
of 2026-10-18): at homopolymer tracts homozygous for a deletion, 926 reads were spelled at the
reference length, and 565 of them carried a one-base deletion touching the tract's edge.

## The change

A deletion touching the edge now counts as one crossing it: its flank bases are recorded beside the
spelling, and the read is complete by the same rule (`MIN_ALIGNED_FLANK_BASES` aligned bases beyond
the tract on each side). Calling then compares the read with each allele over the wider of the two
flank extents (§5.6): against the deletion allele it differs by the misread base alone and is
priced at the substitution rate; against the reference it differs in length.

It is not the variant §4 tried and did not adopt on 2026-10-11, which counted a flank base the read
shows *substituted* as part of the extent. That turned every flank sequencing error into an extra
observation. A substitution error does not produce a deletion in the alignment, so this change does
not.

## What it does

### HG002, tandem-repeat benchmark, one sample, `--defaults`

`main` (outlier weight 0.07) against this branch; correct tract genotypes over the tracts both
built:

| | 5x | 10x | 30x | 50x |
|---|---|---|---|---|
| `main` | 2,857 | 4,686 | 7,176 | 8,904 |
| this branch | 2,873 | 4,709 | 7,225 | 8,953 |
| change | +16 | +23 | +49 | +49 |
| tracts gained / lost | 27 / 11 | 46 / 23 | 75 / 26 | 78 / 29 |

By truth genotype at 50x: homozygous non-reference +46, heterozygous +19, homozygous reference
−16. **Some of the reference losses are the tract scorer's, not the caller's.** The scorer counts a
truth record only if it starts inside the tract; at tracts such as chr5:7,613,161 the branch now
spells 19 reads with a one-base deletion beside the tract and writes it in the widened record, which
a truth record starting just outside the tract would match but the scorer does not see.

Haplotype level, RTG vcfeval 3.13 against the benchmark's truth in its regions (`None` row):

| | 5x | 10x | 30x | 50x |
|---|---|---|---|---|
| true positives, truth side | +15 | −12 | +30 | +32 |
| false positives | −3 | +7 | −25 | −29 |
| F-measure, `main` → branch | 0.4469 → 0.4474 | 0.5856 → 0.5852 | 0.6629 → 0.6638 | 0.6743 → 0.6754 |

### Tomato oracle (four tomato2 samples, macOS host)

`main` at 5c741b11 reproduced the recorded baseline first; every line that moves is this change's.

- **The designed swap, 14 times on the defaults route.** A one-base generic deletion beside a tract
  stops being its own record and is written inside a widened STR record — `AT → A` at
  SL4.0ch01:11,805,609, 1/1 in all four samples, becomes `TAAAAAAAA → AAAAAAAA` at 11,805,610, 1/1
  in all four. Where a tract's alleles re-split, 3 sample genotypes change (heterozygous to
  homozygous). 9 STR records move in QUAL; non-STR records move only in the duplication filter's
  numbers (at most 0.000045 in posterior). Records: 7,509 → 7,508.
- **The fitted stutter loses an inflation.** Reads with a flank deletion are left out of the stutter
  census (§5.3), so touching reads at deletion alleles no longer count as slips: homopolymer slip
  rates at 16–21 repeats fall 15–24% (0.0690 → 0.0522 at 18 repeats); the median slip row moves
  0.5%. The fitted route then swaps records the same way (15 only on the branch, 18 only on `main`),
  3 sample genotypes change and 259 records move in QUAL.
- **The window lines move with the records**: the table has one line per written record per sample,
  so 12 new records add 12 lines a sample, and 7 records near tracts have their depth counted over
  slightly different positions.

### Tests

`cargo test`: 5,244 passed. The one failure,
`a_frequency_that_is_not_a_probability_does_not_reach_the_caller`, fails identically on `main`
(a debug assertion in the regularised incomplete Beta function, from the parameter-fit performance
merge). The cross-platform checksums do not move: their fixture has no deletion beside a tract.

## Not measured

**The stuttering benchmark**, which §8 of the spec asks for beside HG002, could not be run: its
16 simulated runs lived in a worktree that has since been deleted, and regenerating them means
rerunning the simulation and its nine callers.
