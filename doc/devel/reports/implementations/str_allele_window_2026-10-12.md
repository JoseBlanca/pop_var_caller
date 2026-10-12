# The allele window at ±9: what was built and what it measured

**2026-10-12, branch `stutter-multiplier`, commits `38b1fe56` and `ef4b0546`.** Spec:
[`str_allele_window.md`](../../ng/spec/str_allele_window.md). Plan:
[`str_allele_window.md`](../../ng/impl_plan/str_allele_window.md).

## What changed

- **The repeat-tract fit lets an allele sit up to 9 repeats from the reference length**, one past
  the ±8 the census records. It was 6. The window is written as the recorded range plus one and
  checked when the crate compiles, so the two cannot drift apart again. Both fits use it: the usual
  fit, which climbs to the best-scoring numbers, and the counting fit.
- **The two outermost classes mean "this far or further"**, and calling now shares an edge class's
  weight equally among the candidate lengths at or past it. Before, a candidate past the window took
  a floor of 1 in 10¹².
- **The counting fit reads which way a slip went, and how far, from the slip model rather than from
  the read's bucket.** An allele past the recorded range shares its end bucket with its own shorter
  slips and its one-step longer slip, and measuring from the bucket lost them. The code review
  found this (`reports/reviews/str_allele_window_review_2026-10-12.md`).

## A drawn class: the defect, and the window recovering it

A class drawn with a fifth of its chromosomes 7 or 8 repeats short, slip rate 0.08, 20 samples at
10 reads each, reads recorded to ±8.

| fit | window | slip rate | far shares |
|---|---|---|---|
| usual (500 tracts) | ±6 | **0.28** | 0.18 piled in the −6 class |
| usual (500 tracts) | ±9 | **0.081** | −9: 0.000, −8: 0.193, −7: 0.045 |
| counting (2,000 tracts) | ±9 | 0.079 | −9 and −8 together 0.147, against 0.15 drawn |

The counting fit splits the −8 alleles between −8 and −9 (0.03 of the 0.15 sits at −9 on this draw).
That is expected: the two differ only in reads that slipped inwards. The tests check the two
classes together.

## tomato2: 59 accessions, 32 Mb

**Slip rates, counting fit, ±6 against ±9** (`tmp/window/tomato2_counted_span9_fixed.tsv`, against
`tmp/route2/tomato2_counted_span6.tsv`):

- **Short classes barely move.** Homopolymers of 8–16 repeats and dinucleotides of 6–15 repeats
  fall by 3.3% or less. Trinucleotides of 6–11 repeats fall by up to 7%, and of 12 and 13 repeats
  (10 and 9 tracts) by 14% and 24%.
- **Long classes fall.** Homopolymers of 17–21 repeats fall by 5–25%. Dinucleotides of 16–20
  repeats fall by 18–31%. Dinucleotides of 21, 23, 24, 25 and 27 repeats fall by 2.3 to 5.5 times;
  there, 0.23–0.42 of each class's alleles sit 7 or more repeats from the reference, and the ±6
  window had forced their reads to count as stutter. Each of those classes has 11–21 tracts. The
  22-repeat class, with 0.03 of its alleles that far, falls 9%.
- **The counting-fit fix moves the slip rates by at most 0.013.** It moves the shorter share in five
  dinucleotide classes of 20 repeats or more, each with 300 reads or fewer. One of them (23 repeats,
  124 reads) went from 0 to 1, which says the class cannot place a direction at all.

**The usual fit puts mass on the far classes that the counting fit does not.** On the 13-repeat
homopolymer class it puts 0.033 of its alleles 7 or more repeats away, where the counting fit puts
0.001. This is the same flat tail the earlier measurements traced to the usual fit's average over
tract frequencies. With the window at ±9, calling hands that tail to any far candidate.

## What ±9 costs

| fit | cohort | ±6 | ±9 | factor |
|---|---|---|---|---|
| counting | tomato2, all 103 classes | 19.6 s | 70 s | 3.6 |
| usual | tomato2, five classes of 25–274 tracts | 349 s | 1,163 s | 3.3 (1.06 to 5.4 by class) |
| usual | the drawn class above | 42 s | 468 s | 11 |
| usual | the oracle's four accessions, 20 regions | 11 min 04 s | 12 min 28 s | at most 1.13 |

The spec estimated 2.1, from the table width alone. The usual fit needs no more rounds at ±9: over
the five tomato2 classes it took 243 rounds against 247 at ±6, with one walk of 15 running out of
its 40 rounds at each width. So the extra time is the cost of each round. A whole tomato2 run of
the usual fit took 6.7 hours at ±6; at a factor of 3.3 it would take about 22 hours. That is an
estimate from five of 45 classes, not a run. The oracle's figure is an upper bound, because other
jobs were running during it.

## The oracle: four tomato2 accessions over 20 regions

The default-parameter outputs do not move. The two fitted outputs do, and the baseline is
re-recorded in `scripts/promote_ng_oracle.baseline`:

- **The fitted file:** 53 lines move, all repeat-tract rows. The 37 slip rates move by a median of
  1.2%, between 0.953 and 1.044 of their old values. These four accessions' classes are short.
- **The calls made with the fit:** QUAL moves at 302 repeat-tract records, by at most 3.7.
  - One record appears, a heterozygote one repeat short at QUAL 6.8.
  - One record loses its +1 allele: that sample's 2/2 call becomes 0/0.
  - One sample's genotype flips from 1/1 to 0/0, on 2 reference reads against 1.
  - Each genotype that changed was called at a genotype quality of 4 or less, before and after.

## What this settles, and what is still open

**For spec §5.1, whether the counting fit replaces the usual fit first:** at ±9 the usual fit costs
about 3.3 times as much on real classes, and its flat tails now reach candidates through the edge
split. The counting fit costs about a minute on tomato2.

Still open, as in the spec: the bias against long alleles from reads that must cross the whole
tract (§5.3), and whether ±8 is itself too narrow for the longest tracts (§5.4). tomato2's
dinucleotides of 25 and 27 repeats still hold 0.42 of their alleles at the edge.

The kimura measurement (2,169 accessions) has not been rerun at ±9.
