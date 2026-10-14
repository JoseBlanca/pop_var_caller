# A partial read is compared with the end of the tract it showed

*Implementation report, 2026-10-12, branch `str-partial-phase` (from main at a87ce7b6). The defect
was found in review on the STR-bundle branch
(`pop_var_caller-str-bundles/doc/devel/reports/reviews/str_bundles_partial_reads_2026-10-12.md`,
finding M4). Specs changed: [`read_likelihoods.md`](../../ng/spec/read_likelihoods.md) §5.2 and
[`locus_generation_ssr.md`](../../ng/spec/locus_generation_ssr.md) §3, each with a dated note. Every
benchmark run here used the simulator's true stutter, one group per library; no fit was run.
Scratch, not committed: `tmp/` in the branch's worktree.*

---

## 1. The defect

A **partial read** at an STR segment is one that held one flank and ran off its own end inside the
repeat. It shows a lower bound on the tract's length, and the caller scores it by asking how likely
the candidate allele is to be at least that long, times how well the letters match.

The generator stores a partial read's bases in read order. A read that held the **left** flank shows
the start of its carrier's tract; one that held the **right** flank shows the tract's **end**. The
scoring compared every partial read with the **start** of the candidate. For a right-held read whose
length is not a whole number of repeats, that puts every base out of phase: a `CA` tract read from
the right showing `ACACA` was compared with `CACAC`, five mismatches, about 4 in 10¹⁸ against every
candidate. The read then fell to the outlier term.

Measured in review on one simulated sample (PCR, 30×): 14,175 of 22,661 right-held partial reads at
segments were compared out of phase — none at period 1, half at period 2, three in four at periods 4
to 6 — which was 27 in 100 of all partial observations at segments. At a pure tract the read scored
the same tiny number against every candidate, so the loss was evidence thrown away rather than
evidence for a wrong allele.

## 2. The fix

- **The emission takes the side.** `SsrEmissionModel::censored_emission` gains a `HeldFlank` argument
  (`Left` or `Right`, a new type in `locus_generation`). A left-held read is scored as before. A
  right-held read is scored as the mirror image: the read, the candidate's tract and the motif are
  each written back to front and handed to the same left-held scoring. So the read is compared with
  the candidate's last bases, a part-repeat length change is resized at the tract's start instead of
  its end, and a whole-repeat slip is placed in each run as found from the right. A reversed pure
  tract is a pure tract of the reversed motif (`CACACA` under `CA` is `ACACAC` under `AC`), and the
  stutter model reads only the period and the repeat count, so nothing else changes.
- **The side comes from the witness.** The witness is the set of tract positions a read covered, laid
  from the border it anchored. Flush with the left border only means left-held; flush with the right
  border only means right-held (`HeldFlank::of_witness`). The row needs the tract's length for that,
  so `SsrSampleEvidence` now carries it (`tract_len`, taken from the record's region, which is the
  tract).
- **A witness that covers the whole tract cannot say.** Such a read ran past the far border without
  anchoring there, from one side or the other, and two such reads of the same bases are one
  observation in the psp. The row scores it from both sides, each weighted one half. On a pure tract
  of period 1 the two halves are equal.

The generator was not changed, so the psps are the same files; both binaries below called the same
psps.

## 3. Stuttering benchmark

Six runs, main (a87ce7b6) against this branch. A **cell** is one sample at one STR segment. For
variant cells: *right* is the genotype called right, *wrong* a different genotype (including the right
lengths with the wrong sequence), *missed* no call or a reference call. *False* counts reference cells
given a variant call. The 262 compound hexanucleotide segments (`AAAAAC`, `AAAAAG`, `AAAAAT`) are
left out, as in every earlier report on this benchmark.

| run | main right / wrong / missed / false | fixed right / wrong / missed / false |
|---|---|---|
| PCR, 1 sample, 5× | 508 / 198 / 419 / 14 | 508 / 198 / 419 / 14 |
| PCR, 1 sample, 15× | 1,031 / 47 / 47 / 17 | 1,033 / 45 / 47 / 17 |
| PCR, 1 sample, 30× | 1,116 / 3 / 6 / 3 | 1,116 / 3 / 6 / 3 |
| PCR, 50 samples, 5× | 40,074 / 8,257 / 8,290 / 1,441 | 40,196 / 8,179 / 8,246 / 1,406 |
| PCR, 50 samples, 30× | 55,978 / 334 / 309 / 76 | 55,990 / 325 / 306 / 74 |
| PCR-free, 50 samples, 30× | 56,381 / 114 / 126 / 11 | 56,401 / 101 / 119 / 11 |

By period, right variant cells, main → fixed:

| run | period 1 | period 2 | periods 3–6 |
|---|---|---|---|
| PCR, 1 sample, 15× | 204 → 204 | 192 → 192 | 635 → 637 |
| PCR, 50 samples, 5× | 7,374 → 7,379 | 6,682 → 6,682 | 26,018 → 26,135 |
| PCR, 50 samples, 30× | 10,502 → 10,501 | 9,789 → 9,792 | 35,687 → 35,697 |
| PCR-free, 50 samples, 30× | 10,672 → 10,673 | 9,983 → 9,983 | 35,726 → 35,745 |

The 1-sample 5× and 30× runs did not move a cell. False calls fell by 35 at 50 samples, 5×, and by 2
at 30×; no false call is gained anywhere.

**Cells whose verdict changed**, each way:

| run | variant cells to right | variant cells from right | other variant changes | false calls removed / added |
|---|---|---|---|---|
| PCR, 1 sample, 15× | 2 | 0 | 0 | 0 / 0 |
| PCR, 50 samples, 5× | 353 | 231 | 120 (between wrong and missed) | 55 / 20 |
| PCR, 50 samples, 30× | 27 | 15 | 0 | 3 / 1 |
| PCR-free, 50 samples, 30× | 30 | 10 | 1 | 0 / 0 |

The gain is small because, on a pure tract, the misread evidence scored the same against every
candidate and so mostly cancelled; what the fix recovers is the information those reads carry. At
50 samples and 5×, where a sample has the fewest spanning reads, the partial reads matter most, and
that is where the most cells move in both directions.

Many more records move than verdicts: at 50 samples, 30×, 2,270 of 5,708 records differ in some
field (`QUAL`, `GQ`, allele frequency), and at 1 sample, 30×, 571 of 2,409.

## 4. HG002

GIAB HG002's tandem-repeat benchmark, the STR segments of its tier (`TIER_KIND=ssr_locus`), one
sample, `--defaults`. Right genotypes out of the tracts both arms called:

| depth | class | main | fixed |
|---|---|---|---|
| 30× | homopolymer | 3,503 of 3,864 | 3,504 of 3,865 |
| 30× | period 2 and above | 2,566 of 2,865 | 2,567 of 2,864 |
| 50× | homopolymer | 3,595 of 3,902 | 3,594 of 3,902 |
| 50× | period 2 and above | 2,613 of 2,885 | 2,619 of 2,887 |

19 of 6,975 tracts changed verdict at 30× (10 to right, 8 from right, 1 between wrong kinds) and 10 at
50× (7 to right, 2 from right, 1 between wrong kinds). Homopolymers move because an interrupted run
of one letter is not its own mirror image. A pure run of one letter cannot move: a read reversed has
the same number of mismatches against every candidate, so its two sides score alike. *(Corrected in
review, 2026-10-12: this sentence also blamed whole-tract reads scored from both sides, which is
false for a pure run. On the stuttering benchmark, 18 of the 20 period-1 cells that changed verdict
at 50 samples and 5× have an interrupted true allele, as do the one cell that moved in each 30× run.)*

## 5. Tests

In `calling/likelihood/ssr_emission.rs`:

- `a_right_held_read_is_compared_with_the_candidates_end` — the `ACACA` case: from the right it scores
  exactly what `CACAC` from the left does, and read from the start it is below 10⁻¹⁵.
- `left_and_right_reads_of_a_pure_tract_score_alike` — at periods 1 to 6, three repeat counts and
  every read length, a tract's last `k` bases from the right score bit for bit what its first `k`
  from the left do.
- `a_period_one_read_scores_the_same_from_either_side` — a run of one letter, with and without a
  sequencing error: right-held scores equal main's.
- `mirror_image_reads_of_a_symmetric_interrupted_candidate_score_alike` — `AAAATAAAA`, the interrupted
  route's exact sum, a read and its mirror image score the same.
- `the_complete_read_emission_is_unchanged` — the complete-read emission pinned to its formula.

In `calling/likelihood/ssr.rs`: `a_partial_read_is_scored_from_the_flank_its_witness_held` — the row
asks for the left side, the right side, or one half of each for a whole-tract witness.

`cargo test --release --no-fail-fast --lib -- calling:: locus_generation:: run::callers`: 1,668
passed, 0 failed. Every test fixture's partial read is left-held; they now pass a tract length
(`TEST_TRACT_LEN`, 1,000 positions) long enough that their witnesses read as left-held, so none of
their expectations moved.

## 6. Cross-platform checksums

Two VCF checksums moved, and nothing else: in the two-sample cohort with sequencing errors the repeat
tract's `AF` moved in its sixth decimal (0.249923 → 0.249924); in the 21-sample cohort its `QUAL` went
37.2 → 37.3 and its `AF` 0.023717 → 0.023718. The parameters files did not move. With the row forced
back to reading every partial from the left, all four digest tests passed on this branch, so this
change alone moves them. Re-recorded in `src/cli/cross_platform_digests.rs` (Linux container, arm64,
glibc only).

## 7. Assumptions and what is left

- **A whole-tract witness is scored one half from each side.** Keeping the side in the observation
  would remove the guess, but it is a psp format change, and the STR-bundle branch is adding exactly
  that field (`held_flank`, psp 2.1) for bundles. Once it lands, a segment can read it too.
- **Partial reads and flank-deletion alleles.** A partial read records no flank deletion, and the
  censored term compares it with the allele's spelling as stored, as it did before; a right-held read
  against an allele that deletes right-flank bases is compared with the end of that spelling. Not
  measured here.
- `SsrSampleEvidence::new` takes the tract length as a third argument; `LocusLen::from_positions` is
  now `const`.
