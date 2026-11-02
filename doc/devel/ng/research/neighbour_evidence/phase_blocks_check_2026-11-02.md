# Writing phase in blocks: the prototype check (2026-11-02)

*The check that [`neighbour_evidence_phase_output.md`](../../spec/neighbour_evidence_phase_output.md)
§7.1 asks for before the Rust build. Python prototype of the pass (`phaser.py`, 5 rounds, 16
haplotype templates one per sample, rounds 2–5 averaged), founder-simulator cohorts, scored against
the simulator's truth. Scripts: `phase_blocks.py` (block rules), `phase_joins.py` (joining windows),
`phaser.py --extra-draws 8 --region`. Dumps in `tmp/phase_check/` of the `neighbour-pass` worktree.*

## Summary

**The rule passes the 1% ceiling on six of the nine cohorts and fails on all three outbred ones.**
Five of the six are the inbred and mixed-depth cohorts of 50 and 200 samples, at 16 and 64 founders.
The sixth, 10 inbred samples, passes by writing no phase at all, which is the intended behaviour.
The outbred cohorts:

- 50 samples from 16 founders: 1.2% of the links inside blocks are wrong, just over the ceiling;
- 50 samples from 64 founders: 7.8%;
- 10 samples from 64 founders: 36%.

Neither more draws nor a stricter rule fixes the last two: the model there is confident and wrong.

- A *flip* is a wrong link: two consecutive phased heterozygous records of a sample whose relative
  phase is not the truth's.
- **What passes** is a stricter form of the adopted rule. Besides the round's own path, each round
  draws 8 more paths from the same forward probabilities. A link is kept only if all 36 paths (4
  rounds × 9) give the two records the consensus's relative phase, and every path leans by at least
  0.5 at both records.
- **The rule the owner adopted** (every round's path agrees and leans) misses the ceiling on five
  of the nine cohorts.
- **Joining windows works.** 53 joins were made on 4 cohorts, and none was wrong.
- **The share of paths that disagree at a link is not yet an error rate.** It ranks links within a
  cohort, but the links where all paths agree are wrong 0.06% of the time on one cohort, 12% on
  another and 37% (23 of 62) on 10 outbred samples. So the ceiling cannot yet be a user's
  parameter.

## The cohorts and the rules

Records are phased where they are truly heterozygous (the stand-in for a written heterozygote) and
the consensus of the rounds leans by at least 0.5. Links longer than 300 kb are broken, and blocks of
one record are written unphased. "Phased" counts only records in blocks of two or more.

| cohort | true hets phased: round rule / 36-path rule | flips inside blocks: round rule | flips inside blocks: 36-path rule | block N50 in records (kb), 36-path rule |
|---|---|---|---|---|
| 50 samples (10 at 30×, 40 at 3×), inbreeding 0.5, 16 founders | 70% / 57% | 16 of 4,827 (0.33%) | 2 of 3,757 (0.05%) | 13 (78) |
| same, inbreeding 0.8 | 73% / 61% | 19 of 2,244 (0.85%) | 2 of 1,811 (0.11%) | 17 (90) |
| same, inbreeding 0.8, 64 founders | 43% / 28% | 48 of 1,361 (3.53%) | 3 of 859 (0.35%) | 14 (80) |
| 200 samples (40 at 30×), inbreeding 0.8, 64 founders | 58% / 46% | 120 of 8,271 (1.45%) | 18 of 6,419 (0.28%) | 18 (113) |
| 50 samples at 3×, inbreeding 0.5, 16 founders | 66% / 48% | 30 of 4,396 (0.68%) | 3 of 3,025 (0.10%) | 9 (61) |
| 10 samples at 3×, inbreeding 0.95, 64 founders | 2% / 0% | 3 of 3 | 0 of 0 | — |
| 50 samples at 3×, outbred, 16 founders | 34% / 17% | 114 of 4,297 (2.65%) | 24 of 2,054 (1.17%) | 4 (34) |
| 50 samples at 3×, outbred, 64 founders | 15% / 7% | 263 of 1,930 (13.6%) | 65 of 836 (7.8%) | 4 (31) |
| 10 samples at 3×, outbred, 64 founders | 2% / 2% | 25 of 66 (38%) | 15 of 42 (36%) | 13 (105) |

On the cohorts that pass, deep and shallow samples both pass. The 200-sample cohort has 4 flips in
1,513 links between deep samples and 13 in 4,564 between shallow ones. Links that touch a record
filled without reads are no worse: 1 flip in 342 links on the same cohort.

**How many extra draws.** A sweep of 0, 2, 4 and 8 extra draws per round (`EXTRA=` in
`phase_blocks.py`) gives these flips inside blocks:

| cohort | 0 | 2 | 4 | 8 |
|---|---|---|---|---|
| 16 founders, inbreeding 0.5, mixed depth | 0.33% | 0.11% | 0.05% | 0.05% |
| 64 founders, inbreeding 0.8, 50 samples | 3.53% | 1.45% | 1.43% | 0.35% |
| 64 founders, inbreeding 0.8, 200 samples | 1.45% | 0.55% | 0.35% | 0.28% |

The 64-founder, 50-sample cohort needs all 8. Its counts are small (16, 14 and 3 flips), so the step
from 4 to 8 is uncertain.

**A stricter lean** (0.6 or 0.8 instead of 0.5) is in the `.out` files. It does not change which
cohorts pass.

## Why the outbred cohorts fail

At 50 outbred samples from 64 founders, 1,977 links have every one of the 36 paths agreeing, and
244 of them (12%) are wrong. At 10 outbred samples, 62 such links have 23 wrong (37%). The draws are
samples of the model's own uncertainty given its templates. Here the templates do not hold the
sample's true haplotypes: 64 founders give 128 founder haplotypes, and 50 samples carry most of them
once or twice. The model then copies a wrong template with confidence. No count of draws or rounds
measures that, because every draw comes from the same model.

The prototype also fixes the copying length at 300 kb. The Rust pass fits it per sample, choosing
30 kb or 300 kb, and on an outbred cohort it may choose 30 kb, which breaks links over 30 kb apart.
Whether that rescues the outbred cohorts is unmeasured. It cannot be measured until Rust writes
phase.

## The share of disagreeing paths, as an error rate

These are links within 300 kb, before any rule, binned by the share of the 36 paths that disagree
with the consensus at the link:

| cohort | none disagree | ≤ 3% | 3–10% | 10–25% | > 25% |
|---|---|---|---|---|---|
| 16 founders, inbreeding 0.5, mixed depth | 0.06% | 0.77% | 2.3% | 4.4% | 32% |
| 64 founders, inbreeding 0.8, 200 samples | 0.62% | 4.4% | 6.1% | 14% | 37% |
| 64 founders, outbred, 50 samples | 12% | 21% | 22% | 25% | 39% |

The share orders the links within every cohort, but the level moves 200-fold between cohorts. A
threshold on it would not mean one error rate, so the ceiling cannot yet be the user's parameter.

## Joining windows

Each contig (1 Mb) was phased as two windows, 0–600 kb and 400 kb–1 Mb. The first window writes
below 500 kb, and the second's margin is 400–500 kb. Joins used the 3 records of the last written
block nearest the boundary (spec §2.4), under the 36-path rule:

| cohort | joins tried | joined | wrong joins |
|---|---|---|---|
| 16 founders, inbreeding 0.5, mixed depth | 43 | 20 | 0 |
| 64 founders, inbreeding 0.8, 50 samples | 12 | 2 | 0 |
| 50 at 3×, outbred, 16 founders | 86 | 1 | 0 |
| 64 founders, inbreeding 0.8, 200 samples | 83 | 30 | 0 |

No join had the two windows disagreeing on orientation at the shared records. Joins fail because
the last block does not reach the boundary with 3 records in a 100-kb margin. The real margin is
1 Mb, so more will hold. With 2 records, joins rose to 23, 3, 3 and 35, still with none wrong. The
spec's 3 stays.

## What follows

- **Adopt the 36-path rule** as the spec's link rule. It is the adopted rule made stricter, and it
  is what passes. In Rust the 8 extra draws could share one pass per record: the cumulative
  weights of a record's states built once, each draw then picking from them. Cost and memory are
  unmeasured. (The Rust build draws each separately; its cost is in
  `doc/devel/reports/implementations/phase_output_2026-11-02.md`.) Storing each path's orientation bit at every record a sample could phase is about 5
  bytes a record, or 0.6 GB at 3,000 samples and 39,000 records.
- **The outbred cohorts need a check that does not come from the model's own draws.** That is the
  owner's decision (see the chat report of 2026-11-02).
