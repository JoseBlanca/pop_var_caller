# Starting the SNP/indel fit from an earlier run — built and measured, 2026-11-07

*Branch `fit-warm-start`. Design: [`../ng/spec/parameter_prepass_joint_fit.md`](../ng/spec/parameter_prepass_joint_fit.md)
§3.3.3 and the format-10 amendment of [`../ng/spec/parameters_file.md`](../ng/spec/parameters_file.md).
Measured on kimura, 400 tomato accessions, 48 threads. The kimura psps are format 1.0, which main
cannot read, so the runs used an experimental build carrying the same fit code (tag
`fit-subset-experiment`).*

## 1. The answer

**When a few samples join a cohort that was already fitted, the SNP/indel fit takes less than half
the time if it starts from the earlier run's numbers, and it reaches the same answer.** On the 400
kimura samples, adding 20 samples that the earlier fit had never seen took the fit 29 minutes,
against 66 minutes from scratch.

| SNP/indel fit on 400 kimura samples | rounds over every position | fit time | whole `estimate-parameters` run |
|---|---|---|---|
| from scratch (three starting points, then the winner finished) | 132 in the finish | 66 min | 1 h 20 min |
| from an earlier fit of the same 400 | 9 | 4 min | 18 min |
| from an earlier fit of the first 380 by name: 20 new samples, 21 new read groups | 64 | 29 min | 43 min |

Against the fit from scratch, the run that added 20 new samples ends at:

- **log-likelihood** −2.543396e8, against −2.543397e8;
- **mismapped share** 5.844%, against 5.841%, about an eighth of its standard error (±0.020%);
- **error rate at ordinary positions**, per read group: at most 0.04 standard errors away;
- **homozygote excess**, per sample: at most 0.33 standard errors away;
- **error rate at mismapped positions**, per read group: 0.03 standard errors away at the median,
  1.3 at worst. A fit restarted on the same 400 samples moves these rates as far (1.1 at worst), so
  the spread comes from where the stopping rule lets the fit stop, not from the warm start.

**What the saving depends on.** The 380-sample fit already put the mismapped share at 5.88%, so the
cohort-wide numbers barely had to move. Most of the 64 rounds were spent in the slow finish, where
each Newton step is accepted only at 1/8 or 1/4 of its length. An addition that moves the
cohort-wide numbers further — 300 samples growing to 400, say — will cost closer to a fit from
scratch. That slow finish is now the larger cost in both kinds of run.

## 2. What was built

- **The parameters file records the fit's own numbers** (format 10, table `[snp_indel_fit]`): the
  mismapped share, the frequency density, and per read group the two error rates and per sample the
  homozygote excess. The mismapped share and the error rate at mismapped positions were not written
  anywhere before. No calling run reads the table; calling runs carry it through unchanged.
- **`estimate-parameters --start-from <parameters file>`** starts the fit from that table. Samples
  and read groups are matched by sample, @RG ID and library; anything the earlier file lacks starts
  at the median of what it holds. When the earlier fit holds at least half the samples, the one run
  from its numbers replaces the starting points. Below half, every starting point runs as before.
  That floor was chosen, not measured. A file fitted against another reference, or one without the
  table, is refused.

## 3. The idea it replaced: fitting the cohort-wide numbers on a subset of samples

Before this, the same wall time was attacked by fitting the cohort-wide numbers on a subset of the
samples, then the per-sample numbers on all of them. The cohort-wide numbers do not settle with
subset size, so the answer depends on how many samples were used:

| cohort-wide numbers fitted on | mismapped share | SNP/indel fit time |
|---|---|---|
| all 400 | 5.84% | 66 min |
| 200 | 4.28% | 60 min |
| 100 | 2.41% | 12 min |

Using the 200-sample answer only as the starting point for a fit of all 400 reached the same answer
as the fit from scratch, but took 145 rounds and 2 hours, so it saved nothing. Not adopted.

## 4. A side observation

The reported standard error of the mismapped share differed fourfold between runs that reached the
same maximum: ±0.020% from scratch and in the 380-to-400 run, ±0.040% in the 200-then-400 run,
±0.076% in the two runs started from the full 400-sample fit. That points at the same badly
conditioned curvature that makes the finish slow, and belongs with that work.
