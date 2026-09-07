# B1 — measured on the real whole-genome CRAM

**Date:** 2026-09-07
**Branch:** `ng-cram-window`, at `bd789ab7`, against its merge base `52b7b787`
**Plan:** [cram_reference_window.md](../../ng/impl_plan/cram_reference_window.md), Milestone B
**File:** `benchmarks/tomato_big_cram/DRR000741.p1.cram` — one tomato accession, whole genome,
49 GB, coordinate-sorted, CRAM 3.0 from samtools 1.21, 112,140 blocks

## Verdict

**A whole run over 10 Mb of chromosome 1 holds 129 MB less and writes the same VCF, byte for
byte.** Peak resident goes from **277.8 MB to 148.8 MB** at the best of six runs a side, 41 % less
at the median. **No wall-time difference is measurable**: the two minimums are 0.9 s apart on runs
of about 56 s, inside a run-to-run spread of 20 s caused by other work on the machine.

The decode alone, measured on its own, holds **185.9 MB against 11.4 MB** — a factor of 16 — and is
if anything slightly quicker.

## What the whole run costs, either way

`pop_var_caller_exp call-from-alignments --threads 1 --defaults`, one sample, over
`SL4.0ch01:1-10,000,000`. Six runs a side, alternated branch/base so drift in machine load falls on
both. Peak resident is `getrusage(RUSAGE_CHILDREN).ru_maxrss` of the child — the same counter
`/usr/bin/time -v` reports, which is not in the dev container.

| | branch, min | median | max | base, min | median | max |
|---|---:|---:|---:|---:|---:|---:|
| peak resident, MB | **148.8** | 164.2 | 176.1 | **277.8** | 277.8 | 278.4 |
| wall, s | 55.9 | 57.5 | 76.0 | 55.0 | 65.6 | 78.0 |

**The memory figure is solid and the timing figure is not, and the table says which is which.**
The base's peak sits in a 0.6 MB band across all six runs — it is a chromosome, and a chromosome is
the same size every time. The branch's varies over 27 MB, because what is left is the calling
machinery rather than one big fixed allocation. Every branch run is below every base run by at
least 100 MB.

Wall time is a different story: **two other sessions were building on this machine throughout**
(load average 9.2 falling to 5.7), and both sides show a 20-second spread with their slow runs
landing in different rounds. On the minimums the branch is 0.9 s slower, which is 1.6 % and far
inside that noise. **The honest statement is that this measurement cannot see a wall-time
difference**, not that there is none; the clean timing evidence is the decode measured on its own,
below, where the windowed path is 5 % *faster*.

## The VCF is byte-identical

All three of the first three pairs, compared line for line with only `##commandline` and
`##parametersFile` excluded: **0 differing lines**, over 23,838 lines. Two runs of the branch
against each other likewise differ nowhere, so the comparison is not hiding a run-to-run wobble.
The run calls 2,862 repeat tracts and keeps 11,170,190 reads of 13,902,133.

## The decode on its own

`ng_cram_decode_layers`, 60 containers of `SL4.0ch01`, 600,000 records, fastest of seven, each pass
run in its own process so the memory figure is that pass's rather than the union of all of them.

| | whole contig resident | per-block window |
|---|---:|---:|
| peak resident, MB | 185.9 | **11.4** |
| seconds for 60 containers | 0.502 | **0.477** |

**16.3× less resident, and 5 % quicker.** The windowed pass fetches 0.5 Mb of bases across those 60
blocks where the other holds all 90.9 Mb of the chromosome. This is the figure the branch exists
for, and it is measured on the real file rather than on a fixture — the research note's §2 showed
the benchmark CRAMs are a different workload on this path, their blocks spanning megabases.

## The digest the plan named does not exist, and never did

**The plan required `a127ecf5083f535a002a5f461f150ede` from both passes. The measured value is
`c0bdbb0e464a920b966ad487fc3ca678`.** That is not a regression, and the check that establishes it
is worth recording because "the oracle changed" is exactly what a regression looks like at first
glance.

Three builds, same file, same contig, same 60 containers, same 600,000 records:

| build | digest |
|---|---|
| this branch, `bd789ab7` | `c0bdbb0e464a920b966ad487fc3ca678` |
| merge base, `52b7b787` | `c0bdbb0e464a920b966ad487fc3ca678` |
| `48a93310` — **the commit that introduced the harness, which is where the research note's number came from** | `c0bdbb0e464a920b966ad487fc3ca678` |

And within each build, all four passes agree: through `RecordBuf`, read directly, tags discarded,
and decoded against a per-block window. So the decode has produced the same 600,000 records at
every point in this history, including at the commit whose own report quoted a different constant.

`hash_record` is byte-identical between `48a93310` and now, and so is `container_offsets`, so
neither the field set nor the block selection moved. **The number in
[`cram_read_path_2026-09-04.md`](../../ng/research/cram_read_path_2026-09-04.md) §1 was wrong when
it was written** — recalled rather than read off a run — and the plan copied it. A second sign it
could not be right: that sentence claims the same digest covers *both* files the note measured, and
two different CRAMs cannot hash their records to one value.

Both documents are corrected to the measured value, with the command that produces it.

## What was not measured

- **A fragmented reference.** Milestone A′ made a block spanning several chromosomes decodable, and
  this file has none — 0 in its 112,140 blocks. What the two-pass path costs on a draft assembly of
  many short contigs is unmeasured, and would need such a CRAM.
- **A cohort.** One sample, one thread. The memory this saves scales with open files (each holds
  its own window instead of sharing one chromosome cache), so 63 samples is where the shape of the
  saving changes, and it is not measured here.
- **Wall time cleanly.** Worth re-running the six rounds on an idle machine if the 1.6 % matters;
  nothing about the change predicts a cost, and the isolated decode is faster.
- **A real whole-genome human CRAM's block spans.** Still the open question spec §10 point 3 names,
  and still the thing that would decide whether the shared cohort window in §12 is ever wanted.
