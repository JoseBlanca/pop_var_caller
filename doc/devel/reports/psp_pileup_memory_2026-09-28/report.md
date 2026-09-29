# Calling 2,169 tomato samples runs out of memory at read pile-ups

*2026-09-28. Measured on the server `kimura` (96 cores, 125 GB of memory, 23 GB of swap) over
the tomato25 cohort: 2,169 accessions at about three reads a position, stored as psps (one file of
stored observations per sample, 8.6 TB in all) on one USB spinning disk.*

## In short

`call-from-psps` is killed by the kernel's out-of-memory killer at the same place on every run:
**SL4.0ch00, just after position 1,504,000**, about 14 minutes into the calling. Up to that point
it holds a steady 43 GB.

**The cause is one copy per record of a list that can be thousands of entries long.** For every
record it holds, a sample's reader keeps the identifiers of all the reads that span that record's
position — the *live set*, which the reader needs later to decode the record's body. It keeps a
**full copy of that list for each record**. Almost everywhere the list is a handful of reads, so
the copies cost nothing. At a pile-up they do not: from ch00:1,504,001 a median of **5,669 reads**
span each position in each sample, against about three elsewhere, so every record carries thousands
of 8-byte identifiers. One sample's reader was measured holding **20 million identifiers for 8,816
records**, about 160 MB, and the process grew by about 145 MB for every sample that loaded the
stretch — some 300 GB had all 2,169 been allowed to finish.

A position with thousands of reads in a three-reads-a-position sample is almost certainly a
collapsed repeat: many copies of a sequence in the genome, mapped onto one copy in the reference.
The unplaced scaffolds in ch00 are where such stretches are expected.

**The owner's decision (2026-09-28): fix it with a depth ceiling** — drop, as it is read, a record
whose depth is far above what its sample normally has, and report how much ground was dropped.
Where the "normal depth" comes from is the open question; see [the fix](#the-fix).

## What happened

- **10:37** — the pipeline started `call-from-psps` with `--cohort-locus-builder-regions-len 16000`.
- **10:51:56** — killed. The tmux session it ran in was killed with it. The session's systemd log:
  `Failed with result 'oom-kill' … 123.6G memory peak, 23.1G memory swap peak`.
- The last call written was at ch00:1,503,980. Every later run died with the same last call.

The psp prefetch merged earlier the same day (commit `3925440e`) was running, and **is not the
cause**: what it reads is the kernel's file cache, which the kernel takes back before it kills
anything, and a run without it showed the same climb at an earlier, smaller pile-up (next section).

## How it was found

Each step was a run of an instrumented build on the real cohort under a 100 GB memory cap
(`systemd-run --user --scope -p MemoryMax=100G -p MemorySwapMax=0`), so that only the test was
killed. The instrumentation is kept as [`diagnostics.patch`](diagnostics.patch).

What was ruled out, and on what evidence:

| suspected | evidence against it |
|---|---|
| overlapping variants chaining across a long stretch (the merge loads until no observation reaches past what it holds, and nothing bounds that — spec `cohort_merge.md` §4.1) | the longest chain reached about 600 bases past its round |
| individual oversized records | none over 1 MB |
| a record claiming to reach far ahead | no record moved the reach by more than 50 kb |
| a single sample holding too many records | none held more than 100,000 |
| a locus with too many candidate alleles | every call up to ch00:1,503,986 had two or three alleles, the longest 11 bases, memory flat at 42.4 GB |
| the writing of calls | it holds nothing between records |

What pinned it down: a watchdog thread printing memory every half second, with the stage the round
loop was in. Memory climbed **during the loading of the round that starts at 1,504,001**, sample
by sample — 43 GB with no sample loaded, 107 GB after 443 of the 2,169 — about 145 MB per sample.
Counting each sample's live-set arena then gave the numbers above.

**The same effect, smaller, at ch00:0.43 Mb.** There the live sets are 400 to 1,000 reads, and the
earlier eight-minute timing runs went from 16 GB to peaks of 59 GB (without the prefetch) and 67 GB
(with it) crossing that stretch. That is the same mechanism, not a separate problem.

## Where it is in the code

`PspSummarySource::next_summary` in `src/run/psp_source.rs`:

```rust
let ids = self.walk.live_reads().ids();
let live = LiveSpan { start: …, len: ids.len() as u32 };
self.live_ids.extend_from_slice(ids);
```

Every record appends the whole live set to the sample's `live_ids` arena. The arena is released
with the records that point into it (`release_before`), so it is bounded by what the merge holds —
about one round, 16,000 bases here — but within that it grows as **records × reads spanning each
record**.

**Why no earlier run showed it.** The cost is records × depth × samples. On the 63-accession
benchmark the same pile-up costs 63 × 160 MB, about 10 GB at worst, and the benchmark's regions
may not include ch00's pile-ups at all. At 2,169 samples it is 300 GB.

## The fix

### Decided: a depth ceiling

A record whose depth is far above its sample's normal depth is dropped by the reader as it is read,
and the run reports how much ground was dropped that way per sample. This is the maximum-depth
filter other callers apply, and it bounds memory wherever pile-ups occur.

**Proposed ceiling**: 20 times the sample's mean depth, never below 100 reads — so a 3× sample is cut
at 100 reads and a 30× sample at 600. It has to be relative to each sample for the caller to work
across its whole depth range: an absolute number low enough to protect a 3× cohort would cut
ordinary positions in a 300× sample.

**Open: where the sample's normal depth comes from.** The owner's constraints:

- it must work in both modes — from psps, and straight from alignment files
  (`call-from-alignments`), where nothing about the sample has been measured beforehand;
- reading a whole BAM or CRAM to measure it is too costly;
- reading only the start of a sorted file is not representative of the genome.

Not yet checked: whether a psp's header or census, or the parameters file, already records a mean
depth per sample or per read group. If it does, psp mode has its number for free, and the question
is only what direct mode does. One candidate that needs no reading ahead is a **running estimate**
kept as the records go past — the psp source already counts how deep the records it drew were
(`StoredSampleTallies`). Its weakness is the start of each run, before enough ground has gone by,
which is precisely the owner's objection to reading the start of the file; it would need a floor
that protects the first stretch.

### Later: stop copying the list for every record

Store only what changes in the live set from one record to the next, plus a full copy every so many
records, and rebuild the set when a body is decoded. The output is byte-identical and it helps at
every depth, including high-coverage cohorts where hundreds of reads legitimately span every
position. **On its own it does not fit this cohort into 80 GB**: at 5,669 reads spanning a
position and reads of about 150 bases, some 76 reads start or end at each base, so the stored
changes still come to an estimated 15–25 MB per sample in this stretch — 30–50 GB on top of the
43 GB baseline. It is the structural fix; the ceiling is what makes the run fit.

## Testing the fix on the cohort

Apply [`diagnostics.patch`](diagnostics.patch) to a build for testing only. With it,
`PVC_DIAG_FROM=<position>` skips building and calling every round on ch00 that ends before
`<position>` — the files are still read up to there, but nothing is called. From 1,480,000 it
reaches the failure point in about 8 minutes, against 15 for a run from the start:

```
systemd-run --user --scope -q -E PVC_DIAG_FROM=1480000 -p MemoryMax=100G -p MemorySwapMax=0 \
  <build>/pop_var_caller call-from-psps \
  --reference ~/tomato25/genome/S_lycopersicum_chromosomes.4.00.fa \
  --catalog ~/tomato25/genome/S_lycopersicum_chromosomes.4.00.fa.repeats.parquet \
  --psp /media/tomato25_vcfs/snv_calling/psps \
  --output tmp/diag/variants.vcf.gz \
  --parameters /media/tomato25_vcfs/snv_calling/parameters.toml \
  --ploidy 2 --cohort-locus-builder-regions-len 16000
```

The fix works if the run passes ch00:1,520,000 with memory staying near 43 GB. A test far past
that is worth one run from the start, since ch00:0.43 Mb has a pile-up of its own.

## Other facts from the same session

- **Round shape.** With no setting, a cohort above about 1,000 samples gets one region in flight
  (`round_shape_for`, `src/cli/calling_run.rs`), so building and genotyping run on one core. The
  pipeline now passes `cohort_locus_builder_regions_len = 16000` (32 regions of 500 bases). The
  default is worth revisiting.
- **The psp prefetch** (`--psp-prefetch-bytes`, default 8 GiB, commit `3925440e`): with it, the
  same stretch of genome took 80 s instead of 130 s once its look-ahead had filled, the disk read at
  about 200 MB/s instead of 85, and the first 2.8 GB of calls were byte-identical to a run without
  it. It fills the file cache further than its budget at the start of a run (about 27 GB), because
  the disk's own read-ahead (`read_ahead_kb`, set to 8192 on `sdb` that morning and not persistent)
  reads beyond each of its requests.
- **The owner's memory limit** for the caller on this server is 80 GB.
