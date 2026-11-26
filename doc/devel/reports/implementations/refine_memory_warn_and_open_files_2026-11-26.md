# Refinement's memory check only warns; the open-file check counts what pvc holds (2026-11-26)

Two changes to what `pvc` checks before a run, both from the 2026-11-25 runs on rick. Branch
`mem-warn`, from main at 4b29ce57. Measured on the macOS laptop; nothing ran on kimura or rick.

1. **The refinement memory check never stops a run; it warns.** It stopped above 90% of the
   memory the run may use (`refine_peak_guard_2026-11-25.md`). On rick, at 100 samples and 32
   threads, the estimate before calling was 30.7 GB while refinement used at most 3.5 GB
   (`purge_delay_rick_2026-11-25.md`): an estimate that far above the use can stop a run that
   would have fitted, and a stop at a window comes after what may be days of calling. A run the
   system kills for memory after calling has ended keeps its scratch directory with the resume
   point, so `pvc resume-calling` with a smaller budget finishes it anyway (owner's decision).
2. **The open-file check, which already existed, counted a third of what pvc holds.** 157 CRAMs
   on rick passed it (817 needed, by its count) and failed part-way with "Too many open files"
   under the default soft limit of 1,024. It now counts from measurement, says how to raise the
   limit, and the two psp-reading commands check too.

## 1. The memory warning

**When it warns.** Where the expected peak is above **85%** of the memory the run may use (the
smaller of the machine's memory and the control group's limit, as before). The estimate's formula
and constants did not change; `STOP_SHARE` (90%) became `WARN_SHARE` (85%).

- **Before calling**: one warning, saying that the estimate is an upper bound — no window is
  planned yet, so each is taken to be as large as `--refine-window-bytes` and the longest contig
  allow — and that each window is estimated again.
- **At a window**: only where the window's estimate is higher than every window's before it.
  Every window over the share would print a line for each of thousands of windows; only the first
  would hide a later, worse one. Windows filled to the budget count nearly the same, so a new
  highest is rare: of n windows in no particular order, about ln n are (7 of 1,000). The stage
  test confirms fewer warnings than windows when every window is over the limit.

**What it says**: the estimate and its parts, the limit and where it came from, the largest budget
within 85% and, where fewer threads would fit at this budget, the thread count, and how to finish
a killed run. The scratch directory is named in the command: the check now runs just after the
scratch directory is made, still before calling, and the stage carries it to every window. With
the hidden-duplication filter off there is no scratch directory, refinement runs while calling
does, and the warning says a killed run has to be started again.

Before calling, tomato25's cohort (2,176 samples, 48 threads) at a 100-GB budget under a 118-GB
control-group limit (from the unit test, exact):

```text
refinement memory, before calling: about 131.8 GB expected at its peak, 112% of the 118.0 GB this
run may use (the memory.max of the control group pvc runs in), above the 85% at which pvc warns;
the run goes on. This is an upper bound: no window is planned yet, so each is taken to be as large
as --refine-window-bytes and the longest contig allow; each window is estimated again once it is
planned. The estimate is 100.0 GB for a window as --refine-window-bytes counts it (100.0 GB
allowed), 13.6 GB for 48 threads' working arrays (270,328 walked records at 1.05 kB each a thread),
and 16% on top for the memory allocator. Within 85% would be --refine-window-bytes 76000000000
(76.0 GB) or less. If the system kills the run for memory after calling has ended, its scratch
directory keeps the calls, and `pvc resume-calling --scratch-dir /scratch/pvc_tmp_ab12
--refine-window-bytes 76000000000` finishes it.
```

At a window, a fixture of 100 samples on 96 threads whose arrays dominate (from the unit test,
exact):

```text
refinement memory, window from contig 0 position 1: about 46.7 GB expected at its peak, 117% of
the 40.0 GB this run may use (the memory.max of the control group pvc runs in), above the 85% at
which pvc warns; the run goes on. The estimate is 10.0 GB for a window as --refine-window-bytes
counts it (20.0 GB allowed), 30.2 GB for 96 threads' working arrays (300,000 walked records at
1.05 kB each a thread), and 16% on top for the memory allocator. Within 85% would be
--refine-window-bytes 14000000000 (14.0 GB) or less, or 61 threads or fewer (--threads). If the
system kills the run for memory after calling has ended, its scratch directory keeps the calls,
and `pvc resume-calling --scratch-dir /scratch/pvc_tmp_ab12 --refine-window-bytes 14000000000`
finishes it. Later windows warn again only where their estimate is higher still.
```

At the 80-GB budget the tomato25 replay ran at (103.0 GB measured peak, glibc), the estimate is
105.4 GB, 89% of 118 GB: that run is now warned, and goes on, where under the 90% stop it ran
silently. Within 85% the largest budget is 76 GB.

**What was removed**: `--refine-ignore-memory-limit` (all three commands; added 2026-11-25,
unreleased), `StageStop::OverMemory`, `VcfWriteError::RefinementOverMemory`,
`CallingRunError::RefinementOverMemory`. `PeakTooClose` became the warning `PeakNearLimit`;
`MemoryGuard` became `MemoryWatch`, which carries the scratch directory instead of the ignore flag.

**What stays**: the stop for a window over `--refine-window-bytes` even at its smallest size —
that one is certain — and its scratch directory kept whatever `--keep-tmp-on-error` says.
`refinement_stopped_on_memory` is renamed `refinement_stopped_over_budget`, and the kept
directory's log line now reads "because a refinement window did not fit --refine-window-bytes and
calling had finished".

## 2. Open files

### 2.1 What pvc holds open

Measured by sampling `lsof -p` every 0.1 to 0.2 s through whole runs on tomato1 (the BAQ-free
CRAMs, 20 regions of 100 kb), and accounted for in the code:

| Run | Files | Samples | Threads | On the alignments | On the FASTA | Other | Peak |
|---|---|---|---|---|---|---|---|
| call-from-alignments, CRAM | 10 | 9 | 1, 2, 4, 8 | 40 | 108 | 5 | 153 |
| call-from-alignments, CRAM | 30 | 29 | 1, 2, 4, 8 | 120 | 328 | 5 | 453 |
| call-from-alignments, CRAM | 63 | 57 | 1, 32 | 252 | 671 | 5 | 928 |
| call-from-alignments, BAM (the first 10 CRAMs converted) | 10 | 9 | 4 | 10 | 68 | 5 | 83 |
| call-from-psps | 9 psps | 9 | 1, 8 | 9 | 2 | 5 | 16 |
| call-from-psps | 29 psps | 29 | 1, 8 | 29 | 2 | 5 | 36 |
| estimate-parameters | 9 psps | 9 | 1, 8 | 9 | 1 | 4 | 14 |
| generate-psps (one sample at a time) | 10 or 30 | — | default | 8 | 19 | 5 | 32 |

Some files share a sample (63 CRAMs are 57 samples). The peak did not change with the thread count
in any row. "Other" is the standard streams, the log and the scratch directory.

**Where it comes from** (code read; `src/run/callers.rs` has the citations in its constants'
docs). Direct mode builds one walker for every sample before it reads, and the walkers stay open
for the whole run; the threads share them. A walker holds three locus generators — the pileup, the
repeat tract and the repeat bundle — and each keeps one cursor on every file of its sample.

- **A CRAM file: 11.** 4 CRAM readers: one a cursor, and the file's container prefetcher. 7 FASTA
  readers: two a cursor (the mismatch filter's, and the decode's, since a CRAM stores reads as
  differences from the reference) and one for the prefetcher.
- **A BAM file: 4.** 1 BAM descriptor, which every cursor reads through a shared block cache, and
  one FASTA reader a cursor, for the mismatch filter.
- **A sample: 4 more**, the FASTA accessors its generators hold, whatever its file count.
- **Indexes** (`.crai`, `.bai`, `.csi`, `.fai`) are read into memory and hold nothing.

The formula fits the table exactly on the FASTA column: 63 CRAMs × 7 + 57 samples × 4 = 669,
measured 671 with the run's padding-bases reader and one more; 10 × 7 + 9 × 4 = 106, measured 108;
BAM 10 × 3 + 9 × 4 = 66, measured 68.

### 2.2 The formula

```text
needed = 11 × CRAM files + 4 × BAM files + 4 × samples + 32
```

The per-file and per-sample terms are what the code opens, not estimates. **The 32 is the margin**:
5 to 7 were open in the runs above, so 25 or more are slack for the output, the catalog and
whatever the runtime opens on its own. At one CRAM a sample (15 each), Linux's default soft limit
of 1,024 is reached at 66 samples and macOS's 256 at 14; at one BAM a sample (8 each), at 124 and
28. rick's 157 CRAMs need 2,387. psp commands need one a psp and 32.

Checked against the operating system on the 63-CRAM cohort, refinement on, 8 threads: with
`ulimit -n 953` (the formula's need) the run finished; with 952 it was refused before reading.

The old constants were 3 a CRAM, 2 a BAM and 2 a sample: arithmetic from a time when a sample had
one generator. They undercounted a CRAM sample threefold.

### 2.3 The check

At start-up, before any data is read: in `call-from-alignments` (the existing check, before any
alignment file is opened), and now in `call-from-psps` and `estimate-parameters` (once the psps are
listed). `resume-calling` opens only its spill and `generate-psps` one sample at a time, about 32
descriptors whatever the cohort, so neither checks.

The limit is read with `getrlimit(RLIMIT_NOFILE)` through `rustix`, already a dependency, with no
`unsafe` (`src/open_files.rs`). An unlimited soft limit has nothing to check against. pvc never
raises the limit itself. Where the need is above the soft limit the run stops, naming the need, its
parts, the soft limit, and either:

- **the command**: `ulimit -n` a quarter above the need, rounded up to a power of two, at least
  1,024 and at most the hard limit; or
- **where the hard limit is below the need too**: the administrator (`nofile` in
  `/etc/security/limits.conf`, or `LimitNOFILE=` in a systemd unit), and the most samples that fit
  under the hard limit at this cohort's average cost a sample.

It says the thread count does not change the need, so it does not suggest fewer threads.

The 63-CRAM cohort under `ulimit -n 952` (real run):

```text
error: the run stopped: this run needs 953 files open at once and this process may open 952 (its
soft limit, `ulimit -n`): 0 BAMs at 4 each and 63 CRAMs at 11 each, 57 samples at 4 more each for
the reference bases their walks read, and 32 for the log, the output, the scratch files and slack.
Every sample's files stay open for the whole run, so the thread count does not change this. Raise
the soft limit in the shell that starts pvc, then run it again: `ulimit -n 2048` (the hard limit,
unlimited, allows it). In a systemd unit, set LimitNOFILE=2048.
```

The same with the hard limit at 900 as well (real run; the remedy sentence differs):

```text
... The hard limit, 900, is below the need too, so the soft limit cannot be raised far enough
without the administrator: ask them to raise the hard limit to at least 953 (`nofile` in
/etc/security/limits.conf, or LimitNOFILE= in a systemd unit), or call at most 51 samples at once.
```

29 psps under `ulimit -n 40` (real run):

```text
error: this run reads 29 psp files and needs 61 files open at once, one a psp and 32 besides, and
this process may open 40 (its soft limit, `ulimit -n`). Every psp stays open for the whole run, so
the thread count does not change this. Raise the soft limit in the shell that starts pvc, then run
it again: `ulimit -n 1024` (the hard limit, unlimited, allows it). In a systemd unit, set
LimitNOFILE=1024.
```

## 3. Files changed

- `src/neighbour_evidence/peak.rs`, `stage.rs`, `src/cli/calling_run.rs`, the three commands,
  `src/vcf/writer.rs`, `src/run/paralog_filter/resume.rs`, the two refinement examples: §1.
- `src/open_files.rs` (new), `src/run/callers.rs`, `src/run/mod.rs`, `src/cli/call_from_psps.rs`,
  `src/cli/estimate_parameters.rs`: §2.
- `doc/user/neighbour_evidence.md`, `calling.md`, `inputs.md`; the specs
  `neighbour_evidence_pass.md` §6, `command_line.md` §7.5, `run_streaming.md` §7.1a,
  `alignment_cursor.md`.

## 4. Tests

Unit tests: the warning before calling, whole, and at a window; no limit; the budget and threads
that fit at 85%; the before-calling judgement at 70 and 100 GB, carrying the limit and the scratch
directory; a stage under a 1-MB limit writes the same bytes as with no limit and warns fewer times
than it has windows. For open files: the shortfall, the suggested `ulimit -n` (rick's case gives
4,096), a hard limit between the need and that value, a hard limit below the need with the samples
that fit, the psp refusal's whole message, and the existing direct-mode tests with the new counts.
`cargo test --tests` (debug) passed; clippy reports nothing in the touched files.

## 5. For the owner

- **Fewer threads cannot fix an open-file shortfall**, contrary to the brief: the need does not
  depend on threads. Where the hard limit is too low the message names the administrator or a
  smaller cohort instead.
- **psp mode needs a fifteenth of direct mode's descriptors** for a CRAM cohort (1 against 15 a
  sample). The refusal for direct mode does not mention it; `doc/user/inputs.md` does.
- **The three generators each open their own cursor** on every file, which is where most of the 15
  come from. Sharing one reader between them, or opening the prefetcher's FASTA reader lazily,
  would lower the count; not attempted here.
