# Performance Review: repeat-scanner
**Date:** 2026-10-30
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** `pvc repeat-catalog` — the tandem-repeat scan, row building and catalog write — and the catalog read path, on branch `str-adjacent-runs` at c237378b
**Verdict:** Apply the listed wins (applied: H1–H6, L1–L3, L6; open: L5)
**Hot-path evidence:** sampling profile (`sample`) of the whole tomato build; per-stage instruction counts; per-change instruction, peak-RSS and wall measurements

---

## 1. Scope and constraints

- **Reviewed:** the catalog build and its read path, at c237378b (branch `str-adjacent-runs`, which
  had just made the scanner split abutting repeats and emit a joined tract with its pieces).
- **Targets:** wall time and peak memory of the build, which runs once per reference (tomato 782 Mb,
  GRCh38 3.1 Gb), and wall time of the read path, which runs once per command — so per sample in a
  cohort. The owner's rule: a large wall-time gain for a little memory is a win; a change that
  makes the code clearer is accepted even for a small gain; a change that complicates it for a
  small gain is rejected. The catalog must stay byte-identical.
- **Hardware:** Apple M-series laptop (6 fast cores), macOS, native release build
  (fat LTO, one codegen unit, `target-cpu=apple-m1`, mimalloc).
- **Hot-path evidence:** a `sample` profile of the one-thread tomato build, the methodology agent's
  per-stage instruction counts, and before/after counts for every change below.
- **In scope:** `src/tandem_repeat.rs`, `src/repeat_catalog/{builder,row,parquet_file,segments}.rs`,
  `src/reference_info.rs`, `src/region_typing/segment_criteria.rs`.
- **Out of scope:** the caller proper (calling, psp, parameter estimation); vendored crates.
- **Categories dispatched**, one agent each in its own worktree:
  - methodology — always;
  - hot_loops — the per-base scan and per-interval row building;
  - allocations, data_layout — peak memory was about 32 bytes a reference base;
  - concurrency — the per-contig threads;
  - io_and_syscalls — the FASTA pass with its two MD5s, and the parquet writer.

  float_portability was not dispatched: no float reaches the catalog except a purity
  comparison these changes leave alone.

Agent reports, diffs and the shared brief are in `tmp/perf_review_2026-10-30_repeat-scanner/`.
Timings taken by different agents ran concurrently and are not comparable across agents; every
figure in sections 3 and 5 was re-measured by the orchestrator on one quiet machine, step by step.

## 2. Verdict

**Apply the listed wins.** Eight changes were applied in sequence, each measured against the step
before it, and each left the catalog byte-identical. Together, on the whole tomato genome at one
thread:

| | c237378b | after this review |
|---|---:|---:|
| wall | 23.9 s | 13.9 s |
| instructions | 335 G | 162.5 G |
| peak RSS | 3,836 MiB | 404 MiB |

On GRCh38: one thread went from 90.9 s and 7.67 GiB to 52.2 s and 0.97 GiB. Six threads went from
33.5 s and 24.7 GiB to 18.5 s and 4.36 GiB. The GRCh38 baseline is main, not c237378b, so it also
includes the branch's own scanner change. The new build's catalog is identical at one and six
threads.

The owner then accepted H6 and L6, which were applied in follow-up commits. L5 is still open.
With all of it, on the whole tomato genome:

| threads | wall | peak RSS |
|---:|---:|---:|
| 1 | 11.3–11.7 s | 540 MiB |
| 2 | 6.0 s | 922 MiB |
| 4 | 3.7 s | 1,256–1,344 MiB |

GRCh38 at one thread: 41.9 s and 1.36 GiB. At six threads: 11.4 s and 4.68 GiB.

## 3. Measurement plan

The gate every change passed, run by `tmp/perf_review_2026-10-30_repeat-scanner/measure.sh`:

1. `pvc repeat-catalog --threads 1` on tomato chromosome 1 (90.9 Mb) and on its first 5 Mb, under
   `/usr/bin/time -l`.
2. **Instructions retired** decide speed. They repeat to within 4 in 1,000 run to run.
3. **Maximum resident set** on chromosome 1 decides memory. The 5 Mb file never reaches the real
   peak.
4. **Wall time** is a secondary check, interleaved A/B where it disagreed with instructions.
5. `cmp` against the c237378b catalog of the same file.

Owed: the peak-RSS figures are macOS ones. On macOS, maximum RSS also counts buffers a growing
`Vec` outgrew and freed, which is why several memory savings here are larger than the live heap
explains. Re-measure the build on Linux (`/usr/bin/time -v`, dev container) before quoting them as
general; `repeat_catalog.md` §2.3 says so where it quotes them.

Step by step, on chromosome 1:

| step | change | instructions | peak RSS | wall | catalog |
|---|---|---:|---:|---:|---|
| 0 | c237378b | 40.32 G | 2,943 MiB | 2.86 s | — |
| 1 | H3 stack finalised online | 40.99 G | 2,613 MiB | 2.82 s | identical |
| 2 | H2 copy floor before the copy | 34.70 G | 2,613 MiB | 2.60 s | identical |
| 3 | H1 scan returns only what can matter | 28.55 G | 489 MiB | 1.90 s | identical |
| 4 | H4 short segments skip the split | 22.51 G | 489 MiB | 1.69 s | identical |
| 5 | L1 duplicate pass in place | 22.42 G | 421 MiB | 1.71 s | identical |
| 6 | L2 rows not pre-sized | 22.41 G | 387 MiB | 1.72 s | identical |
| 7 | H5 FASTA bases taken in runs | 19.09 G | 387 MiB | 1.64 s | identical |

L3 (bases buffer reuse) acts across contigs, so it was measured on the whole genome: peak RSS went
from 494 to 404 MiB, instructions from 162.7 to 162.5 G, and the catalog stayed identical.

## 4. Build / toolchain configuration

No change. The release profile already has fat LTO, `codegen-units = 1`, `panic = "abort"` and
`target-cpu=apple-m1`. The methodology agent measured mimalloc against the system allocator on
chromosome 1. Without mimalloc, peak RSS went from 3.09 to 6.49 GB and instructions from 40.6 to
about 52 G, so the allocator choice is right. There is no `target-cpu` for Linux aarch64; this was
not measured.

## 5. Code-level findings

### Hot-path

**H1: [tandem_repeat.rs](../../../../src/tandem_repeat.rs) — the scan returned 23.6 M tracts a 91 Mb chromosome when 774,585 become rows** — *applied*
- **Confidence:** High (allocations agent; re-measured, step 3).
- **Mechanism:** the builder drops tracts below the catalog's copy floor. It cannot drop them
  before the period-multiple pass, because a tract below the floor can still remove one above it.
  For example, a 9-base period-2 tract covers half of a 16-base period-4 tract. Filtering at the
  plain floor was tried, and it writes 774,631 rows instead of 774,585.
- **What is safe** is a derived floor, worked out from the longest period down. A tract is needed
  if it clears its own floor, or if it is at least half as long as the shortest needed tract at a
  multiple of its period. That gives 4, 8, 9, 16, 20 and 18 bases for periods 1–6 at the catalog's
  floors.
- **Why the derived floor keeps every needed eliminator, chains included:** a tract that removes a
  needed tract X is at least half of X, and X is at least its period's needed length. So the
  remover is at least half that length, which is no shorter than the remover's own needed length.
- `find_tandem_repeats_for_copy_floor` takes the consumer's floor. `find_tandem_repeats` calls it
  with no floor, so every other caller is unchanged.
- A 300-sequence property test against five floor tables pins the result. It fails if the derived
  floor is replaced by the plain one.
- **Measured:** peak RSS 2,613 → 489 MiB, instructions 34.70 → 28.55 G.
- **Complexity:** one function, one derived table, one test.
- **Side effect for the owner:** the log's `below_copy_floor` now counts only the sub-floor tracts
  the scan returns. On chromosome 1 that is 1,492,417, where it was 22,844,074. The catalog file
  does not carry this number. `BuildTally::below_copy_floor` says so.

**H2: [row.rs](../../../../src/repeat_catalog/row.rs) — `row_for_interval` copied every tract before testing the copy floor** — *applied*
- **Confidence:** High (hot_loops, allocations and data_layout agents all found it; step 2).
- **Mechanism:** an uppercased copy of the tract and of the motif was made for each of the 23.6 M
  intervals, and 22.8 M of them then failed the floor. The checks now run first, in their old
  order, so the tally is unchanged; the motif is borrowed instead of copied.
- **Measured:** instructions 40.99 → 34.70 G. Peak RSS unchanged; the log's tally lines unchanged.
- **Complexity:** none added.

**H3: [tandem_repeat.rs](../../../../src/tandem_repeat.rs) `maximal_scoring_subsequences` — the Ruzzo–Tompa stack held every segment until the contig's end** — *applied*
- **Confidence:** High (three agents; step 1).
- **Mechanism:** the doc comment called the growth pathological and deferred the fix, but on real
  sequence it is the normal case: about 18.4 M segments a chromosome, 1.34 GB at period 1. At a
  Rule-2 push every stacked segment has an `l` at least the new one's, so no later search can
  reach them. They are final, and are emitted there, in the same order. The brute-force property
  test still passes.
- **Measured:** peak RSS 2,943 → 2,613 MiB, instructions +1.7%.
- **Complexity:** three lines; removes a deferral note.

**H4: [tandem_repeat.rs](../../../../src/tandem_repeat.rs) — every segment went through the motif-change split, and the split's callback was no longer inlined** — *applied*
- **Confidence:** High (hot_loops agent; step 4).
- **Mechanism:** 84 M of the 92 M segments at periods 2–6 are shorter than anything `emit` keeps,
  and every piece of a split lies inside its segment. So those segments now return before the
  split, and segments with no mismatch are emitted without it.
- The split is `#[inline(never)]`. Inlined, it made the scan's per-segment closure too big to
  inline into the scan loop. Four interleaved runs each: out of line 1.69 s and 22.5 G
  instructions, inlined 1.77 s and 24.6 G.
- **Measured:** instructions 28.55 → 22.51 G, wall 1.90 → 1.69 s.
- **Complexity:** none net; the no-mismatch test moved from the split into its caller.

**H5: [reference_info.rs](../../../../src/reference_info.rs) — the FASTA pass handled every base byte one call at a time** — *applied*
- **Confidence:** High (io agent; step 7).
- **Mechanism:** `push_window` takes each run of base bytes inside a sequence line in one go.
  Every other byte still goes through `push_byte`, where all the checks and errors live.
- I checked the run path against `push_byte`'s sequence arm:
  - a run never contains a line end;
  - a `>` at a line start still takes the byte path;
  - the MD5 batch flushes at the same points.
- This pass also runs on every catalog read, as the reference check.
- **Measured:** instructions 22.42 → 19.09 G, wall 1.72 → 1.64 s; the 59 `reference_info` tests pass.
- **Complexity:** about 45 lines, plus the invariant that `push_bases` matches `push_byte`.

**H6: [builder.rs](../../../../src/repeat_catalog/builder.rs) `flush_pending` — reading and scanning never overlap** — *applied after the owner's decision (2026-10-30)*
- **Confidence:** High (concurrency agent, with timers around the hand-offs).
- **Mechanism:** contigs are scanned in batches of `--threads`; the batch waits for its slowest
  contig and its write, and only then does reading resume. At two threads on c237378b, the
  scanners sat idle 6.1 s of 17.8 s. Handing each contig to the next free scanner and writing in
  reference order cut wall time from 17.2–18.4 s to 11.8–13.0 s, with the catalog identical.
- **After this review it is still worth it:** whole tomato takes 13.7 s at one thread, 9.35 s at
  two and 6.65 s at four, against perfect overlap at 6.9 s and 3.4 s.
- **Complexity:** about 110 lines — long-lived threads, a hand-off channel, a reorder buffer and
  panic propagation. It would also have to hand the bases buffer back to keep L3.
- **The decision:** with overlap, `--threads 1` uses two cores, a reader and a scanner. The
  owner accepted that.
- **Applied:** rewritten on `crossbeam-channel`, so the scanners share a receiver without a mutex.
  - Each scanner returns the contig's bases buffer with its rows, so L3's reuse survives.
  - Results are collected when a contig starts as well as when it ends, so a returned buffer
    serves the very next contig.
- **Measured on whole tomato, against the build before it, two rounds each:**

  | threads | wall | peak RSS |
  |---:|---|---|
  | 1 | 13.8 → 11.3–11.7 s | 403 → 540 MiB |
  | 2 | 9.1 → 6.0 s | 677–728 → 922 MiB |
  | 4 | 6.3 → 3.7 s | 1,209–1,225 → 1,256–1,344 MiB |

  The catalog is byte-identical at every count.
- **GRCh38:** one thread 52.2 → 41.9 s, six threads 18.5 → 11.4 s.

### Likely

**L1: [tandem_repeat.rs](../../../../src/tandem_repeat.rs) period-multiple pass — the tracts were held three times** — *applied*
- **Confidence:** High (data_layout agent; step 5).
- **The change:**
  - survivors are compacted in place with `retain` instead of copied into a second vector;
  - the search lists hold 16-byte spans, only for periods that have a multiple in range;
  - each list is sized once from the sorted output.
- **Measured:** peak RSS 489 → 421 MiB, instructions flat.
- **Ablation:** without the sizing, peak RSS is 452 MiB and wall 1.90 s, so the sizing stays.
- **Complexity:** the pass is shorter than before.

**L2: [builder.rs](../../../../src/repeat_catalog/builder.rs) `scan_contig` — rows reserved for every interval** — *applied*
- **Measured:** peak RSS 421 → 387 MiB (step 6). One line.

**L3: [builder.rs](../../../../src/repeat_catalog/builder.rs) `flush_pending` — every contig's bases grew from an empty buffer** — *applied*
- **Mechanism:** a 91 Mb chromosome passes through buffers of every power of two below 128 MiB.
  The largest finished buffer is now kept for the next contig.
- **Measured:** whole genome, peak RSS 494 → 404 MiB; catalog identical.
- **Complexity:** six lines. The spec's "nothing accumulates across contigs" now names the one
  buffer that does.

**L4: [tandem_repeat.rs](../../../../src/tandem_repeat.rs) period-multiple pass — a forward-only cursor instead of a binary search per tract** — *experiment shows no gain — closing*
- Measured by the hot_loops agent at c237378b, when the period-1 list held 18.4 M tracts.
- After H1 it holds about a tenth of that. Re-measured on top of step 5, three interleaved runs:
  wall 1.71–1.72 s for both, instructions 22.42 → 22.54 G. More code for nothing; rejected.

**L5: [reference_info.rs](../../../../src/reference_info.rs) `flush_md5` — the two MD5s run as two passes over every base** — *not applied: owner decision*
- **Measured (io agent probe):** one interleaved two-lane MD5 takes 0.12 s against 0.22 s for two
  `md-5` passes on chromosome 1. Both digests are checked against `md-5`.
- **Why it matters:** this pass runs on every command that reads the catalog — once per sample.
- **Complexity:** a hand-written MD5 of about 60 lines, plus per-lane buffering. It would delete
  the two-thread `split_digests` path.
- Diff: `io_md5-two-lanes-probe.diff`.

**L6: [parquet_file.rs](../../../../src/repeat_catalog/parquet_file.rs) — dictionary encoding on the four coordinate columns** — *applied after the owner's decision (2026-10-30), commit 587a4860*
- **Measured (io agent), delta encoding instead, on chromosome 1:**
  - write 0.09 → 0.06 s;
  - file 6.73 → 4.25 MB (−37%);
  - whole-catalog read 0.035 → 0.024 s.
- It changes the catalog's bytes, so it breaks byte-identity with existing catalogs by design.
  The owner accepted that, on condition that the caller's output did not change. Catalogs are
  rebuilt.
- **Re-measured on whole tomato:**
  - file 60.7 → 35.7 MB;
  - build 162.5 → 155.8 G instructions, about 2% fewer cycles at the same clock;
  - `dev type-regions` 0.55 → 0.47 s;
  - typed regions identical.
- Old catalogs still read, because parquet records each column's encoding.
- Spec §3 already said the coordinates delta-encode; the writer had not asked for it.

### Speculative

**S1: [builder.rs](../../../../src/repeat_catalog/builder.rs) — scan one contig's periods in parallel** — *experiment shows no gain — closing*
- Measured by the concurrency agent:
  - at two threads, 18.51 and 19.12 s against 18.35 and 18.17 s for the base;
  - peak RSS higher, because the per-period lists coexist.

### Note

- **Narrowing `RepeatInterval` to `u32` coordinates:** rejected. It would break spec §4's u64
  coordinate rule on contigs over 4.29 Gb.
- **Struct-of-arrays:** does not fit; the loops read start, end and period together.
- **The `canonical` copy (91 MB):** stays. Doing it in place would change rows that contain an `N`.
- **Reserving each contig's bases from the `.fai`:** cut 469 MB in a hard-coded experiment at
  c237378b. L3 now takes most of what remains, and the build is not given a `.fai` today. Not
  pursued.
- **`Vec<RowRejection>`:** collects the rejections only to count them. After H1 it holds about 1.5 M
  bytes a chromosome. Not worth changing.

## 6. Out-of-scope observations

- **Every command that reads the catalog re-hashes the whole reference**, to prove it still matches
  its `.fai`. On chromosome 1 that is 0.19 s of a 0.27 s `pvc dev type-regions`, against 0.03 s to
  open the catalog and read all its rows. `generate-psps` pays it once per sample. H5 and L5 make
  the pass cheaper; caching the proof per reference is a design question for a separate issue.
- **`profiling_environment.md` said instructions cannot be counted on this machine.**
  `/usr/bin/time -l` on the host counts them for the whole process. Corrected in this change.

## 7. What's already good

- **The scan reads one canonicalised copy of the contig:** `find_tandem_repeats` builds it once and
  shares it across all six periods, rather than canonicalising per period
  ([tandem_repeat.rs](../../../../src/tandem_repeat.rs)).
- **The FASTA pass batches both digests per 64 KiB window**, rather than per byte, and the parquet
  writer buffers internally. The I/O agent found no syscall-level waste
  ([reference_info.rs](../../../../src/reference_info.rs)).
- **Byte-identity across thread counts is designed in:** rows travel with their contig's order, so
  completion order cannot reach the file. That is what made every change here checkable with a
  plain `cmp` ([builder.rs](../../../../src/repeat_catalog/builder.rs)).
