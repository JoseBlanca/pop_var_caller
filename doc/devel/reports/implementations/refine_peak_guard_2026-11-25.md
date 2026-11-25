# Refinement's expected peak memory, checked against what the run may use (2026-11-25)

Refinement (the neighbour-evidence pass, off by default) sizes its windows to a byte budget,
`--refine-window-bytes`. The process holds more than that budget counts: on 2,176 tomato samples at
48 threads it peaked at 103.0 GB under an 80-GB budget
(`refine_memory_overhead_2026-11-25.md`). A run whose budget looked safe could therefore be killed
by the system days into a job. `pvc` now estimates refinement's peak, compares it with the memory
the run may use, and stops with an explanation when the estimate is too close to that limit.

Branch `peak-guard`, from `mem-overhead` (af6ba479). Nothing ran on kimura or rick.

## 1. What it does

1. **Before calling starts**, every command that can run refinement (`call-from-alignments`,
   `call-from-psps`, `resume-calling`) estimates refinement's peak from the budget, the sample
   count, the thread count and the longest contig, and compares it with the memory limit.
2. **At every window**, once the window is planned and before it is solved, it estimates the peak
   again from the window's own count and the records its walks read, both then known exactly.
3. If an estimate is above **90% of the limit**, the run stops with an error naming the estimate,
   its parts, the limit and where the limit came from, and what to change. Given
   `--refine-ignore-memory-limit`, it logs the same as a warning and goes on.
4. Otherwise it logs the estimate: once before calling, and for the first window. Each later window
   whose estimate is higher than every one before it goes to the debug log.

Where refinement does not run (the default, one sample, ploidy other than two), nothing is
estimated or logged. The check changes no window and no output: the thread count enters only the
estimate, never the plan (pass spec §2.3).

## 2. The estimate

```text
peak = (counted + 1.05 kB × walking threads × walked records) × 1.16
```

- **counted**: what the window budget counts for the window.
- **walking threads**: the run's threads, or its samples if fewer, since a thread walks one sample
  at a time.
- **walked records**: the records the window's walks read.
- **× 1.16**: the allocator's slack.

| Constant | Value | Where it was measured |
|---|---|---|
| Walks' arrays, a walked record a thread | 1,050 bytes | 1.04 to 1.07 kB, drawn cohort of 300 samples, tomato25's mix of 2.8 alleles a record, 1 to 6 threads, 100- to 400-kb cores (overhead report §3) |
| Allocator slack | 16% | tomato25's highest window, 103.0 GB measured against about 89 GB held (glibc, 48 threads); mimalloc at its defaults held 12 to 25% on the drawn cohort (Linux, 1 to 6 threads) |
| Budget count, a sample a walked record | 170 bytes | 169 on the drawn cohort (2,314 MB, 45,600 walked records, 300 samples); at most 216 on tomato25's 80-GB window |
| Budget count, a sample a record held | 250 bytes | about 240 measured on 400 tomato samples (pass spec §6) |
| A record's own cost, whatever the cohort | 1,000 bytes | **not measured**: an allowance so that a two-sample estimate is not near zero |
| Densest records a kb | 150 | the drawn cohort's density; tomato25's densest stretch holds about 90 |

**Before calling**, nothing about the windows is known yet, so both terms are bounded:

- counted = the budget, or what a window of the longest contig could hold at 150 records a kb,
  whichever is less;
- walked records = the budget over 170 bytes a sample, or what such a window could walk (the core
  and 50 kb either side at 150 a kb, the far margins at `--refine-margin-density`), whichever is
  less.

At a large cohort the budget sets both, and the estimate follows the measurements (§4). At a small
cohort the window's length sets both, and the estimate is an upper bound that stays small (§5).
170 bytes is the lower of the two figures measured, so the before-calling estimate errs high.

**At a window**, counted is the window plan's own count and walked records is the number of records
its walks read (with thinned margins, the walked ones; without, every site). For the budget to
suggest, both are scaled in proportion to the budget.

## 3. The limit, and the margin

**The limit** is the smaller of:

- **the machine's memory**: `MemTotal` from `/proc/meminfo` on Linux; on macOS `hw.memsize`, read
  by running `/usr/sbin/sysctl -n hw.memsize`, because asking the kernel directly needs `unsafe`,
  which the crate forbids (`[lints.rust] unsafe_code = "forbid"`);
- **the control group's limit** (Linux): cgroup v2's `memory.max`, or v1's `memory.limit_in_bytes`,
  for the group `/proc/self/cgroup` names and every group above it, since a parent's limit binds
  its children. A systemd unit's `MemoryMax=` is such a limit. `max`, and v1's value for no limit
  (9223372036854771712), mean no limit.

The code is `src/memory_limit.rs`; it only reads files and runs `sysctl`. When nothing can be read
(another operating system, no `/proc`), the run logs "not checked: the memory this run may use
could not be read", with the estimate, and goes on. Not read: `RLIMIT_AS` (`ulimit -v`), which
limits address space rather than resident memory and which mimalloc's reservations would exceed
long before the resident peak did.

**The margin: stop above 90% of the limit.** The estimate is the middle of what was measured, not
an upper bound. The spread around it:

- the allocator's slack was 12 to 25% under mimalloc's defaults, around the 16% used: the peak up
  to about 8% above the estimate;
- the rest of what the program holds was within 5% of the count, either way (overhead report §3):
  up to about 4% of the peak.

Stopping at 90% keeps about 11% of headroom (1 / 0.9) above the estimate for that spread. A lower
share would refuse budgets that were measured to fit: at 85%, 80 GB under 118 GB would stop,
although the replay at 80 GB peaked at 103.0 GB, 87% of 118.

**Other stages.** With the hidden-duplication filter on (the default), refinement runs in the
filter's writing pass, after calling has ended and let its memory go, so it is the run's peak at
any budget near the machine's size. With `--hidden-dup-fdr 0` refinement runs while calling still
runs, and calling's own memory comes on top of the estimate; the estimate does not include it.
Calling's own peak at 2,176 samples was not measured here (at 100 tomato samples in psp mode it was
10 GB, `doc/user/performance.md`).

## 4. What it says on tomato25

2,176 samples, 48 threads, default windows (2-Mb core, 1-Mb margins, margins thinned to 5 a kb),
chromosome 1 about 90 Mb. Limit 118 GB (118 × 10⁹ bytes):

| Budget | Counted | Walked records | Arrays | Expected peak | Share of 118 GB | Measured peak (glibc) | Outcome |
|---|---|---|---|---|---|---|---|
| 65 GB | 65.0 GB | 175,713 | 8.9 GB | 85.7 GB | 73% | 88.7 GB | runs |
| 70 GB | 70.0 GB | 189,230 | 9.5 GB | 92.3 GB | 78% | — | runs |
| 80 GB | 80.0 GB | 216,262 | 10.9 GB | 105.4 GB | 89% | 103.0 GB | runs, 0.8 GB below the stop |
| 100 GB | 100.0 GB | 270,328 | 13.6 GB | 131.8 GB | 112% | — | stops; suggests 80 GB |

The estimates lie within 4% of both measured peaks: 2% above at 80 GB, 3% below at 65 GB.

**If the unit's limit is `MemoryMax=118G`**, systemd reads that as 118 GiB, 126.7 × 10⁹ bytes, and
the shares are 68%, 73% and 83%; the stop is then at 114.0 GB and the largest budget that fits is
86 GB.

The line the 70-GB run logs before calling:

```text
refinement memory, before calling: about 92.3 GB at its peak, 78% of the 118.0 GB this run may use
(the memory.max of the control group pvc runs in); 70.0 GB for a window as --refine-window-bytes
counts it (70.0 GB allowed), 9.5 GB for 48 threads' working arrays (189,230 walked records at
1.05 kB each a thread), and 16% on top for the memory allocator
```

The error a 100-GB run stops with, before calling:

```text
genotype refinement is expected to need about 131.8 GB, 112% of the 118.0 GB this run may use (the
memory.max of the control group pvc runs in), and pvc stops above 90% so that the system does not
kill the run part of the way through. The estimate is 100.0 GB for a window as --refine-window-bytes
counts it (100.0 GB allowed), 13.6 GB for 48 threads' working arrays (270,328 walked records at
1.05 kB each a thread), and 16% on top for the memory allocator. To go on, lower
--refine-window-bytes to 80000000000 (80.0 GB) or less, or give --refine-ignore-memory-limit to
run anyway, if this run's peak is known to be lower.
```

At a window, the error names the window ("in its window from contig 0 position 43350003"), and,
where fewer threads would fit at the same budget, adds "run on N threads or fewer (--threads)". It
ends: "A run started with --keep-tmp-on-error and the hidden-duplication filter on keeps its calls,
and `pvc resume-calling` finishes it with the new options."

## 5. Small cohorts

- **One sample**: refinement never runs (`decide`: "one sample: no other sample to copy from"), so
  nothing is estimated or checked.
- **Two samples** with `--refine-genotypes on`, default 32-GB budget, any thread count: the
  window's length sets it. Counted at most 0.9 GB (600,000 records at 1.5 kB), 324,500 walked
  records on 2 threads, 0.7 GB of arrays: **1.8 GB expected**.
- **Ten samples on 48 threads**: 2.1 GB counted, 10 walking threads, 3.4 GB of arrays: **6.4 GB**.

These are upper bounds, set by 150 records a kb; a sparser genome holds less. A short genome is
bounded by its longest contig: a 300-kb contig gives a window of 45,000 records at most. The
largest the walks' term can be before calling, at any cohort, is about 0.34 GB a walking thread
(324,500 records × 1.05 kB), so only a machine with less than about 0.4 GB a thread could be
stopped by it on a small cohort.

## 6. The option

**`--refine-ignore-memory-limit`**, under *Advanced*, on `call-from-alignments`, `call-from-psps`
and `resume-calling`. Named after the `--refine-*` family the other refinement options belong to.
`resume-calling` does not inherit it from the stopped run: the resume point does not store it, and
the resumed run reads its own limit and must be given the option again.

## 7. What changed

- `src/memory_limit.rs` (new): the limit, as §3.
- `src/neighbour_evidence/peak.rs` (new): the constants, `PeakEstimate` (before calling, at a
  window), `MemoryGuard`, `check`, the error `PeakTooClose` and its message, and the suggestions.
- `src/neighbour_evidence/stage.rs`: `NeighbourEvidenceDecision::memory_guard`; each window checked
  before it is solved; `StageStop::OverMemory`. The vector-handing `push` and `finish` now return
  `StageStop<Infallible>` (`VectorStageStop`), so they can report the new stop rather than panic.
- `src/vcf/writer.rs`: `VcfWriteError::RefinementOverMemory`.
- `src/cli/calling_run.rs`: `check_refinement_memory`, called in all three commands right after
  the refinement decision; `CallingRunError::RefinementOverMemory`.
- `src/run/paralog_filter/resume.rs`: a resumed decision carries no guard until the command sets it.
- Help of `--refine-window-bytes` (both calling commands), the pass spec §6, and two short additions
  to `doc/user/` (`neighbour_evidence.md`, `calling.md`).
- `examples/refine_replay.rs`, `examples/refine_mem_probe.rs`: the new stop's match arm.

## 8. Tests

Unit tests: the limit parsed from sample `/proc/self/cgroup`, `memory.max`, `memory.limit_in_bytes`
and `/proc/meminfo` text (v2 unit, parent tighter, `max`, container root, v1 and its no-limit value);
the estimate at tomato25's budgets and against the two measured peaks; the small-cohort, short-contig
and unthinned bounds; the stop at 90%, the forced run, no limit; the suggested budget and threads;
the whole message. In the stage, a window over a 1-MB limit stops before it is solved, naming
the window, and the same run with `--refine-ignore-memory-limit` writes the same bytes as a run
with no limit. In `calling_run`, the before-calling judgement at 70 and 100 GB, forced and with no
limit.

## 9. For the owner to check

- **The slack constant is 16%, mimalloc's middle at its defaults.** When the change that sets
  mimalloc's purge delay to 0 lands, the slack measured on the drawn cohort falls to 6 to 7%; it has
  not been measured at 48 threads. The constant should then be lowered once the kimura run of the
  overhead report §7 confirms it.
- **Every window is checked, not only the first.** A later window that would go over stops the run
  then, after calling. With the filter on, that discards the spill unless `--keep-tmp-on-error` was
  given — the same as a window over budget does today. An OOM kill, by contrast, leaves the scratch
  directory behind. Keeping the scratch directory on this stop, whatever the flag, would make it
  always resumable.
- **`doc/user/neighbour_evidence.md` still says** "the process as a whole has used up to about a
  fifth more than the budget"; the measurements say about 1.3 times at 48 threads. I added a bullet
  and left the owner's sentence as it was.
- **The 1,000-byte cost of a record and the ceiling of 150 records a kb** are not measured bounds;
  they matter only where the window's length, not the budget, sets the estimate.
