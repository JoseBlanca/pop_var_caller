# Refinement memory the window budget does not count (2026-11-25)

The window budget (`--refine-window-bytes`) counts what refinement holds for every sample at every
record of a window. On the 2,176-sample tomato25 replay at 48 threads, the highest window's real
memory peaked at 103.0 GB under an 80-GB budget and 88.7 GB under 65 GB
(`refine_replay_tomato25_2026-11-25.md`). A straight line through those two points reads as a
fixed 27 GB on top of the budget. This report measures what that extra memory is made of.

Branch `mem-overhead`, from main cd65e4e2. Every measurement here was made on the laptop (macOS,
arm64, six threads at most), on a drawn cohort, natively and in the Linux development container.
Nothing ran on kimura; §7 gives the commands for one confirming run there.

## 1. What the extra memory is

The extra memory has two parts. Neither is fixed: both grow with the window, so with the budget.

1. **The walks' working arrays, one set a thread.** Each thread solving a sample holds about
   1.05 kB for every record the window walks: the templates' allele shares, the sample's
   likelihoods, the posteriors and the drawn paths. The budget leaves them out on purpose, so that
   the windows, and with them the output, do not depend on the thread count (pass spec §2.3, §6).
   Measured: 33 to 78 MB a thread on the drawn windows (100- to 400-kb cores), a straight line in
   the thread count.
   *Projected* for tomato25's highest window, about 170,000 walked records at 48 threads:
   **about 9 GB**.
2. **The allocator's slack**: memory the process holds resident beyond what the program has asked
   for and not freed. Measured on Linux at 2 to 6 threads: 8 to 21% of the bytes asked for under
   glibc's malloc, and 12 to 25% under mimalloc, `pvc`'s allocator, at its default settings.
   Under mimalloc, most of it is memory freed in the last second and not yet returned. Telling
   mimalloc to return freed memory at once (`MIMALLOC_PURGE_DELAY=0`) brings it down to 6–7%.
   *Inferred* for the tomato25 window, as the 103.0 GB less the rest: **about 14 GB**, 16% of
   what was asked for.

What the budget counts agrees with what the program holds, apart from the walks' arrays: within
5% of the count at every window size tried, above it or below.

So the overhead is a share of the budget, not a constant: about 29% at 48 threads on tomato25.
The two points fit `1.16 × (budget + 0.11 × budget) ≈ 1.29 × budget`, which gives 103 GB at an
80-GB budget, as measured, and 84 GB at 65 GB, against 88.7 measured.

**The replay ran on glibc's malloc**, `pvc` runs on mimalloc (cd65e4e2 records this), and at six
threads mimalloc's slack was larger than glibc's. So a `pvc` run at the same budget would peak at
least as high as the replay did, until mimalloc is told to return freed memory at once.

## 2. How it was measured

Three quantities, per window, each from the end of the window before to the end of this one:

- **Counted**: what the window budget counted for the window's plan
  (`NeighbourEvidenceStage::last_window`).
- **Requested peak**: the most bytes the program held at once, asked for and not yet freed. A
  wrapper around the global allocator counts every allocation and free
  (`tools/counting_allocator`). The figure does not depend on which allocator runs underneath:
  every arm below gave the same requested peaks to the megabyte.
- **Resident peak**: what the kernel says the process holds in memory, read every 20 ms; on Linux
  also the kernel's own peak (`VmHWM`), which gave the same values. It is anonymous memory
  throughout: mapped files were under 4 MB.

**Requested minus counted** is live data the budget does not count. **Resident minus requested**
is the allocator's slack.

**The cohort** is `examples/refine_mem_probe.rs`'s: 300 drawn samples over 1.2 Mb at 150 records
a kb, tomato25's mix of 2 to 6 alleles, 3 fragments covering a base, 90 in 100 samples inbred.
Windows have fixed 200-kb cores and 250-kb margins, thinned to 5 records a kb beyond 50 kb of the
core, as `pvc`'s default. Windows 2 to 4 have the same records (100,000, of which about 45,600
are walked) and are compared below. Window 5 is the run's last two windows together and is left
out. Window 1 has no left margin.

300 samples at 6 threads keeps the per-thread part's share near tomato25's: 48 threads for 2,176
samples is 2.2 threads a hundred samples, 6 for 300 is 2.0.

**Arms**: mimalloc (`pvc`'s, version 3.3) and the system allocator (`--no-default-features`;
glibc 2.36 in the Debian 12 container; kimura's version was not checked), at 1, 2, 4 and 6
threads. mimalloc ran in the container and natively on macOS. Every arm wrote the same output (hash `c79d12db1e6d8797`).

## 3. Live data the budget does not count

Requested minus counted, MB, window 4, by threads:

| Threads | glibc | mimalloc (Linux) | mimalloc (macOS) |
|---|---|---|---|
| 1 | 59 | 59 | 59 |
| 2 | 109 | 109 | 109 |
| 4 | 207 | 206 | 207 |
| 6 | 304 | 304 | 304 |

It rises 49 MB a thread, in a straight line, and is 10 MB at none. 49 MB over 45,600 walked records
is **1.07 kB a walked record a thread**, what the memory build of 2026-11-24 computed for the walks'
arrays (1.10 kB). Nothing else of size goes uncounted at this window.

At other window sizes, mimalloc on Linux, window 4:

| Core | Walked records | Counted | Uncounted at 1 thread | at 6 threads | A thread |
|---|---|---|---|---|---|
| 100 kb | about 31,300 | 1,571 MB | 103 MB | 266 MB | 33 MB, 1.04 kB a walked record |
| 200 kb | about 45,600 | 2,314 MB | 59 MB | 304 MB | 49 MB, 1.07 kB |
| 400 kb | about 74,000 | 3,875 MB | −130 MB | 262 MB | 78 MB, 1.06 kB |

The part that does not grow with threads (the one-thread figure less one thread's arrays) is
+70 MB at 100-kb cores, +10 MB at 200 kb and −208 MB at 400 kb: within 5% of the count either way.
It is a net figure: what the budget leaves out, such as the batch of lines waiting to be read,
less what it counts that the window does not hold at its peak. This report does not break it down
further.

### What the walks' arrays hold

Per thread and walked record, from `solve_round` and `walk`, at tomato25's allele mix (on average
2.8 alleles and 6.1 genotypes a record, 16 templates and the background):

| Array | Bytes |
|---|---|
| The templates' allele shares, `f64`, 17 × 2.8 | about 380 |
| The posteriors and the neighbours' predictions, `f64` a genotype, with their offsets | about 110 |
| The likelihoods, `f64` a genotype | about 50 |
| The definite alleles drawn, a byte each of the shares | about 50 |
| Each site's spans and its `WalkSite` | about 100 |
| Drawn haplotypes, paths, the forward checkpoints and the rest | about 150 to 400 |
| **Measured** | **1,040 to 1,070** |

The shares are a third of it or more.

## 4. The allocator's slack

Resident minus requested, MB, windows 2, 3 and 4 (the requested peak was about 2,550 MB):

| Threads | glibc (Linux) | mimalloc (Linux) | mimalloc (macOS) |
|---|---|---|---|
| 1 | 14, 31, 40 | 426, 453, 483 | 392, 408, 464 |
| 2 | 432, 473, 416 | 504, 515, 548 | 688, 714, 749 |
| 4 | 326, 521, 503 | 351, 310, 461 | 741, 836, 928 |
| 6 | 195, 277, 375 | 441, 463, 627 | 692, 823, 969 |

Settings, Linux, 6 threads:

| Arm | Slack, windows 2–4 | Share of requested | Time |
|---|---|---|---|
| mimalloc, as `pvc` runs it | 441, 463, 627; repeated 475, 622, 541 | 17–24% | 278 s, repeated 279 s |
| mimalloc, `MIMALLOC_PURGE_DELAY=10` | 203, 201, 209 | 8% | 284 s |
| mimalloc, `MIMALLOC_PURGE_DELAY=0` | 152, 157, 171 | 6–7% | 289 s |
| glibc | 195, 277, 375 | 8–14% | 289 s |
| glibc, `MALLOC_ARENA_MAX=2` | 208, 283, 270 | 8–11% | 329 s |

At 400-kb cores, 6 threads, mimalloc's slack grew window after window: 613, 909 and 1,003 MB
(15 to 24%); with the purge delay at 0 it held at 278, 268 and 270 MB (6.5%), in 416 s against 405.

**Where mimalloc's slack comes from.** mimalloc 3 keeps memory it has been given back for a second
(`purge_delay`, 1,000 ms by default; it was 10 ms in mimalloc 2) before returning it to the
system. Two moments show it in the trace of a 6-thread window (`MEMORY_TRACE`):

- *During the rounds*, resident sits 350 to 450 MB above requested, 13 to 17%: each round frees
  and asks again for every thread's arrays and the next round's haplotypes.
- *When the rounds end*, the haplotype store is freed and the core's results are made at once.
  The store's memory is not yet returned when the results ask for theirs, so for about a second
  the process holds both. In the 400-kb window this is the window's resident peak: 5,045 MB against
  4,079 requested, 970 MB of slack against about 450 during the rounds. The budget's rule, that the
  store and the results are never held together, holds for the program's data but not for its
  resident memory.

glibc shows no such moment in its traces: resident memory falls as the store is freed. Its slack is
small at one thread (14 to 40 MB) and 195 to 521 MB from two threads on. Limiting it to two memory
pools (`MALLOC_ARENA_MAX=2`) left the slack the same and cost 14% in time.

**macOS is not a guide for kimura.** The same mimalloc held 1.5 to 2 times the Linux slack there,
growing every window, for reasons not looked into. The Linux figures are the ones to read.

## 5. Projected to tomato25, 2,176 samples at 48 threads

The highest window of the 80-GB replay: core 43,350,003–44,399,984, 148,946 records.

| Part | GB | How obtained |
|---|---|---|
| Counted by the budget | ≤ 80 | The plan's limit; the replay did not log the count (it does now). |
| The walks' arrays | about 9 | **Estimated**: 1.05 kB × 48 threads × about 170,000 walked records (the core's 148,946, about 14,000 in the 50 kb either side read in full, up to 9,500 in the far margins at 5 a kb). The 1.05 kB is measured on the drawn cohort of 300; the arrays are sized by the walked records and the 16 templates, not by the cohort, so it should carry over. |
| Other live data, net of what the count overstates | −4 to +4 | **Measured** on the drawn cohort at within 5% of the count, either way, depending on the window's shape; taken as 0 below. |
| The allocator's slack (glibc) | about 14 (10 to 18) | **Inferred**: 103.0 measured less the lines above. About 16% of what was asked for, inside the 8 to 21% glibc showed at 2 to 6 threads. |
| **Resident peak** | **103.0** | **Measured** (kernel peak, glibc, 2026-11-25). |

How the same window would fare otherwise, *estimated* from the shares in §4:

| Allocator | Slack | Peak |
|---|---|---|
| glibc (the replay) | about 14 GB | 103 GB, measured |
| mimalloc as `pvc` runs it | 12–25% of about 89 GB | about 100 to 111 GB |
| mimalloc, purge delay 0 | 6–7% | about 95 GB |

At other cohort sizes the walks' arrays matter more as the threads per sample rise. Their share
of the budget is about `threads × 1.05 kB / (170 B × samples)` at tomato-like windows: 14% at 48
threads and 2,176 samples; a whole budget's worth at 48 threads and 300 samples. At small cohorts
the window is bounded by its core length, not the budget, and the arrays stay small in absolute
terms: at 57 samples, 48 threads and a 2-Mb core at 26 records a kb, about 64,000 walked records,
about 3 GB.

## 6. Options

### Count it

- **C1. A stated total, not a new count** (recommended). The windows must not depend on the
  thread count (§2.3), so the budget cannot charge the walks' arrays without breaking that. What
  it can do is say, once the first window is planned, what the run will hold: the count, the walks'
  arrays at this run's threads (1.05 kB × threads × walked records, from the plan), and a margin
  for the allocator; and warn when that exceeds the machine's memory. No output change, low risk.
  It turns the budget advice below into a number the run prints.
- **C2. Charge the walks' arrays at a fixed nominal thread count**, 48 say. The windows stay
  independent of the threads, but a small run is charged for threads it does not have, and a run
  above 48 threads is still undercounted. Not recommended.
- **C3. A total memory limit that sets the window budget from the threads.** The output would then
  depend on the thread count, which the owner's rule forbids. Not recommended.

### Cut it

| Option | Expected on tomato25's window, 48 threads | Risk |
|---|---|---|
| **K1. mimalloc's purge delay at 0** in `pvc` (or 10 ms, mimalloc 2's default) | About −5 to −16 GB against mimalloc as it runs now (slack 12–25% → 6–7% of about 89 GB); about −9 GB against the replay's glibc | Low. Output unchanged. On the drawn cohort it cost 3 to 4% in time (278 s twice against 289 s; at 400-kb cores 405 against 416 s); 10 ms cost 2% (284 s) for 8% slack. Needs either `MIMALLOC_PURGE_DELAY=0` in the run's environment, or a call to `mi_option_set` at start-up from a small crate that may use `unsafe`, as the counting allocator does. |
| **K2. The templates' shares computed where they are read**, not held for the whole window | About −3 GB (the shares, about 380 B of the 1.05 kB) | Medium. Output unchanged if computed in `f64` the same way; costs time, since a walk reads each share in the forward sweep, the backward sweep and the draws. |
| **K3. The shares held in `f32`** | About −1.5 GB | Medium. Output moves: the walks' sums are no longer the same bits. |
| **K4. Fewer threads walking at once** | Proportional | High cost in time; not recommended. |

### Budget advice for the tomato25 rerun

At 48 threads the peak is about 1.29 times the budget with the replay's glibc, and should be
expected at least that with `pvc`'s mimalloc as it runs now: **70 GB gives about 90 GB, 75 GB
about 97 GB.** With K1 the factor should fall to about 1.2, so 80 GB would give about 96 GB. The
confirming run of §7 measures all three on kimura.

## 7. Confirming on kimura (prepared, not run)

One replay of the short slice (chromosome 1 at 44.0–45.5 Mb, 130,402 entries) at 48 threads,
fixed 500-kb cores and margins, in three arms: mimalloc as `pvc` runs it, mimalloc with the purge
delay at 0, and glibc. Each logs, per window, the count, the requested peak and the resident peak,
and writes a memory trace. About 15 minutes an arm (821 s for these windows on 2026-11-23). It
builds from a bundle of this branch into a clone under `pvc_work`, so nothing is pushed and
nothing outside `/media/tomato25_vcfs/pvc_work/` is written. The commands:

```bash
# NOT RUN. One confirming replay on kimura, for the owner to launch when the machine is free.
# Three arms on the short slice (chr01 44.0-45.5 Mb, 130,402 entries), same fixed windows, 48
# threads, each counting the bytes asked for: mimalloc (pvc's allocator), mimalloc with
# MIMALLOC_PURGE_DELAY=0, and the system allocator (glibc). About 15 min an arm (821 s on
# 2026-11-23 at these windows).

# --- on the laptop: hand the branch over without pushing it ---
git -C /Users/jose/devel/pop_var_caller-mem-overhead bundle create \
  /Users/jose/devel/pop_var_caller-mem-overhead/tmp/mem-overhead.bundle 894d3a13..mem-overhead
scp /Users/jose/devel/pop_var_caller-mem-overhead/tmp/mem-overhead.bundle \
  joxi@kimura.comav.upv.es:/media/tomato25_vcfs/pvc_work/mem-overhead.bundle

# --- on kimura: everything under /media/tomato25_vcfs/pvc_work ---
W=/media/tomato25_vcfs/pvc_work
git clone -q ~/pop_var_caller $W/mem-overhead-src
cd $W/mem-overhead-src && git fetch -q $W/mem-overhead.bundle mem-overhead && git checkout -q FETCH_HEAD
CARGO_TARGET_DIR=$W/target-mem ~/.cargo/bin/cargo build --release --example refine_replay
CARGO_TARGET_DIR=$W/target-mem-system ~/.cargo/bin/cargo build --release --no-default-features --example refine_replay
cp $W/target-mem/release/examples/refine_replay        $W/refine_scale/refine_replay_count_mimalloc
cp $W/target-mem-system/release/examples/refine_replay $W/refine_scale/refine_replay_count_glibc

cd $W/refine_scale
# 500-kb cores and margins, as the thinning report's runs on this slice: three windows, comparable.
# Three arms: pvc's allocator as it is, the same with freed memory returned at once, and glibc.
for A in mimalloc mimalloc_purge0 glibc; do
  T=count_${A}_455; B=${A%_purge0}
  EXTRA=""; [ "$A" = mimalloc_purge0 ] && EXTRA="--setenv=MIMALLOC_PURGE_DELAY=0"
  MEMMAX=118G THREADS=48 FIXED=500000,500000 \
    ENVS="--setenv=REPLAY_MARGIN_DENSITY=5 --setenv=MEMORY_TRACE=$W/refine_scale/$T/trace.tsv $EXTRA" \
    ./run_replay.sh $T refine_replay_count_$B slice_ch01_44_455.spill 130402 80000000000
  md5sum $T/out.vcf > $T/out.md5   # every arm must write the same bytes
  rm -f $T/refine_replay $T/out.vcf
done
grep -E 'allocator|window core|last windows|done' count_*_455/run.err
```

What to read: per window, `requested_peak − counted` should be close to 48 × 1.05 kB × the
window's walked records (the walks' arrays), and `kernel_peak − requested_peak` is the slack under
each allocator. The three `out.md5` must agree.

## 8. What changed on the branch

- `tools/counting_allocator`: a global-allocator wrapper that counts the bytes asked for and their
  peak, and reads the resident memory (Linux `statm` and `VmHWM`, macOS `proc_pidinfo`), sampled
  on a thread of its own. A crate of its own because it needs `unsafe`, which the main crate
  forbids; a dev-dependency only.
- `examples/common/window_memory.rs`: the allocator declaration both probes share — mimalloc,
  wrapped, by default; the system allocator, wrapped, with `--no-default-features`; none under
  `dhat-heap` — and the per-window peaks. `MEMORY_TRACE` writes every reading to a file.
- `examples/refine_replay.rs`: **now runs on mimalloc**, as `pvc` does; it ran on glibc before.
  Each window's line adds the count (in all and outside the window, records, core and margin
  lengths), the requested peak and the resident peak. The old `peak` field is kept, now the
  kernel's peak.
- `examples/refine_mem_probe.rs`: the same per-window fields.

The replay's new lines were built and checked with clippy under the default, system-allocator and
all-features builds, but not run: there is no spill on the laptop. The kimura run of §7 is their
first.

The measurements' logs are in the worktree's `tmp/mem/`, one `.err` (per-window lines) and one
`.trace` (readings every 20 ms) per arm.
