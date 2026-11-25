# mimalloc's purge delay on a whole `pvc` run (2026-11-25)

`pvc` allocates memory through mimalloc. mimalloc 3 keeps memory the program has freed for one
second (its *purge delay*, 1,000 ms by default) before returning it to the system. The refinement
report of the same day (`refine_memory_overhead_2026-11-25.md` §4) measured what that costs inside
refinement alone: resident memory 17 to 24% above the bytes in use at the default, 6 to 7% with the
delay at 0. This report measures a whole run, and records the change that followed: **`pvc` now
sets the delay to 0 at start-up.**

Branch `purge-delay`, from `mem-overhead` (af6ba479). Everything ran on the laptop (macOS,
arm64), natively and in the Linux development container; nothing ran on kimura or rick.

## 1. Result

On Linux, returning freed memory at once lowered the run's peak resident memory by about 2 GB,
from 6.8–7.0 GB to 4.9 GB (−29%). It cost about 17 s of CPU time in 565 (+3%), most of it the
kernel's work of taking pages back, and about 2.6 s of wall time in 151 (+2%). The VCF was the same,
byte for byte, in every run of every arm, on both systems.

| Linux, 6 threads | Wall, s | CPU (user + sys), s | of which sys, s | Peak resident, GB |
|---|---|---|---|---|
| Delay 1,000 ms (mimalloc's default) | 150.1–152.0 | 560.4–573.2 | 16.8–17.6 | 6.82–7.04 |
| Delay 10 ms (mimalloc 2's default) | 152.7–152.8 | 575.0–575.4 | 24.8 | 5.24–5.26 |
| Delay 0, built into `pvc` | 152.9–154.3 | 578.9–586.5 | 30.1–30.5 | 4.91–4.92 |

The ranges are over three repeats each; §4 has every run.

**Decision: set the delay to 0 in `pvc`.** Memory is what limits the large cohorts (on tomato25 the
refinement's window budget is chosen to keep the run inside kimura's memory); 3% more CPU is a small price for
about 30% less peak memory on this input. 10 ms keeps most of the saving (5.25 GB) for about half
the CPU cost, and remains available to a user as `MIMALLOC_PURGE_DELAY=10`.

## 2. Input and setup

- **Cohort**: tomato1, the 63 tomato accessions of `benchmarks/tomato1/`, in their BAQ-free CRAMs
  (`crams_nobaq/`), over the 80 regions of 100 kb in `regions.bed` (8.0 Mb). About three reads a
  position. Chosen because it runs the whole pipeline — reading, calling and refinement — on a
  real cohort in about 2.5 minutes, so every arm could be repeated three times and interleaved.
- **Command**: `pvc call-from-alignments --defaults --ploidy 2 --threads 6 --refine-genotypes on`,
  with a repeat catalog rebuilt from the reference for this run (the one in
  `benchmarks/tomato1/crams/` is from an older reference digest and is refused).
- **Stages** (from the log, Linux, default delay): calling 2 min 1 s, refinement 29–30 s. The
  refinement windows here are small, 0.01 to 0.04 GB of a 32-GB budget each, so the run's memory
  peak, and the 2 GB saved, are almost certainly in the calling stage, not in refinement. The extra
  time is there too: calling took 2 min 3–4 s with the delay at 0, refinement still 29–30 s.
- **Arms** were run in turn (one of each, then the next repeat), so that load from other work on
  the laptop falls on every arm alike.
- **Measured with**: on macOS, `/usr/bin/time -l`; in the container, the clock for wall time,
  `/proc/<pid>/stat` for CPU and the kernel's peak (`VmHWM`) for resident memory, read every 0.2 s.
- **Builds**: the macOS arms used two binaries, the branch before the change (`base`, the delay set
  only through the environment) and after it (`builtin`). The Linux arms used the built-in binary
  only, with `MIMALLOC_PURGE_DELAY=1000` for mimalloc's default; macOS shows that this matches the
  binary without the change (§3).

## 3. macOS

On macOS the figure that moved is the *peak memory footprint*, what macOS charges the process for:
7.85–7.96 GB at the default, 6.81–6.86 GB at 0 (−1.1 GB, −14%). The *maximum resident set size*
did not move: 9.74–9.90 GB in every arm. Why macOS's resident count keeps the returned pages was not
looked into; the Linux figures are the ones to read for kimura.

The first two repeats ran while another job was compiling on the laptop, which shows in the
spread (wall 132 to 158 s within one arm). The third repeat ran on a quiet machine:

| macOS, repeat 3 | Wall, s | CPU, s | of which sys, s | Peak footprint, GB |
|---|---|---|---|---|
| `base`, default | 132.4 | 543.4 | 10.7 | 7.93 |
| `base`, `MIMALLOC_PURGE_DELAY=0` | 133.5 | 555.2 | 21.7 | 6.82 |
| `base`, `MIMALLOC_PURGE_DELAY=10` | 133.4 | 550.3 | 15.2 | 6.91 |
| `builtin` | 134.5 | 564.7 | 21.5 | 6.85 |
| `builtin`, `MIMALLOC_PURGE_DELAY=1000` | 133.3 | 553.7 | 10.8 | 7.95 |

**The environment variable still wins over the built-in setting**: `builtin` with
`MIMALLOC_PURGE_DELAY=1000` gave the default's footprint (7.81–7.95 GB against 7.85–7.96) and
system time (10.8–13.1 s against 10.7–12.5), in all three repeats.

## 4. Every run

macOS (`/usr/bin/time -l`; GB are 10⁹ bytes):

| Arm | Repeat | Wall, s | CPU, s | sys, s | Max RSS, GB | Peak footprint, GB |
|---|---|---|---|---|---|---|
| `base`, default | 1 | 138.2 | 568.5 | 11.4 | 9.86 | 7.96 |
| `base`, delay 0 | 1 | 149.0 | 619.6 | 24.2 | 9.90 | 6.86 |
| `base`, delay 10 | 1 | 150.0 | 619.7 | 17.6 | 9.79 | 6.87 |
| `builtin` | 1 | 149.0 | 625.4 | 24.7 | 9.80 | 6.84 |
| `builtin`, delay 1000 | 1 | 158.2 | 649.8 | 13.1 | 9.88 | 7.89 |
| `base`, default | 2 | 146.7 | 605.2 | 12.5 | 9.90 | 7.85 |
| `base`, delay 0 | 2 | 143.7 | 596.1 | 23.9 | 9.87 | 6.81 |
| `base`, delay 10 | 2 | 132.1 | 542.8 | 15.2 | 9.83 | 7.02 |
| `builtin` | 2 | 133.0 | 557.3 | 21.2 | 9.89 | 6.82 |
| `builtin`, delay 1000 | 2 | 132.4 | 548.4 | 10.8 | 9.82 | 7.81 |
| `base`, default | 3 | 132.4 | 543.4 | 10.7 | 9.88 | 7.93 |
| `base`, delay 0 | 3 | 133.5 | 555.2 | 21.7 | 9.82 | 6.82 |
| `base`, delay 10 | 3 | 133.4 | 550.3 | 15.2 | 9.78 | 6.91 |
| `builtin` | 3 | 134.5 | 564.7 | 21.5 | 9.74 | 6.85 |
| `builtin`, delay 1000 | 3 | 133.3 | 553.7 | 10.8 | 9.79 | 7.95 |

Linux, in the development container (Debian 12, 6 CPUs, 24 GB), built-in binary:

| Arm | Repeat | Wall, s | CPU, s | sys, s | Peak resident, GB |
|---|---|---|---|---|---|
| delay 1000 | 1 | 152.0 | 573.2 | 17.6 | 7.02 |
| built-in (0) | 1 | 154.3 | 586.5 | 30.5 | 4.90 |
| delay 10 | 1 | 152.8 | 575.3 | 24.8 | 5.26 |
| delay 1000 | 2 | 150.5 | 561.6 | 17.1 | 7.04 |
| built-in (0) | 2 | 152.9 | 578.9 | 30.3 | 4.91 |
| delay 10 | 2 | 152.8 | 575.4 | 24.8 | 5.26 |
| delay 1000 | 3 | 150.1 | 560.4 | 16.8 | 6.82 |
| built-in (0) | 3 | 153.3 | 580.3 | 30.1 | 4.92 |
| delay 10 | 3 | 152.7 | 575.0 | 24.8 | 5.24 |

**Output.** All 24 VCFs, macOS and Linux, have the same MD5 once the `##commandline=` header line,
which names the binary's path, is left out: `5edd694b9227ac31e2f93da33d633489`. Whole files
written by the same binary are identical too (`c819cb86…` for every `base` run, `72199c4b…` for
every `builtin` run on macOS).

## 5. What changed

- `tools/mimalloc_settings`: a crate of its own, because calling mimalloc's C option functions
  needs `unsafe`, which the main crate forbids. `return_freed_memory_at_once()` sets the purge
  delay's *default* to 0 (`mi_option_set_default`). mimalloc reads its environment variables when
  the process starts, before `main`, and marks an option it found there as set; setting the
  default leaves such an option alone, which is how `MIMALLOC_PURGE_DELAY` keeps priority. A test
  checks that the call moves mimalloc 3's 1,000 ms to 0, and leaves a value from the environment
  alone.
- `src/main.rs`: calls it first thing in `main`, under the same features as the mimalloc
  declaration (so not in a `dhat-heap` build).
- `Cargo.toml`: the crate is an optional dependency, switched on by the `alloc-mimalloc` feature
  with `mimalloc` itself.

Not changed: the memory probes in `examples/` (`refine_replay`, `refine_mem_probe`) still run
mimalloc at its default unless the environment says otherwise. Their arms in the refinement
report set the delay through the environment, which still works.

## 6. What this does not cover

- **One input, one corner of the range**: 63 samples at about three reads a position, 6 threads.
  The saving is freed memory held for a second, so it should grow with how fast a run frees and
  asks again — more threads, larger windows — but that was not measured here. The kimura replay
  prepared in the refinement report (§7 there) would measure it at 48 threads on tomato25.
- **CPU at more threads.** The extra cost is system time, the kernel taking back pages: +13 s on
  this run. At 48 threads, with more pages returned concurrently, it could be a larger share; the
  same kimura run would show it.
