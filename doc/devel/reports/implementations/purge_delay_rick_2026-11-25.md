# mimalloc's purge delay on rick: 100 tomato samples, 32 threads (2026-11-25)

`pvc` now tells mimalloc to return freed memory to the system at once (purge delay 0) instead of
holding it for mimalloc's default 1,000 ms (`mimalloc_purge_delay_2026-11-25.md`). That change was
measured on the laptop, on 63 samples at 6 threads. This report repeats the measurement on rick, on
a larger cohort and at 32 threads, with main at 4b29ce57.

## 1. Result

At 32 threads, returning freed memory at once lowered the run's peak resident memory by 3.3 GB,
from 10.6–10.8 GB to 7.4 GB (−31%). It cost about 100 s of CPU time in 6,130 (+1.6%), almost all
of it system time (the kernel taking pages back: 306–312 s against 213–216 s). Wall time did not
change: 802–806 s in both arms, which is the 1-s resolution of the measurement. The VCF was the same
in all six runs.

| Arm | Wall, s | CPU (user + sys), s | of which sys, s | Peak resident, GB |
|---|---|---|---|---|
| A: `MIMALLOC_PURGE_DELAY=1000` (mimalloc's default) | 803–806 | 6,130–6,133 | 213–216 | 10.58–10.78 |
| B: no variable (`pvc`'s built-in 0) | 802–804 | 6,232–6,246 | 306–312 | 7.36–7.39 |

Ranges are over three repeats each, run A B A B A B; §5 has every run. GB are 10⁹ bytes.

**Next to the laptop.** The Linux container on the laptop (63 samples, 8 Mb in 80 regions,
6 threads) saved 2.0 GB of 6.9 (−29%) for +3% CPU and +2% wall. Here the saving is the same share
of the peak and larger in bytes, and the CPU cost is half the share, with no wall-time cost. The
expectation recorded in the laptop report (§6 there), that the system-time cost could grow at more
threads, did not hold at 32: +95 s of system time on a 13-minute run.

## 2. Input and setup

- **Machine**: rick, 32 cores, 125 GiB, Debian 12 (Linux 6.1, glibc 2.36), idle apart from these
  runs. One run at a time.
- **Cohort**: the 157 tomato CRAMs of `~/tomato25/bench_crams/PRJNA454805/*.p1.cram`. Their read
  groups name **100 samples** (several CRAMs share an `SM`), so the VCF has 100 sample columns.
- **Reference**: `~/tomato25/genome/S_lycopersicum_chromosomes.4.00.fa`, with a repeat catalog
  built for these runs by `pvc repeat-catalog` (6,441,172 repeats, 5 s).
- **Region**: `SL4.0ch01:10,000,000–15,000,000` (5 Mb), sized from a 2-Mb pilot (6 min 45 s) so a
  run takes about 13 minutes. 216,955 records called, 122,955 written.
- **Command**: `pvc call-from-alignments --defaults --ploidy 2 --threads 32 --refine-genotypes on
  --regions region5.bed`, with the hidden-duplication filter on (the default), so refinement runs
  after calling has finished and let its memory go.
- **Build**: `cargo build --release` of main 4b29ce57, run natively (rick has no container
  runtime). The same binary in every arm.
- **Measured with** a small wrapper (`tmp/purge_rick/measure2.py`), because rick has no
  `/usr/bin/time`: wall from the clock, user and system CPU and the peak resident set from the
  kernel's resource usage for the finished process (`wait4`, the same source `/usr/bin/time -v`
  reads), and `VmRSS` from `/proc/<pid>/status` every second, to tell the stages apart.
- **Open files**: the default limit of 1,024 is too low for 157 CRAMs at 32 threads (the first
  pilot stopped with "Too many open files"); the runs raised it with `ulimit -n 65536`.

## 3. Which stage holds the peak

Calling holds it, in both arms. Calling took 758–761 s of the 802–806 and refinement 41–43 s. The
peak during calling was 10.40–10.69 GB with the default delay and 7.20–7.32 GB at 0. When calling
ends the resident memory drops to about 0.2–0.4 GB in both arms, and refinement then rises to its
own peak, 3.3–3.5 GB at the default and 2.7–3.0 GB at 0.

So almost all of the 3.3 GB saved is in calling, as on the laptop. Refinement, which the
expected-peak estimate is about, is not the run's peak at 100 samples over 5 Mb.

## 4. Refinement's expected peak against the measured one

The run logs refinement's expected peak twice: before calling, from the budget and the longest
contig, and at each window, from the window's own records.

**Before calling**, at the default 32-GB budget: 30.7 GB. Nothing about the windows is known yet,
so this is an upper bound set by a window of chromosome 1 at 150 records a kb (15.6 GB counted,
10.9 GB of working arrays for 32 threads). The real windows here held about 1.2 GB each, so the
bound was eight times what refinement used. It is far below the 90% stop (23% of the 135 GB
`MemTotal`), so it only matters as a log line.

**At the window**: the highest window (core 12.0–14.0 Mb, 95,023 records, 62,537 walked) was
logged as **3.8 GB** (1.2 GB counted, 2.1 GB of arrays for 32 threads, plus 16% for the
allocator). The resident memory during refinement peaked at:

| Arm | Refinement's peak resident, GB | Estimate, GB | Estimate over measured |
|---|---|---|---|
| A, delay 1,000 ms | 3.29, 3.33, 3.45 | 3.8 | +10 to +15% |
| B, delay 0 | 2.73, 2.92, 2.95 | 3.8 | +29 to +39% |
| C, delay 0, 1.2-GB budget | 2.71 | 3.8 | +40% |

The estimate is above every measured peak, so it errs safe here. Two cautions on reading the gap:

- The resident memory was read once a second, and the windows take 5 to 18 s, so a short spike
  could be missed. The kernel's own peak counter covers the whole run, where calling dominates, so
  it cannot settle this.
- Even at the default delay the measured peak is below counted + arrays (3.3 GB) without any
  allocator slack. Either the arrays term (1.05 kB a walked record a thread, measured at 300
  samples and 6 threads) is high at 100 samples, or not all 32 threads hold their arrays at once in
  a window of 100 samples. This run cannot tell which.

With the delay at 0 the allocator's 16% is too much here, as the guard report expected (§9 there):
the laptop measured 6–7% at delay 0. This does not yet support a new constant, since refinement was
small here; the tomato25 replay on kimura is the case that sets it.

**Arm C** set `--refine-window-bytes 1200000000` to make the budget bind. It bound only on the second
window, whose core shrank from 12.00–14.00 Mb to 12.00–13.95 Mb (1.20 GB of 1.20 GB). Its
before-calling estimate fell to 4.1 GB, now set by the budget, close to the window's 3.8 GB. Its VCF
differs from arms A and B, as expected when a window boundary moves.

## 5. Every run

| Run | Arm | Wall, s | CPU, s | sys, s | Peak resident, GB | Calling, s | Refinement, s | Calling's peak, GB | Refinement's peak, GB |
|---|---|---|---|---|---|---|---|---|---|
| A1 | delay 1,000 | 803.3 | 6,133.4 | 216.2 | 10.78 | 759 | 42 | 10.67 | 3.33 |
| B1 | delay 0 | 803.3 | 6,246.4 | 312.3 | 7.39 | 758 | 43 | 7.32 | 2.95 |
| A2 | delay 1,000 | 803.3 | 6,129.5 | 213.4 | 10.58 | 759 | 41 | 10.40 | 3.45 |
| B2 | delay 0 | 802.3 | 6,231.5 | 310.1 | 7.36 | 758 | 42 | 7.20 | 2.92 |
| A3 | delay 1,000 | 806.3 | 6,133.2 | 212.9 | 10.69 | 761 | 42 | 10.69 | 3.29 |
| B3 | delay 0 | 804.3 | 6,232.9 | 305.7 | 7.37 | 759 | 42 | 7.24 | 2.73 |
| C1 | delay 0, 1.2-GB budget | 803.3 | 6,221.5 | 311.2 | 7.47 | 759 | 42 | 7.29 | 2.71 |

Peak resident is the kernel's figure for the whole run. The two stage peaks are the highest of the
once-a-second samples: calling up to its `calling: done` log line, refinement from 3 s before its
first window's log line to the end. Stage times are from the log's timestamps.

**Output.** With the `##commandline=` and `##pvcLog=` header lines left out, the six A and B VCFs
have the same MD5, `04fab3844c9a72293f3d4dc36fe4f258`. C1's is `33dae2c11ecdab1f288ebf11adb2ef3a`.

Raw logs, the per-second memory traces and the tables are in `tmp/purge_rick/` (not committed).

## 6. What this does not cover

- **Refinement as the run's peak.** At 100 samples over 5 Mb, refinement used under half of
  calling's memory, so this run says how the purge delay acts on calling at 32 threads, not on a
  refinement window of tens of GB. The tomato25 replay on kimura (2,176 samples, 48 threads) is
  still the measurement for that, and for the allocator constant in the estimate.
- **One cohort at one depth.** These tomato samples are low coverage; a deep cohort holds more
  reads a position in calling, and its saving was not measured.
