# Refinement at cohort scale: memory and the walk's speed (2026-11-13)

Branch `refine-memory-budget`, from main a4954e76. Spec changed: `neighbour_evidence_pass.md` §2.2,
§6. Every change but the budget leaves the written records byte for byte as main writes them.

## 1. What was wrong

On 2026-11-13 `pvc call-from-psps` with `--refine-genotypes on` ran on the tomato25 cohort — 2,176
samples at about three reads a position, 48 threads — over the first 3 Mb of chromosome 0, under a
60 GB cap. Calling finished in 16 minutes. Refinement never finished its first window: memory
reached the cap as the solve began, the solve ran at the cap on 48 cores for 10 minutes, and the
kernel killed the run.

Two things were wrong, and the spec's §6 had predicted neither:

- **Memory.** A window was sized to 120 million sample-records — one sample at one record — at an
  assumed 40 bytes each. Measured on 100 to 800 samples, a window held about 500 bytes a
  sample-record, and each thread's walk held about 6.5 KB a window record besides.
- **Time.** §6's arithmetic gave 33 core-hours for 3,000 samples over a 900-Mb genome. Measured,
  the cost grew in proportion to the samples (773, 1,415, 2,871 and 5,281 CPU-seconds for 100, 200,
  400 and 800 samples over the same 1 Mb), and a sample-record cost about 130 µs in a window after
  the fitting ones: about 150 times §6's figure. Two factors make it up: the five phasing rounds,
  the last four drawing nine paths each with reads, about eight walks' work; and each walk running
  about 20 times slower than the operation count, because it waited on memory — perf counted 0.6
  instructions a cycle over the refinement phase.

The samples-squared template choice was not the problem: 3 CPU-seconds a window at 400 samples.

## 2. The changes

| commit | what | written records |
|---|---|---|
| cd865006 | read alleles held packed, 7 bytes instead of 12, each sample named once a record (`PackedRecordFragments`); the averaged rounds' sums in one flat array a sample instead of two vectors a core record | same bytes |
| 20e8c071 | held VCF lines zstd-compressed, 1,024 at a time on the pool (5.4 times smaller at level 1) | same bytes |
| 1803ded5 | the walk keeps two sites' emissions, not the window's, and computes each twice | same bytes |
| f424325b | the window budget in bytes, from what each record holds (`--refine-window-bytes`, 32 GB); a window over it at the 500-kb floor stops the run with an error | windows change only where the budget binds |
| 1a47369f | the walk's posterior loop without two integer divisions a state; a path's haplotypes and its joints in one array each, not three vectors a site | same bytes |
| b3e20d8b | the forward probabilities kept at the first site of every block of 64 and computed again, block by block, in the backward walk; the paths drawn in the backward walk's sweep | same bytes |

Byte identity was checked by comparing the records of whole runs against main's: 200 samples over
1 Mb (24,000 to 26,000 records) and 400 samples over 4 Mb (93,011 records), chromosome 0. The
neighbour-evidence lib tests pass, among them the walk against the Python reference to 1e-9 and the
stage's same bytes at one thread and four.

**Why recomputing is faster here.** The walk read and wrote every site's forward probabilities and
emissions — about 5 KB a site at 17 templates — three times a walk, over windows of tens of
thousands of sites, on 48 threads at once. Holding less and computing more: the two-site emissions
added a fifth to the instructions and took a third off the cycles; the blocks keep a block's arrays,
about 330 KB, in a core's 2-MB L2 cache.

## 3. What it measured

Refinement phase only (from the end of calling), counted by `perf stat` on the run's process; 200
samples over the first 1 Mb of chromosome 0, `--defaults`, 48 threads:

| build | instructions | cycles | instructions a cycle | peak memory |
|---|---|---|---|---|
| main | 3.39 T | 5.54 T | 0.61 | 11.9 GB |
| two-site emissions (1803ded5) | 4.03 T | 3.50 T | 1.15 | 8.4 GB |
| + no division, flat paths (1a47369f) | 3.63 T | 3.26 T | 1.11 | 8.5 GB |
| + forward blocks (b3e20d8b) | 3.81 T | 1.80 T | 2.12 | 4.7 GB |

Cycles rather than seconds because another job shared the machine during some of these runs.
400 samples over 4 Mb: refinement's wall time 284 s on main and 163 s at b3e20d8b, its peak memory
37.7 GB and 16.0 GB.

**The full cohort.** 2,176 samples over the first 3 Mb of chromosome 0, the cohort's own parameters
file, 60 GB cap:

| build | refinement, wall (48 threads) | refinement, CPU | peak memory after calling | outcome |
|---|---|---|---|---|
| main | — | more than 8 core-hours, unfinished | 60 GB, the cap | killed |
| budget in bytes (f424325b) | 848 s | 8.8 core-hours | 43.2 GB | 50,708 records |
| forward blocks (b3e20d8b) | 483 s | 4.2 core-hours | 37.7 GB | the same 50,708, byte for byte |

At 4.2 core-hours for these 3 Mb, the genome's 783 Mb would take about 1,100 core-hours, 23 hours
on 48 cores, if chromosome 0's density holds. That is an upper bound in one respect: these 3 Mb
are the run's first windows, where every sample also fits its walk's two lengths with twelve more
walks (§3.5), about half of such a window's cost; later windows skip them.

## 4. What is left

The profile is now flat: the walk's backward sweep and posterior about 39% of the refinement phase,
the emissions 17% and the forward step 12% (each computed twice). The largest remaining factors
change what is written and are the owner's to rule on:

- the five phasing rounds, each after the first drawing nine paths;
- the pair states, ordered — each pair of templates held twice, `(k₁, k₂)` and `(k₂, k₁)`, though
  the model is symmetric in the two copies; summing them once would halve most of the walk's
  arithmetic and change the last bits of the posteriors;
- a reference panel: all rounds on a panel of samples, the rest solved once against it.

**Not measured here.** Every figure is from chromosome 0, the unplaced-scaffold bin, where a sample
carried about 22 read alleles a record — likely collapsed repeats; elsewhere the read alleles, the
largest per-sample cost, may be fewer. Calling on chromosome 1 with `PROF_REGION` ran out of memory
in seconds at 400 samples, a calling-side problem not looked into.
