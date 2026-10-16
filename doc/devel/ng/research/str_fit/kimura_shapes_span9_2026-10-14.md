# STR counting fit and shapes at ±9 on the tomato25 cohort (kimura, 2026-10-13/14)

Branch `str-fit-c` at **`e03017c9`** ("feat(research): allele-length shapes with the two edge classes free; tomato2 at ±9; kimura prompt").
The cohort is tomato25: 2,169 samples, 2,651 read groups, about 3 reads a position.

**The shape analysis was stopped on the owner's decision on 2026-10-14 at 08:41.** 12 of its 102 single-shape jobs had finished, covering one-base runs of 8 to 15 repeats. None of them used the eight-number `ref_step_power6_edges`: no job fitting it finished. So questions (a) and (b) have partial answers for short one-base runs only, and no answer for the edge shape or for runs of 17 repeats and longer.

## (c) The counting fit at ±9

| | ±9 (this run) | ±6 (2026-10-11) | change |
|---|---|---|---|
| wall time | **52 min 58 s** | 15 min 38 s | 3.4 times longer |
| peak memory | **45.6 GiB** (47,812,836 kB) | 31.5 GiB (33,034,756 kB) | 1.45 times more, +14.1 GiB |
| classes fitted | 34 of 141 | 34 of 141 | same |
| spread (`concentration`) at a bound (0.001 or 10,000) | **none**; range 0.158 (3-base, 6 repeats) to 4.09 (2-base, 16 repeats) | none; range 0.155 to 3.97 | |

- **Same evidence at both windows:**
  - 86,688 tracts in 141 classes;
  - 45,432 tracts left out by the guard;
  - 66,203,224 (tract, sample) pairs with reads, out of the 89,484,264 the kept tracts make (74.0 in every 100);
  - 19.4 GiB once arranged.
- **Where the extra memory and time come from:** the evidence stage costs the same at ±9 as at ±6, so the 14 GiB and 37 minutes extra are spent in the per-class fits.
- **Fewer rounds:** 18 of the 34 classes stopped after 2 rounds of the fit. At ±6, 6 classes stopped after 2 rounds and 2 after 1. The rest took 3 (9 classes), 4, 5, 6 or 10 rounds.
- **Classes not fitted:** the other 107 classes have fewer than 8 tracts with reads.

**Genotype rows at ±9:**
- **Time and memory:** 8 min 35 s, with a 31.7 GiB peak (33,193,200 kB). At ±6 this step took about 5 minutes.
- **Size on disk:** 95 GB for 34 classes, against 46 GB at ±6.
- **Row counts:** the same as at ±6, for example 9,685 rows for the 3-base class at 6 repeats.

## (a) and (b): the classes that finished

All finished classes have 100 or more tracts, capped at 300 for the fit. No class with fewer than 100 tracts finished, so that group has no result. Each class's draws are seeded by the class, so jobs for different shapes of one class share the same tracts, resamples, free fit and noise.

**Single-sample calls changed per 10,000 rows.** The noise column is the free fit swapped for one of its own resampled refits. "—" marks a job that did not finish.

| class | rows | noise | `ref_step_power6` (6 numbers) | `ref_power4_edges` (4 + 2 edges) | `ref_step_power6_edges` (6 + 2 edges) |
|---|---|---|---|---|---|
| 1-base, 8 | 514,608 | 3.8 | — | 8.6 | — |
| 1-base, 9 | 519,907 | 4.3 | 4.0 | 30.3 | — |
| 1-base, 10 | 514,675 | 9.1 | 7.0 | 33.0 | — |
| 1-base, 11 | 495,783 | 10.9 | 12.0 | 65.2 | — |
| 1-base, 12 | 491,930 | 65.0 | 63.5 | 37.7 | — |
| 1-base, 13 | 472,081 | 30.9 | 40.4 | 111.6 | — |
| 1-base, 15 | 445,683 | 30.8 | 102.0 | — | — |

The same calls per 10,000 summed over the finished classes:

| shape | classes | rows | shape | noise on the same classes |
|---|---|---|---|---|
| `ref_step_power6` | 6 (9–13, 15) | 2,940,059 | **36.5** | 24.7 |
| `ref_power4_edges` | 6 (8–13) | 3,008,984 | **46.8** | 20.2 |

What the finished classes show:

- **`ref_step_power6` at 9 to 12 repeats:** it changes as many calls as the noise, 0.8 to 1.1 times, as it did at ±6.
- **`ref_step_power6` at 13 and 15 repeats:** it changes more, 1.3 times the noise at 13 repeats and 3.3 times at 15 repeats (102 against 31). At ±6, 15 repeats was 32 against 37.
- **`ref_power4_edges`:** it changes more calls than the noise in 5 of its 6 classes, 2.3 to 7.0 times. The exception is 12 repeats, where the noise is high (65).
- **Log-likelihood loss,** the free fit's composite log-likelihood minus the shape's: 909 to 3,159 for `ref_step_power6`, and 335 to 8,052 for `ref_power4_edges`. It depends heavily on which 300 tracts are drawn.

**The free fit's edges at ±9 hold little for these lengths.** The share at "9 or more shorter" is 0.000 to 0.001 for 8 to 13 repeats, and 0.006 at 15 repeats. The share at "9 or more longer" is 0.000 to 0.001. These are the classes below where the ±6 bump appeared, so they don't show whether the wider window removes it.

**A smaller bump appears one step inside the longer edge.** In every finished class, the free share at +8 exceeds the share at +7:

| repeats | 9 | 10 | 11 | 12 | 13 | 15 |
|---|---|---|---|---|---|---|
| share at +8 | 0.003 | 0.003 | 0.005 | 0.005 | 0.008 | 0.011 |
| share at +7 | 0.000 | 0.000 | 0.002 | 0.002 | 0.002 | 0.004 |

Reads are recorded only to ±8, which may be the cause. I didn't check this.

**(b) The numbers of `ref_step_power6` at ±9:**

| repeats | ref share | shorter split | −1 of shorter side | +1 of longer side | 2-step weight below | power below | 2-step weight above | power above |
|---|---|---|---|---|---|---|---|---|
| 9 | 0.851 | 0.430 | 0.897 | 0.794 | 0.076 | 3.16 | 0.127 | 2.10 |
| 10 | 0.774 | 0.449 | 0.770 | 0.661 | 0.208 | 3.40 | 0.277 | 2.36 |
| 11 | 0.616 | 0.538 | 0.717 | 0.634 | 0.283 | 3.57 | 0.300 | 2.25 |
| 12 | 0.594 | 0.532 | 0.564 | 0.647 | 0.518 | 3.21 | 0.229 | 1.70 |
| 13 | 0.547 | 0.504 | 0.541 | 0.618 | 0.488 | 2.57 | 0.278 | 1.87 |
| 15 | 0.516 | 0.512 | 0.438 | 0.601 | 0.602 | 1.97 | 0.262 | 1.57 |

How the numbers move from 9 to 15 repeats:
- **Ref share:** falls from 0.85 to 0.52.
- **Share one repeat shorter, of the shorter side:** falls from 0.90 to 0.44.
- **Two-step weight below:** rises from 0.08 to 0.60.
- **Power below:** rises to 3.6 at 11 repeats, then falls to 2.0 at 15. At ±6 it was 3.9 at 11 and 1.46 at 15.
- **Shorter split, share one repeat longer, two-step weight above:** flat, at 0.43–0.54, 0.60–0.79 and 0.13–0.30.
- **Power above:** wanders between 1.6 and 2.4.

The question that mattered was whether the power below still falls to 0 from 19 repeats up once the edge has its own share. **It is not answered:** no class above 15 repeats finished.

## The cost of the shape step

| job | took |
|---|---|
| `ref_step_power6`, one class with 300 tracts | 6 h 41 min to 9 h 42 min (6 jobs) |
| `ref_power4_edges`, one class with 300 tracts | 7 h 2 min to 10 h 1 min (6 jobs) |
| `ref_step_power6_edges`, one class with 300 tracts | **none finished**; the longest ran 11 h 36 min (1-base, 8 repeats) |
| the timing run: all three shapes on 1-base, 8 repeats, one process | **stopped at 17 h 25 min, unfinished** (started 15:16) |

The script prints a class only when the class is done, so how far the stopped jobs had got is not known.

For comparison, at ±6 one process fitted five shapes on that same class in 3 h 34 min. The rows here are twice as wide (190 genotypes against 91).

How the run went:
- **Split:** after the timing run passed 6 hours, the owner chose to split the work into one job per class and shape, 102 jobs.
- **Concurrency:** first 23 jobs at once, then 40, once each job was seen to settle at about 2.2 GB after loading. Loading peaks at about 13 GB.
- **What ran:** 35 jobs started, and 12 finished, all with exit 0. The other 23 were stopped mid-fit, and 67 never started.

## How this differs from the task as written

- **Worktree:** `tmp/pvc_str_fit_c`, not `../`, because the project rules forbid writing outside the repository.
- **No catalog,** as in the ±6 run, which used none either.
- **Parameters file:** `--fitted` is a copy of the 10-06 parameters file that this branch can read, `tmp/fit_precision_kimura_2026-10-06.for_span9.parameters.toml`. `str-fit-c` refused the original with "unknown field `length_spectrum_by_period`" and "unknown field `length_spectrum_by_stratum`", because it renamed those fields to `allele_length_variation_*`. The copy differs from the original in three places:
  - the two renamed lists are empty, under their new names;
  - the fallback concentration's warrant is `supplied` instead of `fitted_here`.

  The examples read only `[inbreeding]`, each sample's homozygote excess, and it is unchanged: 2,169 of 2,169 samples matched. The two refused attempts' logs are kept as `*.attempt1.log` and `*.attempt2.log`.
- **Shapes split:** one process per class and shape, instead of `--jobs 16` per class. This gives the same numbers, because each class's draws are seeded by the class itself.
- Nothing committed; no committed file changed.

## Paths (all under `/home/joxi/pop_var_caller/`)

- Worktree: `tmp/pvc_str_fit_c` at `e03017c9`
- Launchers: `tmp/str_fit_c_kimura_span9.sh` (build, counting fit, rows, timing run), `tmp/str_fit_c_kimura_span9_split.sh`, `tmp/str_fit_c_kimura_span9_rest.sh`, `tmp/str_fit_c_kimura_span9_rest2.sh` (the split shapes run and its two queues)
- Status: `tmp/str_fit_c_kimura_span9.status.log`
- Parameters copy: `tmp/fit_precision_kimura_2026-10-06.for_span9.parameters.toml`
- Counting fit: `tmp/str_counted_kimura_span9.tsv`, log `tmp/str_counted_kimura_span9.log` (failed attempts `tmp/str_counted_kimura_span9.attempt{1,2}.log`)
- Genotype rows: `tmp/str_rows_kimura_span9/` (95 GB), log `tmp/str_rows_kimura_span9.log`; per-class link directories `tmp/str_rows_kimura_span9_byclass/`
- Timing run (stopped, empty): `tmp/str_shapes_kimura_span9_1_8.{txt,log}`
- Shapes, per class and shape: `tmp/str_shapes_kimura_span9_split/<class>.<shape>.{txt,log,numbers.tsv}`; job lists `jobs.txt`, `remaining.txt`, `remaining2.txt`; `done.txt`, which misses some finished jobs, so count finished jobs by `Exit status: 0` in the logs instead
