# Performance Review: estimate-parameters
**Date:** 2026-10-15
**Reviewer:** rust-performance-review skill (single reviewer covering every category; no sub-agents, by the lead's instruction)
**Scope:** `pvc estimate-parameters` — mainly the repeat-tract (STR) penalised fit (`ssr_fit/penalised_fit.rs`, `ssr_fit/expected_counts.rs`), secondarily the SNP/indel fit's passes and convergence (`joint/fit.rs`)
**Verdict:** Apply the listed wins (STR: H1, H2, measured and bit-identical); run experiments for the rest; three SNP/indel convergence changes are owner decisions and carry the largest gain on real data
**Hot-path evidence:** two sampling profiles supplied by the lead; four sampling profiles, two allocation profiles (dhat), a per-E-step count trace and interleaved instruction counts taken in this review

---

## 1. Scope and constraints

**Reviewed against** branch `str-penalised-fit`. The review started at `53a3e4f9`; during it the lead
committed `ae239e26` ("a round's read-length logs are taken once"), which changes `RoundModel::score`.
**Every measurement and patch below that is quoted as final was re-done on `ae239e26`**; the
`53a3e4f9` numbers are kept only where they are the evidence for a count that the commit did not move
(rounds, E-steps, rows). Line numbers cite `ae239e26`.

**What the targets are.** Wall time and peak memory, with the fitted numbers bit-identical — checked
by the reference fixture (`scripts/research/str_fit_fixture.py compare --require exact`) — or moved
by an owner decision with the move measured. The owner prefers simpler code to small gains
(`tmp/handoff_perf_then_D_2026-10-15.md` §2). The range is one sample to thousands and three reads to
hundreds (CLAUDE.md); kimura (2,169 samples, 66.2 million (locus, sample) pairs with reads, 19.4 GiB
of evidence, spec `str_penalised_fit.md` §3) is the scale milestone D3 must time.

**Hardware.** Apple M-series, 18 cores: 6 performance and 12 efficiency. All runs in this review used
`RAYON_NUM_THREADS=6`, on a machine shared with other work (load 2.6 to 17 during the runs), so **wall
times here are indicative only**. The primary measure is `instructions retired` from `/usr/bin/time -l`,
which repeated to within 0.03% for the same binary (1,253,186,783,518 and 1,252,830,563,594 at
1,000 × 30); memory is measured by dhat's peak live heap, which is deterministic.

**What one run of the STR bench means.** `ng_str_fit_bench --period 1 --truth pcr --seed 1
--loci-at-most 100 --methods D_tract_prior`: 14 classes (repeat counts) of one-base repeats, 971 loci,
fitted with production settings including standard errors. "1,000 × 30" is 1,000 samples at 30 reads
a locus; "1,000 × 3" the same at 3 reads.

**Evidence available.**

| what | where | taken by |
|---|---|---|
| 59 tomato2 samples, `pvc estimate-parameters`, 18 threads: 61 min wall, SNP/indel fit 60 min, repeat-tract fit 15 s, peak footprint 3.95 GB | `tmp/perf_ep/t2_59.err`, profile `tmp/perf_ep/sample_snp.txt` | lead |
| STR bench 2,000 × 30, 18 threads, at `53a3e4f9`: 21.3 s wall, 1,062 MiB peak RSS | `tmp/perf_ep/sample_str_p1.txt`, `base_p1_2000x30.*` | lead |
| per-E-step counts (rows, distinct rows, chunks, table bytes) at 1,000 × 30 and 1,000 × 3 | `tmp/perf_review/runs/trace_p1_1000x*.err`, summariser `analyse_trace.py` | this review |
| sampling profiles at 1,000 × 30 and 1,000 × 3, before and after the patches, on `ae239e26` | `tmp/perf_review/prof/` | this review |
| dhat at 500 × 30, before and after, on `ae239e26` | `tmp/perf_review/dhat_base/`, `dhat_v123/`, summariser `dhat_top.py` | this review |
| fixture grid, 144 runs below 2,000 samples, for each patch set | `tmp/perf_review/grid_v1`, `grid_v3`, `grid_new_v123` | this review |

Runs at 2,000 samples and the 59-sample run were not repeated, by instruction.

**In scope:** [penalised_fit.rs](../../../../src/parameter_estimation/joint/ssr_fit/penalised_fit.rs),
[expected_counts.rs](../../../../src/parameter_estimation/joint/ssr_fit/expected_counts.rs),
[fit.rs](../../../../src/parameter_estimation/joint/fit.rs) (`maximise`, `expectation_pass`,
`one_position`, `EvidenceCursor`).

**Out of scope:** calling; the census read; the evidence gathering (`gather_strata`), whose memory
the spec already names as the largest (§3) and leaves to a later step.

**Categories covered:** methodology (always); hot_loops and allocations (the E-step's per-row work);
concurrency (rayon bursts, the idle third); data_layout (per-row tables); float_portability (every
proposal states its effect on bits). io_and_syscalls not applicable.

---

## 2. Verdict

**Apply the listed wins.** Two STR changes are measured, bit-identical on all 144 fixture runs below
2,000 samples, and small in code:

- **H1 — the standard errors stop re-scoring every row for steps that cannot change the scores.**
  At 1,000 × 30 the errors took 6.44 s and now take 2.87 s; the whole fit went from 12.6 s to 8.8 s,
  and the instructions the process retires fell by 32% (1,253 to 853 billion). At 1,000 × 3 the gain
  is 4% of instructions, because there the errors are a small part of the fit.
- **H2 — the table of distinct rows is built without a heap copy per row.** At 1,000 × 3 the process
  retires 19% fewer instructions (503 to 406 billion) and the fit took 4.98 s against 5.63 s; at
  1,000 × 30, 3.6% fewer.

**What decides `estimate-parameters` on real data is not the STR fit.** On 59 tomato2 samples the
repeat-tract fit took 15 s of a 61-minute run; the SNP/indel fit took 60 minutes because each of its
three starts ran to the 200-pass limit. Three ways to stop it sooner (L5a–L5c in §5, all owner
decisions because they change numbers or the stopping rule) are the largest gain available in this
command. **My recommendation is L5a first**: let a later start stop when it heads where an earlier
start is heading, even when that start hit the limit. The three starts ended within 2 log-likelihood
units of each other (−7,985,053, −7,985,051, −7,985,053), so up to two thirds of the 60 minutes were
spent re-deriving one answer.

**Why the STR threads idle a third of the time** is answered in H3: it is serial work on the main
thread between parallel bursts, not imbalance inside the bursts.

---

## 3. Measurement plan

1. **Gate every STR change on the fixture and on instruction counts, not wall time.**
   `RAYON_NUM_THREADS=6 uv run --no-project python scripts/research/str_fit_fixture.py run --out <dir>
   --only <the 144 names below 2,000 samples> --jobs 1 --bench <binary>` then `compare --require exact
   --partial`. CPU: `/usr/bin/time -l <bench> … 2> err; grep 'instructions retired' err`, two runs per
   variant interleaved A B B A. A change is worth keeping at a fall of 3% or more in instructions on
   either the 1,000 × 30 or the 1,000 × 3 cell.
2. **Gate memory on dhat's peak, not RSS.** Build the bench with the dhat hook (§4) under
   `--profile profiling --features bench-fixtures,dhat-heap`, run 500 × 30, read `At t-gmax`. RSS and
   footprint are not repeatable enough here: the same binary at 1,000 × 30 gave footprints of 354,
   431, 438 and 467 MB in four runs.
3. **Add one bench cell with a large class at depth 30** (for example `--samples 500 --depth 30
   --loci-at-most 2000`), because the grid caps classes at 100 loci and the cost that grows with the
   square of a class's rows (L1) is invisible below that.
4. **For the SNP/indel decisions L5a–L5c**, measure on the oracle's four tomato2 samples (about 16 min at
   4 threads) and, once one is chosen, on a 16-sample subset of tomato2: passes per start, each start's
   ending (converged / limit / agreed), final log-likelihood against a 1,000-pass reference, and the
   furthest parameter's distance in standard errors. Merge when the SNP/indel fit's passes fall by a
   third or more with the final log-likelihood no lower.
5. **Before D3, a quiet 18-thread run of the STR bench at 2,000 × 30 with H1+H2**, to see how much of
   the 33% idle remains at the full core mix (H3 predicts most of it).

---

## 4. Build / toolchain configuration

`[profile.release]` (`lto = "fat"`, `codegen-units = 1`, `panic = "abort"`) needs no change. Three
methodology findings:

- **[bench] The STR bench does not run the allocator `pvc` ships.** `pvc` installs mimalloc
  (`src/main.rs:20`, feature `alloc-mimalloc`, on by default); `examples/ng_str_fit_bench.rs` installs
  none, so it runs macOS's system allocator. Every malloc/free figure from the bench (about 6% of thread
  samples, and 8.6% of the main thread's time spent in `_xzm_free_main` at 1,000 × 30) and every RSS
  figure describes an allocator production does not use. Fix: give the bench the same
  `#[cfg(feature = "alloc-mimalloc")] #[global_allocator]` line as `src/main.rs`, and the dhat hook
  this review added to its private copy (`tmp/perf_review/src_new/examples/ng_str_fit_bench.rs`,
  four lines at `main`). Cost: six lines, no new dependency.
- **[profile] The quoted 1,062 MiB is resident set size, not memory held.** The same 2,000 × 30 run's
  `/usr/bin/time -l` reports a peak footprint of 643 MB (613 MiB); the remainder is pages the system
  allocator had freed and not yet returned (its `mach_vm_reclaim` frames are in the profile). dhat at
  500 × 30 puts the live peak at 173.6 MB. Use dhat's peak for memory decisions (plan item 2).
- **[profile] The lead's STR profile predates `ae239e26`.** That commit cut the bench's CPU at
  1,000 × 30 from 171–212 s user to 54–58 s (6 threads, two runs each), so `score` is no longer half of all
  samples (it is now 40% at 1,000 × 30 and 15% at 1,000 × 3) and the serial main-thread share of the
  wall is now larger, not smaller.

---

## 5. Code-level findings

### What one E-step does, and how many there are

An **E-step** is one pass of `expected_counts_by_distinct_rows` over one class: build or take the
class's **table of distinct rows** (each (locus, sample) pair's reads pooled into one row of bucket
counts, identical rows held once), **score** every distinct row under all 190 genotypes, **add each
sample** to its row's weights in fixed chunks of rows, then turn each distinct row into **counts**.
That is four rayon bursts, with serial merges on the main thread between them.

Counted with an instrumented copy (`analyse_trace.py`; the round counts are unchanged by every patch
here because the bits are):

| cell | rounds | class E-steps in rounds | rows in them | distinct rows in them | E-steps in the errors |
|---|---|---|---|---|---|
| 1,000 × 30 | 68 | 474 | 34.2 million | 9.2 million | 588 for curvature (14 classes × 42 steps) over 40.8 million rows, and 971 one-locus steps |
| 1,000 × 3 | 249 | 2,337 | 160.9 million | 3.15 million | the same 588 + 971 |
| 2,000 × 30 (lead's run) | 66 | 455 | — | — | the same structure; errors 10.5 s of the 20.3 s fit |

So at depth 30 **the standard errors scan more rows than all the rounds together**, and at depth 3
the rounds dominate and almost every row is a repeat of another (2% distinct).

### Hot-path

- [expected_counts.rs:582](../../../../src/parameter_estimation/joint/ssr_fit/expected_counts.rs#L582), [penalised_fit.rs:1531](../../../../src/parameter_estimation/joint/ssr_fit/penalised_fit.rs#L1531) — **[Hot-path]** H1: the standard errors re-score every distinct row for 36 of their 42 steps although those steps cannot change a score
- **Confidence:** High
- **Hot-path evidence:** lead's 2,000 × 30 profile: `Period::standard_errors` is 47.7% of the main thread's samples (the wall), with 31.7% of it waiting on the `score` burst; fixture `seconds.errors` 10.5 of 20.3 s. Measured here on `ae239e26`, 1,000 × 30, 6 threads: errors 6.44 s → 2.87 s.
- **Pattern matched:** hot_loops — redundant executions (the only kind of change the August review found to pay in this module).
- **Mechanism:** the curvature takes a step up and down in each of the class's 21 or 22 unbounded numbers (`penalised_fit.rs:1531-1547`). A step in F or in an allele-length share changes only each genotype's *prior* (`one_draw`, `two_draws`); the stutter numbers, so `per_allele`, `ln_read_chances` and every row's genotype likelihoods, are unchanged to the bit. Yet each step builds a `RoundModel` and scores every distinct row afresh. The patch scores the class once at its central stutter, keeps only the scaled likelihoods (`ScoredRow { largest, scaled }`), computes each step's two prior-weighted totals from them, and multiplies by the prior where `counts_of` reads it — the same products in the same order. The share and F columns run first and the table is dropped before the three stutter columns score their own, so at most one table is alive. Each curvature entry receives exactly two additions to zero, which commute, so running the columns in a different order changes no bit.
- **Measurement plan (run):** fixture grid, 144 runs below 2,000 samples: `144 exact … 0 fail --require exact` (`tmp/perf_review/grid_new_v123.compare`, with H2 also applied; H1 alone was also exact at `53a3e4f9`, `grid_v1`). Interleaved, two runs each:

  | cell | variant | instructions retired | user s | wall s | `seconds.errors` |
  |---|---|---|---|---|---|
  | 1,000 × 30 | `ae239e26` | 1,253.2 / 1,252.8 billion | 53.7 / 58.1 | 12.1 / 12.7 | 6.44 |
  | 1,000 × 30 | + H1 | 852.9 / 853.0 billion | 38.5 / 39.7 | 8.9 / 8.9 | 2.87 |
  | 1,000 × 3 | `ae239e26` | 525.2 / 524.9 billion | 24.4 / 25.5 | 5.8 / 6.0 | 0.82 |
  | 1,000 × 3 | + H1 | 502.8 / 502.6 billion | 23.8 / 24.0 | 5.6 / 5.7 | 0.60 |

  The gain grows with the share of rows that are distinct, so with depth: largest at deep reads, about 4% at 3 reads.
- **Complexity cost:** moderate and mostly structural. `expected_counts_by_distinct_rows` splits into "score the rows" and "count from scored rows" (two new `pub(super)` functions, `scored_rows` and `expected_counts_from_scored`); `add_sample` takes the row's two totals as an argument; the errors' column loop becomes two loops. It also halves each scored row (see L2). One extra parallel burst per E-step (the totals) that can be folded into the `score` map in the ordinary path.
- **Suggested experiment / fix:** the patch, against `ae239e26`: `tmp/perf_review/patches/v12_errors_reuse_scored_rows.diff` (251 lines including context; the errors part is the second file in it).

- [expected_counts.rs:473-526](../../../../src/parameter_estimation/joint/ssr_fit/expected_counts.rs#L473-L526), [expected_counts.rs:1090](../../../../src/parameter_estimation/joint/ssr_fit/expected_counts.rs#L1090) — **[Hot-path]** H2: building the distinct-row table allocates and SipHashes a fresh vector for every (locus, sample) pair, and the merge copies every row twice
- **Confidence:** High
- **Hot-path evidence:** lead's 2,000 × 30 profile: `DistinctRows::of` 13.5% of the main thread, its serial merge 10.6%, `SipHash` and `RandomState` frames inside it. This review, `ae239e26`, 1,000 × 3: `in_blocks_of` is 22.4% of all thread samples; dhat 500 × 30: `pooled_reads` is 26.9 million of the run's 96.7 million allocations.
- **Pattern matched:** allocations — allocation in a per-row loop; a clone used only as a lookup key; hot_loops — a slow hasher for internal integer keys.
- **Mechanism:** each block calls `pooled_reads`, which allocates a `Vec<u32>` per pair, then hashes it with SipHash; the serial merge on the main thread calls `entry(pooled.clone())` for every row of every block (a copy even when the row is already known) and `reads.push(pooled.clone())` for a new one. The patch pools into one scratch row per block, looks it up by slice, copies only a new row, keys the merge by borrowed slices, and uses `ahash` (already a dependency). Row numbering is still first-met order and no map is iterated, so the table is the same.
- **Measurement plan (run):** exact on the 144-run grid (with H1). Interleaved against H1 alone:

  | cell | H1 | H1 + H2 |
  |---|---|---|
  | 1,000 × 3 instructions | 502.8 / 502.6 billion | **406.4 / 406.1 billion (−19%)** |
  | 1,000 × 3 `seconds.expectation` | 4.69 | 4.04 |
  | 1,000 × 30 instructions | 852.9 / 853.0 billion | 822.5 / 822.4 billion (−3.6%) |

  dhat 500 × 30, blocks allocated over the run: 96.7 million → 49.8 million with H1 + H2.
- **Complexity cost:** small: a 9-line `pool_reads_into`, a lookup-then-insert in place of `entry`, a borrowed key type in the merge. No new dependency.
- **Suggested experiment / fix:** `tmp/perf_review/patches/v3_distinct_rows_build.diff` (73 lines, applies on top of H1).

- [penalised_fit.rs:1223-1245](../../../../src/parameter_estimation/joint/ssr_fit/penalised_fit.rs#L1223-L1245), [expected_counts.rs:787-845](../../../../src/parameter_estimation/joint/ssr_fit/expected_counts.rs#L787-L845) — **[Hot-path]** H3: the threads idle because the main thread works alone between bursts — that is the "idle third"
- **Confidence:** High for the diagnosis; Medium for how much each fix recovers
- **Hot-path evidence:** lead's 2,000 × 30 profile (18 threads): the 18 workers sit in `__psynch_cvwait` for about 31% of their samples, and the main thread is on the CPU for 29% of the wall (5,450 of 18,819 samples) while it waits on a rayon latch for the other 71%. Since workers only run inside bursts the main thread starts, **the workers' idle share is essentially the main thread's serial share**. This review, `ae239e26`, 6 threads, 1,000 × 30: main thread on CPU 22.3% of the wall (1,823 of 8,189 samples), workers idle 31.5%; with H1 + H2, 19.1% (1,148 of 6,014) and 29%. At 1,000 × 3 the workers also spend 7–10% of their samples in `swtch_pri`/`cthread_yield` — rayon spinning between short bursts.
- **Pattern matched:** concurrency — serial sections between fork-join bursts (Amdahl); allocations — cross-thread frees concentrated on one thread.
- **Mechanism:** what the main thread does alone, at 1,000 × 30 on `ae239e26`: **frees** of the per-row vectors the workers allocated (`_xzm_free_main`, `_free`, `free`: 846 samples, 10.3% of the wall) — the `scored` table's two vectors per distinct row and `per_row`'s one; the **distinct-row merge** (about 4%, removed in part by H2); the **merge loops** — the partial weights summed over every chunk × every distinct row and `counts.add` over every row (4.3% after H1); plus the M-step, which is under 1%. None of these depend on the thread count, so at 18 threads, where the bursts are three times shorter, they become the larger share — which is why the lead's 18-thread profile shows 29% where 6 threads show 19–22%. Within-burst imbalance from the 6 performance + 12 efficiency cores accounts for the remaining idle (about 10 points at 6 threads); the 4,096-row floor on a chunk (`expected_counts.rs:803`) leaves a class of 176,000 rows only 43 chunks for 18 threads, but that burst (adding samples) is cheap.
- **Measurement plan:** count, per E-step, the main thread's on-CPU samples in a `sample` profile (`tmp/perf_review/sample_run.sh`, `selftime.py --thread main-thread`); the merge gate is the main thread's on-CPU share falling, with the fixture exact. Then the 18-thread run of plan item 5.
- **Complexity cost / fixes, in order:**
  1. H1 already halves the frees (one vector per distinct row instead of two). Measured: frees 10.3% → 6.6% of the wall.
  2. Store the scored table and the per-row counts as two flat buffers (`distinct × genotypes` and `distinct × (classes + 9)` `f64`s): one allocation per E-step instead of two per distinct row, so the main thread frees two blocks. Bit-identical (the same values, summed in the same order). Small–moderate: index arithmetic in `counts_of` and the `counts.add` loop.
  3. Under mimalloc (what `pvc` runs) the frees may already be cheaper; measure with the bench on mimalloc (§4) before doing 2.
  4. The serial merges cannot be parallelised without changing the order of float additions; at 4% they are not worth an owner decision.

### Likely

- [expected_counts.rs:803-845](../../../../src/parameter_estimation/joint/ssr_fit/expected_counts.rs#L803-L845) — **[Likely]** L1: each chunk's partial weights span every distinct row of the class, so their memory grows with the square of the class's rows
- **Confidence:** High on the arithmetic; the large-class case is not measured (plan item 3)
- **Hot-path evidence:** dhat 500 × 30: `vec![RowWeights::default(); width]` at line 808 holds 7.6 MB (4.4% of the peak). Trace at 1,000 × 30: 25.7 MB beside a scored table of 151.6 MB in the largest E-step (88,000 rows, 48,702 distinct, 22 chunks).
- **Pattern matched:** allocations — buffer sized to the whole problem per parallel task; data_layout.
- **Mechanism:** the chunk is `max(rows / 256, 4,096)` rows, so a class has up to 256 chunks, and each holds a `RowWeights` (24 bytes) for *every* distinct row: 256 × 24 = **6.1 kB per distinct row once a class passes about a million rows — twice the scored row itself** (3.0 kB on `ae239e26`, 1.5 kB after H1). The bench never reaches that (its largest class is 100 loci). A deep cohort does: a class of 5,000 loci at 2,000 samples × 30 reads is 10 million rows, most of them distinct, and its partials alone would be tens of GB. At three reads it does not matter: distinct rows are a few thousand.
- **Measurement plan:** dhat peak at the large-class cell of plan item 3, before and after. Gate: the line-808 site's share of the peak falls to under 1%, fixture exact.
- **Complexity cost:** moderate. Each chunk keeps only the rows it touched: sort its `(index, excess)` pairs stably by index, accumulate runs in row order, return `Vec<(u32, RowWeights)>`; the merge adds those in chunk order. **Bit-identical**: today's merge adds `+0.0` for every untouched row, and adding `+0.0` to a non-negative sum never changes it, so dropping those additions leaves each total's sequence of non-zero additions unchanged.
- **Suggested experiment / fix:** as above; about 25 lines in `counts_from_scored`.

- [expected_counts.rs:1506-1531](../../../../src/parameter_estimation/joint/ssr_fit/expected_counts.rs#L1506-L1531), [expected_counts.rs:582](../../../../src/parameter_estimation/joint/ssr_fit/expected_counts.rs#L582) — **[Likely]** L2: halving the E-step's table does not lower the peak until the pairs' table is halved too
- **Confidence:** High (dhat, deterministic)
- **Hot-path evidence:** dhat at 500 × 30. On `ae239e26` the peak is 173.6 MB, of which the E-step's scored rows are 87.9 MB (28,910 distinct rows × two vectors × 190 genotypes × 8 bytes) and the drawn evidence 64.4 MB. With H1 + H2 the scored rows are half that, but the peak only falls to 169.4 MB, **because the peak moves to F's pairs**: `ordered_row` holds 83.5 MB there (28,910 rows × 19 × 19 × 8 bytes).
- **Pattern matched:** data_layout — a symmetric matrix stored whole.
- **Mechanism:** `ordered_row` stores the 19 × 19 matrix of allele-pair likelihoods with every off-diagonal value twice (`ordered[first * classes + second]` and its mirror, lines 1523-1524); 190 values carry all of it. Storing the triangle (or the 190-genotype vector H1 already keeps) halves it. The two tables are alive at different moments, so the peak is the evidence plus the larger of them: with both halved, about 130 MB at 500 × 30 against 174 MB (64 MB of evidence plus about 44 MB of table and their companions) — a quarter off, estimated from the dhat sites, not measured.
- **Measurement plan:** dhat `At t-gmax` at 500 × 30; gate a fall of 20% or more with the fixture exact.
- **Complexity cost:** moderate: `PairSums::of` (not read in full here) must index a triangle or the genotype list instead of a square. The values are the same, so the sums are too, if `PairSums` keeps its order of additions — check before claiming bits.

- [penalised_fit.rs:1234](../../../../src/parameter_estimation/joint/ssr_fit/penalised_fit.rs#L1234) — **[Likely]** L3: the distinct-row table is rebuilt every round although it cannot change (owner decision: the memory rule)
- **Confidence:** Medium
- **Hot-path evidence:** at 1,000 × 3, 2,337 builds over 160.9 million rows where 14 would do; after H2, `in_blocks_of` is still 19.4% of all thread samples (5,125 of 26,401) on `ae239e26`. At kimura's depth almost every row repeats another (tomato2's 1-base class at 8 repeats: 268,034 pairs, 620 distinct rows), so building is the main per-row cost of a round there.
- **Pattern matched:** hot_loops — recomputing an invariant.
- **Mechanism:** the table depends on the reads and the homozygote excesses, neither of which a round moves (the errors already build it once per class for that reason, `penalised_fit.rs:1502-1504`). Keeping each class's table for the whole fit removes all but one build per class. **What it costs is memory the arch forbids** ("at most one class's distinct-row table exists", arch `str_penalised_fit.md` contract): at kimura the per-pair part is 16 bytes a pair (`(u32, f64)`), 1.06 GB for 66.2 million pairs — 5.5% of the 19.4 GiB evidence — or 265 MB (1.3%) if only the `u32` is kept and the excess is read through the sample index each round. The distinct rows themselves are small at low depth but approach the evidence's size at 30 reads.
- **Measurement plan:** implement behind the existing `reuse_settled_counts`-style config switch; fixture exact; instructions at 1,000 × 3 (gate: −10% or more) and dhat peak at 1,000 × 30 (report the rise).
- **Complexity cost:** small in code (a `Vec<Option<DistinctRows>>` beside `last_counts`), large in policy. **My recommendation: keep only the `u32` indices for all classes (+1.3% of evidence at kimura)** and rebuild the distinct rows' read vectors only when a class's own table is needed — or decide at D3 with the timed number.

- [penalised_fit.rs:1223-1245](../../../../src/parameter_estimation/joint/ssr_fit/penalised_fit.rs#L1223-L1245) — **[Likely]** L4: at low depth, classes too small to fill the pool are still visited one at a time (owner decision: the memory rule)
- **Confidence:** Medium
- **Hot-path evidence:** 1,000 × 3, `ae239e26`: workers idle 31.6%, and 7.3% of all samples are rayon spinning (`swtch_pri`); an E-step there lasts about a millisecond (about 4,000 E-steps in 5 s), so each of its four bursts is a few hundred microseconds.
- **Pattern matched:** concurrency — fork-join overhead on short work.
- **Mechanism:** the one-class-at-a-time rule bounds memory by the largest class's table; at three reads every table is a few MB. Running classes in batches whose summed distinct rows stay under the largest class's (or under a fixed budget) keeps the same bound and gives each burst several classes' work. Each class's E-step is unchanged, so its counts are bit-identical whatever the batch.
- **Measurement plan:** instructions and workers' `__psynch_cvwait` + `swtch_pri` share at 1,000 × 3 and 60 × 3; gate a fall of 10 points in idle with the fixture exact; dhat peak unchanged at 1,000 × 30.
- **Complexity cost:** moderate: a batching function and a `par_iter` over a batch. It must respect the spec's trap that a thread waiting inside one class can pick up another (spec §3): the budget bounds the batch, not the steal, so the batch's sum is the bound.

- [fit.rs:2394-2509](../../../../src/parameter_estimation/joint/fit.rs#L2394-L2509), [fit.rs:766](../../../../src/parameter_estimation/joint/fit.rs#L766) — **[Likely]** L5: the SNP/indel fit — every start runs to the 200-pass limit; the passes, not the per-pass kernel, are where the hour goes (owner decisions)
- **Confidence:** High that the passes are the cost; Medium on how many each option saves
- **Hot-path evidence:** `tmp/perf_ep/t2_59.err`: three starts × 199 passes, each "stopped at the pass limit", 20 min each; a plain pass 4.5–4.8 s, the final pass, which also sums the information, 8.5–9.2 s (1.76–1.95× a plain pass; the spec measured 1.48–1.80× for the judging passes, `fit.rs:696-699`); a cycle 17–22 s. At the end of each start 22–24 of 185 parameters are more than 0.1 standard error from the maximum, the furthest **the invariant share, 3.9–4.2 standard errors away**. Final log-likelihoods −7,985,053, −7,985,051, −7,985,053. The spec recorded the same on kimura: three starts at the limit, 9 h 25 min (`fit_precision.md` §1).
- **Pattern matched:** methodology — the iteration count, not the instruction count, is the lever.
- **Mechanism and options:**
  - **L5a — agreement against a start that hit the limit.** Spec §3.4 (its second "Trap") says a later start may stop on agreement only with a start that *converged*; none does here, so all three run out. Comparing projected endpoints (each start's value plus its Newton distance, both already computed on every judging pass) lets starts 2 and 3 stop once they head to start 1's endpoint. Up to 40 of the 60 minutes on this run. Numbers: the winner is still the best start, so the answer changes only if a stopped start would have ended higher; here the three ended within 2 log-likelihood units.
  - **L5b — take the Newton step the fit already computes.** Each judging pass sums the information and the scores and solves `I⁻¹g` inside the parameters' box (`newton_step`, `fit/standard_errors.rs:583`), then uses it only to judge. Offering `θ + d` as a candidate, kept only if the next pass's log-likelihood is higher (the same guard the SQUAREM jump has, `fit.rs:2441-2467`), turns the slow linear approach along the invariant-share / allele-frequency-shape ridge into a near-quadratic one where the quadratic model holds. It also changes what is reported: today the returned invariant share is 4 of its own errors from the maximum; a converging fit would return the maximum. This is the larger change in code and in numbers.
  - **L5c — judge every k-th cycle.** Once judging starts (about pass 16 here), every cycle pays an information pass, about 20% of the SNP/indel fit's time on this run (one extra 4.5 s in a 17–22 s cycle). For a start that runs to the limit the parameters do not depend on judging, so they would be bit-identical; for one that converges, it stops up to k − 1 cycles later, closer to the maximum. Smallest change, smallest gain; conflicts with L5b, which wants the information every cycle.
- **Measurement plan:** plan item 4.
- **Complexity cost:** L5a small (the comparison exists; the yardstick changes from "converged" to "best so far, projected"); L5c one counter; L5b moderate (a candidate step, its guard, and its interplay with SQUAREM's `reach` schedule). **Recommendation: L5a, then L5b measured on the oracle; L5c only if neither is taken.**

### Speculative

- [penalised_fit.rs:1550-1574](../../../../src/parameter_estimation/joint/ssr_fit/penalised_fit.rs#L1550-L1574) — **[Speculative]** S1: the errors' per-locus slopes — each locus's slope clones the locus's reads and builds a table and a model from scratch
- **Confidence:** Low
- **Hot-path evidence:** `standard_errors::{closure}` (the per-locus `par_iter`) 2.3% of thread samples at 1,000 × 30 on `ae239e26`; 971 E-steps in the trace.
- **Pattern matched:** allocations — a clone to build a temporary owner; recomputing an invariant (`RoundModel::of`, with 4,750 logs, once per locus).
- **Mechanism / fix:** the per-locus slope at the class's central numbers could read the class's own scored table; to stay bit-identical it must reproduce the one-locus path's order of additions (a locus's rows first-met, `per_row` then `counts.add`). The `tract.clone()` alone can go by giving `expected_counts` a tract slice.
- **Measurement plan / cost:** instructions at 1,000 × 30; not worth more than a borrowed slice unless the share grows at D3's scale.

- [fit.rs:1098-1176](../../../../src/parameter_estimation/joint/fit.rs#L1098-L1176), [fit.rs:1338](../../../../src/parameter_estimation/joint/fit.rs#L1338) — **[Speculative]** S2: the SNP/indel fit's evidence is decoded afresh on each of the 597 passes; the decoding is a fifth of the SNP/indel fit's CPU
- **Confidence:** Low on any fix
- **Hot-path evidence:** `tmp/perf_ep/sample_snp.txt`, self time: `one_position` 387,945 samples (66.0% of all), `expectation_pass` 108,628 (18.5%; the inlined `EvidenceCursor::next_position` and `fill_depth_weights`), idle wait 11.5%. Line-level attribution is not available: 270,000 samples are under unsymbolicated frames (`???`, `0x0`).
- **Pattern matched:** hot_loops — recomputing what the parameters do not touch.
- **Mechanism:** `next_position` re-derives, for every library at every position on every pass, the depth range and the Poisson weights of its range (`fill_depth_weights`, up to 32 divisions); none depends on the parameters. Caching decoded positions is not viable (gigabytes at 2 million positions × 59 libraries), but the weights depend only on (library, fewest reference reads, spread) and could come from a per-library table built once — the same loop, so the same bits.
- **Measurement plan:** first a line-attributed profile of `benches/ng_joint_fit_perf.rs` `generic_by_samples/32` under `--profile profiling`; only if `fill_depth_weights` holds more than a third of the cursor's time, build the table and gate on instructions at −5%.
- **Complexity cost:** small–moderate; the August review's lesson applies (state added to `one_position`'s body measured slower); this sits outside it, in the cursor.

### Note

- **The SNP/indel fit's pool is well used.** In the 59-sample profile the workers wait 6–7% of their samples (11.5% counting the parked main thread); the fixed join tree over 128 chunks (`fit.rs:3237`) suits the core mix. Changing the chunk count moves bits and is not worth it.
- **The 3.95 GB peak of the 59-sample run is not attributed.** The SNP/indel evidence is 120 MiB, contamination's markers 18 MiB, the STR evidence 101 MiB by the run's own log. Before any memory work on `estimate-parameters`, a dhat run of the command on a small cohort (four tomato2 samples) is needed to say where the rest is.
- **H1's extra burst.** In the ordinary E-step the totals can be computed inside the `score` map (`unzip` into rows and totals), so H1 adds no burst there; worth folding in when H1 is applied.

---

## 6. Out-of-scope observations

- `pair_cubics_by_tract` (`expected_counts.rs:1538`) builds its own distinct-row table each time it is called — once per class in each pair alternation (28 builds in the 1,000-sample runs, 2 alternations × 14 classes) and once per class in F's errors; with L3 it would take the cached one.
- `gather_strata`'s evidence (19.4 GiB at kimura, 315 bytes a pair) is, after L1–L2, the dominant memory at every depth; its layout (`Vec<SampleTractReads>` each holding a `Vec<(LibraryIndex, Vec<u32>)>`, two heap blocks per pair) is the place a memory review of D3 should start.

---

## 7. What's already good

- **Determinism by construction**: the add-samples burst uses chunks fixed by the row count and merges them in order (`expected_counts.rs:801-827`), so none of this review's patches had to argue about thread counts — they only had to keep per-row arithmetic.
- **The errors already reuse what they can across steps** (one distinct-row table per class, `penalised_fit.rs:1502-1504`); H1 extends the same idea one level down, to the scores.
- **`ae239e26` took the logs out of the row loop the right way**: a per-round table read in the same order, bit-identical, with the identity check stated in the commit.

---

## Author response convention

Address each finding by its identifier (H1–H3, L1–L5, S1–S2) with one of: `applied in <commit>` /
`experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` / `won't fix because …`.

## Scratch left on disk (the review cannot delete)

`tmp/perf_review/`: `src_copy/` and `src_base/` and `src_new/` (source copies; `src_new` holds H1+H2 on
`ae239e26`), `target/`, `target_new/`, `target_dhat/` (builds, several GB), `bin/` (the measured
binaries), `runs/`, `prof/`, `dhat_base/`, `dhat_v123/`, `grid_v1/`, `grid_v3/`, `grid_new_v123/`
(fixture runs and comparisons), `patches/` (the diffs), `orig/`, scripts `analyse_trace.py`,
`dhat_top.py`, `selftime.py`, `sample_run.sh`, and a stray `pf_copy_check.rs`.
