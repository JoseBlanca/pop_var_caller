# Performance Review: snp-fit
**Date:** 2026-10-18
**Reviewer:** rust-performance-review skill (orchestrator), five category sub-agents, a fresh reviewer for each change
**Scope:** `pvc estimate-parameters`' SNP/indel parameter fit — the passes over the evidence, the information they sum, the contamination fit they call
**Verdict:** Apply the listed wins — applied: main at `d1b6e54f` (ten commits), output bit-identical
**Hot-path evidence:** two sampling profiles of the 59-sample fit taken in this review on the current code, one of them with each callee split out so its share can be read; per-position counts from an instrumented copy; DHAT and a counting allocator for memory; instructions retired and wall time, interleaved

---

## 1. Scope and constraints

**What was reviewed.** The SNP/indel half of step 4's parameter pre-pass: `fit_jointly`, `maximise`,
`expectation_pass`, `one_position` ([fit.rs](../../../../src/parameter_estimation/joint/fit.rs)), the
information scoring ([fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs)),
the Newton step and its box solve ([fit/standard_errors.rs](../../../../src/parameter_estimation/joint/fit/standard_errors.rs)),
and the contamination fit `fit_jointly` calls ([contamination.rs](../../../../src/parameter_estimation/joint/contamination.rs)).
The repeat-tract fit and calling are out of scope.

**Reviewed against** `50115d40` (main on 2026-10-18), on branch `snp-fit-perf`; merged into main as
`d1b6e54f`.

**Targets.** Wall time and peak memory, with every fitted number bit-identical — or moved, measured,
and put to the owner. The fit must work from one sample to several thousand and from three reads a
position to several hundred (`CLAUDE.md`); the cohorts measured here are tomato2 at about 5.5 reads a
position.

**Hardware.** Apple M5 Pro, 18 cores: **6 performance and 12 efficiency**. A thread count above 6
puts work on cores roughly half as fast, so a pass at 18 threads waits for its slowest core. The
machine was shared with other sessions all day; §3 says how each number was protected from that.

**Inputs.**

| name | what it is | positions | passes |
|---|---|---|---|
| 1 sample | SRS1839219, one of the oracle's four | about 2 million | 22, 20, 21 |
| 4 samples | the oracle's four tomato2 accessions, from stored psps | about 2 million | 53, 46, 37 (137 in all) |
| 59 samples | tomato2's 59 accessions | 1,998,212 | 45, 43, 43 (131 in all) |

**Categories dispatched**, each to its own sub-agent in its own worktree: `methodology` (a bench that
reaches the judging passes; build settings), `hot_loops` twice (the information scoring; `one_position`),
`concurrency` (idle workers; what grows with the cohort), `allocations` with `data_layout` (peak memory).
`float_portability` was applied by every agent to every change. `io_and_syscalls` does not apply: the
evidence is read once, before the passes.

---

## 2. Verdict

**Apply the listed wins — done.** Ten commits are on main at `d1b6e54f`. Every one leaves the
parameters file byte for byte as it was, checked on 4 samples at 2, 4, 6 and 18 threads and on 59
samples at 6 and 18 threads, and each was read by a reviewer that had not written it.

| | main before (`50115d40`) | main after (`d1b6e54f`) | change |
|---|---|---|---|
| 59 samples, 18 threads, wall | 1,005 s and 915 s | 917 s and 873 s | **−6.8%** |
| 59 samples, CPU (user) | 15,286 s and 14,592 s | 14,181 s and 13,948 s | −6.0% |
| 59 samples, instructions retired | 245.25 trillion | 225.95 trillion | −7.9% |
| 4 samples, 6 threads, wall | 182 s and 184 s | 161 s and 172 s | −9.0% |
| 4 samples, instructions retired | 29.32 trillion | 25.72 trillion | −12.3% |
| 1 sample, 18 threads, wall | 26.0 s and 26.8 s | 22.7 s and 22.2 s | −15% |
| 4 samples, peak footprint | 152–172 MB | 95–111 MB | about −40% |
| 59 samples, peak footprint | 376–385 MB | 374–387 MB | none measurable |

**Where the gain came from.** Today's profile changed the picture the last review worked from: the
information the Newton step and the standard errors are computed from is now summed on nearly every
pass, and scoring it was 41% of the workers' CPU at 59 samples. Three changes removed work done twice
— in that scoring, and in `one_position` — without moving a bit (§5, H1–H3).

**Where memory went.** At 1 to 4 samples the run's peak was not in the SNP/indel fit at all but in the
contamination fit it calls, which built 58 to 74 MB of arrays for a cohort too small to have a single
marker. That is gone (H4). At 59 samples the peak is elsewhere and did not move.

**Two decisions for the owner.** Cutting the census into 512 chunks instead of 128 (§5, L3): 5.2% less
wall time at 59 samples, the output moved once by at most 4 × 10⁻⁸ of a standard error — recommended,
adopted by the owner, and applied after this report (the oracle's calls do not move). Profile-guided optimisation (§4): 3.4% fewer instructions for a build flow three times as long — not
recommended now.

---

## 3. Measurement plan

What was measured, and how each number is protected from a shared machine:

1. **Profile first.** `sample` on the 59-sample fit for 120 s mid-fit, and again for 90 s on a build
   with `#[inline(never)]` on each of the expectation pass's callees, so that each is a separate symbol
   (`sample` merges addresses, so it cannot attribute lines; the split build can attribute functions).
   Raw files: `prof/sample_judging.txt`, `prof/sample_attr.txt` under
   `tmp/perf_snp-fit_2026-10-18/` in the main checkout (ignored by git).
2. **Gate on instructions retired, not wall time.** `/usr/bin/time -l` reports instructions retired;
   for one binary and one input it repeats to within 0.02%, whatever else the machine is doing. Every
   candidate was gated on it, A B B A, at 2 threads (the owner's rule after the first hour: five agents
   on 6 fast cores made every wall number meaningless). Wall time was then taken once, by the lead,
   alone.
3. **Bit identity.** `cmp` of the parameters file against the baseline's at every step;
   `cli::cross_platform_digests` and `parameter_estimation::joint` (478 tests) at every commit.
4. **Wall time, interleaved, with the load logged.** A queue that starts each run only when the one-minute
   load is below 3. The first 59-sample round ended at loads of 21 to 35 — a Python simulation from
   another session held one core throughout — and showed the new code *slower* (1,280 s against
   1,114 s), which contradicted its 7.9% fewer instructions. It was rerun main, new, new, main with the
   load logged every 30 seconds (`timing2/load.log` there): load 18 to 20 during every run, which
   the fit's own 18 threads account for, except two minutes at 38 to 53 in the first main run. That
   rerun is the 59-sample row of §2.
5. **A bench that sees the hot path** (§4) for the next round.

---

## 4. Build / toolchain configuration

`[profile.release]` (`lto = "fat"`, `codegen-units = 1`, `panic = "abort"`) and `.cargo/config.toml`'s
`-C target-cpu=apple-m1` need no change.

- **[bench] The bench never reached a judging pass, so it never timed 41% of the fit — fixed in
  `d398f0cb`.** The older groups of `benches/ng_joint_fit_perf.rs` run eight passes from one start and
  stop before the fit starts judging. The new group `ng_joint_fit/generic_judging` fits drawn cohorts
  of 1, 8 and 32 samples with judging switched on as early as `maximise` allows; 8 of its 11 passes sum
  the information, which the setup asserts from the fit's own counter (`StartRecord::judging_passes`).
  Its drawn data match tomato2's mix of candidate alleles at 32 samples (67.9 / 26.7 / 5.4% of
  positions with 1 / 2 / 3 candidates, against 67.5 / 27.6 / 4.9%). On the earlier draw, the 32-sample
  case split its CPU between the expectation and the scoring as the real 59-sample fit does, to within
  4 points; the 1- and 8-sample cases spend more of their time between passes than a real fit does and
  are regression checks, not stand-ins. A whole run of the group: 38 s at 2 threads.
  `NG_JOINT_FIT_BENCH_THREADS=2 cargo bench --features bench-fixtures --bench ng_joint_fit_perf -- generic_judging`.
- **[build] Profile-guided optimisation: 3.4% fewer instructions, output bit-identical — owner's
  decision; recommendation: not now.** Trained on the 1-sample fit, measured on the 4-sample fit at 2
  threads: 29,314.6 billion instructions against 28,317.8 billion, two runs each; the parameters file
  identical. It costs a build flow about three times as long (an instrumented build, a 7-minute training
  run, a second build), a profile that goes stale whenever `fit.rs` changes, and an `llvm-profdata`
  matching rustc's LLVM (`rustup component add llvm-tools`). 3.4% of instructions is about half of what
  this review's code changes bought, for a permanent cost to every release; worth revisiting when the
  release process is scripted.
- **[build] `-C target-cpu=native`: no gain.** −0.002% instructions, inside the 0.014% between two
  baseline runs; LLVM reads this chip as `apple-m4`, which adds instructions (`bf16`, `i8mm`) that no
  scalar `f64` code uses, and a `native` binary would not run on an M1. Keep `apple-m1`.

---

## 5. Code-level findings

### What the profile showed

59 samples, 18 threads, mid-fit, on the split build. Shares of the 18 workers' samples (822,000; the
main thread is excluded, because it waits inside each pass for the whole of it):

| what | share |
|---|---|
| `one_position` — each position's likelihood, posteriors and counts | 47.5% |
| scoring each position's slopes for the information (`score_the_segregating_branch` 22.9%, `score_the_duplicated_branch` 14.1%, `add_position` 2.3%, `fill_read_slopes` 1.4%) | 41% |
| workers idle, waiting for the pass's last chunks | 6.6% |
| decoding the evidence (`next_position` 1.15%, `fill_depth_weights` 0.09%) | 1.2% |

The main thread spends all of its samples waiting for the pass's workers, so the work done between
passes is below the profile's resolution at 59 samples.

**Why the scoring is so large now.** Since the Newton step (2026-10-16), judging starts by about pass 13
and every Newton try sums the information too, so nearly every pass pays for the scoring, where before
only the last did.

**The 2026-10-15 review's S2 is refuted.** It put decoding the evidence at about a fifth of the fit's
CPU, from a profile in which the scoring was inlined into the same symbol as the decoding. Split out,
decoding is 1.2%.

### Hot-path

- [fit/information.rs:557](../../../../src/parameter_estimation/joint/fit/information.rs#L557), [fit/information.rs:676](../../../../src/parameter_estimation/joint/fit/information.rs#L676) — **[Hot-path]** H1: the scoring re-read each sample's values at every node — applied in `1d0aaea4`
- **Confidence:** High
- **Hot-path evidence:** the two functions are 22.9% and 14.1% of the workers' samples (table above).
- **Pattern matched:** hot_loops — hoist what the inner loop does not change.
- **Mechanism:** each branch was walked node by node and, within a node, sample by sample, so a
  sample's three likelihoods, its libraries' read slopes, its row of scores and its homozygote excess
  were fetched again at each of the rule's 16 nodes. The walk is now sample by sample and node by node
  within a sample; each node's weight is formed once a candidate into a small scratch vector. Every
  accumulator still receives the same terms in the same order — a sample's row slot over the nodes in
  node order, each node's slope in the frequency over the samples in sample order — so the sums are the
  same bits.
- **Measured:** alone on `50115d40`, 4 samples, 2 threads, A B B A: 29,316.5 / 29,311.5 billion
  instructions against 27,894.0 / 27,893.3 billion (**−4.85%**); a drawn 40-sample cohort, where the
  information is summed in blocks as at 59 samples: −6.9%.
- **Complexity cost:** about 70 lines; three scratch vectors on `PositionScores`; each function
  destructures its output to split the borrows. A `debug_assert` checks the rows and the scratch hold
  the same samples.

- [fit.rs:3814](../../../../src/parameter_estimation/joint/fit.rs#L3814) — **[Hot-path]** H2: the attribution re-summed the branch totals the likelihood had just summed — applied in `7c8bab29`
- **Confidence:** High
- **Hot-path evidence:** `ln_sum_exp` and libm's `log` together were about a fifth of `one_position`'s
  busy samples at 4 samples (8,270 + 7,863 against 41,708).
- **Pattern matched:** hot_loops — redundant executions (the only kind of change the August review found
  to pay in this function).
- **Mechanism:** per noise class, the likelihood block sums each branch over its candidates and nodes
  (`ln_sum_exp`, one `exp` a term) and then the four branches. The attribution then rebuilt the same
  buffers and summed them again, and `score_the_shares` summed the four branches a third time: at one
  candidate, about 74 `exp` and 8 `ln` repeated at every position. The likelihood now keeps the four
  totals in the scratch, and the attribution forms each share's exponent from them in the same order.
- **Measured:** −5.30% of instructions on top of H3 (4 samples, 4 threads); H2 and H3 together −8.06%
  (2 threads, A B B A).
- **Complexity cost:** four `[f64; 2]` fields in the scratch; two buffer rebuilds and five sums deleted
  (+44 −39 lines).

- [fit.rs:3440](../../../../src/parameter_estimation/joint/fit.rs#L3440) — **[Hot-path]** H3: logarithms that depend only on the parameters were taken at every position — applied in `70611a73`
- **Confidence:** High
- **Hot-path evidence:** libm's `log`, 7,863 self samples at 4 samples.
- **Pattern matched:** hot_loops — hoist invariants.
- **Mechanism:** the branch shares, the class shares, `ln 3` and the log of how many alleles each
  candidate stands for were taken afresh at every position — 30 to 50 logarithms a position — and the
  scorer took some of them again. `BranchShares`, which the scorer already built once a pass, moved into
  `PassModel`, which both read. A simplification as well as a saving (+80 −83 lines).
- **Measured:** −2.90% of instructions (4 samples, 4 threads, two runs each side).

- [contamination.rs:755](../../../../src/parameter_estimation/joint/contamination.rs#L755) — **[Hot-path]** H4: a cohort too small to have a marker built the census-length arrays anyway — applied in `ebbd6a2f`
- **Confidence:** High
- **Hot-path evidence:** DHAT at 4 samples: the run's peak, 125.4 MB, is in `fit_contamination_over`,
  73.9 MB of it the marker search's arrays (37 bytes a position); the log says "0 markers of 1998212
  positions".
- **Pattern matched:** allocations — don't build a container you can show is useless.
- **Mechanism:** a marker needs reads from at least 8 samples (`MIN_SAMPLES_WITH_DATA`), and a sample's
  libraries are pooled before samples are counted, so a cohort of 2 to 7 samples can never have one.
  The search now returns no markers before allocating, with the same refusal and log line.
- **Measured (DHAT, 4 samples):** peak 125,402,015 → 71,206,145 bytes; the run's peak moves to the
  repeat-tract fit. Footprint 152–172 MB → 95–111 MB in the timing runs.
- **Complexity cost:** four lines.

### Likely

- [contamination.rs:770](../../../../src/parameter_estimation/joint/contamination.rs#L770) — **[Likely]** L1: the marker search held its choosing totals beside the marker index — applied in `1a4d8b92`
- **Confidence:** High on the mechanism; the saving is arithmetic, because no cohort small enough to
  run under DHAT here reaches 8 samples, and at 59 samples the run's peak is elsewhere.
- **Mechanism:** the per-position count of samples with reads is a `u32` (it counts samples), and the
  three choosing totals are dropped once the markers are chosen, before the position-to-marker index is
  built: 37 → 29 bytes a position at the peak, 76 → 58 MB at two million positions, for any cohort of 8
  or more. The same integer tests on the same totals; the 59-sample parameters file is byte-identical.
- **Complexity cost:** small.

- [fit.rs:3721](../../../../src/parameter_estimation/joint/fit.rs#L3721) — **[Likely]** L2: the reported pass held all 128 chunks' totals before adding any — applied in `9092ad1e`
- **Confidence:** High
- **Mechanism:** each start's last pass, which also keeps every position's posteriors, collected every
  chunk's totals and then added them in order — a second copy of the per-position lists (8 MB at two
  million positions) and every chunk's sums, which grow with the cohort by about a kilobyte a sample.
  It now computes a window of twice the pool's width in parallel and adds each window's chunks in chunk
  order: the same additions, so the same bits at any pool width.
- **Measured:** the SNP/indel fit's own peak at 4 samples 60.3 → 51.8 MB (counting allocator); by
  arithmetic about 0.2 GB → 0.06 GB at 2,000 samples. **The cost:** each window ends with a wait, and
  the last pass took 9.35 s against 9.00 s at 18 threads, 16.45 s against 15.88 s at 6 (59 samples, one
  run each) — about 4% of a pass that runs three times in a fit of about 130.
- **Complexity cost:** six lines.

- [fit.rs:3758](../../../../src/parameter_estimation/joint/fit.rs#L3758) — **[Likely]** L3: the workers' idle time is the tail of each pass; 512 chunks shorten it and move the output once — **adopted by the owner, applied after this report was filed** (the commit that changes `CHUNKS_PER_CENSUS`)
- **Confidence:** High on the mechanism; see the measurement for the gain.
- **Hot-path evidence:** workers idle 6.6% of their samples at 59 samples, 18 threads. Per-chunk
  timestamps on 4 samples (18 threads, 137 passes): 5.77 of the 6.36 points of idle are the tail, from a
  thread's last chunk to the pass's end; 0.45 the head; 0.13 the gaps between one thread's chunks.
- **Mechanism:** the 128 chunks cost the same to within 2% of their mean, and which chunk is slow
  changes from pass to pass (the correlation between consecutive passes is −0.016). At 18 threads 128
  chunks make 7.1 rounds, so most threads finish their seventh and wait while two or three run an eighth.
  No schedule that keeps the bits can fix that: a chunk's total is a running sum over its positions, so
  it cannot be split, and starting chunks longest-first does not help when the lengths are equal and
  unpredictable (simulated over the measured times: 8.4% idle against 8.6%). Only more chunks shorten
  the tail: a model of equal chunks with ±8% noise gives 8.8% idle at 128 chunks and 2.6% at 512 on this
  machine's mix of cores.
- **What it costs in bits:** the sums are split at different positions, so the output moves once — the
  same at every pool width (checked at 2 and 18 threads). On 4 samples the SNP/indel numbers move by at
  most 1.45 × 10⁻⁷ of their own standard error, with the same 137 passes. The repeat-tract fit, fitted
  afterwards from those numbers, moves too: 99 lines, one of them by 7% (§6).
- **Measured:** instructions unchanged (−0.009% at 2 threads, inside the noise). Wall, 59 samples, 18
  threads, 128 then 512 then 512 then 128 chunks, the load logged every 30 s and never above 24:
  880.9 s, 827.8 s, 829.6 s, 866.9 s — **−5.2%** (828.7 s against 873.9 s), with the same 225.96 trillion
  instructions and the same passes (start 1 converged after 45). At 59 samples the output moves less
  than at 4: 329 of 6,509 lines change; the largest relative change of any number is 4.5 × 10⁻⁶ (a
  repeat-tract standard error), the largest move of a fitted value 3.9 × 10⁻⁸ of its own standard error.
- **Complexity cost:** one constant, and its comment re-measured. The oracle's baseline and the
  cross-platform digests would be re-recorded once, with the move stated.
- **Recommendation: adopt it.** A 5% shorter fit on the hardware this caller is developed on, for a
  one-constant change whose effect on the output is four orders of magnitude below anything a standard
  error can see. It costs nothing measurable in instructions, and nothing in memory now that the
  reported pass holds a window rather than every chunk (L2). Measured on this machine only.

### Rejected, with their measurements

| candidate | what it would have done | measured | why rejected |
|---|---|---|---|
| read the genotype prior's slope from the pass's table | the scoring calls `genotype_frequencies_slope_in_frequency` where `quadrature.prior_slopes` holds the same values | −0.4% alone; **+0.15%** on top of H1 (4 samples), +0.3% (40 drawn) | with the sample loop outermost the function is cheaper than three strided table loads |
| hand the scoring each node's slope in the frequency from `one_position` | `one_position` already forms it at the nodes it does not skip | **+2.8%** on top of H1 (40 drawn) | two fills a position and a branch in the innermost loop cost more than the division saved; adds state to `one_position` |
| reuse a sample's read likelihoods across candidates it has no read on | the values are identical | counted: about 0.02% of instructions at 4 samples, about 0.5% at 59 (estimated) | not worth the comparison logic |
| take `tally_of`'s inputs out of the candidate loop | expected reference reads do not depend on the candidate | at most 0.3% at 59 (estimated) | small at 5 reads a position; **worth re-measuring on deep data**: at a few hundred reads a stored depth stands for up to 22 depths and the call loops over them |
| cache whole results for positions where every read is on the reference | 95% of positions at 4 samples | — | at 59 samples the libraries' depth codes almost never repeat, so it does not hold across the range of cohorts |
| start chunks in another order | shorten the tail | simulated: no gain | §L3 |
| add the iterating passes' chunks as they finish | bound the partial totals held | — | changes the order of additions, so the bits; about 0.1 GB at 2,000 samples |
| flatten the per-sample vectors in `Statistics` and `InformationSums` | fewer allocations a chunk | about 6 × samples allocations a chunk; nothing per position | no time at the profile's resolution |

### Note

- **Nothing is allocated per position.** 7,200 allocations a pass at 4 samples against 2 million
  positions; `Statistics::new` 129 times a pass at 128 chunks (once a chunk, plus the total). The memsets the profile shows (0.8% of the workers) are
  per-position resets of buffers the next position needs.
- **What grows with the cohort, measured on drawn cohorts with judging from the start.** The work done
  between passes on the main thread — the maximisation (about 35 µs a sample), the standard errors, the
  Newton step's box solve, the quadrature tables — is 75 ms a pass at 2,000 samples, against a pass of
  about 240 s at 18 threads (extrapolated linearly in positions × samples from 7 s a pass at 59). The
  box solve's known limit (`fit_precision.md` §2) **does fire on realistic data**: an outbred drawn
  cohort of 2,000 samples, half of whose homozygote excesses sit at 0, fixed 940 to 995 parameters at a
  bound in 86 to 137 rounds, 375 ms a judging pass (at most 469 ms) — about 0.2% of a pass. An inbred
  one took 8 ms. Not worth replacing the solve for speed.
- **Memory at 2,000 samples is set elsewhere** (arithmetic from measured per-unit sizes): the census
  evidence about 4.2 GB (2.1 MiB a read group), contamination's per-marker vectors about 2.5 GB, and the
  repeat-tract fit's evidence and table, which `estimation_memory.md` §8 puts at about 23 GB. The
  SNP/indel fit's own structures are about 0.2 GB.

---

## 6. Out-of-scope observations

- **A repeat-tract number moves 7% when its SNP/indel inputs move by one part in 10⁸.** With 512
  chunks (L3) one tract stratum's `between_locus_share` (period 3, 30 repeats) went from 0.002217 to
  0.002372; the file gives it no standard error. Either that number is poorly determined and the file
  should say so, or the repeat-tract fit is sensitive where it should not be. Worth a look in the
  repeat-tract fit's own review.
- **The kept-position list** ([loci.rs:351](../../../../src/parameter_estimation/joint/loci.rs#L351)) holds
  16 bytes a position, 33.5 MB at two million positions, for the whole `estimate-parameters` run; its
  only reader in this command is a digest. At 4 samples it is now 47% of the run's peak.

---

## 7. What's already good

- **The fixed join tree** ([fit.rs:3587](../../../../src/parameter_estimation/joint/fit.rs#L3587)): chunks
  sized from the census and added in an order fixed by their number made every change here a question
  about per-position arithmetic only — none had to argue about the pool.
- **One scratch per chunk** (`Scratch`, `PositionScores`): two million positions a pass, no allocation per
  position, which is why the allocation categories found nothing on the hot path.
- **The pass model's per-pass tables** (`ReadLogs`, the quadrature's priors): H3 extended the same idea to
  the branch shares rather than inventing a new one.

---

## Author response convention

Address each finding by its identifier with one of: `applied in <commit>` / `experiment shows no gain —
closing` / `disputed because …` / `deferred to <issue>` / `won't fix because …`. H1–H4, L1 and L2 are
applied; L3 was adopted by the owner and applied afterwards; profile-guided optimisation waits on the
owner.

## Where the evidence is

In the main checkout, under `tmp/perf_snp-fit_2026-10-18/` (ignored by git, so on that machine only):
`prof/` (the profiles), `runs/`, `timing/`, `timing2/` and `timing3/` (every timed run's parameters file,
log, `/usr/bin/time -l` output and load; `timing3/` is the 128- against 512-chunk comparison),
`oracle_main_fit/` and `oracle_main5490_fit/` (main's fitted files the 512-chunk move was measured
against), `brief_common.md` (what every sub-agent was told), and the scripts `fit.sh`, `queue.sh`,
`lines.py`, `compare_move.py`. The sub-agents' own worktrees, and the one this review worked in, were
deleted once the work was merged.
