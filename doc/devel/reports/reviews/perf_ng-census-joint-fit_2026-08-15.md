# Performance Review: ng-census-joint-fit
**Date:** 2026-08-15
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** the joint parameters fit and the census — `src/ng/parameter_estimation/joint/` plus `generic/depth_bins.rs`, 12,801 lines never before measured
**Verdict:** Apply the listed wins — **seven were applied and measured during the review, taking 70% off the repeat-tract fit** (§2a); one design decision remains, and it is the owner's (§2)
**Hot-path evidence:** a sampling profile and six timed runs, all taken during this review; none existed before it

---

## 1. Scope and constraints

**What was reviewed.** One module: the joint route of step 4's parameter pre-pass, which estimates
error rates, heterozygosity, inbreeding, contamination and repeat-tract slippage once for a whole
cohort, before any variant is called. It is not the variant caller.

**Reviewed against** commit `866a46b3`, branch `ng-census-encoding`, worktree
`/Users/jose/devel/pop_var_caller-ng-census-encoding`.

**Targets.** Per [CLAUDE.md](../../../../CLAUDE.md) §0 the fit must degrade gracefully from **one
sample to several thousand** and from **three reads a position to several hundred**. The census
file exists for a memory bound: peak resident is meant to be *samples × the largest section*, not
the sum of the file. Target hardware is this macOS arm64 host (18 cores, 64 GB) and the Debian
aarch64 dev container that mirrors the production target.

**Hot-path evidence available.** A `sample(1)` profile of one run and six timed runs, all taken for
this review and recorded verbatim in
[tmp/perf_review_2026-08-15_ng-census-joint-fit/_evidence.md](../../../../tmp/perf_review_2026-08-15_ng-census-joint-fit/_evidence.md).
Before this review there was no benchmark, no profile and no phase attribution of any kind for this
module; `benches/` holds nine benches and none names `parameter_estimation`.

**In-scope files.**

| file | lines |
|---|---|
| [src/ng/parameter_estimation/joint/census.rs](../../../../src/ng/parameter_estimation/joint/census.rs) | 3,726 |
| [src/ng/parameter_estimation/joint/fit.rs](../../../../src/ng/parameter_estimation/joint/fit.rs) | 3,302 |
| [src/ng/parameter_estimation/joint/ssr_fit.rs](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs) | 1,618 |
| [src/ng/parameter_estimation/joint/census_file.rs](../../../../src/ng/parameter_estimation/joint/census_file.rs) | 1,459 |
| [src/ng/parameter_estimation/joint/contamination.rs](../../../../src/ng/parameter_estimation/joint/contamination.rs) | 1,362 |
| [src/ng/parameter_estimation/joint/loci.rs](../../../../src/ng/parameter_estimation/joint/loci.rs) | 1,304 |
| [src/ng/parameter_estimation/generic/depth_bins.rs](../../../../src/ng/parameter_estimation/generic/depth_bins.rs) | 763 |
| [examples/ng_joint_records_walk.rs](../../../../examples/ng_joint_records_walk.rs) — as the measurement vehicle only | 1,631 |

**Deliberately out of scope.** `src/psp/`, `src/pileup/`, `src/ssr/`, `src/var_calling/`,
`src/vcf/` — production, frozen, ng does not edit them. The other `examples/ng_joint_*` harnesses
and `examples/ng_depth_term_family.rs` — one-off measurement programs, not code that ships.
Building a census from an existing pileup — unbuilt and blocked, since ng writes no pileup file.
Two gates are red and both pre-date this branch's first commit, and neither was chased: the
aggregate `cargo clippy --all-targets --all-features -- -D warnings`, and `cargo test
--all-targets`, which panics in `benches/psp_writer_perf.rs:386`.

**Categories dispatched.** Six, at most three at a time, all read-and-reason with every timed run
serialised through the orchestrator — two agents timing at once produce numbers that are not
comparable, and a wrong number here looks exactly like a right one.

| category | why |
|---|---|
| `methodology` | always; and nothing here has ever been benchmarked |
| `hot_loops` | the profile puts 44% of samples in two inlined numeric closures |
| `allocations` | nested `Vec` shapes throughout, and a memory bound to check |
| `data_layout` | the innermost loops walk `Vec<Vec<f64>>` |
| `concurrency` | rayon in three files, and an identity oracle that pins its CPU count |
| `io_and_syscalls` | the census file, its directory and its seeking reader |

**Per-category findings**, left in place as an audit trail, are in
[tmp/perf_review_2026-08-15_ng-census-joint-fit/](../../../../tmp/perf_review_2026-08-15_ng-census-joint-fit/).

---

## 2. Verdict

**Apply the listed wins** — eight of them are named in a profile, are value-preserving or nearly
so, and cost between one line and one small table each. But they are the second thing to do, not
the first, because of what the measurements say about size.

**The repeat-tract fit takes about six days on the 63-accession tomato cohort, and micro-wins
cannot rescue that.** The cost is linear in one quantity — how many tracts a fit reads, summed over
the distinct pooled sets `fit_strata` actually fits — at about **0.14 seconds a tract at 8 samples
on four threads**, measured over a 26-fold range (§3). Tomato's selection puts that quantity near
530,000, and the per-sample term takes it to roughly 0.95 s a tract at 63 samples. Every hot-loop
finding in §5 put together is worth a factor of two or three. The gap is a factor of a thousand.

**The design decision, and it is the owner's.** Two knobs already exist and both are set to
"unbounded" on the only path that runs this:

- **The per-stratum cap.** `SelectionTerms::ssr_cap` bounds how many tracts a stratum keeps.
  [examples/ng_joint_records_walk.rs:167](../../../../examples/ng_joint_records_walk.rs#L167) sets
  it to 1,000,000 with the comment "No stratum is capped here: the tomato catalog's largest holds
  far fewer loci than this". That is true — the largest holds 217,812 — and it is why nothing is
  ever capped. The design's own figure is 5,000.
- **The borrowing floor.** A stratum under 1,000 tracts takes in rings of same-period neighbours
  until it reaches the floor. Measured: turning borrowing off cut the same run from **1,036.8 s to
  155.5 s**, because 7,824 tracts fitted became 1,032. Borrowing is what lets a thin stratum carry
  an answer at all, so this is a price, not a defect — but it had never been priced.

**Recommendation: cap the tracts a fit reads, and size the cap by measurement.** Fit each pooled
set on a bounded sample of its tracts, raise the bound until the slippage numbers stop moving, and
set it there. This is not a new trade to authorise — it is the knob the specification already
carries, left at infinity. What it spends is statistical precision on the fat strata, which are
precisely the ones that have precision to spare; what it buys is the difference between a fit that
finishes and one that never has. **The alternative — leaving it uncapped and making the inner loop
two or three times faster — turns six days into two, which is not a different answer.**

Everything in §5 is worth applying regardless of which way that decision goes, and none of it
depends on it.

---

## 2a. What was applied and measured

Seven findings — H2 through H7 and one of S1's siblings — were applied to
[ssr_fit.rs](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs) and
[fit.rs](../../../../src/ng/parameter_estimation/joint/fit.rs) during this review, and an eighth
was tried and refuted. Each was timed on its own, one change at a time, on the same host-native
build, four rayon threads, 8 tomato accessions.

| step | change | 6 spans (299 tracts fitted) | 12 spans (590) |
|---|---|---|---|
| baseline | — | 47.1 s | 79.1 s |
| 1 | H5 + H6: hoist `ln B(a,b)` out of the 60-step bisection, and test the argument swap before building the front factor | **39.1 s** | — |
| 2 | H3 + part of H4: build the genotype prior once with the quadrature, and sum the by-descent term over the 13 homozygous slots instead of all 91 | **23.1 s** | — |
| 3 | H2: carry a rescaled running product, one logarithm a quadrature point instead of one a sample | **21.3 s** | 38.1 s |
| 4 | S1's sibling: stop the quantile bisection on bracket width rather than always at 60 steps | **19.7 s** | — |
| 5 | *(hoist the per-stick shapes and their log-Beta out of the 256-point loop)* | *19.8 s — **no gain**, reverted* | — |
| 6 | H7: the ordinary-position fit's per-pass table of error-rate logarithms | 19.7 s (that phase: 10.4 s → 9.8 s) | — |
| 7 | the rest of H4: four running sums in the 91-element dot product instead of one | **15.8 s** | — |
| 8 | the same loop as `wide::f64x4` lanes instead of four scalars | **14.1 s** | **24.6 s** |
| 9 | bump the pinned compiler 1.95 → 1.97.1, nothing else changed | **13.6 s** | — |
| 10 | flatten `TractLikelihoods::scaled` to one buffer with a stride | 13.5 s — **no time gain**, kept for the allocations | — |
| | **cumulative** | **−71.3%** | **−68.9%** |

The whole 6-span run went from 66.8 s to 32.2 s, and the ordinary-position fit from about 10.4 s
to 9.9 s.

**Step 5 is recorded because it failed.** Hoisting the two Beta shapes and their log-Beta out of
the 256-point loop looked like the same class of win as step 1 and moved nothing — 19.7 s to
19.8 s — because step 1 had already reduced the log-Beta to once a bisection, and a suffix sum
over thirteen classes is twelve additions. It was reverted rather than kept, and the reason is
written at the site so nobody re-derives it.

**Steps 7 and 8 are the answer to "can Rust's faster float arithmetic help here", and the answer
is yes but not by that route.** The dot product is a serial `f64` accumulation, so it ran at the
latency of one addition and could not use the machine's vector lanes — Rust will not let the
compiler reassociate it. Rust 1.98's `algebraic_add`/`algebraic_mul` exist to lift exactly that,
but they are unstable on this repository's pinned 1.95 (`error[E0658]: use of unstable library
feature 'float_algebraic'`) and `rust-toolchain.toml` pins the version deliberately, with a comment
saying autovectorisation decisions shift silently between versions. Splitting the sum in the source
instead needs no toolchain change and is reproducible rather than left to a compiler flag: four
explicit accumulators took 19.7 s to 15.8 s, and writing the same association as `wide::f64x4`
lanes — `wide` is already a dependency — took it to 14.1 s. **Together they are worth 28.4% of the
repeat-tract fit, more than any other single change here.**

**The compiler bump was measured on its own, which is what the pin exists for.** Moving
`rust-toolchain.toml` from 1.95 to 1.97.1 — the newest stable this machine has, and still below the
1.98 that stabilises the algebraic operators — was worth 14.1 s → 13.6 s with no source change, and
left every fitted number where it was. `Containerfile` moves from `rust:1.95-bookworm` to
`rust:1.97-bookworm` to match; **until the image is rebuilt, every ephemeral container run will
download the 1.97.1 toolchain**, so the rebuild is worth scheduling.

**Step 10 was kept despite measuring nothing, and the reason is not the clock.** Flattening the
per-tract likelihood rows from `Vec<Vec<f64>>` into one strided buffer moved 13.6 s to 13.5 s,
inside noise — at eight samples a tract carries only a few rows, so there were few pointer hops to
remove. It is kept because it turns one heap allocation *per sample per tract* into one *per tract*:
on a 5,000-tract stratum at 63 samples that is 315,000 allocations becoming 5,000, which is the axis
the allocations review priced this structure on. Unlike step 5, it is not neutral on every axis —
but it should not be reported as a speed win, because it is not one.

### Where the remaining time is, profiled after all ten steps

`sample(1)`, 94,270 samples over 5 threads, 73,356 of them busy:

| | share of busy CPU |
|---|---|
| the per-tract likelihood loop | 35.2% |
| **the ordinary-position kernel `one_position`** | **30.4%** |
| the incomplete-Beta continued fraction | 11.1% |
| `log` | 8.9% |
| `exp` | 7.3% |
| `Scorer::refresh` | 0.6% |

**The two halves are now level.** The repeat-tract half has had 71% taken off and the
ordinary-position half about 6%, so `one_position` is the largest single thing left. Its inner loop
is a running product over samples with a rescale test, node by node — serial by construction, and
not open to the accumulator trick that paid in the repeat-tract half without changing what the
rescale means. The next honest step there is the benchmark seam of §3, not another edit. **Steps 1 and 2 are exactly value-preserving** —
both keep every float addition in its original order, which is why the fitted rows are unchanged.
**Step 3 is not**: a product of likelihoods replaces a sum of their logarithms, so the result may
differ in the last bits. On both inputs every printed fitted quantity — slippage level, shorter
share, fall-off, concentration, and the ordinary-position fit's log-likelihood of 72042 — agreed
with the pre-review baseline at the precision the harness prints. That is agreement to four
decimals rather than bit-identity, and **the owner has ruled that sufficient** (2026-08-15), so no
re-baseline of the identity oracle is required for this change.

Checks run after all steps: `cargo test --lib ng::parameter_estimation` — **644 passed, 0 failed**,
including the positive control `a_drawn_stratum_returns_the_numbers_it_was_drawn_with`,
`a_drawn_cohort_comes_back_at_the_parameters_it_was_drawn_at`, and
`a_cohort_fitted_from_files_gives_the_parameters_it_gives_from_memory`. `cargo clippy --release
--lib` — clean.

**What this does and does not change about §2.** The per-tract cost falls from about 0.14 s to
about 0.045 s at 8 samples, so the six-day extrapolation becomes roughly two. The decision in §2
stands unchanged: taking two thirds off a number that is three orders of magnitude too large is
worth having and is not an answer. Only bounding how many tracts a fit reads is.

Still unapplied and still worth taking: H1 (the owner's call), and everything under Likely — in
particular L1, which would let the identity oracle stop pinning its CPU count, and the last piece
of H4, flattening `TractLikelihoods::scaled` from `Vec<Vec<f64>>` to one buffer with a stride,
which would also let the dot product read contiguous lanes without the per-sample pointer hop.

---

## 3. Measurement plan

The first deliverable was the measurement, and it was taken. This section records what was run so
it can be re-run, then names what is still missing.

### What was run

Host-native release build (`cargo build --release --example ng_joint_records_walk`, landing in
`target/`, which is the sanctioned exception to "cargo runs in the container" — the container
builds Linux binaries the macOS sampler cannot attach to). `RAYON_NUM_THREADS=4` and
`MAX_PASSES=60` throughout; generic census target 60,000 positions; regions are the first N
100-kb spans of `benchmarks/tomato1/regions.bed`. **These wall numbers are not comparable with the
container oracle's ~88 s** — different binary, different thread count, different load.

| run | spans | samples | tracts kept | tracts fitted | ordinary fit | gather | repeat-tract fit | whole run | peak RSS |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 6 | 8 | 349 | 299 | 11.1 s | 0.0 s | **47.1 s** | 66.8 s | — |
| 2 | 6 | 4 | 349 | 299 | 5.3 s | 0.0 s | **29.9 s** | 40.0 s | 572 MB |
| 3 | 12 | 8 | 669 | 590 | 10.8 s | 0.0 s | **79.1 s** | 103.6 s | 860 MB |
| 4 | 24 | 8 | 1,300 | 7,824 | 10.9 s | 0.0 s | **1036.8 s** | 1070.6 s | 994 MB |
| 5 | 24 | 8 | 1,300 | 1,032 | — | 0.0 s | **155.5 s** | 188.3 s | 994 MB |
| 6 | 6 | 8 | 349 | 299 | 10.8 s from files | — | — | — | 929 MB |

Run 5 is run 4 with `SSR_BORROWING_FLOOR=0`. Run 6 is run 1 with `CENSUS_FILES` set.

**The law.** Seconds a tract fitted, at 8 samples and 4 threads: 47.1/299 = 0.157, 79.1/590 =
0.134, 155.5/1032 = 0.151, 1036.8/7824 = 0.133. At 4 samples, 29.9/299 = 0.100. Linear in tracts
fitted at about 0.14 s a tract, and about `0.043 + 0.014 × samples` seconds a tract across the two
sample counts measured.

**Why tracts fitted is not tracts kept.** `fit_strata` fits each *distinct pooled set* once
([ssr_fit.rs:611-648](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L611-L648)). While a
whole motif length holds fewer than 1,000 tracts, every stratum in it borrows the identical set and
the dedup collapses them to one fit — runs 1 and 3 each did exactly one. Once a motif length passes
1,000 tracts each stratum stops at a different ring, the sets stop being identical, and one fit
becomes one fit a stratum. Run 4's eight distinct sets, from the harness's own "borrowed from"
column: 1127, 1126, 1124, 1121, 1116, 1102, 1032 and 76 tracts.

**The extrapolation, labelled as such.** Tomato's selection is 462,701 tracts in 141 strata. Every
stratum contributes at least its own tracts and a stratum under the floor is topped up to about
1,100, so tracts fitted is at least 462,701 and plausibly near 530,000. At 63 samples the law gives
0.95 s a tract, so **about 500,000 seconds — near six days — at four threads.** This stretches a
two-point sample line eight-fold beyond its range: read "days" as solid and the exact figure as
not. For scale, the harness's own comment
([ng_joint_records_walk.rs:288](../../../../examples/ng_joint_records_walk.rs#L288)) records the
*ordinary-position* half of the same cohort at 883 s.

### The profile

`sample <pid> 100`, 1 ms period, whole process, 5 threads (main + 4 rayon workers), against run 1.
Full output:
[sample_tomato_run1.txt](../../../../tmp/perf_review_2026-08-15_ng-census-joint-fit/sample_tomato_run1.txt).
Totals: main thread 50,598 samples, each worker 48,605, 245,018 in all. Self time, verbatim from
`sample`'s own ranking:

```
        rayon::iter::plumbing::bridge_producer_consumer::helper::hb9e994a96ba58f7a        107067
        __psynch_cvwait  (in libsystem_kernel.dylib)        54146
        log  (in libsystem_m.dylib)        26130
        pop_var_caller::ng::parameter_estimation::joint::fit::one_position        25304
        pop_var_caller::ng::parameter_estimation::joint::ssr_fit::ln_gamma        7715
        pop_var_caller::ng::parameter_estimation::joint::ssr_fit::regularised_incomplete_beta        7705
        exp  (in libsystem_m.dylib)        5858
```

**Reading the flat list alone gets this wrong, and the review's first attempt did.** The 44% top
entry is the fully inlined body of the parallel closures in `ssr_fit.rs`, but `Scorer::score` opens
three separate parallel sites and only the call tree separates them. Summing every subtree under
each:

| call site | what it is | samples | share of `fit_strata` |
|---|---|---|---|
| `ssr_fit.rs:926` | the `par_iter` over tracts → `ln_tract` | 30,422 | **77.5%** |
| `ssr_fit.rs:910` | `dirichlet_points` — the 256-point quadrature rebuild | 8,514 | **21.7%** |
| `ssr_fit.rs:899` | `refresh` → `TractLikelihoods::of` | 338 | **0.9%** |

So the quadrature rebuild is worth a fifth of the repeat-tract fit — `ln_gamma` and
`regularised_incomplete_beta` sit underneath it — and `refresh` is worth nothing. Percentages taken
over all 245,018 samples are also diluted by 54,146 waiting samples, ~48,600 of them the main
thread parked by construction; the honest busy denominator is **190,872**, which puts `log` at
13.7% and `one_position` at 13.3%.

The 22% parked entry splits
48,518 on the main thread (which rayon parks by construction, since it is not a pool worker) and
5,628 across the four workers — **the workers idle 2.9% of their samples, so the pool is not
starved.**

### What is still missing, in the order it unblocks other work

1. **A benchmark seam.** `fit_strata` and `fit_stratum` own 70% of the run and cannot be timed
   without opening CRAMs and a reference — so a wall delta on the harness cannot separate a
   regression in `ln_tract` from one in `noodles_cram::Block::decode`. Both functions are already
   `pub` on a fully `pub mod` path and take only plain data with public fields, so **a criterion
   bench can call the real production function with zero lines changed in `src/`**. The fixture
   already exists: `draw_stratum`
   ([ssr_fit.rs:1399-1452](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L1399-L1452))
   builds a `StratumEvidence` from a seeded draw and is trapped inside `#[cfg(test)] mod tests`.
   Cost: a ~150-line bench file, three lines in `Cargo.toml`, and moving ~60 lines behind a
   `bench-fixtures` feature so the bench and the module's own positive control share one generator.
   The generic half's seam is a ~100-line builder generalising the `two_sample_cohort()` helper at
   [census.rs:2811](../../../../src/ng/parameter_estimation/joint/census.rs#L2811), feeding
   `fit_jointly` an in-memory resident census with no CRAM.
   **Build this before applying anything in §5** — every finding there is unfalsifiable until it
   exists, and this review's own first reading of the profile shows how easily a flat self-time
   list misattributes work. This is the single highest-value missing measurement.
2. **A file-backed run of the repeat-tract half.** It has never happened. `refit_from_files`
   ([ng_joint_records_walk.rs:1043](../../../../examples/ng_joint_records_walk.rs#L1043)) calls
   `fit_jointly` only, so `gather_strata`'s file path is covered by unit tests and by nothing else
   — and that is the path the whole by-section memory bound exists for.
3. **Peak resident of the file path alone**, with the resident cohort dropped before the refit.
   Today the harness holds both, so no whole-process peak can separate them. One run.
4. **A thread-count sweep of `fit_strata`** at 1, 2, 4, 8, 18 on the 12-span input. Decides whether
   the per-tract parallel grain is the ceiling (L6, S2). None exists.
5. **A DHAT run.** No allocation profile has ever been taken. Needed to confirm the footprint
   arithmetic in L2, L3 and L9, all of which are computed from `size_of` rather than measured.

---

## 4. Build / toolchain configuration

`[profile.release]` is already right and needs no change: `lto = "fat"`, `codegen-units = 1`,
`panic = "abort"`, `debug = "line-tables-only"`. `rust-toolchain.toml` pins the channel with an
explicit note about keeping criterion baselines comparable. The nine existing benches all use
`black_box` and `harness = false`.

Two items:

- **`.cargo/config.toml` has no entry for aarch64 Linux** — the dev container and the stated
  production target — where it pins `x86-64-v3` for x86 Linux and `apple-m1` for macOS. So every
  container build compiles at rustc's default `generic`, i.e. baseline Armv8.0-A. **This is a
  reproducibility gap, not a speed one, and should not be sold as a win.** Armv8.0-A already
  mandates NEON and f64 FMA; the two loops the profile names are a serial float accumulation and a
  continued fraction with a data-dependent break, neither of which LLVM may vectorise at any
  `target-cpu`; `log`/`exp` are libm ifunc calls that ignore the flag; and LSE atomics are already
  covered by rustc's default `outline-atomics` for this triple. What is missing is that the
  deployment this repo builds most is the only one whose codegen baseline nobody chose. Fix is four
  lines naming `neoverse-n1` (Armv8.2-A, satisfied by Graviton2-and-newer and every M-series core),
  with a comment saying it buys reproducibility.
- **An allocator swap is not indicated, and the question is now closed.** Summing every
  `libsystem_malloc` symbol in the profile gives **516 of 245,018 samples — 0.21%**, against 11%
  for `log` alone. No allocation change in this module can be sold as a wall-clock win on the
  evidence we have; every allocation finding below is filed on the memory axis instead.

---

## 5. Code-level findings

### Hot-path

**H1: [src/ng/parameter_estimation/joint/ssr_fit.rs:602](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L602) — the repeat-tract fit's work is unbounded, and the two knobs that would bound it are both set to infinity on the only path that runs it**

- **Confidence:** High
- **Hot-path evidence:** runs 1, 3, 4 and 5 of §3, and the 39,116-of-50,598 main-thread samples
  under `fit_strata` at `ssr_fit.rs:648`.
- **Mechanism:** cost is linear in tracts fitted at 0.14 s a tract; tracts fitted is set by the
  per-stratum cap (`ssr_cap`, disabled at 1,000,000 in the harness) and the borrowing floor
  (1,000 tracts, which multiplied run 4's work 6.7-fold). Neither is a hot loop; both are policy.
- **Measurement plan:** with the benchmark seam of §3 item 1, fit one drawn stratum at 500, 1,000,
  2,000, 5,000 and 20,000 tracts and plot every fitted slippage number against tract count. The cap
  goes where the numbers stop moving relative to the spread already recorded in
  `str_stratum_size_sweep_2026-08-13.md`. Threshold: the smallest cap at which no fitted quantity
  moves by more than that spread.
- **Complexity cost:** none in code — `ssr_cap` already exists, is already part of the recording
  terms two runs are compared on, and is already plumbed. The cost is statistical and is stated in
  §2.
- **Fix:** set `ssr_cap` in the harness to the measured value rather than 1,000,000, and correct
  the comment at line 166, which explains why the cap never fires as though that were the intent.

**H2: [ssr_fit.rs:846](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L846) — `ln_tract` takes one logarithm per sample per quadrature point, where the sibling half of the same fit already proved one per point is enough**

- ✅ **APPLIED** — step 3 of §2a, 23.1 s → 21.3 s at 6 spans. Not bit-identical by construction.
- **Confidence:** High
- **Hot-path evidence:** `log` at 26,130 of 245,018 samples (11%), which the profile attributes to
  `sum.ln()` in this loop and to `TractLikelihoods::of`.
- **Mechanism:** the ordinary-position half solved this at
  [fit.rs:1627-1632](../../../../src/ng/parameter_estimation/joint/fit.rs#L1627-L1632) with a
  `RESCALE`/`LN_RESCALE` running product — "one logarithm per node instead of one per sample, which
  is a sixty-fold saving in the innermost loop". `ln_tract`'s rows are scaled so the best genotype
  is exactly 1.0, so the product only underflows downward and one threshold suffices. **The saving
  scales with cohort size**: 8 logs a point become 1 at the harness, 3,000 become 1 at the top of
  the committed range.
- **Measurement plan:** the seam bench, plus runs 1 and 2 of §3 (8 and 4 samples) — a change that
  helps at 8 and not at 4 is one whose benefit is in the sample loop, which is the claim. Merge at
  ≥5% on `fitted in`. This changes the floating-point result: require every printed slippage number
  within 1e-9 relative, `cargo test …::ssr_fit` green, and a new test that a tract with several
  hundred samples does not return `-inf`.
- **Complexity cost:** two constants and an underflow branch. Honest cost: both halves of the fit
  would then carry a hand-rolled rescale, and should share one.

**H3: [ssr_fit.rs:819](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L819) — `ln_tract` rebuilds the genotype prior tables for every tract, when they depend only on the quadrature point**

- ✅ **APPLIED** — step 2 of §2a, with H4: 39.1 s → 23.1 s at 6 spans. Value-preserving.
- **Confidence:** High
- **Hot-path evidence:** the 44% inlined-closure entry.
- **Mechanism:** the block filling `identical` and `independent` over 91 slots at each of 256 points
  reads nothing from the tract, yet runs once per tract. The tables belong on `Quadrature`, built
  once in `dirichlet_points` — 186 kB, L2-resident, shared read-only across workers.
- **Measurement plan:** seam bench at one fat stratum, where the redundancy factor is the tract
  count; merge at ≥5%. Bit-identical output.
- **Complexity cost:** two flat `Vec<f64>` on `Quadrature`, which grows from ~26 kB to ~240 kB and
  stops being per-call scratch — worth a line in its doc comment.

**H4: [ssr_fit.rs:838](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L838) and [:747](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L747) — the innermost dot product runs over 91 genotypes for a vector that is zero in 78 of them, through a `Vec<Vec<f64>>` the compiler cannot reason about**

- ✅ **APPLIED** — the by-descent split landed in step 2 of §2a, and the dot product itself in
  steps 7 and 8 (four accumulators, then `wide::f64x4`): 19.7 s → 14.1 s, the largest single win
  here. **Flattening `scaled` is still open** and would let those lanes read contiguous memory.
- **Confidence:** High
- **Hot-path evidence:** the 44% inlined-closure entry. Two categories filed this independently.
- **Mechanism:** `identical[slot]` is set to zero for every heterozygous pair three lines above, so
  78 of 91 products are `0.0 × row[slot]` — a 13-element dot product wearing a 91-element loop, and
  91×4 multiply-adds where 91×2 + 13×2 would do, about 43% fewer. Separately, `scaled` is
  `Vec<Vec<f64>>`, so each sample's row costs a dependent load whose length LLVM cannot relate to
  `genotypes.len()`, leaving a bounds check inside the loop that blocks widening. Flattening to one
  buffer with a stride of 91 fixes the second.
- **Measurement plan:** seam bench, merge at ≥5%; and `cargo asm` on the inlined body, checking for
  `fmla` over `v`-registers rather than scalar `fmadd`. **Splitting the loop reorders the additions**,
  so require the printed numbers within 1e-9 relative.
- **Complexity cost:** a 13-entry list of homozygous slots built with the quadrature, and a stride
  invariant on `scaled`. The loop reads slightly less like the formula and wants a one-line comment
  saying the by-descent prior is zero off the diagonal.

**H5: [ssr_fit.rs:1103](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L1103) — `beta_quantile` recomputes the log-beta constant on all sixty bisection steps, when its two shapes never move**

- ✅ **APPLIED** — step 1 of §2a, with H6: 47.1 s → 39.1 s at 6 spans. Value-preserving.
- **Confidence:** High
- **Hot-path evidence:** the call tree puts `dirichlet_points` at 8,514 samples — **21.7% of the
  repeat-tract fit**, which is itself 70% of the run — and `ln_gamma` (7,715) plus
  `regularised_incomplete_beta` (7,705) are what sits under it. **This is the second-largest single
  win in the report**, and unlike H1 it costs no precision at all.
- **Mechanism:** each bisection step recomputes `ln_gamma(a+b) − ln_gamma(a) − ln_gamma(b)` from
  scratch, and `(a, b)` is fixed for the whole bisection and for all 256 points of a stick. That is
  `256 × 12 × 60 × 3 = 552,960` `ln_gamma` calls a quadrature build where twelve would do. The
  reflection branch swaps `a` and `b`, and log-beta is symmetric, so one constant serves both sides.
  The doc comment at `ssr_fit.rs:874-876` already records the quadrature build at 19 ms; this is
  where most of it goes.
- **Measurement plan:** seam bench; merge at ≥10%. Values change only by float association, so
  require 1e-9 relative agreement.
- **Complexity cost:** one private function gains an argument; the current signature stays as a thin
  wrapper so no other caller moves.

**H6: [ssr_fit.rs:1067](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L1067) — `regularised_incomplete_beta` computes three `ln_gamma`, two `ln` and an `exp`, then throws them away on the symmetric branch**

- ✅ **APPLIED** — step 1 of §2a, with H5. Value-preserving. The bisection above it also gained a
  bracket-width exit (step 4, 21.3 s → 19.7 s).
- **Confidence:** High
- **Hot-path evidence:** the same 21.7% subtree as H5.
- **Mechanism:** `front` is computed before the argument-swap test and is dead on the branch that
  swaps. The condition reads only `a`, `b`, `x`, so it can be tested first at no cost.
- **Measurement plan:** seam bench; merge at ≥3%. Exactly value-preserving, so tests must pass
  **unchanged**, not merely within tolerance.
- **Complexity cost:** none. Two statements swap order.

**H7: [fit.rs:717](../../../../src/ng/parameter_estimation/joint/fit.rs#L717) — `ln_reads_given_genotype` recomputes three logarithms per genotype per candidate per sample per position, of probabilities that are constant for a whole pass**

- ✅ **APPLIED** — step 6 of §2a: that phase went 10.4 s → 9.8 s. Bit-identical; the
  log-likelihood of 72042 is unchanged.
- **Confidence:** High
- **Hot-path evidence:** `one_position` at 25,304 of 245,018 samples (10.3%), against a phase clock
  of 11.1 s for 60 passes over 59,900 positions.
- **Mechanism:** the three probabilities derive from `alt_copies`, `ploidy` and `error_rate` only,
  and the first two are fixed for the run while the third is fixed for a pass. A table of
  `2 classes × read groups × 3 copies` built once in `expectation_pass` replaces up to 18 calls a
  sample a position. Only `total.ln()`, which depends on the sample's own depth spread, must stay.
  This is the same argument the file already makes for `BetaQuadrature` one level in
  ([fit.rs:2384-2386](../../../../src/ng/parameter_estimation/joint/fit.rs#L2384-L2386): "**Computed
  once for the whole pass**, because it depends on nothing that varies from position to position,
  and it sits in the innermost loop").
- **Measurement plan:** run 1's ordinary-position phase clock; merge at ≥8%. Bit-identical, so
  require the printed log-likelihood — 72042 in the baseline — to match exactly.
- **Complexity cost:** a small struct threaded from `expectation_pass` into `one_position`, which
  already carries an `#[allow(clippy::too_many_arguments)]`. The table must be rebuilt whenever the
  rates move, which is once a pass, beside the quadrature that is already rebuilt there.

### Likely

**L1: [fit.rs:1428](../../../../src/ng/parameter_estimation/joint/fit.rs#L1428) — the parallel chunk size is a function of the thread count, and this is the one line the identity oracle's CPU pin exists for**

- **Confidence:** High on the mechanism; the wall consequence is unmeasured.
- **Mechanism:** `POSITIONS_PER_CHUNK.min(positions.div_ceil(rayon::current_num_threads()))` moves
  the chunk boundaries with the thread count whenever `positions < 16,384 × threads` — at 59,900
  positions, 4 chunks of 14,975 at four threads against 6 of 9,984 at six — and the `.reduce` that
  recombines them joins in whatever order the pool finishes, following a split tree also sized from
  the thread count. The float additions are therefore reordered twice over by the CPU count, which
  is exactly the eighth-digit trajectory difference `tmp/run_oracle.sh` pins CPUs to hide. The
  sibling branch four lines above already collects and joins in index order, and says why.
- **Measurement plan:** run the harness at `RAYON_NUM_THREADS` 1, 2, 4, 6, 8 and diff the printed
  fit. Before: differences at the eighth digit. After: byte-identical at every thread count. That
  equality is the threshold. Also watch peak RSS, since `collect` holds every chunk's `Statistics`
  where `reduce` folds them away — at 3,000 samples and 2,000,000 positions that delta is ~145 MB
  and wants an ordered streaming fold instead.
- **Complexity cost:** one constant replaces one expression and the two branches collapse into one
  shape. **This is a deliberate behaviour change requiring a one-time oracle re-baseline** — do not
  merge it as a silent optimisation. What it buys is that the oracle stops needing a CPU pin at all.

**L2: [ssr_fit.rs:1233](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L1233) with [census.rs:1744-1766](../../../../src/ng/parameter_estimation/joint/census.rs#L1744-L1766) — `gather_strata` asks for every stratum at once, so a file-backed run decodes the whole repeat-tract half of every sample simultaneously**

- **Confidence:** High on the call shape; the cost is computed, not measured.
- **Mechanism:** `with_strata` fills every sample's sections before borrowing any, so peak is
  `samples × Σ(sections in the band)` — and the band is everything. At tomato's selection that is
  8.33 GB at a thousand samples for the lent sections alone. The lending *contract* holds — no
  caller can retain a section — but the *bound* the design argues for does not. The code comment at
  `ssr_fit.rs:1228-1232` is honest that this is deliberate and unresolved; this finding is what it
  costs. **Borrowing never crosses motif length**, so the natural band is one period, not all 141
  strata, and looping `with_strata` per period would restore the bound while leaving the borrowing
  rule exactly as it is.
- **Measurement plan:** §3 items 2 and 3 — a file-backed repeat-tract run, and peak RSS with the
  resident cohort dropped. Threshold: if the high-water mark of decoded bytes grows with the total
  file rather than plateauing at the largest band, the bound is not held.
- **Complexity cost:** moderate and structural — `gather_strata` returns after several scoped calls
  rather than one, and `fit_strata` is called per period. No new types.

**L3: [ssr_fit.rs:745-798](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L745-L798) — `Prepared` is the largest structure in the fit at 760 bytes a tract-sample, and it is built three times over**

- **Confidence:** High on the arithmetic; footprint only, not time (see §4's 0.21%).
- **Mechanism:** each row is a `Vec<f64>` over 91 genotype pairs = 728 B heap + 32 B overhead.
  `refresh` collects a new `Vec<TractLikelihoods>` before dropping the old, so peak doubles at the
  moment of assignment; and `fit_pooled` builds a fresh `Scorer` inside the starting-point loop, so
  the whole thing is thrown away and rebuilt for each of three starts. At a 5,000-tract cap and 63
  samples that is 239 MB; on tomato's largest stratum uncapped it is 10.4 GB, and 20.8 GB with the
  refresh transient. **This is the number that decides whether the uncapped selection is runnable at
  all**, and it is a second, independent argument for H1.
- **Measurement plan:** print `tracts × samples × 91 × 8` from `refresh` and compare against peak
  RSS with `MAX_PASSES=0`; or a unit test fitting a drawn 20,000-tract, 100-sample stratum under a
  DHAT ceiling. Threshold: reuse the buffers if `Prepared` exceeds 25% of peak RSS anywhere.
- **Complexity cost:** hoisting the `Scorer` out of the starting-point loop is a one-line move with
  no new invariant and removes two of the three builds. `collect_into_vec` removes the doubling.
  Flattening `scaled` is the same change as H4.

**L4: [ssr_fit.rs:669-725](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L669-L725) and [:228](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L228) — borrowing deep-clones the whole stratum even when it borrows nothing, and re-clones the growing pool once per ring**

- **Confidence:** High
- **Mechanism:** `let mut pooled = evidence.clone()` runs *before* the early return for a stratum
  that needs no borrowing. And `pooled = pooled.pooled_with(neighbour)` starts with `self.clone()`
  inside the ring loop, so taking *k* rings copies the accumulated pool *k* times — O(k²) bytes
  where O(k) suffices. `tracts_with_reads()` is an O(tracts) rescan called on every ring, and again
  from the neighbour filter for every candidate stratum, making that pass O(strata² × tracts).
- **Measurement plan:** a counter, not a timer — total tracts copied, printed after `fit_strata`,
  on the 12-span input. Threshold: fix if it exceeds 3× the tracts in the selection.
- **Complexity cost:** near zero — reorder two statements, return `Cow`, carry a running count, and
  take `self` by value so the extend is in place.

**L5: [ssr_fit.rs:623-648](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L623-L648) — the pooled set is built before the dedup map that exists to avoid building it is consulted**

- **Confidence:** High
- **Mechanism:** `borrow_up_to_the_floor` pools and clones at line 623; `done.entry(key)` is
  consulted at 646. On a hit the whole pooled `StratumEvidence` is discarded unused. The key needs
  only each candidate's `tracts_with_reads()`, not its tracts, so splitting neighbour choice from
  pooling makes the key computable before any copy — and the function's own doc comment says the
  hit case is the common one: "on tomato that is 71 fits over 4,164 tracts collapsing to a handful".
- **Measurement plan:** the same copied-tracts counter as L4; merge if it falls by more than 2×.
- **Complexity cost:** low, and it composes with L4. Needs `tracts_with_reads()` memoised once per
  stratum, which also fixes L4's rescan.

**L6: [ssr_fit.rs:602-665](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L602-L665) — the fit opens roughly 9,200 to 20,500 rayon fork-joins per pooled fit from a thread rayon must wake with a condvar each time**

- **Confidence:** Medium — the fork-join count is exact from source; its cost is not separable in a
  wall-clock sampler.
- **Hot-path evidence:** 48,518 of the main thread's 50,598 samples are `__psynch_cvwait` under
  `Registry::in_worker_cold` → `LockLatch::wait_and_reset`. Only 2,080 are anywhere else.
- **Mechanism:** counted from source, `fit_pooled` makes 4,608 `Scorer::score` calls at one slippage
  group and 10,278 at eight, each opening a `par_iter`, plus a second in `refresh` and a third in
  `dirichlet_points`. Because the caller is the main thread and not a pool worker, each is a full
  latch round trip — inject, sleep on a condvar, be signalled, wake. A pool worker calling the same
  `par_iter` recurses through `join` and executes one half inline instead. Running `fit_strata`
  inside `rayon::scope` removes every round trip. **Nothing about the arithmetic changes**: every
  `par_iter` here collects in index order and the only reduction is a serial sum.
- **Measurement plan:** `fitted in` at 1, 2, 4, 8 threads, best of three, before and after; plus a
  fresh `sample` capture where the main thread's parked count should fall from 48,518 towards the
  worker level of ~1,400. Merge if wall falls beyond noise **and** the fit is byte-identical.
- **Complexity cost:** three lines, no new type. One invariant: `fit_strata` must not then be called
  from inside another rayon closure.

**L7: [census_file.rs:427-456](../../../../src/ng/parameter_estimation/joint/census_file.rs#L427-L456) — opening a census reads a fixed mebibyte that the byte counter never sees, and the same constant is a correctness cliff**

- **Confidence:** High
- **Hot-path evidence:** measured. Run 6's eight censuses were 0.068 MB apiece and `HEAD_READ_BYTES`
  is `1 << 20`, with `read_as_much_as_there_is` stopping at end of file — so all eight were read
  whole before the fit. Real traffic was about 1.85× what the files hold; `bytes_read` reported 0.85×.
- **Mechanism:** the counter that the records spec's re-read argument rests on is blind to the larger
  half of the traffic at small file sizes. On a real tomato census, megabytes a sample, the fixed
  mebibyte is a small share — so this is a blind counter rather than a slow path. Separately, a
  directory past one mebibyte (~40,000 sections) makes a well-formed census decode as `Malformed`.
- **Measurement plan:** count the head read into `bytes_read` and re-run run 6; the printed ratio
  should become ~1.85×. For the cliff, a unit test with a directory of 40,001 sections.
- **Complexity cost:** a 64 KiB probe that grows only when the read did not stop at end of file —
  one loop. **The cliff half is correctness-adjacent and should be routed through a correctness
  review rather than merged as a performance change.**

**L8: [census.rs:1262-1272](../../../../src/ng/parameter_estimation/joint/census.rs#L1262-L1272) — one `lseek` and one `read` per section, where a band's sections abut on disk**

- **Confidence:** Medium
- **Mechanism:** per scoped call per sample: one `open`, then a seek and a read per section, raw and
  unbuffered, each section in a single read — which is the right shape for one section, since the
  buffer is sized to it and a `BufReader` would only add a copy. But the sections of one read group
  are contiguous and asked for in offset order, so a band could be one `pread`. At a thousand
  samples that is 141,000 seeks and 141,000 reads for the tract half against 1,000 reads coalesced.
  Also `buffer.resize(len, 0)` zeroes bytes that `read_exact` immediately overwrites.
- **Measurement plan:** the file-backed repeat-tract run of §3 item 2, with `strace`/`dtruss` counts
  or a syscall counter; merge if coalescing cuts wall time on a cold page cache by ≥10%.
- **Complexity cost:** low — group adjacent extents before reading. The zeroing fix folds into it.

**L9: [contamination.rs:430-433](../../../../src/ng/parameter_estimation/joint/contamination.rs#L430-L433) — four dense `samples × positions` arrays are built where a `samples`-length scratch would do**

- **Confidence:** High on the arithmetic; the site is invisible at 8 samples.
- **Mechanism:** 16 bytes per sample-position, allocated while the lent generic sections are still
  alive. At 63 samples and the production 2,000,000-position target that is 2.02 GB; at a thousand
  samples, 32 GB. Three of the four are read one position at a time in an ordered pass, so the
  `positions` dimension buys nothing — `EvidenceCursor` on the fit's own path already shows the
  single-pass shape.
- **Measurement plan:** raise the harness's generic target from 60,000 to 1,000,000 at 8 samples and
  watch peak RSS; the predicted delta from these four arrays alone is 128 MB. Threshold: restructure
  if peak rises by more than 100 MB.
- **Complexity cost:** moderate — the `major` allele must still be known first, so the totals pass
  stays and only the four arrays collapse.

**L10: [contamination.rs:543](../../../../src/ng/parameter_estimation/joint/contamination.rs#L543) — the ancestry gram matrix is a scalar indexed triangular loop, serial and quadratic in cohort size**

- **Confidence:** Medium
- **Mechanism:** the inner loop is an axpy written with a computed index into a flat `Vec`, so every
  iteration carries a bounds check the compiler cannot hoist. The `samples²` growth is inherent to a
  gram matrix and is not the finding; the constant being several times what it needs to be is. At
  the committed several-thousand end it is 4.5 million multiply-adds a marker on one thread while
  four to eighteen sit idle. **Parallelising over the row index preserves every entry's addition
  order exactly and is bit-identical; parallelising over markers would not be and must not be used.**
- **Measurement plan:** `examples/ng_joint_sample_count_sweep.rs` at 8, 63 and 500 samples. Merge if
  500 improves ≥20% and 8 does not regress; confirm the widening with `cargo asm`.
- **Complexity cost:** one hoisted local and a slice split for the scalar half; one `into_par_iter`
  for the parallel half, at the cost of recomputing the centring per row.

**L11: [.cargo/config.toml](../../../../.cargo/config.toml) — the aarch64-Linux target has no `target-cpu` entry.** Detail and the reason it is a reproducibility fix rather than a speed one are in §4.

**L12: [loci.rs:837](../../../../src/ng/parameter_estimation/joint/loci.rs#L837) — the repeat catalog is reopened and its parquet footer and page index re-parsed once per contig per pass, and the selection makes two passes**

- **Confidence:** Medium
- **Hot-path evidence:** the selection phase is 3.0-3.3 s in every run, against 66.8 s to 1070.6 s
  totals — so this is small today and grows with contig count.
- **Measurement plan:** time the selection phase at 13 contigs and at a reference with hundreds;
  merge if it grows with contig count rather than with bases.
- **Complexity cost:** low — open once and reuse the parsed metadata.

**L13: [census.rs:194](../../../../src/ng/parameter_estimation/joint/census.rs#L194) — `never_walked` writes every depth code through a bit-shifting setter when the sentinel is all-ones and one `memset` says the same thing**

- **Confidence:** High on the pattern; the site is cold-ish (once per read group per sample).
- **Mechanism:** `NEVER_WALKED_CODE` is `0b11111`, so the buffer is `0xff` bytes. At the two million
  kept positions the type's own comment mentions, that is two million iterations of assert, two
  divisions and a masked read-modify-write, against one `memset`.
- **Measurement plan:** the walk phase clock at 12 spans; if it does not move, take it anyway as a
  one-line simplification and file at Note.
- **Complexity cost:** none, but the equality between the sentinel and the all-ones mask becomes
  load-bearing and needs a `const` assertion.

**L14: [fit.rs:918](../../../../src/ng/parameter_estimation/joint/fit.rs#L918) — `Statistics::genotypes` is `Vec<Vec<[f64; 3]>>`, written in the innermost sample loop right beside a flat array with the opposite index order.** Flatten it node-major, sample-minor to match `BetaQuadrature::priors`. Measurement: the ordinary-position phase clock, merge at ≥3%.

**L15: [contamination.rs:193-220](../../../../src/ng/parameter_estimation/joint/contamination.rs#L193-L220) — `Marker` is five per-marker heap vectors and `fit_alpha` reads one element of each, hundreds of times a sample.** Project the one sample's columns into a flat cache before the golden-section search. Measurement: the contamination phase at 63 accessions and 2,000,000 positions; merge at ≥10% of that phase.

### Speculative

**S1: [ssr_fit.rs:1131](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L1131) — `climb_scalar` always spends sixteen golden-section evaluations on a coordinate, with no bracket-width exit.** Each evaluation is a full pass over every tract, so the multiplier is real; sixteen steps shrink a bracket by 0.0007, and after the first round most coordinates are already settled. **This is an estimator question before it is a performance one** — measure accuracy first, and merge only at ≥15%, because a change to the search must earn more than a mechanical hoist.

**S2: [ssr_fit.rs:602-665](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L602-L665) — parallelise across distinct pooled sets instead of within one stratum's tracts.** The independent unit of work is the pooled set, and per-fit cost varies enormously; the current grain spreads at most a stratum's tracts. **The profile does not support this**: workers idle 2.9%, so there is no idle time to reclaim on the measured input, and no thread-count sweep exists. It changes no float sum — each score collects in tract order and sums serially — but the inner iterators must go serial or the outer one re-enters the same pool. Do L6 and the thread sweep first.

**S3: [census.rs:1260-1272](../../../../src/ng/parameter_estimation/joint/census.rs#L1260-L1272) — `mmap` the sections instead of open/seek/read.** Worth a look only after L8; brings its own failure modes on a truncated file.

**S4: [census_file.rs:573-580](../../../../src/ng/parameter_estimation/joint/census_file.rs#L573-L580) — the directory's no-overlap check is O(n²) in section count, paid once per file opened.** At 141 strata it is ~20,000 comparisons; at a thousand samples it is paid a thousand times. Sort by offset and check neighbours.

**S5: [census_file.rs:316](../../../../src/ng/parameter_estimation/joint/census_file.rs#L316) — `encode_ssr`/`decode_ssr` move nine offset counts one `u16` at a time into a `Vec` that never reserves.** Both the sizing and the batching are two lines; the final size is exactly computable, and a `debug_assert_eq!` pins the expression to the encoder.

**S6: [ssr_fit.rs:521](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L521) — `climb_one_round` clones the whole `Parameters` on each of ~10,000 objective evaluations a fit.** Expect no wall-clock change (§4's 0.21%); the save-and-restore alternative adds an invariant to a numerical routine, which is exactly where a silent bug would not show in the output. Do not take it without a bit-for-bit test on a drawn stratum.

**S7: [fit.rs:484](../../../../src/ng/parameter_estimation/joint/fit.rs#L484) — `PositionEvidence::depth_weights` reserves 32 slots a sample and uses one at three reads a position.** Size the stride to the run's own widest bin.

**S8: [census_file.rs:202](../../../../src/ng/parameter_estimation/joint/census_file.rs#L202) — `write_census` encodes every section into memory before writing a byte**, because the directory needs every length and offset first. A second whole copy of the census at peak. Land S5's size expressions first; they are exactly what a streaming version needs.

**S9: [ssr_fit.rs:775](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L775) — `TractLikelihoods::of` takes a logarithm per genotype per bucket per sample per tract, of a value that depends on none of those four.** `probability` is a function of `(group, genotype, bucket)` alone, and the whole table is `groups × 91 × 9` f64 — 6.5 kB a group, L1-resident — buildable in `refresh` beside the input it derives from. **This review's first pass filed it as Hot-path and was wrong**: the flat self-time list attributes `log` to this site and to `ln_tract` together, but the call tree puts the `refresh` subtree at 338 samples — **0.9% of the repeat-tract fit**, because `refresh` only reruns when slippage actually moves, which the climb does far less often than it moves the spectrum. The hoist is still correct and still bit-identical; it is simply worth almost nothing at this cohort size, and it would grow only with the number of slippage groups. Take it if H3 and H4 are being done anyway, since it touches the same loop nest; do not schedule it on its own.

### Note

- **The guard-list scan is not quadratic in practice, and the measurement says so.**
  [ssr_fit.rs:1291-1299](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L1291-L1299)
  scans a stratum's whole guard list twice per tract per sample per read group, which is
  `O(tracts × samples × groups × guard_len)` — three categories filed it independently, and it was
  the orchestrator's own leading hypothesis for why the full cohort never finishes. **On real tomato
  reads the guard list holds 1 entry across the whole 8-sample run and 7 across the 24-span one**,
  and `gathered in 0.0 s` in every run. The shape is real and cheap to fix (bucket the guard by
  locus in one pass); the constant is nil. Revisit only if a cohort ever shows a large guard list.
- **The ordinary-position half already contains the pattern the repeat-tract half needs.**
  [fit.rs:1541](../../../../src/ng/parameter_estimation/joint/fit.rs#L1541)'s `Scratch` is sized
  once and reused — "so a two-million-position pass allocates nothing" — and `EvidenceCursor` reads
  the lent sections in one ordered pass with a per-sample cursor rather than a binary search per
  position. H2, H3, L3 and S9 are all asking `ssr_fit.rs` to do what `fit.rs` already does. **The
  fix is a transplant, not an invention.**
- **`generic/depth_bins.rs` has no findings, and three categories agree.** The ladder is built once
  per run; `bin_for` is a `partition_point` over a 30-element slice and `depth_range` is index
  arithmetic. It is not a hot-loop surface, allocates nothing on the hot path, and performs no I/O.
- **`SampleAtPosition` ([fit.rs:432](../../../../src/ng/parameter_estimation/joint/fit.rs#L432)) is
  exactly 64 bytes and the array-of-structs choice around it is right — do not "fix" it.**
- **The repeat-tract half has never been run against a file-backed cohort.** `refit_from_files`
  calls `fit_jointly` only, so `gather_strata`'s file path is exercised by unit tests and nothing
  else. This is a test-coverage gap with performance consequences, and it is the path the whole
  by-section design exists for.
- **Every sample keeps its own full `RecordingTerms`** after `CohortCensusEvidence::new` has proved
  they are all equal — at a thousand samples, a thousand identical copies of a per-cohort table.

---

## 6. Out-of-scope observations

- [src/ng/reference_info.rs:518-520](../../../../src/ng/reference_info.rs#L518-L520) — the FASTA
  pass walks the reference **one byte at a time** through a `?`-returning call on a `&mut dyn`
  observer, for every base. This sits inside the 3.0-3.3 s selection phase that every run in §3
  pays. The I/O around it is fine (a 64 KiB buffer, no `BufReader`, which is the right shape); the
  per-byte dynamic dispatch is not. Separate PR, and it belongs to whoever owns `reference_info`.
- `benches/psp_writer_perf.rs:386` panics under `cargo test --all-targets`, and the aggregate clippy
  run is red on five errors — three in `ssr_fit.rs`'s test module, two in
  `examples/ng_duplicated_class_harness.rs`. Both pre-date this branch and both are already recorded
  as standing items in `PROJECT_STATUS.md`. Noted only so a reader of this report does not mistake
  them for something it found.

---

## 7. What's already good

- **The ordinary-position fit allocates nothing per position.** `Scratch`
  ([fit.rs:1541](../../../../src/ng/parameter_estimation/joint/fit.rs#L1541)) is sized once and
  reused across a whole pass, and `EvidenceCursor` reads each sample's sparse list with a cursor
  rather than a binary search per position — which is why `one_position` shows up as arithmetic and
  not as allocator traffic.
- **The read likelihoods are computed once per candidate and reused across every quadrature node**
  ([fit.rs:1548-1552](../../../../src/ng/parameter_estimation/joint/fit.rs#L1548-L1552)), with a
  comment saying this is what made the program runnable at all. It is the same insight H3 and S9
  ask the repeat-tract half to apply.
- **The parallel sum is deliberately kept in index order** at
  [ssr_fit.rs:919-927](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L919-L927) — "a
  parallel float sum reorders the additions run to run, and a fit's whole output is a difference
  between one set of parameters and another". That discipline is why L1's fix is available at all:
  the repeat-tract half is already thread-count-independent, and only the generic half is not.

---

## Author response convention

Address each finding by its identifier (H1, L2, …) with one of: `applied in <commit>` /
`experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` /
`won't fix because …`. The "no gain" path is expected and welcome — that is what the measurement
plans are for.
