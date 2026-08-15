# Performance Review: ng-joint-fit, round two
**Date:** 2026-08-15
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** the joint parameters fit and the census — `src/ng/parameter_estimation/joint/` plus `generic/depth_bins.rs`, reviewed a second time after round one applied ten changes to it
**Verdict:** Apply the listed wins — but the wins are two pieces of measuring equipment and one settled decision, not code in the hot loop. **Every code-level candidate this round measured a regression.**
**Hot-path evidence:** a criterion benchmark and an allocation profile, both built during this round and neither of which existed before it; plus round one's sampling profile, re-attributed line by line

---

## 1. Scope and constraints

**What was reviewed.** The same module as round one: the joint route of step 4's parameter
pre-pass, which estimates error rates, heterozygosity, inbreeding, contamination and repeat-tract
slippage once for a whole cohort, before any variant is called. It is not the variant caller.

**Reviewed against** commit `07e552c3`, branch `ng-joint-fit-perf2` cut from `main`, worktree
`/Users/jose/devel/pop_var_caller-ng-joint-fit`.

**Targets.** Per [CLAUDE.md](../../../../CLAUDE.md) §0 the fit must degrade gracefully from **one
sample to several thousand** and from **three reads a position to several hundred**.

**Target hardware, and a correction to round one.** This host is an 18-core Apple M-series with
**6 performance cores and 12 efficiency cores**, not 18 equal ones (`sysctl hw.perflevel0.logicalcpu`
= 6, `hw.perflevel1.logicalcpu` = 12). Every measurement in this report and in round one was taken
at `RAYON_NUM_THREADS=4`, which fits inside the performance cores; above six, added workers land on
efficiency cores and a wall-clock number stops being a statement about the code. **Round one's
proposed thread sweep to 18 would have measured the core mix, and any future one must report the
knee at six as hardware.**

**Deliberately out of scope.** `src/psp/`, `src/pileup/`, `src/ssr/`, `src/var_calling/`,
`src/vcf/` — production, frozen. The other `examples/ng_joint_*` harnesses and
`examples/ng_depth_term_family.rs`.

**Out of scope by the owner's instruction, mid-review: everything about parallelism.** The
architecture of the parallelisation is being decided separately. Two working patches and a
thread-count sweep were therefore dropped after being written; their evidence is recorded in §6 for
whoever takes that decision, and neither patch is applied.

---

## 2. Verdict

**Apply the listed wins.** There are three, and none of them is a faster loop.

**1. The fit can now be benchmarked without opening a CRAM file.** `benches/ng_joint_fit_perf.rs`
drives `fit_stratum`, `fit_strata` and `fit_jointly` on drawn in-memory evidence in about two
minutes, along both axes [CLAUDE.md](../../../../CLAUDE.md) §0 commits to — one sample to 32, and
32 tracts to 128 — with **zero lines changed in the body of any production function**. Round one
named this the single highest-value missing measurement and said to build it before applying
anything else. It is built, committed, and it is what refuted every code candidate below.

**2. The module has an allocation profile for the first time, and it says the memory is not where
the design put it.** 86% of the ordinary-position half's peak resident heap is the contamination
arm, and 69% of the whole peak is four arrays on four lines. The estimator's own reusable scratch
holds 191 kB at that instant.

**3. The tract cap was never an open decision — it was decided on 2026-08-13 and the harness never
applied it.** Round one filed this as the owner's call. It is not: the measurement exists, the spec
records the question closed, and the only defect is one line in a measurement harness. Details in
§2a.

**And the finding that shapes the next round: the ordinary-position kernel did not respond to any
change tried against it.** Six candidates were implemented, each verified in the disassembled
binary to do what it claimed — vector loads appearing where there were none, divisions removed,
bounds checks removed, dependent pointer loads removed — and every one measured **slower**. That
result and its control are §5.

---

## 2a. The tract cap: closed, and closed before this review began

Round one's H1 said the repeat-tract fit's work is unbounded and that bounding it was a design
decision for the owner. **Both halves of that are wrong, and the correction matters because H1 was
ranked above every other finding in the report.**

**The measurement exists.**
[str_stratum_size_sweep_2026-08-13.md](../../ng/reports/str_stratum_size_sweep_2026-08-13.md) drew
repeat tracts at a known truth and fitted them at 50, 100, 250, 1,000, 5,000 and 20,000 tracts, at
20 samples, at three and at six reads a site, with twelve draws at the small counts. **At three
reads a site — tomato's depth — the floor is 5,000 tracts**; at six it is 1,000.
[parameter_prepass_joint_loci.md](../../ng/spec/parameter_prepass_joint_loci.md) §6 question 1
records the question as closed on that date, and its §4.5 carries the consequence: at a cap of
5,000 tomato keeps 86,688 of its 462,701 tracts and **8 of its 141 strata are capped at all**.

**What breaks first is not the number anyone would guess.** The slippage level is carried by every
read at every tract and is still within 4% of truth at 250 tracts. What goes first is the fall-off
— how fast two-repeat slips fall away against one-repeat slips — and the concentration, which says
how alike a stratum's tracts are. **The concentration is counted in tracts, not in reads**, so a
cohort sequenced at 30× needs the same tract count as one at 3×; depth cannot buy it. That is why
the floor is set by the shallow case.

**This round adds the half that report says it lacks: real reads.** Its own §6 states that no real
reads were in it. A sweep was run here on 8 tomato accessions over 24 spans of 100 kb, with
borrowing off so the cap is the only thing that moves between runs, and with nested subsets — the
100-tract run's tracts are the first 100 of the 200-tract run's.

| motif | ref. repeats | tracts | reads | level | shorter share | fall-off | concentration |
|---|---|---|---|---|---|---|---|
| 1 | 8 | 97 | 5,197 | 0.0042 | 0.860 | 0.709 | 1.574 |
| 1 | 8 | 192 | 10,397 | 0.0047 | 0.729 | 0.640 | 0.976 |
| 1 | 8 | 388 | 23,282 | 0.0033 | 0.652 | 0.629 | 0.550 |
| 1 | 8 | **625 (all)** | 38,090 | **0.0026** | **0.654** | **0.676** | **1.667** |
| 1 | 9 | 94 | 4,574 | 0.0046 | 0.574 | 0.528 | 5.235 |
| 1 | 9 | 190 | 10,427 | 0.0043 | 0.553 | 0.575 | 2.955 |
| 1 | 9 | **302 (all)** | 16,960 | **0.0038** | **0.579** | **0.668** | **2.440** |
| 1 | 10 | 94 | 4,719 | 0.0044 | 0.903 | 0.671 | 8.832 |
| 1 | 10 | **105 (all)** | 5,222 | **0.0046** | **0.873** | **0.708** | **7.849** |

**It agrees with the drawn sweep.** 625 tracts is an eighth of the measured floor, so numbers still
moving there is what that report predicts. The shorter share — the one of the four that report says
is learned fastest — is the one that settles here too: 0.652 at 388 tracts against 0.654 at 625.
The concentration, which that report says is learned slowest, swings threefold across the four caps
with no direction.

**The cap of 50 produced nothing at all**, and that is worth recording as a floor on the floor: the
fit refuses a stratum holding fewer than 50 tracts *with reads*, and capping at 50 leaves 47 or 48
once the tracts no read crossed drop out. The lowest usable cap sits a little above the refusal
floor.

**The cost law, refined.** Round one modelled the repeat-tract fit as 0.045 seconds a tract. These
four runs, which hold the number of distinct fits at exactly three and move only the tracts, say it
is **3.9 seconds a fit plus 0.036 seconds a tract**, at 8 samples on four threads (21.8 s at 285
tracts, 28.4 at 487, 40.8 at 795, 48.5 at 1,032; the straight line through those four points misses
none by more than 0.8 s). The fixed part is the 256-point quadrature build and the coordinate
climb's own overhead, and it does not fall with the cap.

**So temper what the cap buys.** 462,701 tracts become 86,688 — a factor of 5.3, not the hundredfold
round one's §2 was reaching for. And 141 strata × 3.9 s is about ten minutes of fixed cost that no
cap touches. **After the cap, borrowing is the larger lever**, and unlike the cap it is not free: 68
of tomato's 141 strata hold fewer than a hundred tracts each, are far under any cap already, and
reach a fittable size only by borrowing from neighbouring repeat counts.

**Applied:** [ng_joint_records_walk.rs:167](../../../../examples/ng_joint_records_walk.rs#L167) now
sets `ssr_cap: 5_000` with the measurement and its consequence written at the site. The previous
value was 1,000,000, above the largest stratum tomato has, with a comment explaining that the cap
therefore never fires — which read as intent and was a defect.

---

## 3. Measurement plan

Round one's §3 listed five missing measurements. Three were taken this round, one was dropped by
the owner, and one remains.

### Taken: the benchmark seam (round one's item 1, "the single highest-value missing measurement")

```sh
cargo bench --features bench-fixtures --bench ng_joint_fit_perf -- --save-baseline <name>
cargo bench --features bench-fixtures --bench ng_joint_fit_perf -- --baseline <name> '<filter>'
```

Five groups, ten cases, about two minutes a whole run at `sample_size(10)`. Baseline on an idle
machine at `RAYON_NUM_THREADS=4`:

| group | case | time |
|---|---|---|
| `stratum_by_tracts` | 32 | 437.56 ms |
| | 128 | 670.79 ms |
| `stratum_by_samples` | 1 | 374.16 ms |
| | 8 | 437.91 ms |
| | 32 | 665.26 ms |
| `strata` | `stands_alone` | 2.3793 s |
| | `borrows_and_shares` | 599.46 ms |
| `generic_by_positions` | 5,000 | 47.204 ms |
| | 20,000 | 185.13 ms |
| `generic_by_samples` | 1 | 17.261 ms |
| | 8 | 47.382 ms |
| | 32 | 160.48 ms |

Every confidence interval above is within 0.5% of its own point estimate.

**A run of this bench on a busy machine is worthless, and the difference is not subtle.** The first
attempt at this baseline was taken while two review sub-agents were running test suites at about
900% CPU each; `strata/stands_alone` came back as `[7.5385 s 10.737 s 13.606 s]` against the
`[2.3781 s 2.3793 s 2.3807 s]` above, and `generic_by_positions/5000` at 159.80 ms against 47.20 ms.
It was discarded. **The check before timing must look for the compiled test binary, not for
`cargo` or `rustc`** — a `cargo test` that has finished compiling shows neither name in the process
list, which is how the contaminated run got started.

### Taken: the allocation profile (round one's item 5)

```sh
cargo build --profile profiling --features dhat-heap --example dhat_ng_joint_fit
target/profiling/examples/dhat_ng_joint_fit <all|generic|tracts|gather|prepared>
```

`examples/dhat_ng_joint_fit.rs` needs no CRAM and no reference, changes no line of `src/`, and
reaches both halves of the fit. Built under `--profile profiling` rather than `release` because
release's `lto = "fat"` collapses the frames dhat attributes by. Results in §5.

### Taken: the tract-cap sweep on real reads

`tmp/perf_tract_cap_sweep_v2.sh`, driving `tmp/perf_scaling_tomato.sh` at a fixed four rayon
threads with `SSR_BORROWING_FLOOR=0` and `SSR_TRACT_CAP` moving. Results in §2a.

### Dropped by the owner: the thread-count sweep (round one's item 4)

The parallel architecture is being decided separately. Not run, and the design point recorded above
about six performance cores against twelve efficiency ones is what any future sweep must account
for.

### Still missing: the repeat-tract half has never been run from census files (round one's items 2 and 3)

`refit_from_files` ([ng_joint_records_walk.rs:1025](../../../../examples/ng_joint_records_walk.rs#L1025))
calls `fit_jointly` only, so `gather_strata`'s file path is exercised by unit tests and by nothing
else — and that is the path the whole by-section memory bound exists for. Peak resident of the file
path alone, with the resident cohort dropped first, is the same run. **This is the top of the next
round's list**, and the allocation profile below sharpens why: `gather_strata` builds a structure
between 2.4 and 5.7 times the size of the sections it derives it from, and holds both at once.

---

## 4. Build / toolchain configuration

`[profile.release]` needs no change and round one's reading of it stands. Three items, one of them
new and important:

- **`cargo asm` reports false negatives on this repository, and every future run of this review
  must know it.** `[profile.release] lto = "fat"` makes rustc emit each crate through the
  *pre-link* optimisation pipeline, which defers LLVM's loop vectoriser to link time. So `cargo asm`
  shows scalar code for loops that end up fully vectorised. This was caught with a control — a
  textbook `out[i] = a[i] * b[i]` vectorises standalone under `rustc -O` and is scalar inside this
  crate under `cargo asm` — and every codegen claim in §5 is therefore read from `objdump -d` on the
  **linked** binary instead. A second, smaller trap: Apple's aarch64 syntax writes vector
  instructions as `fmul.2d v0, v1, v2`, not `fmul v0.2d, …`, so a grep for the latter silently
  returns nothing. **Both belong in
  [profiling_environment.md](../../../../ai/skills/rust-performance-review/performance_review/profiling_environment.md).**
- **That file's toolchain pin is out of date.** It says `rust-toolchain.toml` pins 1.95; round one
  moved it to 1.97.1. Every sub-agent reads that file before writing a measurement plan.
- **The aarch64-Linux `target-cpu` gap stands**, unchanged from round one: `.cargo/config.toml`
  names x86-64 Linux and macOS aarch64 and not the container's own target. It remains a
  reproducibility gap rather than a speed win, and the benchmark seam now makes it directly
  testable.
- **The allocator question stays closed on the time axis** (0.21% of samples in round one's
  profile) and is now open on the memory axis, where §5 files it.

---

## 6. Out-of-scope observations

### The parallelism evidence, recorded and not acted on

The owner ruled mid-review that the parallel architecture is being decided separately. Two patches
were already written and verified by then; both are left unapplied and their evidence is set down
here so that decision does not have to re-derive it. Diffs are in
[tmp/perf_review_2026-08-15_ng-census-joint-fit_v2/](../../../../tmp/perf_review_2026-08-15_ng-census-joint-fit_v2/).

- **The fit's answer depends on how many cores the machine has, and that is why the identity
  oracle pins its CPU count.** The ordinary-position pass sizes its chunks as
  `POSITIONS_PER_CHUNK.min(positions.div_ceil(rayon::current_num_threads()))`
  ([fit.rs:1489](../../../../src/ng/parameter_estimation/joint/fit.rs#L1489)) and recombines them
  with `.reduce`, so both the split and the join order follow the thread count and the float
  additions are reordered twice over. Measured on the same 6-span, 8-accession input, the pass-1
  `largest move` column reads **615262.063054 at four threads, 615262.050187 at six, 615262.047243
  at eight**. A patch that sizes the chunk from the position count alone and folds with a fixed
  binary tree over position order makes all three **615262.062660**, with all 205 printed lines
  agreeing at every thread count. Two further facts for whoever takes this: rayon re-seeds the
  reduce splitter on *every steal*
  (`rayon-1.12.0/src/iter/plumbing/mod.rs:273`), so the join tree was never even a pure function of
  the thread count; and the unpatched binary run twice at a pinned eight threads gives identical
  output, so the pin does currently hide the problem rather than merely narrowing it.
- **Moving the repeat-tract fit onto a pool thread is a regression, and this one is measured.**
  Round one argued that `fit_strata` opens thousands of parallel loops from a thread rayon must
  wake with a condition variable each time, and that `rayon::scope` would turn each handover into an
  inline `join`. On the benchmark seam it is **slower**: +5.1% at 128 tracts, +2.8% at 32 samples,
  +1.1% at 8 samples, no change at 1. The ceiling was small to begin with — summing the round-two
  profile per thread, the four workers idle 2.9% of their samples, and the main thread's parked
  time is not reclaimable CPU — and the change costs more than the ceiling is worth at four
  threads. **Whether it changes sign on a wider pool is unmeasured**, and is one of the things a
  thread sweep would settle.

### From round one, still standing

- [src/ng/reference_info.rs:518-520](../../../../src/ng/reference_info.rs#L518-L520) — the FASTA
  pass walks the reference **one byte at a time** through a `?`-returning call on a `&mut dyn`
  observer. It sits inside the 3.0–3.3 s selection phase every harness run pays. Separate PR, and
  it belongs to whoever owns `reference_info`.
- `cargo clippy --all-targets --all-features -- -D warnings` and the `-D warnings` rustdoc build are
  **already red at this commit** for reasons that pre-date this branch — `src/ssr/cohort/sim.rs`,
  four examples, two `useless_vec` sites in `ssr_fit.rs`'s own test module, and 25 rustdoc link
  errors. Noted only so that a reader does not mistake them for something this review introduced,
  and so that nobody reads a green aggregate gate as evidence a patch is clean. `cargo clippy
  --release --lib` is clean.

---

## 7. What's already good

- **The ordinary-position fit's peak heap is 191 kB while it holds 30,000 positions of 24 samples.**
  `Scratch` ([fit.rs:1603](../../../../src/ng/parameter_estimation/joint/fit.rs#L1603)) is sized
  once and reused across a whole pass, and the allocation profile now confirms the claim its comment
  makes rather than merely repeating it — the peak of that half of the run is the contamination arm,
  not the estimator.
- **The parallel sum is deliberately kept in index order** at
  [ssr_fit.rs:919-927](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L919-L927), with the
  reason written at the site: "a parallel float sum reorders the additions run to run, and a fit's
  whole output is a difference between one set of parameters and another". That discipline is why
  the repeat-tract half is already thread-count-independent while the generic half is not.
- **`generic/depth_bins.rs` has no findings, and now four categories over two rounds agree.** The
  ladder is built once per run, `bin_for` is a `partition_point` over a 30-element slice with an
  exact-region fast path already in place, and the file appears nowhere in the profile.

---

## 5. Code-level findings

### 5a. The measured result, first — six candidates against the ordinary-position kernel, one win

Round one left the ordinary-position kernel `one_position` at **30.4% of busy CPU**, the largest
single thing in the run and the half that had had only one fix applied to it. This round the
`hot_loops` category attributed the profile line by line — parsing the whole call tree rather than
reading the flat self-time list, and reconstructing all 94,270 samples exactly as the check — and
found `one_position` and everything under it at **28,469 samples: 30.2% of all samples, 38.8% of
the busy ones**, distributed like this:

| lines | what the block does | share of `one_position` |
|---|---|---|
| 2034–2060 | segregating branch: the per-node, per-sample genotype weights | **27.0%** |
| 2097–2137 | duplicated branch: the same shape over carriers | **14.1%** |
| 1753–1787 | the read likelihoods, per class × candidate × sample × genotype | **11.9%** |
| 1790–1817 | the integral over the position's own allele frequency | 9.3% |
| 1826–1861 | the integral over the carrier frequency | 8.3% |
| 2061–2091 | segregating branch: the read tallies | 7.5% |
| 1864–1919 | the four branches and the two classes | 7.2% |
| 2138–2168 | duplicated branch: the read tallies | 5.3% |

Six candidates were implemented against those blocks, each verified in the **linked** binary to do
what it claimed. Each was then timed on its own, one patch in the tree at a time, against the same
saved criterion baseline, at four rayon threads on an idle machine.

| candidate | what it does | codegen change, verified | 5,000 pos. | 20,000 pos. | 1 sample | 8 samples | 32 samples |
|---|---|---|---|---|---|---|---|
| **hoist the reference-read term** | turns the loop nest inside out so a term that does not depend on the candidate allele is computed 3 times a sample instead of `1 + 3 × candidates` | bounds checks 109 → 104 | **−3.5%** | **−2.6%** | **−1.1%** | **−1.7%** | **−2.5%** |
| flatten the genotype counts node-major | `Vec<Vec<[f64;3]>>` → one buffer, removing a dependent pointer load and two bounds checks per element | two integer loads and one branch removed per element | +2.1% | +3.3% | −0.4% | +10.5% | +1.6% |
| the reference-probability table | stops recomputing a probability, including a division, that the pass already has in a table | scalar `fdiv` 21 → 12, instructions 4,655 → 4,539 | +3.7% | +5.0% | +2.9% | +4.1% | −0.02% |
| the vectorised term buffer | splits the per-sample factor out of the serial product so the factor pass can widen | `ld3.2d` 0 → **16**, `fmul.2d` 12 → **76**, `fadd.2d` 9 → **41** | +2.7% | +4.8% | +6.0% | +4.8% | +1.8% |
| one reciprocal, no branch | three divisions and a data-dependent skip become one reciprocal and a select | scalar `fdiv` 21 → 18 | −0.01% | +1.3% | +0.7% | +1.9% | +2.1% |

Negative is faster. Every cell has p < 0.05 against the baseline except the four reading within
half a percent of zero.

**Only one of the six is a win, and it is the one that removes *executions* rather than
instructions.** `ln_reference_reads` reads the sample, its depth weights and the copy count — not
the candidate allele — so with the candidate loop outermost it and `SampleAtPosition::non_reference`
both run ten times per sample per noise class where three and one would do. Turning the nest inside
out changes no float's order, so the result is bit-identical.

**Four candidates got the codegen they were designed to get and were slower anyway.** The starkest
is the vectorised term buffer: it took the frequency integral from zero vector loads to sixteen and
from twelve vector multiplies to seventy-six, with the emitted loop unrolled four ways producing
eight samples' factors an iteration — and it is **slower at every size measured, from +1.8% at 32
samples to +6.0% at one.** The reference-probability table strictly removes work — nine scalar
divisions and 116 instructions — and is 5.0% slower at 20,000 positions.

**What that means for the next round, and it is the most useful thing this one found.**
`one_position` is a single large function that the compiler inlines whole and register-allocates as
one body; the four refuted candidates all add state to it — a scratch buffer, a changed layout, a
hoisted table, a restructured index — and each apparently costs more in that body's register
pressure and code layout than it saves in the loop it targets. **The one that won removed nothing
from the body and instead ran the body fewer times.** A third round on this kernel should be
looking for redundant *executions*, not for better instructions, and should not trust codegen
inspection as a proxy for anything: on this evidence the correlation between "the assembly improved"
and "the program got faster" is zero at best.

**A methodology point that generalises past this module.** Round one recorded that reading `sample`'s
flat self-time list alone produced a wrong finding. This round records the sequel: reading the
*disassembly* alone produces wrong findings too, and produced four of them. Both are inputs to a
hypothesis; only the benchmark is an answer. That is precisely why round one ranked the benchmark
seam above every code change, and this round is the demonstration.

---
