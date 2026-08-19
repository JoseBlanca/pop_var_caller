# Performance Review: ng cohort merge
**Date:** 2026-08-19
**Reviewer:** rust-performance-review skill (orchestrator), with five per-category reviewers in isolated worktrees
**Scope:** `src/ng/run/cohort_merge/` — all three drivers, the observation cache, the k-way walk and the builder
**Verdict:** Apply the listed wins — seven were applied and measured; eight were measured and closed
**Hot-path evidence:** two sampling profiles (macOS `sample`, container `perf record -e cpu-clock`), plus interleaved A/B wall-clock on real reads and deterministic instruction counts from callgrind

---

## 1. Scope and constraints

**What was reviewed.** The module `src/ng/run/cohort_merge/` — `mod.rs`, `serial.rs`, `parallel.rs`,
`observation_cache.rs`, `close.rs`, `build.rs`, `organise.rs` — against
[cohort_merge.md](../../ng/spec/cohort_merge.md) and [arch/cohort_merge.md](../../ng/arch/cohort_merge.md).
The module turns k samples' per-locus observations, each stream in coordinate order, into one stream
of cohort observations ready to be called, through three drivers that must produce byte-identical
output: the single-threaded oracle, the same merge read through the observation cache, and the
builders working a round of regions at a time on a rayon pool.

**Reviewed against** commit `f8b69436` on branch `ng-merge-perf`, which is `b0875294` (the
2026-08-19 per-sample keep rule) plus a probe-only change. Every measurement quoted below that
says *before* was taken on `f8b69436`.

**Targets.** One sample to several thousand, and three reads a position to several hundred
([design_principles.md](../../specs/design_principles.md) §0). Target hardware is the dev
container's aarch64 release build on 8 CPUs. There is no latency target; the merge is one stage of a
run whose upstream generator costs fourteen to twenty-three times as much
([cohort_merge.md](../../ng/spec/cohort_merge.md) §6.2), so what this review is worth to a whole
run is bounded by that ratio — see §6.

**Hot-path evidence available.** Both a CPU sampling profile and an off-CPU phase split, on both
platforms:

- macOS host, `sample <pid> 25 1` against a natively built `profiling` binary, attached once the
  reads were walked and only the merge remained. This is what gave the phase split.
- Linux container, `perf record -e cpu-clock -F 499 -g`, on the production architecture.
- Interleaved wall clock on the real tomato panel, both binaries kept side by side and alternated,
  through `examples/ng_cohort_merge_real_cost.rs`.
- `valgrind --tool=callgrind` instruction counts where an effect was smaller than this machine's
  wall-clock swing, which two reviewers measured at 13–26% between runs of one unchanged binary.

**What could not be measured, and it constrains one category.** This machine exposes no hardware
counters — the PMU is not virtualised in the Apple-`container` VM, and the macOS equivalent is
GUI-only. So there is no cache-miss count, no branch-miss count, no `perf c2c` and no `perf sched`
anywhere in this review. Layout findings were settled by wall-clock A/B and callgrind instead, and
one of them (§5, closed C3) could be shown not to pay without being able to say why.

**In-scope files.** `src/ng/run/cohort_merge/{mod,serial,parallel,observation_cache,close,build,organise}.rs`,
plus `src/ng/locus_generation/mod.rs` for the two accessors the walk calls at every record.

**Deliberately out of scope.** The generic locus generator itself; `SampleLocusObservations`'s
four-allocation shape, which is the largest remaining lever and belongs to the generator's spec
(§6); everything downstream of `CohortObservation`.

**Categories dispatched**, one reviewer each, every one in its own git worktree so that no two
benchmarks shared a target directory:

| category | why |
|---|---|
| `methodology` | always; and the allocator A/B and the benchmark's own honesty are its |
| `allocations` | the profile's largest single symbol group is the allocator |
| `data_layout` | the walk iterates 6.09 M records and the tournament chases pointers |
| `concurrency` | rayon, a round barrier, and eight threads buying 1.26× |
| `hot_loops` | the k-way walk is three quarters of the single-threaded merge |

`io_and_syscalls` was **not** dispatched: the merge performs no I/O. Its source is an in-memory
iterator, and the psp path that will read from a file is not built.

---

## 2. Verdict

**Apply the listed wins.** Seven changes were applied, each measured before and after on the real
tomato panel with the runs interleaved, each committed on its own with its numbers, and all 236 of
the module's tests green after every one — including the three that pin the drivers to byte-identical
output. Eight further candidates were built and measured and did **not** pay; they are recorded in
§5 with the numbers that closed them, so the next reviewer does not re-derive them.

**The two questions the review was asked, answered.**

*Why is the parallel driver on 8 threads slower than the single-threaded oracle?* It was not, once
the two are asked to do the same work. The oracle borrows the cohort's observations and frees
nothing while it is timed; the cached drivers own theirs and return 6.4 million blocks to the
allocator inside the clock, against the oracle's 130 thousand — counted with dhat, not estimated.
Timed like for like at 63 accessions the eight-thread driver was already the fastest of the three
(373 ms against the owning oracle's 430 and the cached serial driver's 466). What was true is that
threads bought only 1.26×, and the reason is in the next answer.

*Where does the parallel driver's time go?* Onto one thread. Per round of the merge, from the host
phase split: **182 ms freeing the records eviction dropped, 33 ms drawing the readers forward — both
on the organiser's thread with every builder idle — and 267 ms of processor time in the builders,
finishing in about 40 ms of wall.** The builders scale at about 6.7× on 8 threads; 84% of the merge
was never handed to them. The organiser's own work — resolving overlaps, releasing loci — is 4
samples of 23,204 and is not worth looking at.

**What the six changes are worth, end to end**, on the tomato panel, 63 accessions over 100 kb of
SL4.0, container release build, the merge timed alone (median of 12 rounds, two runs of each):

| driver | at the branch point | after the six changes |
|---|---|---|
| single-threaded oracle | 317.7, 315.5 ms | **155.3, 155.9 ms** (−51%) |
| one reader per sample, 1 thread | 497.2, 502.6 ms | **269.1, 272.1 ms** (−46%) |
| builders on 8 threads | 377.2, 376.2 ms | **253.5, 255.9 ms** (−32%) |

And **peak resident memory fell with it** — 5.66 GB to 5.12 GB on 8 threads, 5.45 to 5.09 on one —
because the allocator change is smaller as well as faster. Nothing in this review trades memory for
speed.

**The seventh change came after that table and changes what the table means**, so it is stated
apart from it (§5, H6). Every number above was taken with the cohort's records copied *before* the
clock started, which charges the merge for none of the cost of producing them. Timed with
production inside — which is what the psp path is, since there a record is decoded from a file
rather than handed over ready-made — the eight-thread driver is **2.2 times** the single-threaded
cached one, 223 ms against 499, where with production hoisted out it was 6%. **The threading earns
its keep on the path that matters and barely earns it on the one measured here.**

**Across the committed range**, which is the test a gain measured at 63 accessions and 11 reads a
sample has to pass:

| corner | before | after |
|---|---|---|
| **one accession**, tomato, oracle | 2.06, 2.07 ms | **1.44, 1.44 ms** (−30%) |
| **one accession**, tomato, 8 threads | 12.61, 12.67 ms | **9.95, 9.81 ms** (−22%) |
| GIAB trio, 3 samples at **313 compared reads**, oracle | 345.7, 348.7 ms | **81.2, 81.4 ms** (−77%) |
| GIAB trio, 8 threads | 105.8, 107.7 ms | **66.0, 65.9 ms** (−38%) |
| 1,000 samples, fabricated dense ground, 8 threads | 45.3 ms | **36.4 ms** (−20%) |
| 3,000 samples, fabricated dense ground, 8 threads | 165.9 ms | **129.1 ms** (−22%) |
| 1,000 samples, fabricated sparse ground, 16 in flight | 25.3 ms | **18.3 ms** (−28%) |
| 3,000 samples, fabricated sparse ground, 16 in flight | 88.3 ms | **62.5 ms** (−29%) |

The high-coverage corner gains most, and that was predicted rather than discovered: two of the six
changes are per-record costs that grow with the reads at a position, and GIAB carries 313 compared
reads a sample against tomato's 11. **The single low-coverage sample — the corner
[design_principles.md](../../specs/design_principles.md) §0 names as the hardest, and the one an
optimisation is most likely to quietly cost — gains 22 to 30%**, because the two changes that could
have hurt it were guarded: the deferred free is skipped when the pool has one thread, and the
parallel sweep splits nothing when the cohort is one sample. **No cell measured anywhere got
worse.**

---

## 3. Measurement plan

Most of this plan was executed during the review. What follows is what it left in the tree, and
what is still owed.

### What the probe now does, and why each part was needed

`examples/ng_cohort_merge_real_cost.rs` gained four things, all committed:

1. **`NG_REAL_ONLY=oracle|oracle_owned|cache|parallel` runs one driver by itself**, with
   `NG_REAL_ROUNDS`, `NG_REAL_WIDTH`, `NG_REAL_THREADS`. A sampling profiler attributes whatever the
   process was doing, and a run that walks 63 CRAMs and then sweeps five widths across three pool
   sizes gives a tree in which no driver is a majority of anything. It prints `# profile-start:`
   when the reads are behind it, so a wrapper can wait for that line and attach.
2. **The merge is timed apart from the copy that precedes it.** A cache consumes its readers, so
   every round after the first needs a fresh copy of the cohort; timing that with the merge charged
   the cached drivers for work the oracle is never asked to do.
3. **`oracle_owned`** merges the same loci from a copy the round owns and drops inside the clock.
   **It, and not `oracle`, is the baseline any cached-driver number should be quoted against.**
4. **The median with its fastest and slowest, and peak resident memory**, instead of a bare mean.
   Two runs of one unchanged binary differ by 13–26% here, which is larger than most of what this
   probe is asked to judge.

### The gate that settles a comparison, and it is not the clock

Under `--features dhat-heap` the probe reports blocks allocated, bytes allocated and **blocks freed
inside the merge**. These are identical on every run of the same code on the same input. Two
drivers whose freed counts differ by millions are not doing the same work whatever their
milliseconds say — which is exactly how the oracle/cache comparison had been wrong. Gate on the
count; use wall time as the check.

### What is still owed

1. **A cohort of thousands on real observations.** Nothing real exceeds 63 accessions, and the
   fabricated ground is a hundred times thinner than what the generator emits. The cheapest
   honest instrument is drawing a 250/1,000/3,000-sample cohort by replicating the 63 real
   accessions — real density, synthetic count.
2. **The cache's peak, measured rather than counted.** [cohort_merge.md](../../ng/spec/cohort_merge.md)
   §8 states 33 records held per sample at 1,000 samples and 16 regions in flight, and says what one
   record costs in bytes is unmeasured. The record is 120 bytes plus its four heap blocks (measured
   here), and at real density the held count is about 3,100 rather than 33. §8's table should be
   re-derived on real observations.
3. **The three-reads-a-position end.** Tomato at 11 compared reads a sample is the thinnest real
   data available; down-sampled CRAMs would reach 3.
4. **The whole run, end to end.** §6.2 puts the generator at fourteen to twenty-three times the
   merge using the cached driver's inflated time, so that ratio is a floor. What this review is
   worth to a run cannot be stated until the run exists.

---

## 4. Build / toolchain configuration

**Applied — the allocator, and it is the largest single change in this review.**
`alloc-mimalloc` is now a default feature (`1c15a80f`). The merge frees far more than it allocates —
6.4 M blocks a round against 216 k — because the observations it walks were allocated by the stage
upstream and are released as the merge passes them; glibc takes an arena lock per free and those
frees are concentrated on one thread. Measured at 63 accessions: 8 threads 333.4 and 335.0 ms →
249.4 and 258.0, one thread 450.6 and 452.9 → 270.1 and 269.8, with peak resident falling 5.66 →
5.12 GB and 5.45 → 5.09. The 9% off peak resident is the same figure the production `var-calling`
path recorded at 50 samples, so this is a second measurement agreeing with one already in the tree.
Cost: a vendored C allocator in every build, and a heap profile is now
`--no-default-features --features dhat-heap`.

**Closed — `target-cpu`.** `.cargo/config.toml` has no `aarch64 + linux` entry, so the production
build does not get ARMv8.1 atomics (`+lse`). Three alternating pairs measured −4.7%, −5.2% and
**+1.3%** — the sign flips and every value is inside this machine's swing. Not worth acting on. If
it is ever revisited, name `+lse` explicitly and never `native`.

**Unchanged and correct.** `[profile.release]` is `lto = "fat"`, `codegen-units = 1`,
`panic = "abort"`; `[profile.profiling]` exists precisely so a sampling profile can resolve
functions that fat LTO would otherwise inline away, and it is what both profiles in this review
used. `rust-toolchain.toml` pins 1.97.1 — note that
[profiling_environment.md](../../../../ai/skills/rust-performance-review/performance_review/profiling_environment.md)
still says 1.95.

---

## 5. Code-level findings

### Hot-path

#### H1: [observation_cache.rs:224](../../../../src/ng/run/cohort_merge/observation_cache.rs#L224), [parallel.rs:139](../../../../src/ng/run/cohort_merge/parallel.rs#L139) — eviction frees a round's dead records on the organiser's thread — **applied in `ab0c0fcd`**

- **Confidence:** High
- **Hot-path evidence:** host phase split, per round of the 8-thread merge: `evict_before` 182 ms of
  a 255 ms merge, nearly all of it `drop_glue::<SampleLocusObservations>`. Container `perf`: the
  glibc allocator family is ~39% of the process. Confirmed independently by two reviewers.
- **Mechanism:** deciding what to drop is one walk of a window's prefix; what costs is returning
  each record's four heap blocks. Those records are unreachable from every window the round hands
  out, so freeing them is disjoint from what the builders read — the one piece of the organiser's
  work that can run beside them. `evict_before_into` moves them to a buffer the driver owns and
  `rayon::join` gives that buffer to a worker to empty while the builders build.
- **Measured:** 63 accessions on 8 threads 367.5 / 371.3 / 372.2 ms → 358.5 / 362.9 / 360.3 (−3%);
  1 accession on 8 threads 12.44 / 12.35 → 10.42 / 10.56 (−15%); 63 accessions on 1 thread
  unchanged. Every overlapped run beat every non-overlapped one.
- **Why 3% and not the 16% the phase split predicts** — the free does not stop costing, it moves,
  and a worker freeing while the builders allocate contends with them for the same arena lock. This
  is the finding whose value changes most under H0's allocator; it should be re-measured there.
- **Complexity cost:** `S: Send` on the parallel driver, one buffer, and a branch on the pool size.
  On one thread `join`'s second closure runs inline, so the move would be pure loss — measured at
  +6% by the allocations reviewer, which is why the driver asks `rayon::current_num_threads()`.
- **Memory:** one round's evicted set stays resident until the worker finishes — peak graveyard
  25,595 records against peak cache 25,611, so the merge holds about twice what it did.
  `the_cache_holds_a_whole_round_and_not_one_region` still passes unchanged, but what it pins is now
  what the *cache* holds rather than what the *merge* holds.

#### H2: [close.rs:290](../../../../src/ng/run/cohort_merge/close.rs#L290) — a tournament node held an index into a second array — **applied in `ae8e9ac9`**

- **Confidence:** High
- **Hot-path evidence:** host profile puts `PendingHeads::replay_from` at 5,837 of 20,312 samples
  inside the oracle's merge (29%); container `perf` at 7.02% of the process; callgrind attributes
  273.7 M of the module's 656.8 M instructions to it (42%).
- **Mechanism:** the match at each level read `tree[node]`, a leaf number, and then `keys[sitting]`,
  the key that leaf named. The second load could not issue until the first returned, so a six-level
  climb at 63 accessions was twelve loads in six dependent pairs, half of them scattered. The node
  now holds the key. The key is one `u128` — contig in the top 32 bits, position in the middle 64,
  sample in the bottom 32 — so genome order is integer order and a match is one comparison instead
  of up to three.
- **Measured:** oracle at 63 accessions 315.2 / 316.2 / 315.8 ms → 262.3 / 260.6 / 261.6 (−17%);
  8 threads −1%, which is the expected shape rather than a disappointment, since the parallel
  driver's remaining time is the organiser's. Cohort axis on the walk-alone probe: −10% at one
  sample, −28% at 63, −26% at 1,000, −32% at 3,000.
- **Two reviewers reached this independently by different routes** — one packing to `u128`, one
  keeping three fields and replacing the sample with the leaf. The packed form was taken because it
  also collapses the comparison.
- **Complexity cost:** the tree is 16 bytes a node rather than 4 — 48 kB at 3,000 samples, L2-resident
  at every cohort size this caller commits to.

#### H3: [observation_cache.rs:335](../../../../src/ng/run/cohort_merge/observation_cache.rs#L335) — the readers were drawn forward by one thread — **applied in `f4e6a025`**

- **Confidence:** High
- **Hot-path evidence:** host phase split, 33 ms of a 255 ms merge, on the organiser's thread with
  every builder idle. The second of the two things the organiser does alone.
- **Mechanism:** covering a round sweeps the cohort repeatedly until a whole sweep moves the chain's
  reach no further. The sweep now runs the samples concurrently and takes the widest reach reported.
  **The fixpoint and the held window are identical**: drawing a sample is monotone in the reach it
  is given, the reach only grows, so both schedules climb to the same least fixpoint, after which
  every sample's last draw is against that same reach. What differs is the sweep count — the serial
  form lets sample *j* see the widening sample *i* &lt; *j* made inside the same sweep.
- **Measured:** 63 accessions on 8 threads 360.7 / 356.5 / 358.8 ms → 346.9 / 345.1 / 344.5 (−4%);
  1 accession unchanged. Every parallel-sweep run beat every serial one.
- **Complexity cost:** `E: Send`, and a schedule whose sweep count is data-dependent in a second way.

#### H4: [locus_generation/mod.rs:397](../../../../src/ng/locus_generation/mod.rs#L397) — one base was compared through a `memcmp` call — **applied in `da604df3`**

- **Confidence:** High
- **Hot-path evidence:** host profile puts `memcmp` at 1,614 of 20,312 samples inside the merge (8%);
  callgrind puts `bcmp` plus the inlined slice-compare frames at 6.6% of the module's instructions,
  and 1.6 M of the region's 8.75 M indirect branches.
- **Mechanism:** slice equality on `[u8]` compiles to a length test and a call whatever the length,
  and the generic mint writes one record per covered position — 1.00 reference bases a record on this
  panel — so the call's prologue cost more than the comparison. Naming the one-byte shape in a
  `match` lets the compiler compare two bytes in place; every other length falls through unchanged.
- **Measured with H5** (both are per-record costs in the same loop): oracle at 63 accessions
  261.4 / 261.0 / 261.8 ms → 248.3 / 249.0 / 248.9 (−5%). Callgrind attributes −1.6% of module
  instructions to this half alone.
- **One number cuts the other way:** callgrind's simulated conditional mispredicts rose 6.37 M →
  8.43 M, the in-place comparison being data-dependent where the call was not, and fell back to
  6.42 M once H5 halved the number of comparisons. This machine has no branch-miss counter, so that
  is comparative evidence between two variants and not an absolute rate.

#### H5: [close.rs:658](../../../../src/ng/run/cohort_merge/close.rs#L658) — a record's two read counts were asked separately — **applied in `da604df3`**

- **Confidence:** High
- **Hot-path evidence:** `LocusCloser::next` is 35% of the host profile and both accessors are inside
  it; callgrind puts the change at −2.24% of module instructions on its own.
- **Mechanism:** `non_reference_reads()` and `reads_compared_with_reference()` each built
  `complete_observations()`, filtered on the same witness and read the same `num_obs`. They are the
  numerator and the denominator of one question. `non_reference_and_compared_reads` answers both in
  one walk.
- **The orchestrator's framing of this was wrong about the size and right about the mechanism**, and
  the correction is the useful part: at 1.03 observations a record there is no iteration to halve.
  **What it removes is the second setup, and its value grows with depth rather than with cohort
  size** — which the GIAB corner then confirmed, where the whole stack is worth −77% against −51% on
  tomato.

#### H6: [observation_cache.rs](../../../../src/ng/run/cohort_merge/observation_cache.rs) — the merge freed records nobody needed freed — **applied in `05d9debb`**

- **Confidence:** High
- **Hot-path evidence:** the decisive experiment is an ablation rather than a profile. Making the
  merge **leak** the records it evicts instead of freeing them: 245 ms → 90 on 8 threads at 63
  accessions, and 265 → 204 on one. So freeing is **63% of the eight-thread merge** and 23% of the
  single-threaded one — after the allocator change, after the deferred free, after everything else
  in this review.
- **Mechanism:** the records were allocated by the stage upstream and are released as the merge
  passes them. A source is now a trait rather than an `Iterator`, and `next_observation` takes a
  record the merge has finished with; a source that mints records fills that one instead of
  allocating. The offer is not an obligation, and a blanket implementation over every `Iterator`
  ignores it, so nothing that exists today had to change.
- **Measured:** with a source that refills against one that allocates, 8 threads
  241.7 / 240.9 / 240.2 ms → 224.3 / 223.3 / 223.4 (−7%); one reader per sample
  531.6 / 544.6 / 538.0 → 497.8 / 500.2 / 498.9 (−7%). Peak resident unchanged.
- **Seven percent for 84% fewer allocator calls, and that ratio is the finding.** Counted with
  dhat at 3 accessions: 1,309,793 heap blocks allocated per round become 212,486. **A record costs
  4.5 heap blocks** — measured, not assumed. If calling the allocator were the cost this would have
  been worth far more; the leak says the class is worth 63%. So what costs is **touching four
  scattered pieces of memory per record, six million times**, which leasing keeps because it still
  copies each record's content into the buffers it reuses. §6 says what would remove it.
- **Complexity cost:** a public bound changes from `Iterator` to `ObservationSource`, and the
  window's fuse becomes a flag of the cache's own. Both are carried by the blanket implementation
  for every existing caller.

### Likely

#### L1: [observation_cache.rs:440](../../../../src/ng/run/cohort_merge/observation_cache.rs#L440) — the window's left edge was found by walking — **applied in `1a9419ba`**

- **Confidence:** High on the mechanism; the size depends entirely on the driver
- **Hot-path evidence:** container `perf` puts `first_reaching_index` at 4.15% of the parallel
  process and **absent from the serial one**. Reproduced by two reviewers.
- **Mechanism:** a sample's records are disjoint and ascending, so reach is monotone across the
  window and the predicate is false over a prefix — `partition_point`'s shape. The split between
  drivers is the eviction schedule, not the search: the cached serial driver evicts immediately
  before each cover so its window starts at the region's left edge, while the parallel driver evicts
  once per round, so a region late in the round walked past every earlier region's records.
- **Measured:** 8 threads at 63 accessions 344.4 / 344.5 / 345.7 ms → 336.4 / 337.2 / 336.5 (−2%);
  cached serial 462.7 / 466.2 → 452.2 / 452.4 (−3%). Every halving run beat every walking one.
- **Complexity cost, and it is the only one in this review that is not free:** the walk would have
  given the right answer on an unordered window, only slowly; the bisection gives a wrong one. No
  new assumption — `build_region` already refuses a sample whose records overlap — but the
  precondition now decides an answer rather than a speed, so it is stated where the search is.

#### L2: [close.rs:479](../../../../src/ng/run/cohort_merge/close.rs#L479) — eight cohort-sized allocations per building region — **open**

- **Confidence:** Medium
- **Hot-path evidence:** below the noise floor at 63 accessions; pattern plus arithmetic above it.
- **Mechanism:** `LocusCloser::over` allocates five vectors the length of the cohort,
  `PendingHeads::over` three more, and `ObservationCache::with_observations` one per builder call —
  once per building region, 500 of them per 100 kb at the default width. At 3,000 samples that is
  about 96 MB of zeroing per 100 kb of genome, whatever the region holds.
- **Measurement plan:** the fabricated cohort probe at 1,000 and 3,000 samples, with the scratch
  hoisted into a per-builder arena reused across regions. Threshold: 5% at 3,000 samples.
- **Complexity cost:** a scratch type owned by the builder and threaded through `build_region`,
  which is a signature change on the module's one public entry point.

### Speculative

#### S1: [close.rs:357](../../../../src/ng/run/cohort_merge/close.rs#L357) — replacing the tournament with a per-locus scan over the cohort

The tournament costs `log k` per record. A linear pass over samples per **locus** — not per
observation, which is what the old shape did and what lost — would cost one comparison per sample
per locus, and on this data almost every sample contributes to almost every locus, which is the
corner where a scan wins. It must close byte-identical loci, and the reason the walk is one pass is
that a deletion widens the locus after later samples were visited
([cohort_merge.md](../../ng/spec/cohort_merge.md) §11's first trap), so the honest shape is a hybrid
that keeps the tournament for sparse cohorts. Left unbuilt: H2 took 13.5% off the module without
touching the algorithm, and this is the one candidate in the review that could change the answer.

#### S2: [observation_cache.rs:328](../../../../src/ng/run/cohort_merge/observation_cache.rs#L328) — `draw_to` restarts at index 0 on every sweep of one cover

Within one cover the window only grows and the reach only grows, so re-reading a record already
inside the reach cannot move it; carrying the scan position across the sweeps of *one* cover is
provably redundant work removed, and is **not** the persistent mark the doc argues against. Built and
measured: a few ms in 115, under the fixture's noise floor. Reverted, recorded here because the
argument is sound and the cost may surface at a cohort size nothing real reaches.

#### S3: [close.rs:595](../../../../src/ng/run/cohort_merge/close.rs#L595) — `cursors_at_open` copies the whole cohort's cursors once per locus

97,408 loci × 63 samples on this fixture, and it grows with the cohort. Not separable from the
per-locus scan beside it without changing what a closed locus carries.

### Closed with a number — measured, did not pay

These were built and timed. They are recorded so that the next review does not spend the same days.

| # | candidate | what closed it |
|---|---|---|
| C1 | Spread eviction's frees across the pool (`par_iter_mut`) | **3× slower**: 138 ms → 402 / 407 / 422. `perf` says why: `__aarch64_cas8_rel` 0.72% → 32.9%, `_int_free` 5.3% → 36.1%. glibc keeps a chunk in the arena it was allocated in, so eight threads freeing one thread's records take one lock. Measured independently by two reviewers. **This is why H1 keeps the free on one thread and only moves which one.** |
| C2 | Skip the member list for the loci the verdict discards | 94,751 of 98,726 closed loci get a `Vec` filled and dropped unread, and removing it is worth −7.5% of module instructions and −2.3% wall *before* H1 — and **nothing after** (99.6 ms against 100.3). It also makes `ClosedLocus::members` mean something different per verdict, breaking four of `close.rs`'s own tests and the probe. Not worth a contract change for zero. |
| C3 | A dense `(start, reach)` side array per held window | `size_of::<SampleLocusObservations>()` is **120 bytes**, of which the cache's scans read 24 — the mechanism is real. Callgrind on the cache driver went **up** 0.4% (1,279 M → 1,285 M) and wall clock showed nothing at 200 or 1,000-base regions. Why it does not pay cannot be said here: it needs a cache-miss counter this machine does not have. |
| C4 | A branchless tournament match | That branch does carry 31% of the region's simulated mispredicts, but `Ord::min`/`max` cost **+34.6%** instructions and a `u128` select **+5.7%**, and the conditional-branch count barely moved in either — the compiler kept the branch both times. |
| C5 | Carry `draw_to`'s scan position across a cover's sweeps | S2 above: under the noise floor. |
| C6 | `target-cpu` with `+lse` | −4.7%, −5.2%, **+1.3%** — the sign flips. |
| C7 | Regions in flight at twice the thread count | Worth −14% on the unchanged code (122.0 → 100.5 ms at 8 threads), but flat once H1 overlaps the free, and it costs the same doubled memory H1 costs. An **alternative** to H1, not an addition. |
| C8 | Widen `rayon::join` to cover the next round too | `merge_ms` 122.6 → 50.5, but the enclosing phase rose 1.14 s → 1.78 s and the total moved 2.37 → 2.28 s. Mostly displacement into the next allocating phase, not a gain. Also close to [cohort_merge.md](../../ng/spec/cohort_merge.md) §6.2's already-rejected owned-windows arrangement. |

### Note

- [build.rs:102](../../../../src/ng/run/cohort_merge/build.rs#L102) — `LocusReferenceBases::over`
  gathers the locus reference **byte by byte with an assertion per byte**, so two members overlapping
  cannot disagree silently. It is the shape that blocks vectorisation, and it is 0.71% of the
  container run because it runs only for the 4% of loci that are built. **Do not turn it into a
  `copy_from_slice`** — the assertion's message names the failure it exists to stop, samples called
  against different references. If it ever becomes hot, copy first and compare the written span in a
  second pass, keeping the check and losing the per-byte branch.
- [build.rs:1400](../../../../src/ng/run/cohort_merge/build.rs#L1400) — the `ReadSighting` sort is
  987 of 20,312 host samples and 0.4–0.7% of the callgrind profiles. It is already `sort_unstable`
  on a 16-byte `Copy` key laid out in the order the algorithm after it consumes, and it is reached
  only where a sample has more than one record inside a built locus. Nothing to do.
- [organise.rs:22](../../../../src/ng/run/cohort_merge/organise.rs#L22) — the reorder `BTreeMap` is
  bounded by regions in flight rather than by cohort size, and the organiser is 4 samples of 23,204.
  Cold, and the cohort axis does not change it.
- **No false sharing to look for.** The builders hold no shared mutable state — each returns an owned
  `RegionOutcome` that rayon collects — so there are no atomics to contend on. Worth stating because
  it is the question this machine's missing `perf c2c` could not have answered.

---

## 6. Out-of-scope observations

- **A record costs 4.5 heap blocks, and that is the biggest lever left anywhere near this module —
  now with the number that sizes it.** `SampleLocusObservations` carries a boxed reference sequence
  and a `Vec<SequenceObservation>`, and each observation carries a boxed sequence and a
  `Vec<ChainId>`; dhat counts 1,309,793 blocks allocated for 289,581 records. **Of a 241 ms round,
  about 90 ms is the merge's own work and about 150 ms is making these records and taking them
  apart** — the 90 measured by leaking, the rest by difference. H6 shows that removing the
  allocator calls recovers only a seventh of it, so the cost is the scattered memory itself. Two
  shapes would remove it and both belong to the generator, not here:
  - **A flatter record**, with its bytes inside it and a heap fallback only for the rare wide one.
    Almost every record here is one base, one sequence and a few read identifiers — a few dozen
    bytes of content inside four allocations.
  - **A record covering a run of positions.** The generator writes one per covered base — 96,605
    per sample over 100,000 — and the merge discards 85,375 of the 97,408 loci it closes as too
    quiet. A record saying *"1,000 to 1,240 are all reference at depth 11"* replaces 240 with one.
    The keep rule needs each position's compared-read count, so a run must end where depth changes;
    depth changes slowly. **Unmeasured** — this is arithmetic from the record counts, not an
    experiment.
- **`examples/dhat_ng_merge.rs` profiles a different merge** — `SampleReads`' file merge, with zero
  mentions of `cohort_merge`. The cohort merge had no heap profile at all until this review wired
  counting into `ng_cohort_merge_real_cost.rs`. Either rename that example or point it at this
  module.
- **Three statements in the module's documents inherit the unfair baseline** and should be reworded
  where they are next touched: [cohort_merge.md](../../ng/spec/cohort_merge.md) §14 question 1's
  "62.6 ms against 16.6" is conservative rather than wrong and needs a clause; §6.2's "fourteen to
  twenty-three times the merge" is a **floor**, not an estimate; and §6.2's phase split was measured
  on fabricated ground where cover leads eviction, where the real profile puts eviction 5.5× ahead of
  cover. **Not** affected, and worth saying so because a blanket correction would sweep them in:
  `observation_cache.rs`'s module doc and the C1/C2 review numbers compare the oracle against
  *itself* at different region counts, with both sides borrowing.
- [profiling_environment.md](../../../../ai/skills/rust-performance-review/performance_review/profiling_environment.md)
  has drifted: `hyperfine` **is** on the host, the toolchain pin is 1.97.1 and not 1.95, `nproc`
  inside the container reports 9 on `--cpus 8`, and nothing in it warns that a reviewer working in a
  worktree has no `benchmarks/tomato1` CRAMs because they are gitignored.

---

## 7. What's already good

- **Three drivers that must agree byte for byte, pinned by property tests over random layouts at
  every region width and thread count.** That equivalence is what made this module optimisable at
  all: every change above was checked by running the same 236 tests, and a change that sped one
  driver while breaking the others would have failed within seconds rather than in a run months
  later. `serial.rs` and `parallel.rs` both carry it.
- **Scratch buffers that live for the whole walk, not the locus.** `LocusCloser`'s
  `alt_reads_per_sample`, `compared_reads_per_sample` and `cursors_at_open` are filled and returned
  to zero over exactly the samples that contributed, so a locus is never revisited and nothing is
  allocated per locus. `build.rs`'s `ReadAlleleScratch` does the same for the per-read allele walk.
- **Release-level assertions where a wrong answer would otherwise be silent**, each with a comment
  saying which silent failure it stops — the coordinate-order check in the cache, the
  never-mixed-kinds check in the walk, the disjointness check in the builder. Two of them were
  candidates for removal on speed grounds during this review and both survived reading their own
  justification, which is what those comments are for.

---

## Author response convention

Address each finding by its identifier — H1, L2, C4 — with one of: `applied in <commit>` /
`experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` /
`won't fix because …`. C1 through C8 are already closed with their numbers and need no response.
