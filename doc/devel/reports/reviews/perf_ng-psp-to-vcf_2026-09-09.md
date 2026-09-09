# Performance Review: ng psp → VCF, round two

**Date:** 2026-09-09
**Reviewer:** rust-performance-review skill (orchestrator), with four per-category reviewers in isolated worktrees
**Scope:** `call-from-psps` — a cohort of stored `.psp` files to a VCF, in one process
**Verdict:** Apply the listed wins — five were applied and measured; the whole command is 3.8× faster on about a tenth more memory
**Hot-path evidence:** a macOS `sample` profile of the real command on the real cohort broken down per thread, DHAT allocation counts, instrumented call counts, and interleaved A/B on every change

---

## 1. Scope and constraints

**What was reviewed.** The same path the 2026-09-07 review took:
[call_from_psps.rs](../../../../src/pop_var_caller_exp/call_from_psps.rs),
[calling_run.rs](../../../../src/pop_var_caller_exp/calling_run.rs),
[psp_caller.rs](../../../../src/ng/run/psp_caller.rs),
[psp_source.rs](../../../../src/ng/run/psp_source.rs),
[callers.rs](../../../../src/ng/run/callers.rs),
[src/ng/run/cohort_merge/](../../../../src/ng/run/cohort_merge/),
[src/ng/psp/](../../../../src/ng/psp/),
[src/ng/calling/](../../../../src/ng/calling/),
[src/ng/run/paralog_filter/](../../../../src/ng/run/paralog_filter/),
[src/ng/vcf/](../../../../src/ng/vcf/),
[src/ng/window_coverage/](../../../../src/ng/window_coverage/),
[reference_info.rs](../../../../src/ng/reference_info.rs).

**Reviewed against** commit `90551ef1` on `main`. Every change is on branch `ng-psp-vcf-perf-2`.

**Targets.** Wall clock and peak resident, over the range `design_principles.md` §0 commits to:
one sample to several thousand, three reads a position to several hundred. Every number here was
taken on one corner of that — the tomato benchmark's **63 accessions at about 3 to 26 reads a
position**, over the whole 8 Mb of SL4.0 that `benchmarks/tomato1/regions.bed` covers — except
where a cohort-size sweep is given, which runs 1, 4, 8, 16, 32 and 63.

**Hardware, and one thing about it that shaped the measurements.** macOS 15 on an **Apple M5 Pro:
6 "Super" cores and 12 "Performance" cores**, 18 logical, 64 GB. Native host release build, fat
LTO, one codegen unit, `-C target-cpu=apple-m1`, mimalloc, toolchain 1.98.0. **Two core tiers mean
a wall-clock figure depends on which cores a run landed on**, which is a fact about this
development laptop and not about the target: ng is meant for Linux servers with many uniform fast
cores, and **nothing here is tuned to this machine** (owner's instruction, 2026-09-09). Where a
change is small, the number quoted is **instructions retired**, which does not depend on the core
tier and repeats to within 0.03% between runs of one binary.

**Deliberately out of scope:** `src/var_calling/`, `src/pileup/`, `src/ssr/`, `src/psp/` — the
frozen production caller, which ng must not edit; the `generate-psps`, `generate-census` and
`estimate-parameters` commands; [locus_score.rs](../../../../src/ng/paralog/locus_score.rs), a
guarded verbatim copy of a production file that `copy_fidelity.rs` fails on any edit to.

**Categories dispatched**, one reviewer each, every one in its own git worktree:

| category | why |
|---|---|
| `concurrency` | 1.57 of 18 cores busy during a pass that is 90% of the run |
| `allocations` | 21.6 million allocations and 3.06 GB of churn in the calling pass |
| `hot_loops` | `lgamma` and `log` were 18% of the one busy thread |
| `io_and_syscalls` + `data_layout` | 2.5 s of startup on every re-run, and the memory slope against cohort size |
| `methodology` | run by the orchestrator; the build configuration was reviewed on 2026-09-07 and is unchanged |

**The fan-out worked and the isolation earned its place.** All four reviewers were benchmarking at
once, and the same unmodified binary ran the 63-accession cohort in anywhere from 87 s to 251 s
across the session. Every reviewer said so at the top of its file and led with a deterministic
measure. **Two of the five adopted changes could not have been separated on wall clock at all**
and were adopted on counts.

---

## 2. Verdict

**Apply the listed wins.** Five are applied, measured and committed. On 63 accessions over 8 Mb,
interleaved on an idle machine:

| | wall clock | calling pass | instructions retired | peak resident |
|---|---|---|---|---|
| `90551ef1` | 87.34 / 88.37 s | 78.89 / 79.91 s | 3,240.4 / 3,240.2 G | 527 / 530 MB |
| `b7042ddf` | 23.17 / 23.27 s | 14.01 / 14.09 s | 2,929.1 / 2,929.1 G | 607 / 590 MB |

**3.8× the speed on about a tenth more memory, 9.6% fewer instructions, and the VCF body is
identical** — 199,641 records, byte for byte, at 1, 4 and 18 threads.

**The calling pass was 90% of the run and one thread; it is now 67% of a run less than a third as
long.** Cores busy during it went from 1.57 of 18 — the calling thread at 96% and all eighteen
workers together at 0.61 — to 2.44 measured under comparable load in the reviewer's prototype, and
the thread sweep is still climbing at eighteen.

Where the time goes now, at 63 accessions: **calling 14.1 s (67%), the hidden-duplication filter's
scoring 6.4 s (31%), startup 2.3 s, writing 0.45 s.** Scoring was 7% of the run before this branch
and is 31% of it now; startup was 3% and is 10%. Both are priced in §5 and neither is fixed here.

The cohort-size sweep, wall clock, every arm's VCF identical to the streaming driver's:

| accessions | 1 | 4 | 8 | 16 | 32 | 63 |
|---|---|---|---|---|---|---|
| `90551ef1` | 3.52 s | 4.17 s | 6.13 s | 10.91 s | 29.49 s | 83.24 s |
| this branch | 3.46 s | 3.91 s | 4.65 s | 6.16 s | 11.21 s | 24.31 s |
| peak resident | 185→209 MB | 221→236 | 247→277 | 325→355 | 490→510 | 541→577 |

**The speedup grows with the cohort and the memory overhead does not** — it stays between 4% and
13% at every size, because the round is sized to hold the ground one region held before.

**One correctness defect was found and fixed on the way**, in a driver no production caller
reached; it is H2 and it has been routed as a correctness matter rather than folded into a
performance claim.

---

## 3. Measurement plan

The workload and the harnesses are reusable and are in
`tmp/perf_review_2026-09-09_ng-psp-to-vcf/`, described in its `EVIDENCE.md`.

1. **The fixture already exists** — 63 tomato accessions' psps over the whole 8 Mb region set, in
   `tmp/perf_2026-09-09/tomato/ng` (2.6 GB), with a fitted parameters file beside it. Subsets of
   1, 2, 4, 8, 16 and 32 samples are symlink directories `psps<n>` in the review directory.
2. **Wall clock and peak resident** — `run.sh <tag> <binary> [threads]` (63 accessions),
   `run8.sh` (8), `runN.sh <tag> <n> <binary>` (any subset), all through `/usr/bin/time -l`, arms
   interleaved in one sitting, full ranges quoted.
3. **Instructions retired**, from the same `/usr/bin/time -l`. **This is the measure to lead with
   on this machine**: it repeats to 0.03% within an arm, it does not depend on which core tier a
   thread landed on, and it separated two changes that wall clock could not.
4. **Allocation counts** — `cargo run --release --example ng_call_from_psps_cost
   --no-default-features --features dhat-heap <reference> <catalog> <psp-dir>`. Deterministic to
   the block. **It stops at the record handover**, so the VCF encoder and the three paralog-filter
   passes are outside every count it gives.
5. **CPU profile** — macOS `sample <pid> <secs> 1`, attached once the calling pass is under way.
   Wait primitives are leaves, so their counts are exactly the parked samples; subtract to get
   cores busy. **Attach late enough**: at 15 s in on a loaded machine the run was still reading
   the reference, and that profile is about startup rather than the calling pass.
6. **Correctness gate on every change** — the VCF body identical (`diff` after filtering `^##`),
   and not dependent on the thread count.

**What this machine cannot measure**: no hardware counters, so no cache-miss count, no
branch-miss count, no `perf c2c`, no `perf sched`. Layout questions were settled by `size_of`,
peak resident and wall-clock A/B.

**A trap worth recording, because it cost an hour.** The `sample` profiler put **6,996 of the main
thread's 36,867 samples — 19% of the one busy thread — in a `_platform_memmove` whose parent frame
was `LiveSet::from_sorted_slice`'s `to_vec`.** Removing that copy entirely changed instructions
retired by 0.1% and cycles by 0.9%, both inside the run-to-run spread. DHAT then priced the same
site at 1,091,609 allocations and 60 MB over the whole pass, which cannot be seconds. **The
attribution was wrong.** Treat a `sample` leaf in a system library as naming a *cost*, not a
*caller*, and confirm the caller with a second instrument before acting on it.

---

## 4. Build / toolchain configuration

**Unchanged and still sound** — reviewed on 2026-09-07: fat LTO, one codegen unit,
`panic = "abort"`, `debug = "line-tables-only"`, `x86-64-v3` on Linux x86_64 and `apple-m1` on
macOS aarch64, mimalloc. The one gap that review recorded — no `[target]` entry for **aarch64
Linux**, the production target, so the container builds at the generic ARMv8-A baseline — is still
open, and is still an owner's decision about the oldest ARM server the project must run on rather
than a measurement.

---

## 5. Code-level findings

### Hot-path

#### H1: [serial.rs:253](../../../../src/ng/run/cohort_merge/serial.rs#L253), [callers.rs:889](../../../../src/ng/run/callers.rs#L889) — the calling pass genotypes one locus at a time on the calling thread

**Confidence:** High. **Applied in `2137b8b3`.** Found and prototyped by the concurrency reviewer.

A 45-second `sample` over the calling pass at 63 accessions: `__psynch_cvwait` is **642,371 of the
700,473 samples across 19 threads, 91.7%**. Busy CPU is **1.57 cores of 18** — the calling thread
at 96.2% of its own samples, all eighteen rayon workers together at 0.61 between them. The
`--threads` help said as much in as many words: *"what the threads parallelise is the reading, and
everything after — building the loci, calling them, writing the VCF — stays on one thread in
genome order."*

Building a locus and genotyping it are pure functions of that locus.
[rounds.rs](../../../../src/ng/run/cohort_merge/rounds.rs) takes
[parallel.rs](../../../../src/ng/run/cohort_merge/parallel.rs)'s round structure — evict, one cover
for the whole round, the round's regions worked at once — and lets the worker carry each locus past
the cohort observation into the caller's own closure, so **a round never materialises a
`Vec<CohortObservation>`** and only the finished records come back. The sink stays on the calling
thread and still sees genome order, because `par_iter` over a round's regions is *indexed* — so
`collect` restores region order whatever order the workers finished in — and rounds are taken in
genome order.

Of the thirteen things the per-locus closure touched that were not per-locus, nine are per-thread
scratch or order-independent folds, one (`padding_reference`) is inert on today's mint, and three
are genuine ordering constraints and all three live at the sink.

**The memory is a tenth and not a multiple because the round divides the ground rather than
multiplying it.** A round holds `regions_in_flight × cohort_locus_builder_regions_len` bases, so the
two knobs are one lever: `round_shape_for` takes the width one region held before and cuts it up —
**fifteen regions of 529 bases at 63 accessions where the streaming driver held one of 7,936.**
Measured at 63 samples with the width left at its default, each extra region in flight costs
175–190 MB; holding the product constant holds peak resident to within 13% at every cohort size
from 1 to 63.

**Above about a thousand samples the rule hands back one region and today's behaviour**, because
`round_width_for` is already at its 500-base floor there and there is no ground left to divide. A
large cohort keeps its memory and gives up the parallelism rather than the other way round. **That
end is unmeasured** and it is the first thing to establish before this reaches a cohort of that
size; `examples/ng_cohort_merge_parallel_cost.rs` is the existing synthetic that walks it.

**Complexity cost:** a new 420-line driver beside `parallel.rs`, whose round structure it repeats;
one new field on `MergeParameters`; `Source: Sync` and `G: Sync` on the psp caller; the frontier
rule reimplemented in fifteen lines rather than through `Organiser`, because the organiser's
reorder buffer has nothing to do when the iterator is indexed. Four tests of its own, described in
H2.

#### H2: [observation_cache.rs:1132](../../../../src/ng/run/cohort_merge/observation_cache.rs#L1132) — `merge_cohort_in_parallel` cannot produce a correct answer over stored psps

**Confidence:** High. **A correctness defect, in a driver no production caller reaches.** Fixed on
this branch as part of `2137b8b3`; **it should be looked at by a correctness review rather than
signed off here.** Found by the concurrency reviewer.

[`ClosedLocusRanges::resolved_into`](../../../../src/ng/run/cohort_merge/close.rs#L113) numbers a
locus's members by their position in the **window** `with_observations` handed out;
[`build_at`](../../../../src/ng/run/cohort_merge/observation_cache.rs#L1138) indexes the sample's
**whole held list**. The two agree only where that window starts at the head of the list, which is
what a driver evicting at each region's own first base gets for free and one evicting once for a
round of regions does not. So **every region after a round's first would build its members from the
region before it.** It surfaced as an abort naming two positions exactly 16,000 apart — one
building region at that cohort size.

**Why the oracle battery was green on it.** Direct mode never calls `build_at`: when a source hands
its records over whole, `with_observations` gives the window the records themselves. **Every
fixture in `cohort_merge` was that shape**, so the byte-for-byte comparison between the serial and
parallel drivers — including its two hundred random layouts — could not reach the defect.

`fixtures::KeepingSource` is the missing shape: a source that hands out summaries and keeps its
records, the way a psp reader does. `rounds::tests` compares the round driver against the streaming
one through it, at four widths and four round sizes. **Removing the fix fails three of the four
tests**, each with the original abort; the fourth, `one_region_in_flight_is_one_region_at_a_time`,
correctly stays green, because with one region there is no prefix to forget.

**The fix applied is the cheap one** — `window_starts_into`, one `partition_point` per sample per
region, added back to the index. Two better ones exist and both are larger: have
`with_observations` hand the closer absolute indices, or have the window carry body handles rather
than indices. **The third makes the class of bug unrepresentable** and is what a correctness review
should weigh.

#### H3: [close.rs:968](../../../../src/ng/run/cohort_merge/close.rs#L968) — the walk mints a member vector for every locus it closes, and 98 in 100 are thrown away

**Confidence:** High. **Applied in `d53a9c0d`.** Found and measured by the allocations reviewer.

DHAT over the calling pass, 8 accessions, 8 Mb: **7,871,163 of the pass's 21,573,752 allocations
(36.5%) and 1,439,667,624 of its 3,055,153,920 bytes (47.1%)** are one `Vec::with_capacity` — six
times the next site. The walk closes 7,871,163 loci and the builder assembles 127,804 of them, so
**1.6 loci in 100 are ever built** and the rest take a heap block for a vector nothing reads. The
vector is as wide as the covering cohort, so the bytes grow with cohort size.

`LocusCloser` now keeps one spare and lends it; `recycle` takes it back. A caller that does not
hand it back is still correct — `next` clears and reserves either way — so this is a loan and not a
contract.

**Two more from the same reviewer went in beside it**: `resolved_into` no longer collects a vector
of member spans it walks once, since a member contributes exactly its own length and the spans are
the running sum; and a sample's supported-allele rows are reserved from the tally that bounds them
rather than collected from a `Filter`, which reports a lower bound of zero and grows 0 → 4 where
one or two are used. Together:

| | allocations | bytes |
|---|---:|---:|
| before | 21,573,752 | 3.06 GB |
| after | 13,594,358 | 1.46 GB |
| | **−37.0%** | **−52.2%** |

**What it is worth on the clock is under one percent**, and that is the honest figure: interleaved
at 8 accessions on an idle machine, instructions retired 275.4 / 275.7 / 276.1 G against 272.5 /
273.7 / 274.7 G and the calling timer 3.68 / 3.76 / 3.75 s against 3.63 / 3.71 / 3.70 s — every
pair in the same direction, about 1%. At 63 accessions it is 0.28%, because the site's count grows
with the ground while the run's total work grows with the loci found.

**The reviewer's patch was reshaped before it landed.** It needed the caller to hand the vector back
on four exit paths, with an `Option` threaded past the one arm that consumed the locus. That arm
consumed it for no reason — `resolved_against` took the locus by value but only ever read it — so
making it borrow collapsed the protocol to one line at the end of the loop body **and removed four
`.clone()` calls at test call sites**.

**Complexity cost:** one field, one public method, and a `for` that became a `while let`. The
failure mode is benign: `recycle` only ever keeps a buffer, never contents.

#### H4: [quality/mod.rs:418](../../../../src/ng/calling/quality/mod.rs#L418) — the site quality's count prior is rebuilt from `lgamma` at every locus and is the same vector every time

**Confidence:** High. **Applied in `03a3b801`.** Found by the hot-loops reviewer.

Every term of the Beta-Binomial prior over the cohort's allele count is a function of the fitted
spectrum's two concentrations and the cohort's chromosome count, and of nothing the locus carries.
It costs `4 × (2N + 1) + 5` `lgamma` a locus at N samples: **73 at eight, 513 at sixty-three,
24,005 at three thousand** — the top of the committed range. Over 8 accessions and 8 Mb that is
9,329,692 of the calling pass's 46,278,940 `lgamma` calls (instrumented count).

Measured at 63 accessions, two interleaved rounds: **instructions retired 3,243.2 G → 3,138.6 G,
−3.2%; cycles 938.9 G → 908.2 G, −3.3%.** Bit-identical by construction — the same expressions on
the same operands in the same order, evaluated once.

**The reviewer's patch used a `thread_local!` and the shipped form does not.** The memo lives in
`CallingScratch` beside the other four buffers this function already borrows from it, keyed by what
it was built from — the chromosome count and the two concentrations, compared by bit pattern — so a
caller handing the next locus a different spectrum gets a rebuild rather than a stale row.

#### H5: [dirichlet_multinomial.rs:127](../../../../src/ng/calling/genotype_prior/dirichlet_multinomial.rs#L127) — every genotype's prior pays two `lgamma` where a rising product needs one `ln`

**Confidence:** High. **Applied in `b7042ddf`, after an owner's ruling.** Found and measured by the
hot-loops reviewer.

For an integer copy count `k`, `lgamma(α + k) − lgamma(α)` is exactly the logarithm of the rising
product `α(α+1)…(α+k−1)`. **36,949,248 of the calling pass's 46,278,940 `lgamma` calls came from
this one file.** At diploidy the identity is one `ln` of one multiply instead of two `lgamma`, and
an allele the genotype carries once needs no transcendental at all.

Measured at 63 accessions, three interleaved rounds: **instructions retired 3,141.1 / 3,141.3 /
3,141.3 G against 2,929.5 / 2,929.9 / 2,930.4 G — −6.7%**, and additive with H4 to within a third
of a percent, as the disjoint call counts predict.

**The difference form is the less accurate of the two, by an amount that grows with the cohort.**
It subtracts numbers of order `α·ln α`, and `α` is the leave-one-out concentration — the run's seed
plus the cohort's expected allele copies. Against an exact reference: at a concentration of 3 the
two agree; at 6,001 — three thousand diploid samples — the difference form is hundreds to thousands
of units in the last place adrift where the product form is at most one.

**What it cost was a test, and the owner ruled on it (2026-09-09): an approximation is acceptable,
flipped genotypes are not.** `the_port_matches_production_bit_for_bit` asserted that ng's prior
matched the production primitive it was ported from bit for bit; the two spellings differ in their
last bits, so keeping it meant keeping the slower and less accurate form. Two tests replaced it:
one holding ng and production to within **1e-9 nats** over the same grid — an absolute bound,
because what it bounds is production's cancellation, whose size is set by `α·ln α` and not by the
row entry — and one comparing both against an exact reference over the 23 grid cells where the
rising product is a whole number under 2^53, in which **ours is within one unit in the last place
everywhere and is never the further of the two from exact**.

**Nothing moved in the output:** the VCF body is byte-identical at 1, 8, 32 and 63 accessions —
434,954 records in all, and 12,577,383 genotype calls at 63 alone.

### Likely

#### L1: [reference_info.rs:748](../../../../src/ng/reference_info.rs#L748) — the reference verification uses two of eighteen cores where the `.fai` makes it splittable

**Confidence:** High on the size, Medium on whether the split is worth it. **Not applied.** Found
and measured by the io-and-layout reviewer.

Startup, stopwatched between every step: **reading and verifying the reference is 2.2–4.4 s of a
2.4–3.6 s startup — 93% of it.** Opening the psps is 0.6 ms a file (32 files in 7–20 ms) and
`segments_over` is 0.15–0.25 s and flat in cohort size. Nothing else reaches a millisecond.

It is CPU-bound, not I/O-bound: reading the 795 MB file warm costs 56 ms and `md5` over all of it
costs 1.25 s. Split three ways with an env switch: **the scan alone is 1.35 s, the two digests add
0.97 s through the shipped `rayon::join`, and 2.24 s when forced serial** — so the side-by-side
hashing the 2026-09-07 review added is working and is saving about 1.3 s.

**This matters more after H1 than before it.** Startup was 3% of the run and is now 10%. A
contig-parallel scan — every byte range is known, because the `.fai` is read before the thread is
spawned — would take the scan term to roughly the largest chromosome's share and leave the
whole-reference MD5 at 1.23 s as the only serial term: **ideal wall about 1.25 s against today's
2.3 s.**

**Complexity cost, and it is why this is not applied:** `reference_info.rs` is 2,720 lines and the
pass carries geometry validation with line-numbered errors, the `.fai` write, and the duplicate-name
check. A contig-parallel version needs each worker to own a handle and a range, and the errors must
keep naming the right line of the right contig.

**It is also an alternative to, not a complement of, the 2026-09-07 review's L1**, which is still
open and still the owner's: deferring the reference *digest* past the calling pass removes the whole
term at the cost of reporting a wrong reference at the end rather than the start. Re-measured on
today's tree: **20–30% of an 8-accession run, 5% at 32, about 2% at 63** — worth most exactly where
the cohort is smallest.

#### L2: [psp_source.rs:259](../../../../src/ng/run/psp_source.rs#L259) — the per-record live-read snapshot is the largest identified term of peak resident

**Confidence:** High on the size. **Not applied.** Found and measured by the io-and-layout reviewer.

`live_ids` is **94.4 MB of a 494–520 MB peak at 32 samples — 18–19%**, and larger than the other
four terms of the merge window together. Over that run 2,891,498,467 identifiers passed through for
240,979,720 heads; among the records held at peak it is 23.3 a record, which is the top of the 3–26
reads a position this fixture carries.

**It conflicts with `e6db04bc` and must be taken after it, not instead of it.** That commit made the
body decoder read the arena span directly rather than copy it into a `LiveSet`; narrowing or
re-deriving the arena would put a copy back at exactly that call site.

#### L3: [observation_cache.rs:713](../../../../src/ng/run/cohort_merge/observation_cache.rs#L713) — the window lookup binary-searches a 470 kB array from the middle for a query that only moves forward

**Confidence:** Medium. **Not applied.** Found and measured by the hot-loops reviewer.

984,904 calls searching slices of 14,796 entries on average — 13,700,191 probes on 8 accessions.
Bit-identical. Measured at 63 accessions: **cycles 917.5 / 927.0 G against 906.5 / 907.2 G, about
−1.6%, while instructions retired move 0.15%** — instructions flat and cycles down is the signature
of a locality win. Nothing at 8 accessions, where the arrays fit in L2.

**Not applied because the shipping form costs a parameter on `CohortObservation::over` and
`window_coverage_at` and about thirty test call sites, for 1.6%.** It grows with the cohort — the
working set is `samples × 470 kB` — so it is worth re-measuring at the large-cohort end rather than
closing.

#### L4: [record.rs:997](../../../../src/ng/psp/record.rs#L997) — every built record mints four buffers, and both ends of the recycling hook already exist

**Confidence:** Medium. **Not applied.** Found and measured by the allocations reviewer.

After H3, **4,674,412 allocations — 34% of what is left** — are the four buffers
`decode_record_body` mints per built record. The hook is already there on both ends and discarded:
`ObservationSource::next_observation` takes a `spare` that `PspSummarySource::next_drawn` opens by
dropping, and `resolved_into`'s `scratch.clear()` frees every inner buffer.

**What blocks it is a type change that leaves this scope**: `reference_bases` and
`SequenceObservation::bases` are `Box<[u8]>`, which cannot be refilled in place. As `Vec<u8>` they
could be, and `decode_record_body` would write into a record handed back to it.

### Speculative

- **S1: [summarise_condition.rs:1374](../../../../src/ng/calling/inference/summarise_condition.rs#L1374) and [assemble.rs:206](../../../../src/ng/vcf/assemble.rs#L206)** — a called genotype is two four-byte heap blocks, one to mint it and one to clone it into the record: **1,685,175 allocations, 12.4% of what is left after H3.** `Genotype` is `Genotype(Box<[AlleleId]>)`; inline storage for the ploidies this caller meets removes both. The type is in `src/ng/types.rs`, outside the reviewer's scope.
- **S2: [records.rs:178](../../../../src/ng/run/records.rs#L178)** — every run sample gets its own eight-byte zeroed read-count vector for every written record: 686,784 allocations, 5%. One flat vector indexed by stride is one allocation a record. Costs the struct a lifetime or an offset pair.
- **S3: [build.rs:528](../../../../src/ng/run/cohort_merge/build.rs#L528)** — each distinct allele's bytes are boxed twice, once for the table and once for the key that unifies it: 844,875 allocations, 6.2%. The module's own doc already defers this to milestone B3.
- **S4: [rounds.rs](../../../../src/ng/run/cohort_merge/rounds.rs)** — the round driver allocates a fresh `Vec<WindowCoverage>` and a `Box<VcfRecord>` per written record, 199,640 times, where the streaming driver reused one buffer. It did not show against the 3.8×; it is the obvious next allocation item on the new path.

### Note

- **A `sample` leaf can name the wrong caller, and it did — see §3.** 19% of the one busy thread was attributed to a copy that DHAT prices at 60 MB over the whole pass. Recorded so nobody re-files it.
- **The tournament's conditional swap is 57.9% taken over 180 million matches and rewriting it branchless is measurably slower** — 5.28 / 5.17 / 5.26 / 6.19 s against 4.88–5.01 s over four interleaved pairs, instructions unchanged. The key is a `u128`, so the branchless form pays two conditional selects on each half plus an unconditional store against a branch the hardware evidently predicts better than 57.9% suggests.
- **The peak-resident slope against cohort size is two numbers, not one**: 10.0 MB a sample from 8 to 32, then **0.88 MB a sample from 32 to 63**, over a 184 MB floor at one sample. The flattening is `round_width_for` narrowing the builder region as samples are added — held records a sample go 16,210 at n=8 to 8,157 at n=63. **The builder region is the memory governor and it is working.**
- **Nothing in the merge's types wastes padding.** `KeptRecord` is exactly 64 bytes — one cache line; `LocusSummary` exactly 32; `MemberRange` and `GenomeRegion` 24 with no hole. No field-reordering finding.
- **`held_observations` is zero at every cohort size**, which is the psp path's deferred build working: the two largest types are never held, and one body is built for every 8.8 record heads walked at 32 samples.
- **The 184 MB floor is the catalog parquet read's transient decode**, not anything the run keeps — what survives is 0.7 MB of segmentation. An `mi_collect` after the read returns 29 MB at one sample and nothing at 8 or 32, because the merge reuses those pages.
- **Returning the arenas' reserved-but-empty capacity makes things much worse**: `shrink_to_fit` took an 8-sample run from 258 MB to **644 MB**. The growth ratchet is doing useful work; only its frequency has been shown to be a problem.
- **The paralog filter's whole contribution to peak resident is 10 MB at 8 samples and 5 MB at 32**, measured with `--paralog-fdr 0`. It is not a memory problem.
- **The parquet catalog reader is already the right shape** — row groups pruned by the contig column's statistics, three columns projected for the row filter, the region predicate pushed down.
- **Opening the psps is not a startup cost and should not be parallelised**: 0.6 ms a file, and the loop is serial so a refusal names the first file in cohort order, which is worth more than 20 ms.
- **The read-group tables are already dense slices indexed by id**, not `BTreeMap` lookups, with the conversion outside the loop.
- **Hoisting the error-spread `ln` in the genotype loop loses at biallelic diploid** — the loop takes the error branch for one genotype in three, so hoisting two logarithms replaces one. It becomes a win at four alleles and up.
- **`fold_samples_into_allele_counts` is quadratic in the cohort** — about `3N²` multiply-adds a locus, 12,000 at 63 samples and 27 million at 3,000. The reviewer could not get `cargo asm` on it under load, so whether the multiply-add loop autovectorises is unestablished. It is the first thing to check at the large-cohort end.

---

## 6. Out-of-scope observations

- **The scoring pass is now 31% of the run** (6.4 s of 20.97 s after startup), where it was 7%. It was parallelised by the 2026-09-07 review and is not obviously wrong; it simply has not been looked at since it became the second-largest term.
- **`merge_cohort_in_parallel` is still there, still has no production caller, and now has H2 against it.** It is either superseded by `rounds.rs` or it should be fixed and tested through a deferring source; leaving it as a third driver that cannot serve one of the two modes is the worst of the three options.
- **Direct mode (`call-from-alignments`) does not use the round driver.** `MergeParameters` carries the count for both, and only psp mode reads it. Whether the same arrangement helps a run that decodes CRAMs is a separate measurement, and its bottleneck is known to be elsewhere.
- **The integration test `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele` still fails on `main`**, as the 2026-09-07 review recorded, and three examples still do not compile. Neither is this branch's.
- **`cargo test --release --lib` reports 8 failures on an unmodified tree** — all `#[should_panic]` tests waiting on debug assertions a release build does not compile in. 6,653 pass on this branch and the same 8 fail.

---

## 7. What's already good

- **The builder region's width is already chosen from the cohort's size** ([calling_run.rs:401](../../../../src/pop_var_caller_exp/calling_run.rs#L401)), which is what made H1's memory neutrality a matter of dividing an existing budget rather than inventing one. The measured slope — 10.0 MB a sample to 32, then 0.88 to 63 — is that governor working.
- **`par_iter` over a round's regions is indexed**, so ordering the output costs nothing ([parallel.rs](../../../../src/ng/run/cohort_merge/parallel.rs)'s own `in_region_order` says so). The round driver could keep the sink serial and in genome order without a reorder buffer because of it.
- **The run reports where its own time went**, split into calling, fitting, scoring and writing ([finish.rs:187](../../../../src/ng/run/paralog_filter/finish.rs#L187)). On a machine whose wall clock swung 3× under load, that line was the only wall-time measurement worth reading.

---

## Author response convention

`H1` applied in `2137b8b3`. `H2` fixed within `2137b8b3` and **routed to correctness review** —
the fix applied is the cheap one and the reviewer named two better ones. `H3` applied in
`d53a9c0d`. `H4` applied in `03a3b801`. `H5` applied in `b7042ddf` after the owner's ruling on the
bit-parity test.

`L1` open — and it is an alternative to the 2026-09-07 review's L1, which is also still open and
still the owner's. `L2` open, and must follow `e6db04bc` rather than replace it. `L3` open,
re-measure at the large-cohort end. `L4` open, and it is the largest allocation item left.
`S1`–`S4` open, none worth its complexity at this size.

**The one thing this review did not establish, and it is the largest gap:** what H1 does at a
thousand samples and at three thousand. The rule degrades to today's behaviour there by
construction, so the risk is a lost opportunity rather than a regression — but "degrades to today's
behaviour" is arithmetic from the width rule, not a measurement.
