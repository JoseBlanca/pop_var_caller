# Performance Review: ng psp → VCF

**Date:** 2026-09-07, finished 2026-09-08
**Reviewer:** rust-performance-review skill (orchestrator), with four per-category reviewers in isolated worktrees
**Scope:** `call-from-psps` — a cohort of stored `.psp` files to a VCF, in one process
**Verdict:** Apply the listed wins — six were applied and measured; the whole command is 1.96× faster on 57× less memory
**Hot-path evidence:** two macOS `sample` profiles of the real command on the real cohort, DHAT allocation counts, a resident-set trace, and interleaved wall-clock A/B on every change

---

## 1. Scope and constraints

**What was reviewed.** The path `pop_var_caller_exp call-from-psps` takes:
[call_from_psps.rs](../../../../src/pop_var_caller_exp/call_from_psps.rs),
[psp_caller.rs](../../../../src/ng/run/psp_caller.rs),
[psp_source.rs](../../../../src/ng/run/psp_source.rs),
[callers.rs](../../../../src/ng/run/callers.rs),
[src/ng/run/cohort_merge/](../../../../src/ng/run/cohort_merge/),
[src/ng/psp/](../../../../src/ng/psp/),
[src/ng/calling/](../../../../src/ng/calling/),
[src/ng/run/paralog_filter/](../../../../src/ng/run/paralog_filter/),
[src/ng/vcf/](../../../../src/ng/vcf/),
[src/ng/window_coverage/](../../../../src/ng/window_coverage/),
[reference_info.rs](../../../../src/ng/reference_info.rs), and the build configuration.

**Reviewed against** commit `cfad6b71` on `main`. Every change is on branch `ng-psp-vcf-perf`.

**Targets.** ng has to degrade gracefully from one sample to several thousand and from three
reads a position to several hundred
([design_principles.md](../../specs/design_principles.md) §0). Every number below was taken on
one corner of that: the tomato benchmark's **63 accessions at about three reads a position**,
walked over the first 20 intervals of `benchmarks/tomato1/regions.bed` (2 Mb of SL4.0) and, for
the final check, over all 80 (8 Mb). Nothing here is a claim about a thousand samples or about
300× coverage except where it says so. This is also the stage a person re-runs when they change
a parameter, so a fixed startup cost is paid over and over.

**Hardware.** macOS 15 on Apple Silicon, 18 logical cores, 64 GB. Native host release build
(`cargo build --release --bin pop_var_caller_exp`): fat LTO, one codegen unit, `panic = "abort"`,
`-C target-cpu=apple-m1`, mimalloc. **These are not the container's numbers** — the committed
baselines in this repo are taken in the Debian 12 aarch64 container, and a cross-machine
comparison with them is not sound.

**Hot-path evidence, and it is the good kind.** Two sampling profiles of the real command on the
real cohort (macOS `sample`, 63 accessions): one over the calling pass, one over the
hidden-duplication filter's scoring pass. DHAT allocation counts, which are deterministic. A
resident-set-over-time trace. Every applied change was then measured by interleaved A/B in one
sitting, medians and full ranges quoted. Raw material is in
`tmp/perf_review_2026-09-07_ng-psp-to-vcf/`.

**Deliberately out of scope:** `src/var_calling/`, `src/pileup/`, `src/ssr/`, `src/psp/` — the
frozen production caller, which ng must not edit. The `generate-psps`, `generate-census` and
`estimate-parameters` commands. One finding lands in `src/ng/paralog/locus_score.rs`, which is a
guarded verbatim copy of production's file; it is filed and not applied, for that reason among
others.

**Categories dispatched**, one reviewer each, every one in its own git worktree so that no two
benchmarks shared a target directory:

| category | why |
|---|---|
| `allocations` | the allocator is a quarter of the calling pass's busy CPU |
| `hot_loops` | `log1p` is 59% of the scoring pass and a varint decoder's cold arm was third in the calling pass |
| `concurrency` | eighteen threads bought 10% over one |
| `io_and_syscalls` + `data_layout` | 63 open files, a spill read twice, and the types the merge holds by the million |
| `methodology` | run by the orchestrator; the build configuration was already sound |

**One thing went wrong with the fan-out and is worth recording.** All four reviewers were
checked out at `a33ada0f` rather than at the commit under review, and two directories in scope —
`src/ng/run/paralog_filter/` and `src/ng/window_coverage/` — did not exist there. Each reviewer
detected it from the step-0 check in its prompt and re-pointed its own worktree before reading
anything. The check is what saved the review; without it four agents would have benchmarked a
tree that was missing half the scope.

---

## 2. Verdict

**Apply the listed wins.** Six are applied, measured and committed. On 63 tomato accessions over
2 Mb, interleaved on an idle machine, three rounds each:

| | wall clock | peak resident |
|---|---|---|
| `cfad6b71` | 43.55 / 45.25 / 44.06 s | 22.89 / 23.02 / 22.88 GB |
| `2ba1f7a9` | 22.86 / 22.48 / 22.50 s | 0.39 / 0.40 / 0.40 GB |

**1.96× the speed on 1/57th of the memory, and the VCF body is identical** — the two files
differ only in the two header lines that record the command line and the parameters file's name,
both of which name the output path.

**The memory defect was the one that mattered, and it was not a tuning question.** Peak resident
grew with samples × ground and nothing was released until the run ended, so the ceiling was the
machine rather than the design. Over the full 8 Mb region set the same 63 accessions now peak at
**573 MB**; the old binary passed **45.6 GB** after 48 seconds of walking, still rising, and was
stopped rather than allowed to exhaust a 64 GB machine.

Peak resident against cohort size, all over 2 Mb, after the changes: **8 accessions 221 MB, 16
266 MB, 32 389 MB, 63 416 MB** — about 3 MB a sample over a 190 MB floor, and flattening.

Threads now buy something. Before, eighteen against one was 48.93 s → 44.06 s (1.11×); now it is
45.09 s → 22.50 s (**2.0×**), and the VCF is identical at both.

Three further candidates were measured and are **not** applied; §5 gives each a number and a
reason.

---

## 3. Measurement plan

The workload and the harness are reusable; both are described in
`tmp/perf_review_2026-09-07_ng-psp-to-vcf/EVIDENCE.md`. To repeat any of this:

1. **Build the fixture** — `generate-psps` over 63 tomato CRAMs and the first 20 intervals of
   `benchmarks/tomato1/regions.bed`, into `tmp/perf_psp_vcf/psps` (657 MB, ~1.9 M stored loci a
   sample, 2 minutes). An 8-accession subset in `tmp/perf_psp_vcf/psps8` runs in 3 s and is what
   every A/B below used; the 63-accession set is for confirmation.
2. **Wall clock and peak resident** — `/usr/bin/time -l`, arms interleaved in one sitting,
   median and full range quoted. This machine's spread between runs of one unchanged binary
   reached a full second while four agents were benchmarking, and is about 0.05 s when idle;
   quote which.
3. **Where the time went inside the run** — the command's own per-pass timer line, which splits
   calling from scoring from writing and does not move with machine load the way wall clock does.
4. **Allocation counts** — `cargo run --release --example ng_call_from_psps_cost
   --no-default-features --features dhat-heap <reference> <catalog> <psp-dir>`, which now prints
   what the calling pass allocated beside the timings it already printed. Deterministic:
   identical run to run, so it settles a change whose wall-clock effect is under the drift, which
   two of the six adopted changes were.
5. **CPU profile** — macOS `sample <pid> <secs> 1`, attached once the run is under way. Rayon's
   parked workers appear as `__psynch_cvwait` with the full sample count; filter them out.
6. **Correctness gate on every change** — the VCF body must be identical (`diff` with `^##` lines
   filtered), and it must not depend on the thread count.

**What this machine cannot measure**, and it constrained one category: there are no hardware
counters — the PMU is not virtualised — so there is no cache-miss count, no branch-miss count, no
`perf c2c` and no `perf sched`. Layout questions were settled by wall-clock A/B and `size_of`.

---

## 4. Build / toolchain configuration

**Already sound, and nothing here blocks the code-level work.** `[profile.release]` is fat LTO,
one codegen unit, `panic = "abort"`, `debug = "line-tables-only"`; `.cargo/config.toml` lifts the
CPU floor to `x86-64-v3` on Linux x86_64 and `apple-m1` on macOS aarch64.

**One gap, and it is a note rather than a finding.** There is no `[target]` entry for **aarch64
Linux**, which is the production target: the dev container builds at the generic ARMv8-A baseline
while the macOS host builds at `apple-m1`. Every number in this review is a host number, so this
changes nothing here. Adopting a floor there means naming the oldest ARM server the project must
run on, which is an owner's decision and not a measurement.

**mimalloc is the default and it earned its place again, but for a different reason than the
`Cargo.toml` comment gives.** Before the memory fix it *cost* 4.2 GB of peak resident on this
path (24.4 GB against 20.1 GB with the system allocator) for no wall-time gain — the opposite
direction from the merge measurement that comment cites. After the fix the retention it was
holding is gone and the question is moot at this scale. What mimalloc does buy on this path is
making a 51-byte allocation cheap enough that removing 15 million of them does not move the wall
clock, which is why two of the findings below are gated on counts rather than on time.

---

## 5. Code-level findings

### Hot-path

#### H1: [src/ng/run/psp_source.rs:233](../../../../src/ng/run/psp_source.rs#L233) — the source keeps every body and every head the file ever handed over

**Confidence:** High. **Applied in `e768cabc`.**

`PspSummarySource` holds each record's evidence bytes in an append-only arena and one
`KeptRecord` per record beside it. Nothing truncated, drained or cleared either. The merge's
window is bounded and the source behind it was not, so peak resident grew with samples × ground:
a resident-set trace of the 63-accession run rises monotonically from 8 MB to 23.6 GB over the
25 seconds of the calling pass and then sits flat.

This is what `cohort_merge_psp_path.md` §3.2 already asks for — *"the retained window is the raw
bytes of every record between the two cursors, per sample — a FIFO the cover advances and
eviction drains"*. The FIFO was there and the drain was not.

`ObservationSource` gains `release_before`, the other half of `build`: the cache says where its
window's left edge went, and a source that keeps evidence drops the same prefix. The default does
nothing, which is right for every source that builds each record as it is drawn. Body ranges are
measured from the file rather than from the buffer, so a range the merge holds across a release
still names the same bytes.

| | peak resident before | after |
|---|---|---|
| 8 accessions, 400 kb | 782 MB | 127 MB |
| 8 accessions, 2 Mb | 2.61 GB | 215 MB |
| 63 accessions, 2 Mb | 24.4 GB | 397 MB |

Wall clock is unchanged (5.26 s against 5.13 s median on 8 accessions). **Complexity cost:** one
trait method with a no-op default, one call site in `evict_before`, and an offset on the source.
**Nothing in any output shows the difference**, so both halves have tests of their own:
`releasing_what_the_merge_has_passed_leaves_the_rest_building_the_same` and
`a_source_that_keeps_its_evidence_is_told_what_the_eviction_passed`.

#### H2: [src/ng/run/paralog_filter/pass_two.rs:290](../../../../src/ng/run/paralog_filter/pass_two.rs#L290) — the hidden-duplication filter scores every record on one thread while seventeen sleep

**Confidence:** High. **Applied in `cf5b238d`.** Found by the concurrency reviewer.

18.21 s of a 41.65 s calling-and-scoring span at 63 accessions, with the main thread busy for
11,445 of 11,473 profile samples and every worker asleep. `WhereTheTimeWent`'s own doc says spec
§8 deferred parallel scoring *until there was a number saying whether it was worth it*. It is
41%.

Scoring one parked record is a pure function of that record, so the spill is read serially in
batches of 2,048, the batch is scored on the pool, and the ratios are folded back in the batch's
own order — which is what keeps the histogram and the kept vector value-for-value what a serial
walk produced.

**Scoring falls from 18.21 s to 1.57 s at 63 accessions**, and it scales with the record count
rather than the sample count: 769 / 208 / 122 / 67 ms at 1 / 4 / 8 / 18 threads on the 8-accession
fixture. Peak resident unchanged.

**Complexity cost:** one batching loop. The reviewer's patch kept a serial arm beside it; that was
measured and dropped — on a pool of one thread the batched loop costs what the plain `for` cost
(769 ms against 773), so the second copy bought nothing and was a second place for the fold's
order to be got wrong.

#### H3: [src/ng/run/cohort_merge/observation_cache.rs:1420](../../../../src/ng/run/cohort_merge/observation_cache.rs#L1420) — closing a cover measures every sample's coverage on the calling thread

**Confidence:** High. **Applied in `aeee3f58`.** Found by the concurrency reviewer.

The samples were drawn on the pool one statement earlier; folding each one's newly drawn records
into its own accumulator then happened in a serial loop. The bases are one borrowed slice every
sample reads by offset and a sample touches nothing but its own accumulator.

Interleaved on 63 accessions, four rounds each, reading the run's own calling-pass timer: **21.26
/ 21.19 / 20.97 / 20.90 s against 19.79 / 19.79 / 19.80 / 19.74** — medians 21.08 and 19.79, no
overlap. 6.1% off the calling pass, about 4.9% off the whole command.

**Complexity cost, and it is the largest of the six.** The naive form needs `S: Send, E: Send`,
and those bounds cascade into the *serial* public drivers, where 36 test sites build sources on
`Rc` and stop compiling. The shipped form passes the measure loop down as a strategy parameter on
two private methods, which confines the bound to `cover_in_parallel`. A closure parameter for 5%
is a fair trade but not an obvious one; it is the change to revert first if this file becomes hard
to read.

#### H4: [src/ng/run/psp_source.rs:268](../../../../src/ng/run/psp_source.rs#L268) — a record's live reads are their own heap block; the bytes beside them are not

**Confidence:** High. **Applied in `84bfa038`.** Found independently by the allocations and
io-and-layout reviewers.

`KeptRecord` owned a `Vec` of the reads live at its position — the last variable-length field in
the type that still boxed, where the record's bytes have been in a per-sample arena all along. It
now carries a span into a second arena, drained on its own cut.

| | calling pass, 3 interleaved rounds each |
|---|---|
| 8 accessions | 962 / 957 / 957 ms → 884 / 877 / 882 ms (**8.0%**) |
| 63 accessions | 19.12 / 19.48 / 19.60 s → 19.01 / 19.13 / 19.24 s (**1.8%**) |

DHAT over the calling pass on 8 accessions: **20,059,090 allocations against 4,993,031**. Peak
resident is a wash at both sizes.

**Complexity cost:** a third arena with a cut of its own, and that cut is a trap. A record
contributes as many bytes as its body has and as many identifiers as it has reads, so the two
arenas fill at different rates and one cut cannot serve both. A release that used the bodies' cut
hands a surviving record another record's read set, which decodes into a plausible record rather
than a refusal — so the module's release test now builds its fixture with chain identifiers, and
an over-releasing mutant of the cut fails it.

**It also corrects a claim in the same file.** `build`'s doc said the encoder writes no chain ids,
so every live set was empty. It reasoned from `encode_record_body`, which drops them from a
record's *body*; the *head* still writes the live-set changes. Counted over the benchmark's stored
files: **96,728,802 identifiers over 15,074,770 records, 6.42 a record, only 6 records in 10,000
empty.** Two facts in two files, and the comment checked one.

#### H5: [src/ng/run/psp_source.rs:673](../../../../src/ng/run/psp_source.rs#L673) — every stored record was cloned to file it

**Confidence:** High. **Applied in `04877fdb`.** Found by the allocations reviewer.

`next_drawn` cloned each record to push it onto `heads`, then read two plain values out of the
original and dropped the clone. Since a record owns its live-read set, that copy allocated and
freed a heap block on every record walked.

DHAT over the calling pass on 8 accessions: **35,125,430 allocations against 20,059,267** — one
per record walked, 43% of the pass's allocations, and 774 MB of copying. Peak live blocks and
bytes unchanged to the byte: this is churn, not retention.

**Wall clock does not separate the two** (5.15 s against 5.12 s median, on a machine whose spread
was over a second). The argument for it is the deterministic count, three lines of diff, and that
the cost per record is the length of the live set — the fixture is three reads a position and ng
has to hold at three hundred.

#### H6: [src/ng/reference_info.rs:748](../../../../src/ng/reference_info.rs#L748) — a reference's two digests were hashed one after the other

**Confidence:** High. **Applied in `2ba1f7a9`.** The startup cost was found by the io-and-layout
reviewer; the fix is the orchestrator's and is not the one that review proposed (see L1).

Verifying a reference computes two whole MD5s over the same bytes — one per contig, one over the
assembly. On the tomato reference that is about 2.0 s of a 2.7 s pass (`md5` over the same 795 MB
file costs 1.02 s here), and it ran serially while every other core was idle.

Each 64 KiB window is now hashed into both at once. Five interleaved rounds on 8 accessions after
a warm one: **3.55–3.61 s against 3.10–3.12 s**, no overlap. The 0.47 s is the whole of it — the
calling pass does not move — so it is the same 0.47 s on a cohort of any size, and on every other
command that reads a reference.

**A run with one thread does not pay for a fork-join it cannot use.** The pass flushes about
twelve thousand windows, and splitting them cost `--threads 1` 5.04 s against 4.80. The pool's
size is read once when the pass starts; guarded, one thread measures 4.80 against 4.81.

### Likely

#### L1: [src/pop_var_caller_exp/call_from_psps.rs:476](../../../../src/pop_var_caller_exp/call_from_psps.rs#L476) — the reference verification could be joined after the calling pass instead of before it

**Confidence:** High that it works; **not applied**, and the reason is a trade the owner should
make rather than a measurement. Found and measured by the io-and-layout reviewer.

`read_reference_verifying_or_creating_fai` already spawns the FASTA pass on its own thread, and
the handle's own doc says to join it "before the caller commits any output". This command joins it
eight statements later, before the cohort is even opened. Nothing between the join and the end of
the calling pass consumes a digest.

Measured with the whole verification moved past the calls: **5.18 s → 3.81 s median on 8
accessions, 26%**, clean separation, VCF byte-identical. That is 2.7 s where H6 recovers 0.5.

**Why it is not applied.** A wrong catalog, or a reference that no longer matches its index, would
be reported **after** the calling pass instead of before it — on the 63-accession run that is 20
seconds of work thrown away before the operator is told, and on a genome-scale run it is hours.
Only the *digest* half moves; the name, length and order checks still run up front. The patch also
makes the ordering conditional (`--parameters` and `--paralog-fdr 0` must still join early, so
there are two orderings where there was one) and re-derives the run's parameters after the join.

**The recommendation is to leave it.** The 2.2 s H6 does not recover is a fixed cost, so it is 2%
of a 63-accession run over 2 Mb and 0.6% over 8 Mb — it only looks large next to a small run. If
the owner would rather have it, the working patch is at
`tmp/perf_review_2026-09-07_ng-psp-to-vcf/ab/patch_deferred_verify.diff` and the honest cost is
"a wrong reference is reported at the end".

#### L2: [src/ng/paralog/locus_score.rs:495](../../../../src/ng/paralog/locus_score.rs#L495) — `log_sum_exp3` pays two `log1p` where one max-shift needs one `log`

**Confidence:** High that it is faster; **not applied.** Found and measured by the hot-loops
reviewer.

Written as two nested pairwise log-add-exps, it evaluates the stabilising shift twice: two `exp`
and two `log1p`, where a single max-shift over the three terms needs three `exp` and one `log`.
`log1p` was 6,776 of the 11,473 samples in the scoring profile. Applied, the run's own scoring
timer went from 788 ms to 645 ms median on 8 accessions (18%).

**Why it is not applied, and the first reason is the decisive one.** H2 already took that pass
from 782 ms to 67 ms by putting it on the pool, so 18% of the serial cost is now about 12 ms of a
67 ms pass. The two remaining reasons would each have been enough on their own: the file is a
**guarded verbatim copy** of production's `src/paralog/locus_score.rs` and `copy_fidelity.rs`
fails on any edit, so the change lands in both trees or the guard is released; and it is a float
reassociation whose bit-parity test puts the disagreement at **8 units in the last place** on the
likelihood ratio, which is a change to what the filter drops, not merely to how fast it decides.

#### L3: [src/ng/calling/genotype_prior/dirichlet_multinomial.rs:257](../../../../src/ng/calling/genotype_prior/dirichlet_multinomial.rs#L257) — the inbreeding mixture computes two `lgamma` and one `ln` per homozygous genotype that are discarded when F is zero

**Confidence:** Medium. **Not applied.** Found by the hot-loops reviewer.

F = 0 is the shipped default, and the shortcut is byte-identical — it discards **5.56 million
transcendental calls** on the 8-accession fixture. The reviewer could not separate it from the
drift at that size and filed it with a plan for the 63-accession fixture.

**Recommendation: take it, but measure it on 63 accessions first**, where the calling pass is 20 s
rather than 0.9 and `lgamma` was 1,385 profile samples. It is byte-identical and the guard is one
branch, so the only question is whether it is worth anything; if the 63-accession A/B cannot
separate it either, close it.

#### L4: [src/ng/run/cohort_merge/serial.rs:343](../../../../src/ng/run/cohort_merge/serial.rs#L343) — decoding bodies and genotyping stay on the calling thread

**Confidence:** High that the parallelism is unclaimed; **not applied**, and it is the largest
remaining lever. Found by the concurrency reviewer.

**The calling pass is now 83% of the run and it is one thread.** Body decode is 29% of the main
thread's calling-pass work and genotyping 34%; the region-parallel driver that would move them
(`merge_cohort_in_parallel`) is written and tested and no production caller uses it. Inside the
current shape the ceiling is small — the main thread was blocked waiting for the cover for only
872 of 15,942 samples (5.5%) — so this is a topology change, not a scheduling one.

**Recommendation: this is the next piece of work, and it is a piece of work rather than a
finding.** It makes a whole round of regions resident at once, so it wants H1's release in place
first (it is), a memory measurement of its own, and the byte-for-byte agreement the module already
holds its three drivers to.

### Speculative

#### S1: [src/ng/vcf/encode.rs:137](../../../../src/ng/vcf/encode.rs#L137) — a VCF line is built out of about ten short-lived `String`s plus six per sample

Measured at under 0.3% of the profile on this cohort; filed for its shape, which is per sample per
record, rather than for today's share. **Not applied.**

#### S2: [src/ng/run/paralog_filter/pass_one.rs:190](../../../../src/ng/run/paralog_filter/pass_one.rs#L190) — each parked record builds a fresh per-sample `Vec` that is encoded into a reused buffer and dropped

Same size and same reason as S1. **Not applied.**

#### S3: [src/ng/paralog/locus_score.rs:427](../../../../src/ng/paralog/locus_score.rs#L427) — H2's non-carrier term is rebuilt once per carrier configuration although it does not depend on one

Superseded by H2 for the same reason as L2, and it is in the same guarded copy. **Not applied.**

### Note

- **[src/psp/varint.rs:91](../../../../src/psp/varint.rs#L91) — the `#[cold]` multi-byte arm is
  correctly marked.** The evidence file put `decode_u64_leb128_cold` third in the calling pass's
  self-time and implied the annotation was wrong. It is not: the hot-loops reviewer counted the
  widths and **98.6% of 113 million decodes are a single byte**. A wider inline arm measured no
  change. Recorded so nobody re-files it.
- **Swapping `libsystem_m` for the pure-Rust `libm` in `log_add_exp` is 34% slower** (863 ms
  against 1,160). Recorded for the same reason.
- **An exact "the smaller branch cannot move the larger" early-out in `log_add_exp` fires on 8.0%
  of calls and is a wash.**
- **`READ_CHUNK_BYTES` is already at its measured optimum** — 3,791 `read(2)` calls for the whole
  8-accession walk, and the constant's own sweep found 16 KiB best.
- **The filter's two reads of the spill are not an I/O cost** — 12 profile samples in the decoder
  against 6,776 for `log1p`. The 43% was arithmetic, not disk.
- **The cohort-observation build and the whole of `src/ng/calling/` are not allocation-bound** —
  largest allocator sites 47 and 38 samples. Two negative results worth keeping.
- **The three parallel arrays in `SampleWindow` cost what the module claims they do** — fusing
  them would put `draw_to`'s stride at 168 bytes instead of 32.
- **Fork-join granularity in the parallel cover was checked and is not a cost.**

---

## 6. Out-of-scope observations

- **`generate-psps` writes what `call-from-psps` then reads back**, and this review only measured
  the reading. The 63-accession fixture took 2 min 23 s to produce 657 MB of psp for 2 Mb of
  ground; nobody has profiled that direction.
- **`merge_cohort_in_parallel` has no production caller** — it exists, is tested, and is reached
  only from tests and examples. That is either a lever (L4) or dead weight, and it should not stay
  ambiguous.
- **`cargo test --release --lib` cannot be used as a gate without a filter**: 8 of 6,635 tests fail
  on an unmodified tree, all of them `#[should_panic]` tests waiting on debug assertions a release
  build does not compile in.
- **The repository is not `rustfmt`-clean**: 28 files differ from `cargo fmt` output at
  `cfad6b71`, and it is still 28 after this branch. Worth a one-off pass so that a formatting
  failure means something.

---

## 7. What's already good

- **The merge keeps summaries and records in parallel arrays rather than one array of pairs**
  ([observation_cache.rs:413](../../../../src/ng/run/cohort_merge/observation_cache.rs#L413)),
  which is exactly what lets a run over stored files fill the first and leave the second empty
  until a locus survives. That design is why H1 was a missing drain rather than a rewrite.
- **The run reports where its own time went**, split into calling, fitting, scoring and writing
  ([finish.rs:187](../../../../src/ng/run/paralog_filter/finish.rs#L187)). Half of this review's
  A/B work used that line instead of wall clock, because it does not move with machine load.
- **A single-position record's depth is served from its head without building the body**
  ([depth.rs:76](../../../../src/ng/window_coverage/depth.rs#L76)) — all but about one record in
  900. Without it the deferred build would have been undone by the coverage measurement.

---

## Author response convention

`H1` applied in `e768cabc`. `H2` applied in `cf5b238d`. `H3` applied in `aeee3f58`. `H4` applied
in `84bfa038`. `H5` applied in `04877fdb`. `H6` applied in `2ba1f7a9`.

`L1` deferred — the trade is the owner's, patch kept. `L2` closed, superseded by `H2` and blocked
by the copy guard. `L3` open, with the experiment named. `L4` open, and it is the next piece of
work. `S1`–`S3` won't fix at this size.
