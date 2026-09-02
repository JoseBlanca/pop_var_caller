# Performance Review: ng-calling
**Date:** 2026-09-02
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** the `call-from-alignments` path — CRAMs to VCF in one process
**Verdict:** Apply the listed wins — three are applied and measured in this branch; the rest are ranked with measurement plans
**Hot-path evidence:** two `sample(1)` profiles of the real command on the real cohort, a thread sweep, and eleven whole-benchmark timings

---

## 1. Scope and constraints

**What was reviewed.** The path `pop_var_caller_exp call-from-alignments` takes: the
run driver (`src/ng/run/`), locus generation (`src/ng/locus_generation/`), read and
reference input (`src/ng/read/`, `src/ng/ref_seq.rs`, `src/ng/raw_chrom_reader.rs`),
the calling loop (`src/ng/calling/`), the VCF writer (`src/ng/vcf/`), and the build
configuration (`Cargo.toml`, `.cargo/config.toml`).

**Reviewed against** commit `f6d38ea9` on `main`. The changes this review applied are
on branch `ng-perf`.

**Targets.** ng has to degrade gracefully from one sample to several thousand and from
three reads a position to several hundred (`design_principles.md` §0). The workload
every number below was taken on is one corner of that: the tomato benchmark's **63
accessions at about three reads a position**, over the 80 intervals (8 Mb of SL4.0) of
`benchmarks/tomato1/regions.bed`, or over the first four of them (400 kb) where a
faster loop was wanted. Nothing here is a claim about a thousand samples or about 300×
coverage except where it says so.

**Hardware.** macOS 15 on Apple Silicon, 18 logical cores. Native host release build
(`cargo build --release --bin pop_var_caller_exp`): fat LTO, one codegen unit, `panic =
"abort"`, `-C target-cpu=apple-m1`. **These are not the container's numbers** — the
committed baselines in this repo are taken in the Debian 12 aarch64 container, which is
the production target, and a cross-machine comparison with them is not sound.

**Hot-path evidence, and it is the good kind.** Two sampling profiles of the real
command on the real cohort, taken with the macOS `sample` tool:

- 18 threads, 20 BED regions, 45 s of sampling — `tmp/perf_review_2026-09-02_ng-calling/sample_par.txt`
- one thread (`RAYON_NUM_THREADS=1`), 20 BED regions, 60 s — `.../sample_serial.txt`

plus a thread sweep and eleven whole-benchmark runs. The raw numbers and the
verbatim profile extracts are in `.../EVIDENCE.md`.

**In-scope files** are listed with each finding.

**Deliberately out of scope:** `src/var_calling/`, `src/pileup/`, `src/ssr/`,
`src/per_sample_pileup/` — the frozen production caller, which ng must not edit
(`src/ng/mod.rs`'s freeze paragraph). One finding lands in `src/bam/alignment_input.rs`,
which ng *reuses* rather than copies; it is filed with that ownership question attached.

**Categories dispatched**, one sub-agent each, each in its own git worktree:

| category | why |
|---|---|
| `methodology` | always — and it found the largest single defect |
| `allocations` | the allocator is a quarter to a third of busy CPU in both profiles |
| `concurrency` | 1.92× of speedup for 2.22× of CPU on 18 cores |
| `hot_loops` | `PileupGenerator::next_locus` is the largest non-idle leaf in both profiles |
| `io_and_syscalls` | 63 CRAMs held open, and a third of the serial run is inside the CRAM reader |

**One deviation from the skill, and it is the orchestrator's call.** The skill tells
each sub-agent to measure its own candidates in its own worktree. Five agents
benchmarking at once on one machine produce plausible wrong numbers rather than a failed
build, and the orchestrator was running timed whole-benchmark arms throughout, so every
agent was forbidden to build, bench or profile. Every sub-agent finding is therefore a
*proposal with a measurement plan*; every number in section 2 and section 3 was measured
by the orchestrator, serially, on a quiet machine, with the arms alternated.

---

## 2. Verdict — and what is already applied

**Apply the listed wins.** Three landed on `ng-perf` and are measured; two were built,
measured, and reverted for showing nothing; the rest are ranked in section 5.

### The headline

**A calling run of the 63-accession tomato cohort over the whole 8 Mb benchmark went
from 212.0 s to 116.8 s — 1.82× — writing a byte-identical VCF, and holds 1.7× the
memory while it does.**

Every row below is one run of `pop_var_caller_exp call-from-alignments` over
`benchmarks/tomato1/regions.bed` with all 63 accessions, `--defaults`, 18 threads,
`/usr/bin/time -l`:

| build | round width | wall | CPU | peak resident |
|---|---:|---:|---:|---:|
| `f6d38ea9`, as shipped | 500 | **212.0 s** | 899 s | **1,113 MB** |
| this branch, system allocator | 500 | 206.9 s | 908 s | 1,155 MB |
| this branch, system allocator | 8,000 | 125.7 s | 773 s | 1,660 MB |
| this branch, mimalloc | 500 | 193.2 s | 765 s | 1,593 MB |
| this branch, mimalloc | 8,000 | 115.3 s | 685 s | 1,811 MB |
| **this branch, nothing set** | 7,936 (chosen) | **116.8 s** | 688 s | **1,897 MB** |

**The VCF is identical in every row.** Bodies compare equal byte for byte, and so does
every header line except `##commandline` and `##parametersFile`, which record the
binary's path and the output's.

**Every row of that table was taken with the machine otherwise idle, and the ratio is a
fact about that.** The rows came from two back-to-back batches that agree with each
other to within 5 s on the settings they share, so they are mutually comparable.
Partway through the session two virtual machines appeared on the host and took about two
cores; re-run under that load, the same two binaries alternated back to back, twice:

| | wall | CPU | peak resident |
|---|---:|---:|---:|
| `f6d38ea9`, as shipped | 320.7 s, 273.8 s | 1004 s, 950 s | 1,217 MB, 1,233 MB |
| this branch, nothing set | 214.4 s, 204.8 s | 845 s, 903 s | 1,708 MB, 1,815 MB |

**1.50× and 1.34× rather than 1.82×**, on VCFs that again compare identical. That is what
a change buying its time by overlapping reading across threads should do when there are
fewer free cores to overlap on. **Quote the ratio with the machine state; the seconds
alone are not portable between them.**

**The trade is real and it is memory.** Of the 784 MB the last row holds above the
first, about 500 MB is the wider round holding more observations at once, and about
150 MB is mimalloc. Both are dial-able from the command line, and the round width's
default already narrows itself as the cohort grows (below).

### What was applied

**A1 — `src/main_exp.rs` never installed mimalloc, so every calling run this project
has ever timed used the macOS system allocator.** `#[global_allocator]` is resolved per
*binary*, not per crate: `src/main.rs` declares it, `src/main_exp.rs` did not, and the
`alloc-mimalloc` default feature only makes the crate available. Both profiles show it —
every allocator leaf is `_xzm_free`, `_xzm_xzone_malloc_tiny`, `xzm_realloc`, which are
Apple's xzone malloc. The family is **25% of busy CPU on one thread and 34% at
eighteen**. Three lines, copied from `src/main.rs`.

Measured, alternated arms, three rounds each: 400 kb of ground and 63 accessions goes
25.7 s → 21.2 s at one thread and 13.4 s → 12.2 s at eighteen; the whole 8 Mb goes
212.0 s → 193.7 s. It costs peak resident, and on this path — unlike on the merge, where
the feature's own note in `Cargo.toml` records a 9% *saving* — it costs a lot of it:
1,155 MB → 1,593 MB at the shipped round width. `MIMALLOC_PURGE_DELAY=0` does not
recover it (197.9 s / 1,607 MB against 193.2 s / 1,593 MB), so it is not deferred
purging.

**A2 — how much ground one round of locus building covers is now a flag, and its
default is chosen from the cohort's size.** `CohortLocusBuilderRegionsLen` has said
"**A command-line parameter**" in its own documentation since it was written, and
`call-from-alignments` never took one: every run used the compiled-in 500 bases. That
number was measured on the *merge* reading pre-built `.psp` files, where a round costs a
scan over records already in memory. A calling run draws its records from one CRAM per
sample instead, and its rounds are exactly where that reading is overlapped across
threads — so a narrow round pays a fan-out, a barrier and eighteen thread wake-ups for a
few microseconds of work each time, and the waste grows with the cohort. The same doc
comment anticipated this: *"a cohort far from that should sweep it again — it is a run
parameter precisely so that costs no code."*

Swept on the 400 kb ground, 63 accessions, 18 threads, VCF identical at every width:

| width | 500 | 1,000 | 2,000 | 4,000 | 8,000 | 16,000 |
|---|---:|---:|---:|---:|---:|---:|
| wall | 11.81 s | 11.16 s | 10.08 s | 9.19 s | 8.26 s | 7.83 s |
| peak resident | 896 MB | 890 MB | 921 MB | 1,020 MB | 1,170 MB | 1,463 MB |

**What costs memory is the product `width × samples`, not the width** — a round holds
roughly one observation per covered base per sample. Measured across three cohort sizes
on the same 400 kb, the extra resident memory tracks that product: 63 samples at 8,000
holds 292 MB more than at 500 (product 472,500), and 16 samples at 32,000 holds 254 MB
more (product 504,000). So the default now bounds the product rather than the width:
about half a million observations, clamped to [500, 16,000] bases. That gives 7,936
bases at 63 samples, 500 — today's number — at a thousand, and the ceiling below 32
samples, where the gain has saturated anyway (four accessions: 3.29 s at 8,000, 3.18 s
at 32,000, 3.22 s at 64,000).

`--threads` came with it, for consistency with `estimate-contamination` and
`repeat-catalog`, which both have one. It is a CPU-cost knob, not a wall-clock one: at
63 accessions the threads past eight buy 1.0 s of wall for 14.5 s of CPU.

**A3 — three per-locus or per-read allocations removed.** Together they are 1.4% at one
thread and nothing measurable at eighteen; each is also a plain improvement, which is why
they stayed.

- `raw_chrom_reader.rs` cloned the contig name into an owned `String` on the **success**
  path of every reference read, so that an error that does not happen could name its
  contig. Binding the read's result first ends the borrow that forced it, and the clone
  moves to the failure path where every other error in that file already clones. Three
  sites.
- `aligned_reads_reader/container.rs` built a fresh `RecordBuf` per CRAM record — name,
  CIGAR, a quality-score push loop doubling 0 → 4 → … → 256 for a 150-base read, and the
  whole auxiliary-tag table, all from capacity zero and freed on the next iteration. One
  buffer now serves the whole container: noodles-sam 0.85 exposes
  `try_clone_from_alignment_record`, which clears and refills. The file's comment saying
  the allocations "happen either way" was stale.
- `read/input/reference.rs` did an unconditional `swap` on the shared resident-contig
  atomic per region query by every sample's draw. A load first is the same answer — when
  the stored contig already is this one the swap writes back what it read — and it keeps
  the cache line shared instead of taking it exclusive.

### What was tried and thrown away

Both were built, measured against the same tree on the same machine with the arms
alternated, showed no change outside run-to-run noise, and were reverted:

- **Devirtualising `try_ordinary_column`'s `&dyn RefSeq`.** Its only caller is
  `process_position<F: RefSeq>`, which holds a concrete `&F`, so the erasure buys
  nothing and the vtable call plus the un-inlinable fetch body looked worth removing.
  20.57 s → 20.55 s at one thread over three alternated rounds. Reverted: it adds a
  type parameter and a monomorphisation for a change nothing can see.
- **Borrowing rather than taking the region in `GeneratorSet::next_locus`.** The
  function moved a `TypedRegion` out of an `Option` and straight back at every locus for
  a borrow the compiler already allows. Same measurement, same null. Reverted.

Two further candidates were **not** built, and the reason is the same in both: the fix
is bigger than the finding. `SmallVec<[ChainId; 4]>` for the per-observation chain-id
vector is 7.2% of busy CPU at the system allocator, but the type is read in about twenty
files across `calling/`, `parameter_estimation/` and `psp/`, and the whole reason it
looked large is the allocator this branch has now replaced. Reading the reference ahead
in `RawChromReader::append_forward` instead of one byte at a time moves a documented
memory bound that tests assert against.

---

## 3. Measurement plan

The instruments this review used are in the repository and should be the ones the next
one starts from.

**Reproducing the headline.** `tmp/perf_review_2026-09-02_ng-calling/` holds
`run_ng.sh` (the command with the cohort and the reference filled in), `ab.sh`
(alternates two binaries over N rounds, printing wall and peak resident from
`/usr/bin/time -l`), `run_ng_n.sh` and `width_by_cohort.sh` (the same over a cut-down
cohort), and `regions4.bed` / `regions20.bed`. A 400 kb round is 8–12 s a side; the
whole benchmark is 2–4 minutes a side.

**The output gate, and the obvious form of it does not work.** Two runs that agree on
every call still differ by a few bytes, because `##commandline` and `##parametersFile`
record the paths they were given. Compare with those two lines stripped:

```
cmp <(grep -v '^##commandline\|^##parametersFile' a.vcf) \
    <(grep -v '^##commandline\|^##parametersFile' b.vcf)
```

**What to profile with.** `sample <pid> <seconds> -file out.txt` on the macOS host,
against a natively-built release binary. It symbolicates inline and its *Sort by top of
stack* section is a self-time ranking. Rayon's parked workers show as
`__psynch_cvwait` with the full sample count and must be excluded before any share is
computed. In the container, `perf record -e cpu-clock` (the default `cycles` event
returns `<not supported>` — the PMU is not virtualised).

**Three gaps in the measurement set-up, in the order they matter.**

1. **No benchmark or harness in the repository reads a CRAM.** The shared fixture
   (`examples/shared/synthetic_alignment.rs`) writes a BAM. So the single largest cost
   family in the profile — CRAM block decode plus the reference MD5 that goes with it,
   about a third of the serial run — has no A/B instrument at all. Anything aimed at
   H2 or H9 needs one first.
2. **A 295-second run reports counts and not one second of timing.** `merge-timing`
   exists and is wired at forty-odd call sites, but only
   `examples/ng_call_cohort_end_to_end.rs` prints it, and that example writes no VCF.
   Giving the real command a way to print the phase split is what would make a long run
   self-describing. Check first whether the counters survive E1's parallel cover —
   `timing.rs` documents them as non-reentrant and E1 made the sweep concurrent.
3. **Four harnesses that drive this path omit `#[global_allocator]`**, including
   `benches/ng_generic_pileup_perf.rs` and `examples/ng_call_cohort_end_to_end.rs`, so
   they measure a different allocator from the binary. `examples/ng_psp_against_production.rs`
   already records what that omission cost on a sibling harness: 27%.

---

## 4. Build / toolchain configuration

**Nothing is left on the table in `[profile.release]`.** `lto = "fat"`,
`codegen-units = 1`, `panic = "abort"` and `-C target-cpu` are all set;
`debug = "line-tables-only"` costs the running program nothing (debug info is a section
the loader never maps) and is what let both profiles resolve to `cursor.rs:679` and
`slice.rs:123` without a rebuild — keep it.

**Two things absent, neither urgent.**

- **No `target-cpu` for aarch64 Linux**, which is the production container.
  `.cargo/config.toml` has entries for x86-64 Linux and aarch64 macOS only. On aarch64
  NEON is already in the baseline, so expect much less than the x86 case; the more
  defensible reason to care is comparability, since the repo's convention is that the
  container holds the committed baseline. Whether to set a floor is a portability
  decision nobody in this review can make — it turns on what aarch64 hardware this
  ships to. Either measure it with `hyperfine` in the container and set one, or write a
  comment saying the omission is deliberate, so the next reviewer does not re-raise it.
- **No profile-guided optimisation.** Sequence it last, after everything in section 5:
  it is the change that makes every other measurement harder to attribute.

**One thing to know rather than change.** `[profile.bench]` inherits `release` but cargo
silently drops `panic = "abort"` for bench and test profiles, because the harness needs
to unwind. So a criterion number on this path is taken on different codegen from the
shipping binary. That is an argument for building the calling-path A/B harness as an
`examples/` binary rather than a criterion bench — `cargo run --release` does get
`panic = "abort"`.

**And one that was checked so nobody chases it.** MD5 is 11% of busy single-threaded CPU
(H2) and there is no build lever for it: `md-5` 0.11 has no `asm` feature and its backend
selector picks the software implementation on everything but loongarch64. Any reduction
has to come from calling it less.

---

## 5. Code-level findings

Applied findings are in section 2 and not repeated. What follows is open.

### Hot-path

**H1 — decoding and calling never overlap, and this is the structural ceiling.**
[serial.rs:295-352](../../../../src/ng/run/cohort_merge/serial.rs#L295-L352). *Confidence: High for the diagnosis, Medium for the prize.*

Every one of the nineteen threads is parked 73–77% of the run. The merge thread does
8,170 samples of work and waits 27,252, all of them in one place — rayon's
`in_worker_cold` blocking on the fan-out's latch — and while it waits the eighteen
workers deliver under a third of the thread-time available to them. The loop is strictly
sequential: evict, then cover, then build-and-genotype-and-write, with no term in common.
Under Amdahl on E1's own 88/12 split a perfectly scaling cover would give 5.9×; the sweep
measures 1.92× before this branch's changes and 1.71× after them.

The structural answer is to stop *pulling* the draw through a barrier and start
*pushing* it: a pool of producers draws each sample's observations forward into a
bounded queue, and the merge consumes. A sample's record sequence is a function of its
source alone, so the fixpoint and the loci are unchanged; what changes is that region
*n+1* decodes while region *n* is genotyped and written.

*Measure the ceiling before building it.* `merge-timing` already counts both halves of a
cover — `COVER_BUSY_NANOS` sums the samples' own drawing across threads and
`COVER_NANOS` is the cover's wall — so `cover_busy / (threads × cover_wall)` is the
fan-out's efficiency with no new counter, and
`cover_busy/threads + (merge_wall − cover_wall)` is the perfect-overlap floor. Build only
if that floor is well under what the run already reaches.

*Complexity: the largest here, and it should be last.* The walkers move onto producer
threads for the run rather than being borrowed per sweep; the memory bound has to be a
cohort-wide budget of records in flight, not a per-sample depth, or a thousand-sample
run buys a thousand queues; and a producer that fails ahead of the merge must surface
that failure where the merge would have needed the record. It overlaps Milestone G
(H3) and should be designed with it.

**H2 — noodles re-MD5s the reference span of every CRAM slice, once per sample, and
there is no supported way to switch it off.** `noodles-cram-0.93.0` `slice.rs:359-363`,
reached from [container.rs:287-330](../../../../src/ng/read/input/aligned_reads_reader/container.rs#L287-L330). *Confidence: High.*

`md5::compress` is **11.2% of busy single-threaded CPU** — the largest non-allocator,
non-`PileupGenerator` leaf in both profiles, and 31% of the whole CRAM read path. For
every slice whose header carries a non-zero reference MD5, noodles hashes the entire
reference subsequence the slice *spans* and compares.

**How much that is, measured from the CRAM indexes rather than asserted.** Summing
`alignment_span` over a `.crai`'s entries gives exactly the bases hashed:

| file | slices | mean span | total per sample |
|---|---:|---:|---:|
| `benchmarks/tomato1/crams/SRR7279481.p1.bench.cram` | 65 | **8.23 Mb** | 535 Mb |
| `benchmarks/human_genome_bottle/crams/HG002_reads_selected_1000_rg.cram` | 1,081 | 2.48 Mb | 2.7 Gb |
| `benchmarks/tomato_big_cram/DRR000741.p1.cram` (whole file) | 112,140 | **7,075 bp** | 0.8 Gb |

**So most of this is an artefact of how the benchmark fixtures were made, and the report
of §2 should say so.** A CRAM slice holds a fixed number of *records*; on a file sliced
down to 80 scattered 100 kb windows, one slice's records run from one window to the next
and the span covers the megabases of empty reference between them — 1,163× the mean span
of the same species' whole-genome CRAM. Across 63 accessions the tomato run MD5s
**33.7 GB** of reference. A caller reading ordinary whole-genome CRAMs pays a small
fraction of that.

*What ng can do.* Nothing, at noodles-cram 0.93 — the builder has one setter
(`set_reference_sequence_repository`), `external_reference_sequence_is_required` is read
from the file's own compression header, and the repository cannot answer lazily because
noodles takes `&sequence[interval]` as a real slice the record decode also reads. The
options are: check whether 0.99 exposes a toggle; carry a `[patch.crates-io]` fork that
either skips the check or memoises it per span (and note that the `@SQ M5` check at open
is strictly stronger than a per-slice sub-span check, so nothing is lost by skipping);
or **re-slice the benchmark CRAMs with a smaller `seqs_per_slice`**, which changes the
fixtures rather than the caller but would tell you how much of the 11% is the fixture.

**H3 — the merge hands its finished record back for reuse and the walker frees it.**
[walker.rs:226](../../../../src/ng/run/walker.rs#L226). *Confidence: High.*

Of the 12,882 allocator samples under `SampleWindow::draw_to`, **12,647 are frees and
essentially none are allocations** — 8.2% of busy CPU at 18 threads, at the system
allocator. Eviction moves finished records into a per-sample spare list instead of
freeing them; `draw_next` pops one and offers it down; `AlignmentFilesWalker::next_observation`
drops it on its first line. So eviction does not avoid a free, it defers one into the
draw, and the draw then allocates a new record four layers below.

**This is the number Milestone G said it was missing.** `ng_run_driver_g2_2026-09-01.md`
dropped the milestone on an allocation-*count* ceiling of 20.7–23.9% and said in terms:
*"None of this is wall time … Getting the time would take either building G1 and timing
it, or a sampling profiler — and this machine cannot run one."* One has now run. It does
not overturn the owner's ruling, whose stated reason was to get a working caller first —
but the free side alone, at one call site, was 8.2% of busy CPU before this branch
changed the allocator, and it should be re-measured on the new allocator before G1 is
reopened.

*Before building anything, run the null experiment* — two lines, in `evict_before` and
`evict_before_in_parallel`: stop filling the spare list at all, so eviction frees where
it evicts. If wall time does not move, the deferral is free and the whole cost is the
record's own free, which is what leasing removes. If it improves, the spare list is
costing something today for a reuse nobody performs, and that is an immediate win.

**H4 — every observation's chain-id vector starts at capacity zero and grows one push
at a time.** [fast_column.rs:327](../../../../src/ng/locus_generation/pileup/fast_column.rs#L327),
[open_record.rs:903](../../../../src/ng/locus_generation/pileup/open_record.rs#L903). *Confidence: High.*

10,991 samples of `RawVec<u64>::grow_one` under `next_locus` — 7.2% of busy CPU at the
system allocator, the single largest identified allocation site. Growth is 0 → 4 → 8, so
one allocation per observation at three reads a position. `SmallVec<[ChainId; 4]>`
removes it; `smallvec` is already a dependency and `open_record.rs:844` already uses it.

*Not applied, and the reasons are two.* The type is read in about twenty files, so the
change needs `Deref<Target = [ChainId]>` to keep read sites compiling and the psp writer
and reader construct it explicitly. And `SmallVec` adds a branch per element access, so
at 300× coverage it is slightly slower — the sweep has to run at both ends of the depth
range, not just at three. Route it past the `chain_id_dead_weight` note first: the better
question may be whether those vectors should exist at all on observations that agree with
the reference.

**H5 — the fan-out fires at least twice per round and splits down to one sample per
job.** [observation_cache.rs:483-501](../../../../src/ng/run/cohort_merge/observation_cache.rs#L483-L501). *Confidence: High for the cost, Medium for the saving.*

Rayon's own parking, unparking and stealing is 14,814 samples, about 5.2 s of the 31 s
of extra CPU the eighteen-thread run spends. All 3,668 `__psynch_mutexwait` samples are
rayon's per-worker sleep-state lock at `sleep/mod.rs:291` — **not any lock of ours**.
Three costs stack: the cover runs until a whole sweep moves nothing, so **every cover
pays a sweep that changes nothing and exists only to prove the fixpoint**; the
`par_iter_mut` over 63 samples has no `with_min_len` and splits to singletons (a worker
stack shows six nested `join` levels above a ~200-sample leaf); and the Jacobi schedule
costs a sweep per link of a chain the serial form follows inside one sweep.

A2 already divides the number of covers by about sixteen, which is most of what this
finding was worth, and should be re-measured before either fix is attempted. If it is
still worth it: `with_min_len` is one line and no new invariant. Having `draw_to` return
*where* it stopped, so the driver proves the fixpoint without another sweep, is exact and
larger — but it turns the sweep from a `try_reduce` into a `map` into a scratch buffer.

**H6 — every worker hammers one process-global atomic once per ordinary column.**
[fast_column.rs:372](../../../../src/ng/locus_generation/pileup/fast_column.rs#L372),
declared at [pileup/mod.rs:409](../../../../src/ng/locus_generation/pileup/mod.rs#L409). *Confidence: High that it happens; unmeasured in isolation.*

`FAST_COLUMNS.fetch_add(1, Relaxed)` is on the return path of the lane that answers about
eight columns in ten, and unlike every other counter in that module it is not behind the
census gate. A relaxed `fetch_add` still takes the cache line exclusive, so with eighteen
workers each walking its own sample the line ping-pongs once per ordinary column — a
million loci a sample-pass, times 63 samples. This is true sharing, not false sharing.
Its own doc prices it *"against the ~90 reads that column would otherwise fold"*, which
is the 130× figure; at this cohort's three reads a position the column folds three.

The value is a run total nothing reads until the walk ends, so it can be a plain `u64` on
`WalkerState` folded into the static at region and chromosome end. **The cheapest probe
is to delete the increment behind a `cfg` and re-run the thread sweep** — not shippable,
but it measures the finding in one line. The diagnostic that says the mechanism was right
is the user/real ratio falling, not just the wall.

**H7 — loci the region clamp discards are built at full cost and then freed.**
[generator.rs:1125-1129](../../../../src/ng/locus_generation/pileup/generator.rs#L1125-L1129). *Confidence: Medium — the cost is measured, the discard rate is not.*

The allocator work whose nearest enclosing frame is the closing brace of the
`Some(Ok(locus))` arm is 4.0% of busy CPU and is **entirely frees**. The walk builds
records in the halo beyond the region and throws away the ones that start outside it.
Count the discard rate first — `counts.records_outside_region` already exists — because
if it is small the finding is small.

**H8 — the per-locus reference base goes through a vtable, a mutex and a one-byte
`BufReader` refill.** [fast_column.rs:166](../../../../src/ng/locus_generation/pileup/fast_column.rs#L166). *Confidence: High.*

`WindowedRefSeq::fetch_into`, `RawChromReader::fetch` and `read_raw_bases` are 2,028
self samples together. Because the walk moves one base at a time, `fetch` always lands in
the extend-forward branch and refills **exactly one byte** through the whole
`fill_buf` / byte-loop / `consume(1)` stack, per locus, per sample. The devirtualisation
half of this was tried and showed nothing (section 2). The read-ahead half was not: it
moves the resident-reference bound that `resident_reference_bases()` and its tests are
written against, which is more than a one-line change and needs its own measurement.

**H9 — every read is deep-cloned at least once, and a replayed read once per region that
replays it.** [cursor.rs:670](../../../../src/ng/read/input/cursor.rs#L670), [:684](../../../../src/ng/read/input/cursor.rs#L684). *Confidence: High for the clone, unmeasured for the replay factor.*

`AlignedRead` owns four `Vec`s, so each clone is four allocations and four memcpys;
the derived clone's inlined field copies are 1.9% of busy CPU. The module's own comment
says *"a read is returned by about a dozen consecutive regions"* — that is the code's
claim and not a measurement. **Print `CursorCounts::reads_replayed` against
`reads_decoded` first**; that fixes the multiplier for nearly nothing and decides whether
this is Hot-path or Likely. The honest fix is narrower than it looks: the preparer takes
the read by value and mutates its CIGAR, so an `Rc` would deep-copy anyway whenever the
retention deque still holds a reference.

### Likely

**L1 — all 63 samples decode their CRAMs against one shared `RwLock`, and its acquire
costs twenty-five times more CPU per unit of work in parallel than serial.**
[reference.rs:216-227](../../../../src/ng/read/input/reference.rs#L216-L227).

`noodles_fasta::Repository::get` is 742 samples in the parallel profile against 21 in the
serial one, and 696 of the 742 are at the shared-cache read acquire. Work-normalised that
is a ratio of at least 23×, the largest of any frame in either profile — and about 0.25 s
in absolute terms. **The module's own documentation predicted this and asked for exactly
this measurement**: *"Today's callers are single-threaded, so it costs nothing yet; it is
worth re-timing when one is not."* E1 made them multi-threaded on 2026-09-01 and it was
not re-timed. It is a small number and a large ratio, which is what a lock looks like just
before it becomes the problem. The fix is not to unshare the bases — 63 private copies of
tomato chromosome 1 is 90 MB each — but to unshare the *lock* while sharing the bases:
a per-sample `Repository` whose adapter is backed by a shared `OnceLock<Arc<Sequence>>`
per contig.

**L2 — the reference is read by 189 independent file sweeps that duplicate bases the
run's one `fasta::Repository` already holds.** [walker.rs:1554-1610](../../../../src/ng/run/walker.rs#L1554-L1610), [ref_seq.rs:594-624](../../../../src/ng/ref_seq.rs#L594-L624).

Three `WindowedRefSeq` per sample, each its own `File::open`, 64 KiB `BufReader` and
growing buffer, rebuilt at every contig change — 2,268 FASTA opens on a 63-sample,
12-contig run. Worth about 1.5% of busy CPU by the profile's own leaves.
`ResidentRefSeq` already exists and already serves both views out of the repository the
CRAM decode has forced resident anyway. Two caveats: a BAM-only run never opens that
repository, so the switch has to be conditional; and `Repository::get` takes the lock per
fetch (L1), so the accessor needs to cache the `Arc<Sequence>` or this trades one cost for
another.

**L3 — a forward region is always served by reading through the gap, however large.**
[cursor.rs:559-593](../../../../src/ng/read/input/cursor.rs#L559-L593).

The reuse test is `region.start >= last_region_start` with no distance term, so a forward
jump of any size walks every record in between rather than seeking. **It costs nothing on
the measured workload** — the bench CRAMs hold only the BED's reads, so the gaps are
empty — and it is filed for the case `design_principles.md` §0 asks about: a whole-genome
CRAM with a sparse `--regions` BED, one gene panel over a human cohort, where chromosome 1
alone walks 82.5 Mb of coordinate to serve 900 kb. One `u64` threshold and one extra term;
`jump_to` is already the unconditional path for a backward region, so it is known-correct
for any region.

**L4 — the SAM header is re-read and re-parsed once per file per chromosome and thrown
away.** [open_bam.rs:432-447](../../../../src/ng/read/input/open_bam.rs#L432-L447).
Invisible on tomato's 12-contig header; 72,000 parses on a 3,000-sample human cohort.

**L5 — two `lgamma` calls per locus-sample-pass compute a quantity that is one `ln` away
from a value computed five lines above.** [dirichlet_multinomial.rs:277-278](../../../../src/ng/calling/genotype_prior/dirichlet_multinomial.rs#L277-L278).

`lgamma(Σα + m) − lgamma(Σα)` is exactly `Σ_{j<m} ln(Σα + j)`, the identity this same file
uses as its test oracle. At ploidy 2 that is `ln(Σα) + ln(Σα + 1)`, and `ln(Σα)` is
already in hand. `lgamma` and its `log` children are about 1,710 samples. **Not
bit-identical** — but more accurate, since the current spelling subtracts two nearly equal
numbers and the file's own doc measures the cancellation at 1.5e-11. Needs a tolerance
gate rather than a byte-identity gate, which is why it is Likely and not applied.

**L6 — a one-byte heap allocation per locus and another per observation.**
[fast_column.rs:356](../../../../src/ng/locus_generation/pileup/fast_column.rs#L356), [:379](../../../../src/ng/locus_generation/pileup/fast_column.rs#L379).
`vec![ref_base].into_boxed_slice()` for a single base. 3.8% of busy CPU at the system
allocator. An inline representation for the one-base case would remove it and would touch
`SampleLocusObservations`, which is a wide type change.

**L7 — `sort_unstable` on a witness that is one run in every DNA-seq observation.**
[witness.rs:50](../../../../src/ng/locus_generation/witness.rs#L50).

**L8 — `find_overlapping` stays an out-of-line call on a table that is almost always
empty.** [open_record.rs:1468](../../../../src/ng/locus_generation/pileup/open_record.rs#L1468). 1,127 self samples.

### Speculative

**S1 — `ReadEvent` is sized by a variant holding a `Vec<u8>`**, so the general path moves
a wide `SmallVec` per contributor per column. [decompose.rs:15-37](../../../../src/ng/locus_generation/pileup/decompose.rs#L15-L37).

**S2 — the two measurement gates on the per-locus path are `OnceLock` acquire loads**,
where a `static AtomicBool` read once per region would do. [generator.rs:1118](../../../../src/ng/locus_generation/pileup/generator.rs#L1118), [genome_walk.rs:1260](../../../../src/ng/locus_generation/pileup/genome_walk.rs#L1260).

**S3 — every walker `next()` returns a `Result` as wide as `WalkerError`, whose variants
carry `String`s.** [genome_walk.rs:506-528](../../../../src/ng/locus_generation/pileup/genome_walk.rs#L506-L528).

### Note

- **The resolved round width is not recorded in the run's output.** With A2 the width
  depends on the cohort's size, so two runs over the same ground with different cohorts
  choose different widths, and `##commandline` — which records what was typed — says
  nothing about it. The same is already owed for `cohort_locus_builder_regions_in_flight`
  by that field's own doc comment. Neither changes an answer; both change memory, and a
  reader of the output cannot tell.
- **`WindowedRefSeq`'s mutex is not serialising anything**, and the hypothesis should not
  be re-raised: not one `__psynch_mutexwait` sample has a `ref_seq.rs` parent, and each
  sample holds its own accessor. Sixty-three locks, not one.
- **The `merge-timing` stopwatches compile out.** Both types are zero-sized without the
  feature; the calls are not a cost in a shipping build.
- **`f64::algebraic_*` is essentially unavailable on this path.** Every `f64` reduction in
  scope has its summation order documented as load-bearing for byte identity, which is
  exactly what those methods relax — `q_sum`'s accumulation, `one_genotypes_log_prior`'s
  fold, `concentration.sum()`. The reductions are also short (three terms at three reads a
  position). The place worth looking is the read-likelihood accumulation in
  `summarise_condition.rs`, which is per-sample per-genotype and long at depth.
- **The fast lane's linear observation scan is not quadratic in depth** and should not be
  filed as such: the scan is bounded by distinct bases × read groups, not by depth, and
  `open_record.rs`'s own comment records that hash-keying it measured worse.
- **The VCF write path is invisible** — every `ng::vcf` symbol together is 15 samples of
  672,882 — though it does about 300 `format!` allocations a record at 63 samples, which
  will matter at a cohort size this benchmark does not reach.

---

## 6. Out-of-scope observations

- [`src/bam/alignment_input.rs:1067-1084`](../../../../src/bam/alignment_input.rs#L1067-L1084) —
  `read_exceeds_mismatch_fraction` is 1,288 self samples (2.5% of busy CPU at 18 threads).
  Its per-base loop does three bounds-checked `get().unwrap_or()` and two data-dependent
  branches, which blocks vectorisation. **This is production's file, reused by ng rather
  than copied**, so the freeze rule applies and the ownership question has to be settled
  before it is touched.
- [`src/pileup/per_sample/baq_engine.rs:406`](../../../../src/pileup/per_sample/baq_engine.rs#L406) —
  `prepare_passthrough` clones the whole quality string and then frees the identical
  original; `mem::take` is the one-token fix. Also production's file. The ng-side
  alternative costs the field-for-field parity guarantee, so this is a decision rather
  than a change to make.
- `noodles-cram`'s `ExternalDataReaders` is a `HashMap` with the default `RandomState`,
  keyed on integer block content ids; `hash_one` is 534 self samples in the parallel
  profile. A dependency's choice, not ours.

---

## 7. What's already good

- **The expensive question is asked before the expensive work.**
  `decode_container_at` resolves which sample owns a CRAM record from the read-group
  number *before* building anything, so a record belonging to another sample is counted
  and dropped without a byte being copied
  ([container.rs:331-341](../../../../src/ng/read/input/aligned_reads_reader/container.rs#L331-L341)).
- **The two per-position data structures that could have been proportional to depth are
  both gated by an exact `false`.** The ceiling-loss `BinaryHeap` is behind an
  `is_empty()` on a heap that receives no push unless the hold ceiling binds, and the
  mate-overlap sort is behind `ActiveReads::may_have_mate_overlap_at`'s O(1) answer.
  Neither costs anything at a quiet position, which is almost every position.
- **The reference walk's read amplification was already found and fixed once.**
  `read_raw_bases` takes bytes out of the `BufReader`'s own buffer and tracks the logical
  file offset, where the shape it replaced issued a fresh 64 KiB `read(2)` per call and
  discarded everything past the base it wanted
  ([raw_chrom_reader.rs:215-259](../../../../src/ng/raw_chrom_reader.rs#L215-L259)).

---

## Author response convention

Address each finding by its identifier with one of: `applied in <commit>` /
`experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` /
`won't fix because …`. Two are already answered: the devirtualisation of
`try_ordinary_column` and the borrow-not-take in `GeneratorSet::next_locus` are
**experiment shows no gain — closed**, and the per-category files under
`tmp/perf_review_2026-09-02_ng-calling/` are the audit trail.
