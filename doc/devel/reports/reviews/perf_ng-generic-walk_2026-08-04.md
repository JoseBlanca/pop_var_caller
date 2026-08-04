# Performance Review: ng-generic-walk
**Date:** 2026-08-04
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** ng's generic (SNP/indel) locus generator, the new cursor-based read path it drives, the reference accessors, and the typed-region stream
**Verdict:** Apply the listed wins
**Hot-path evidence:** two CPU sampling profiles on two fixtures, a syscall census, four DHAT runs, instruction-count A/Bs on both fixtures, and an orchestrator re-verification of the combined patch on a quiet host

---

## 1. Scope and constraints

**What was reviewed.** `src/ng/locus_generation/pileup/`; the **new** cursor-based ingestion
([cursor.rs](../../../../src/ng/read/input/cursor.rs),
[sample_cursor.rs](../../../../src/ng/read/input/sample_cursor.rs),
[aligned_reads_reader/](../../../../src/ng/read/input/aligned_reads_reader/),
[region_raw_aligned_reads.rs](../../../../src/ng/read/input/region_raw_aligned_reads.rs));
[filtering.rs](../../../../src/ng/read/filtering.rs) and
[reference_free_first_filter.rs](../../../../src/ng/read/reference_free_first_filter.rs); the
reference accessors ([ref_seq.rs](../../../../src/ng/ref_seq.rs),
[raw_chrom_reader.rs](../../../../src/ng/raw_chrom_reader.rs),
[reference_info.rs](../../../../src/ng/reference_info.rs)); and the typed-region stream
([region_typing/](../../../../src/ng/region_typing/),
[tandem_repeat.rs](../../../../src/ng/tandem_repeat.rs)).

**Reviewed against.** `d19d4ab` on branch `ng-generic-perf`, whose tree
(`a3800ec02c527f8636343eef77151a2f9d536a82`) is **byte-identical to `main` at `54e5f5e`** —
verified, not assumed. Every sub-agent detached its own worktree at that SHA and confirmed the
tree hash before measuring.

**Targets and hardware.** Human WGS, one sample, HG002 30×. The generator is single-threaded;
the parallel fan-out is a later plan and out of scope. Wall time and peak RSS both in scope.

Host, re-checked this session rather than carried over:

```
$ sysctl -n hw.model hw.ncpu hw.memsize
Mac17,9
18
68719476736
$ sw_vers | head -2
ProductName:		macOS
ProductVersion:		26.5.2
```

**⚠ The 18 cores are not 18 equal cores** (owner, mid-review; verified):

```
$ sysctl -n hw.nperflevels hw.perflevel0.physicalcpu hw.perflevel1.physicalcpu
2
6
12
```

**6 high-performance cores and 12 low-energy cores.** macOS fills the fast six first and spills
under load. Three consequences run through this whole report: the **usable parallel-worker count
is 6, not 18** (and `available_parallelism()` returns 18 — measured); six agents measuring
concurrently oversubscribe exactly the six cores that matter, which is why almost every
wall-clock A/B in the fan-out came back unresolvable and why §5's numbers are instruction counts;
and it is the mechanism behind the 2026-07-31 review's H8 discrepancy (−2.6 % reported under
five concurrent agents, −1.2 % re-verified quiet).

**⚠ The fixture is not the target workload, and this still governs the report.**
`HG002_TR_v1.0.1_Tier_30x.bam` is tandem-repeat-**targeted**, not WGS: 30× *inside the TR
benchmark regions only*. The probe's own counters give 1,541,788 loci over 240,227,974 generic bp
on chr1 = **0.64 % of positions covered**. No true 30× WGS alignment exists on this host.
**Ratios transfer; microseconds do not.** Every headline finding below was therefore checked on a
second, denser fixture with a different container format — the tomato CRAM `SL4.0ch01`
(1,711,775 loci over 87,236,524 bp ≈ 2.0 %) — and findings that depend on locus density are
marked as such.

**Hot-path evidence available.** A `sample` CPU profile of the chr1 walk (three threads,
9,190 samples each) and one of the tomato CRAM walk (4,626 main-thread samples); a
`[profile.profiling]` + `#[inline(never)]` decomposition of the largest site; a counted
`read(2)`/`lseek(2)` census on three fixtures; four DHAT runs; instruction-count A/Bs on both
fixtures from four agents; and an **orchestrator re-verification of the combined patch, serial,
on a quiet host**. Raw artefacts, per-category files, patches, `BASELINE.md` and `VERIFICATION.md`
are the audit trail in gitignored `tmp/perf_review_2026-08-04_ng-generic-walk/`.

**Deliberately out of scope.** The SSR/STR generator except where it shares the cursor; the
cohort and multi-sample paths; `src/bam/` and `src/var_calling/` (frozen production — ng must not
edit it); designing the parallel fan-out (its taxes were catalogued, it was not designed).
`parity.rs`, `tests.rs`, `copy_fidelity.rs`, `mock_reference.rs` are test-only.
`cigar_cursor.rs`, `decompose.rs`, `chain_id_allocator.rs` are byte-identical production copies
enforced at compile time by `copy_fidelity.rs`, so findings against them cap at **Speculative**
and name the owner decision required.

**Categories dispatched.** All six, each in its own git worktree. `methodology` (always);
`io_and_syscalls` (32–50 % of the main thread in `read`/`lseek`); `hot_loops` (the 34.6 % site);
`allocations` (5.65 M allocations per chr21 walk); `data_layout` (3,376 lines of never-audited
new cursor code); `concurrency` (a fixed verification tail that had grown to 83 % of a chr21 run).

---

## 2. Verdict

**Apply the listed wins.** Four independent changes, one per file, compose without conflict and
are re-verified by the orchestrator serially on a quiet host after all agents had stopped:

| fixture | baseline walk | combined | change |
|---|---:|---:|---:|
| chr21 | 1.836 s | **1.114 s** | **−39.4 %** |
| chr1 | 11.620 s | **6.927 s** | **−40.4 %** |
| tomato CRAM `SL4.0ch01` | 6.867 s | **2.914 s** | **−57.6 %** |
| chr21 peak RSS | 21.32 MB | **19.10 MB** | −10.4 % |

Ranges disjoint on every fixture; all four dumps `cmp`-identical; probe counters exact;
`cargo test --lib` 2,869 passed; clippy `-D warnings` clean. Full tables in `VERIFICATION.md`.

**The last review's largest open lever is closed, and by the work it was deferred to.** H3/H4 —
the per-region BAM query re-decoding the same records 30.3× — went to the alignment cursor. It
landed, and it delivered:

| regions (chr21, `loci=256391` at every grain) | 2026-07-31 | 2026-08-04 |
|---:|---:|---:|
| 116,775 (400 bp) | 4.602 s | **1.043 s** |
| 4,671 (10 kb) | 1.914 s | **0.959 s** |
| 1 (whole contig) | 1.382 s | **0.963 s** |

The region-grain penalty was **3.33×** and is now **1.08×**; BGZF/CRAM decode fell from 43 % of
self time to **0.8 %**. Region grain is no longer a lever, and `region_query.rs` and `merge.rs`
are gone.

**Where the time went instead — the shape changed completely, and neither successor was where
the last review was looking.** The freed time did not spread evenly; it exposed two fixed costs:

1. **The reference FASTA, not the alignment file** — `read` + `__lseek` + `read_raw_bases` were
   **32.4 % of the chr1 main thread and 49.9 % of the tomato CRAM's**. Counted: chr21 issued
   **289,308 reads returning 18.96 GB to deliver 53 MB of bases — 357× amplification**, at 1.23
   reads *per locus*. This is H1, and it is most of the −39 %.
2. **`TypedRegionIterator::next` at 34.6 %** of the chr1 main thread — the single largest CPU
   site, and it still had no benchmark. Decomposed for the first time this round: it *is*
   `find_tandem_repeats`, which is 75 % Ruzzo–Tompa loop. This is H2.

**One finding is now the binding constraint on everything else, and its fix is an owner
decision rather than a win to bank.** The FASTA md5 verification is a fixed
**207 G instructions / ~11.3 s**. On chr21 it was already 83 % of the run. After the combined
patch it no longer overlaps on chr1 either: the walk falls **11.62 → 6.93 s** while the wall
falls only **11.63 → 10.64 s**. **Beyond this point, making the walk faster buys almost no wall
time.** But the measured fix (a persisted digest keyed on `(len, mtime)`) is a real weakening of
a check that exists to catch a FASTA not matching its `.fai` — see H3.

**A methodological result that changes what future numbers mean.** Two consecutive `cargo bench`
invocations of an **identical binary** produced **four of six points as "statistically
significant", three at `p = 0.00`, spanning −9.2 % to +15.9 %**. Raising `sample_size` does not
fix it and makes it worse in principle. **Criterion's `change:` line is not admissible as
cross-commit evidence on this host.** The technique that replaced it — `instructions retired`
from `/usr/bin/time -l`, floor-subtracted with `PVC_PROBE_MAX_LOCI=1` — is near-deterministic
under load and carried four of the six agents' results.

---

## 3. Measurement plan

Ordered by what unblocks what.

1. **Decide H3 (the verification tail) before optimising the walk further.** It is now the
   binding constraint on every fixture. The question is not "is it faster" but "what invalidation
   rule is acceptable" — §5 H3 has the measured fix, its negative controls, and the two ways it
   is weaker. **Threshold:** none; this is an owner decision, not an experiment.
2. **Get a real 30× WGS BAM for HG002 and re-run the baseline.** Everything here is 0.64 %-covered
   data. **Threshold:** if the per-locus reference-fetch rate holds, H1's ratio transfers intact;
   if locus density rises, H1 gets *larger* (it is per-locus) and H2 roughly holds (it is per-base
   over the contig, independent of coverage).
3. **Adopt `instructions retired` + floor subtraction as the project's A/B instrument**, and
   record the discipline (§4 B4). **Threshold:** two consecutive unchanged runs report zero
   significant points, twice in a row, before any criterion `change:` line is quoted again.
4. **Land the typed-region bench** written this round (§4 B3) and use it to gate H2 and anything
   downstream of it. It needs no BAM and has zero I/O in the timed body.
5. **Build and measure H5** (`folded_reads` as a sorted `Vec`). It is priced at −3.7 % by
   instruction calibration and was measured at −3.9 % in 2026-07-31 from an independent
   direction, but it has not been built on today's post-H6/H7 tree. It is the largest remaining
   lever in the pileup half.
6. **Then re-profile.** After H1+H2 the walk is ~40 % shorter and the profile will have moved
   again — as it did this round. Do not queue further code work against the *old* shape.

**The reproducible commands.**

```
REF=$HOME/genomes/h_sapiens/gca_grch38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
BAM=…/benchmarks/ssr_hg002/bam/30x/HG002_TR_v1.0.1_Tier_30x.bam
TREF=$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa
CRAM=…/benchmarks/ssr_tomato1/crams/SRR5079860.p1.bench.cram

/usr/bin/time -l ./target/release/examples/ng_generic_walk_probe "$REF" "$BAM" chr21
# knobs: PVC_PROBE_MAX_LOCI=n  PVC_PROBE_WHOLE_CONTIG=1  PVC_GENERIC_REGION_CHUNK_BP=n
#        PVC_PROBE_MAX_RECORD_SPAN=n
# chr21 walk ≈1.8 s is the iteration fixture; chr1 ≈11.6 s; tomato CRAM ≈6.9 s.
```

**Note the fixture paths are `$HOME`-relative on this host** — there is no `/genomes`. The
prompt's paths and the 2026-07-31 report's both elide it.

**The correctness gate for every experiment**, and it caught nothing this round only because
every agent ran it: all four dumps byte-identical by `cmp` (**not** line counts) —
251,792 / 4,406 / 1,718,914 / 11,945 lines, md5s in `BASELINE.md` — plus probe counters exact
(`loci=236081 observations=251786 reads_admitted=54709`).

---

## 4. Build / toolchain configuration

**The `[profile.*]` audit is still a clean pass** — re-verified, not re-derived. Fat LTO,
codegen-units = 1, panic = abort, `debug = "line-tables-only"`, `[profile.soak]` arming
assertions with `overflow-checks`, `[profile.profiling]`, `[profile.bench]`, per-target
`target-cpu`, `unsafe_code = "forbid"`, and the toolchain pin **in effect** (`rust-toolchain.toml`
says `1.95`; `rustc --version` reports `1.95.0 (59807616e 2026-04-14)`). No change recommended.
PGO was considered and deliberately **not** filed — its precondition is a stable workload, and the
generator is still moving weekly. One caveat worth a comment: `[profile.profiling]` sets
`lto = false, codegen-units = 16`, so CPU self-time taken under it does not transfer to release;
this review's baseline profile was correctly taken on `--release`.

**B1 — `cargo test --release` is red on a clean tree with nine failures in three root causes, and
the 2026-07-31 report both under-counted and misclassified them.** Measured:

```
test result: FAILED. 2860 passed; 9 failed; 5 ignored; 0 measured; 0 filtered out; finished in 7.28s
```

Debug passes 2,869; 2,860 + 9 = 2,869 exactly, so these nine are the *entire* profile difference.
The three classes need three different fixes and one blanket `ignore` would hide a real defect:

- **Five are genuine `debug_assert!` tests** (`note: test did not panic as expected`) — four of
  them never previously listed.
- **Two panic in release for a different reason with a stale `expected`** — `left_align_repeated`
  and `left_align_structured` get `range start index 10 out of range for slice of length 6`
  against an expected substring `reference_offset`. `left_align_structured`'s own doc comment
  claims the release panic is real; the test does not check it.
- **One is not a `debug_assert!` at all.** `ssr_marginal_sequence::…the_epsilon_endpoints_…`
  fails on `got 0.012345679012345678` against an **absolute** `1e-18` tolerance — i.e. ~1e-16
  *relative*, at the f64 epsilon boundary, sensitive to fat-LTO/`target-cpu` rounding.
  **2026-07-31 filed this as a debug_assert test; it is not**, and ignoring it in release would
  hide a genuine profile sensitivity. Widen to a relative tolerance.

The run also **stops at `--lib`**, so the release status of integration tests and examples is
unknown; validate any fix with `--no-fail-fast`. The cheapest fix for the first two classes is to
run the suite under the profile the project already built for it: `cargo test --profile soak --lib`.

**B2 — `benches/psp_writer_perf.rs:386` does not merely abort `--all-targets`; the benchmark has
never produced a number.** `cargo bench --bench psp_writer_perf -- flush_block_one` panics
identically: `index out of bounds: the len is 3300000 but the index is 3300000`. The prime loop
exhausts all records without reaching the 16 MiB flush boundary and indexes one past the end. The
panic is the *lucky* outcome — a `.get()` there would have silently measured a non-flush for
months. Add an assertion that the prime loop reached the boundary, then raise `NUM_RECORDS_SNP`.
(psp is out of scope for edits; the bench is in scope, and it is a hard trap for exactly the
workflow a perf review uses.)

**B3 — the typed-region stream had no benchmark; one was written and is ready to adopt.**
`benches/ng_region_typing_perf.rs` (preserved in the audit trail): in-memory reference, **zero
I/O in the timed body**, `black_box` on both sides, and a verification assertion outside the
timed body pinning both `regions > 0` and `counts.ssr_loci > 0` — the latter because a fixture
whose planted tracts fell below a copy floor would skip classification entirely and report a fast,
meaningless number. It already produced a directional result: one 100 kb window runs at
44.5 MiB/s against 40.7 MiB/s for ten windows, so the 34.6 % is **per-base scan work, not
per-window setup**.

**B4 — criterion's `change:` line is not admissible as cross-commit evidence on this host, and
the project has no written measurement discipline.** The null experiment (§2) is the evidence.
Recommend `doc/devel/measurement_discipline.md` stating: the 6+12 core split;
`available_parallelism()` returns 18 and is **not** the usable worker count; wall-clock A/Bs need
alternating binaries and ≥8 runs a side; report **raw values and min**, not just medians
(bimodality is the signature of core-class migration); the noise floor is ~2 % quiet and ~15 %
under load (both measured); prefer `instructions retired` with a `PVC_PROBE_MAX_LOCI=1` floor
subtraction; and `--features dhat-heap` **always** with `--target-dir target-dhat`.

**B5 — the dhat feature-varying-build hazard is still documented nowhere, and the newly-committed
standing instrument now propagates it.** `grep -rn "target-dhat"` returns **one hit repo-wide**:
the 2026-07-31 report's own unapplied recommendation. All eight `examples/dhat_*.rs` still document
the bare form — and `examples/ng_generic_walk_probe.rs:5` was committed since with the hazardous
invocation, while its own comment ~60 lines below *names the incident it causes*. Add
`--target-dir target-dhat` to all nine usage blocks and `target-dhat/` to `.gitignore`.

**B6 — the standing instrument needs three things it does not have.** Its counters are genuinely
good (dispatcher-sourced, destructured by name so a transposed pair cannot compile). But it has
**no mode isolating the two halves of the walk** (`PVC_PROBE_WHOLE_CONTIG` is a bound, not an
isolation — it changes the loci, 256,391 against 236,081); `render()` **echoes none of its four
knobs**, so pasted output is un-attributable and `PVC_PROBE_MAX_RECORD_SPAN` silently changes what
a long deletion can reach; and it prints **no contention signal**, which on this host is the
difference between a regression and a descheduled run.

---

## 5. Code-level findings

### Hot-path

**H1: [src/ng/raw_chrom_reader.rs:205-235](../../../../src/ng/raw_chrom_reader.rs#L205-L235),
[:364-367](../../../../src/ng/raw_chrom_reader.rs#L364-L367) — the reference is re-read from the
kernel once per locus: 289,308 reads move 18.96 GB to deliver 53 MB of bases, a 357×
amplification.** *Confidence: High.*

- **Hot-path evidence.** Profile: `read` + `__lseek` + `read_raw_bases` = 2,982 of 9,190
  main-thread samples on chr1 (**32.4 %**) and 2,339 of 4,626 on the tomato CRAM (**49.9 %**),
  all of it under `WindowedRefSeq::fetch_into`, none of it the alignment file. Then **counted**,
  not inferred, with an instrumented `CountedFile` newtype:

  | fixture | reads | bytes read | lseeks | bases delivered |
  |---|---:|---:|---:|---:|
  | chr21 | 289,308 | 18,960,089,088 | 288,842 | 53,028,593 |
  | chr1 | 1,895,837 | ~124.2 GB | 1,893,349 | — |
  | tomato | 1,920,794 | ~125.9 GB | 1,919,887 | — |

  **1.23 reads and 1.22 lseeks per emitted locus**, every one of them exactly 65,536 bytes.
- **Mechanism.** `OpenPileupRecordTable::open_or_get`
  ([open_record.rs:1254](../../../../src/ng/locus_generation/pileup/open_record.rs#L1254)) fetches
  the reference per opened record, typically `span` = 1 base. `RawChromReader::fetch` calls
  `append_forward(1)`, which `seek_to`s and then calls `read_raw_bases`, which issues
  `reader.read(read_buf)` into a **64 KiB scratch**, copies the one base wanted, and **discards
  the rest**. The next locus seeks back and re-reads the same chunk. The scratch was never a
  buffer — nothing was retained between calls. **The contrast is internal and decisive: the
  typed-region stream uses the same accessor on the same file in 100 kb windows and costs 103
  samples.** The difference is entirely the request size.
- **Fix, built and measured.** Wrap the `File` in a `BufReader<File>`, read via
  `fill_buf`/`consume`, and track the logical file offset so `seek_to` is a no-op when already
  positioned. Result: chr21 **289,308 → 1,452 reads** and **288,842 → 720 lseeks**; chr1
  1,895,837 → 8,566; tomato 1,920,794 → 3,234. Three independent statistics agree — `sys` time
  chr21 0.72 → 0.20 s, chr1 3.84 → 0.25 s, tomato 5.25 → 0.12 s; instructions retired
  (floor-subtracted) chr21 −22.0 %, tomato −40.8 %, a stable ~9,600 instructions per eliminated
  read+lseek pair on both; post-fix profile `read`+`__lseek` = **0.5 %** of the main thread, was
  31.4 %.
- **Complexity cost.** One `BufReader`, one `file_pos: u64` field, and one new invariant: the
  tracked offset must stay in sync with the reader's real position, so every seek path must go
  through `seek_to`. No `unsafe`, no dependency, no build flag. **But see the owner decision in
  §6:** this file documents itself as a diffable copy of production's
  `ManualEvictChromRefFetcher`, and this is its second behavioural divergence.
- **In the §2 combined patch.**

**H2: [src/ng/region_typing/mod.rs:1055](../../../../src/ng/region_typing/mod.rs#L1055)
(`scan_and_absorb`) — the scanner is asked for intervals its only consumer's copy floor is
guaranteed to discard.** *Confidence: High.*

- **Hot-path evidence.** `TypedRegionIterator::next` is **3,176 of 9,190 main-thread samples
  (34.6 %)** on chr1 and ~24 % on the tomato CRAM — the largest single site in the profile.
  Decomposed for the first time via a `[profile.profiling]` build with `#[inline(never)]`:
  `maximal_scoring_subsequences` 3,905 samples, `perf_post_pass` 618, `is_primitive_motif` 103,
  and `prefilter` + `classify` + `absorb` + `emit_into` **together under 1 %**. So the site *is*
  `find_tandem_repeats`, and that is ~75 % Ruzzo–Tompa loop.
- **Mechanism.** `scan_and_absorb` passes `min_copies = 2` while `prefilter`, the only consumer
  of `detections` **on this path**, immediately drops everything below the per-period floor
  `[6,4,4,3,3,3]`. Hoisting the weakest consumer floor into the detector stops the intervals that
  cannot survive it from ever being materialised or motif-checked.
- **Output-neutrality — I checked the argument against the source rather than accepting it.**
  Removing a below-`F` interval `a` can only matter if it stops eliminating some longer-period
  `iv`. Elimination requires `2·len(a) ≥ len(iv)`
  ([tandem_repeat.rs:584](../../../../src/ng/tandem_repeat.rs#L584)); the post-pass only considers
  proper divisors, so `p ≥ 2q`
  ([:571](../../../../src/ng/tandem_repeat.rs#L571)); and the scanner's floor is
  `copies >= params.min_copies` ([:522](../../../../src/ng/tandem_repeat.rs#L522)). Then
  `copies(a) < F ⟹ len(a) < F·q ⟹ len(iv) < 2F·q ⟹ copies(iv) < F`, so `prefilter` drops the
  survivor anyway. **All three premises hold.**
- **Measured.** Intervals materialised **9,659,219 → 2,694,953**; `is_primitive_motif` calls
  **10,693,464 → 2,833,017**; region output identical
  (`generic=102938 ssr=82747 bundle=19912 satellite=278 repeat_bp_with_no_locus=740144`).
  Instructions retired on the isolated site: **−15.2 %** (chr21) and **−16.9 %** (tomato),
  disjoint. In the whole chr21 walk, floor-subtracted: 25.222 G → 24.311 G = **−3.6 %** (the
  author took the smaller full-walk figure over the isolated one; correct).
- **Complexity cost.** One `effective_scan()` method deriving the floor from the criteria, one
  call-site change. **One correction to the patch as written:** its doc comment says `prefilter`
  is the *only* reader of `ScannedWindow::detections`. That is false in general —
  `RegionScanner` reads them raw with no copy floor
  ([tandem_repeat.rs:1121](../../../../src/ng/tandem_repeat.rs#L1121)) — and is true only on the
  path changed. The author correctly left the other path un-raised, which is what let it serve as
  an independent oracle. **Restate the invariant as path-scoped and pin it with a test**, or a
  later tidy-up that pushes the floor into `scan_window` will silently change output.
- **In the §2 combined patch.**

**H3: [src/ng/reference_info.rs:1018](../../../../src/ng/reference_info.rs#L1018) — the FASTA md5
verification is a fixed ~207 G instructions / ~11.3 s, and after H1 it is the binding constraint
on every fixture.** *Confidence: High.* **The finding is banked; the fix is an owner decision.**

- **Hot-path evidence.** A `PVC_PROBE_MAX_LOCI=1` floor run isolates it exactly: baseline
  `seconds=0.118, real=11.30`, 207.47 G instructions — **99 % of a floor run and 88 % of a full
  chr21 run**. It is CPU-bound (md5 57 %, the per-byte loop 41 %, `read` 1.9 %), so only *not
  recomputing* it helps. On chr21 it was already ~83 % of the wall before this review.
  **The orchestrator's own re-verification is the sharpest statement of it:** with the combined
  patch, chr1's walk falls **11.620 → 6.927 s (−40.4 %)** while its wall falls only
  **11.625 → 10.640 s (−8.5 %)**. It used to overlap; it no longer does.
- **The topology argument, now measured.** Against 0/3/6/12 competing spinners, instructions
  retired span **0.06 %** while `real − seconds` grows **12.27 → 15.91 s**. Same work, more wall:
  the tail is contention for the six performance cores, so "it overlaps, so it is free" fails
  once a fan-out saturates them — which the 2026-07-31 review could not argue, having recorded
  the box as 18 uniform cores.
- **Fix, built and measured.** A persisted `<fasta>.refdigest` sidecar guarded by the same
  `(len, mtime)` tuple `cache_key` already uses in memory. chr21 `real` **12.09–13.02 → 1.94–2.41 s**,
  walk `seconds=` untouched, disjoint; a third arm (same binary, cache disabled) sits on the
  baseline, ruling out binary layout. Negative controls pass: a stale length and a stale mtime
  each cost the full re-read, and six unit tests pin that a re-wrapped FASTA is still caught with
  `FastaFaiMismatch { field: "line_bases" }`.
- **Why this is not banked as a win.** The stat-based rule **is a real weakening across process
  boundaries**: it misses an in-place edit whose mtime is restored and whose length is unchanged.
  The check exists to catch a FASTA that does not match its `.fai`. It also writes a second file
  beside the reference, contradicting spec §3.11 and the function's own doc comment. The author's
  recommendation, which I endorse: **ship it opt-in behind the `PVC_REF_DIGEST_CACHE` knob already
  in the diff** — the edit-measure loop the tail actually hurts is exactly the one a developer
  opts into, which buys most of the value with none of the weakening. **Route the always-on
  variant through a correctness review, not this one.**

**H4: [src/ng/locus_generation/pileup/open_record.rs:928](../../../../src/ng/locus_generation/pileup/open_record.rs#L928)
— `OpenPileupRecordTable` is a `BTreeMap` that never holds more than two entries.**
*Confidence: High.* Instrumentation shows the table peaks at **2 records** on both fixtures while
236,081 (chr21) and 1,711,775 (tomato) records are opened and closed through it — so every
insert, lookup and removal pays `BTreeMap` node machinery to manage a pair. A key-sorted `Vec`
gives the same ordered semantics: **−1.61 % of the chr21 instruction budget and −0.98 % of
tomato's, disjoint on both**, dumps `cmp`-identical, `--lib` 2,869 green. `open_record.rs` is
already released from `copy_fidelity.rs`, so **no owner decision is needed**. *Complexity:* one
container swap behind the existing accessors. **In the §2 combined patch.**

**H5: [src/ng/locus_generation/pileup/open_record.rs:389](../../../../src/ng/locus_generation/pileup/open_record.rs#L389),
`:545`, `:640`, `:1728` — `folded_reads` as a sorted `Vec` instead of a seeded hash map.**
*Confidence: High. Carried from 2026-07-31 H5, re-priced, still not built.* The original
measured **−3.9 % on chr21, 6/6 paired runs, output identical**, and was deferred only for a
rebase against H6/H7 — which have since landed. Re-priced this round from an independent
direction: counters give `folds=6,339,812`, `finalise_calls=363,312`, `finalise_ids=6,326,850`
(mean 17.4 per record), and calibrating **37 instructions per avoided hash+probe** against the
25.22 G chr21 walk prices it at **−3.7 %**. Two methods, two directions, same answer.
**It is the largest remaining lever in the pileup half.** *Complexity:* one `FoldedReads` newtype;
it *retires* an invariant (the two sorts that exist only to recover determinism from a
per-process-seeded map) rather than adding one. **Not built — §3 item 5.** Note the smaller
`entry`-over-`remove`+`insert` change measured this round (−0.94 %) is **subsumed by** H5, not
additive; take H5 and drop it.

### Likely

**L1: [src/ng/read/input/region_raw_aligned_reads.rs:203](../../../../src/ng/read/input/region_raw_aligned_reads.rs#L203)
— the sorted early stop deep-clones a whole `RecordBuf` once per region.** *Confidence: High.*
DHAT names it at **673,585 blocks = 11.9 % of every allocation the walk makes**, 59.7 MB — 7.3
allocations per region, because `RecordBuf::clone` rebuilds the name, CIGAR, sequence, quality
scores **and the whole aux-tag table with one heap string per tag**. It is *not* visible in the
CPU profile (`AlignmentCursor::next_read` self is 17 samples, 0.2 %), which is why it was measured
rather than asserted. Two buffers trading places removes it: **−11.93 % of all blocks**;
instructions retired **−1.01 % (chr21) / −0.62 % (tomato), disjoint on both**. *Complexity:* one
field changes type, one `RecordBuf` parking slot, two `mem::swap`s replace one `clone`. It rests
on the contract the method already documents — *"After `Ok(false)`, `buf` holds an unspecified
record the caller must not read"* — and **I verified the single call site
([cursor.rs:715](../../../../src/ng/read/input/cursor.rs#L715)) returns `None` immediately
without reading**. The honest cost: it promotes that contract from benign-if-violated to
load-bearing, with nothing enforcing it. **In the §2 combined patch** (at −1 % it would not be
worth merging for speed alone; it is worth merging because it is free, deletes 12 % of the walk's
allocations, and is byte-identical).

### Speculative

**S1: [src/ng/locus_generation/pileup/cigar_cursor.rs:60](../../../../src/ng/locus_generation/pileup/cigar_cursor.rs#L60)
— `EventsAt`'s second inline slot is used in 0.036–0.063 % of contributions**, but at 6.4 KB the
buffer is L1-resident either way, so the cache-miss mechanism does not fire. Re-filed from
2026-07-31 S3 with occupancy measured. Capped here because the file is copy-fidelity frozen;
releasing it is an owner decision.

**S2: [src/ng/tandem_repeat.rs:208](../../../../src/ng/tandem_repeat.rs#L208) —
`RepeatInterval` is 24 bytes for coordinates that are window-local** and would fit in less.
Inside the 34.6 % function, but the two direct attacks on that function's memory traffic (N1, N2)
both came back *slower*, so the prior here is poor.

**S3: the Ruzzo–Tompa loop's remaining lever is algorithmic, not a perf edit.** A microbenchmark
puts scoring at 0.49 ns/pos (14 % of the loop) and the random branch at 0.15 ns/pos, leaving
**~86 % in the stack machinery at ~43 cycles per push**, with only 1.011–1.33 hops per search —
so the search is not the problem and the stack is never walked. After H2 cuts what is pushed by
72 %, what remains needs a different algorithm, not a tighter one.

**S4: [src/ng/ref_seq.rs:590](../../../../src/ng/ref_seq.rs#L590) and
[generator.rs](../../../../src/ng/locus_generation/pileup/generator.rs) — `PileupGenerator` is
`!Send` and so is `Arc<WindowedRefSeq>`** (proven by compiler output: the `RefCell` inside).
A fan-out must construct the whole per-worker stack inside the worker; only the `Send + Sync`
triple `(PathBuf, Arc<ContigList>, Arc<fai::Index>)` can cross. The design is right; the signature
does not say so. Worth a constructor contract before a caller assumes otherwise. Tax catalogued
at **N = 6**, not 18.

### Note

**N1 — narrowing `RtSeg` 40 → 32 bytes is byte-identical and retires 0.59 % *more* instructions**;
40 → 20 bytes is **+8…19 % slower**. Instrumenting the search explains it: **1.011 hops per
search**, so the 982 KB stack is never walked and there is no locality to win. Recorded so the
"shrink the struct in the hot function" reflex is not re-run. Caveat stated honestly by its
author: a null *instruction* delta cannot refute a *locality* mechanism on a host with no PMU.

**N2 — run-length-compressing the Ruzzo–Tompa score stream is +25…33 % slower.** Implemented and
refuted.

**N3 — 2026-07-31's L1 (boxing `ReadWitness`) is dead, and now has its missing denominator.**
Sizes are unchanged (32/112/64/120), but the walk holds **one finalised locus at a time, carrying
≤ 15 observations** — so the whole saving is **360 bytes**, or 2.2 KB at the real worker count of
6. It wins in no currency on the single-sample target. Its 21 % lives only in a cohort merge,
multiplied by *sample* count — which is out of scope here.

**N4 — the cursor's residency worry is closed.** `kept` peaks at **203 reads = 104 KB, 0.4 % of
peak RSS**. Every live structure in the data-layout scope sums to **under 1.2 MB** against peaks
of 25 MB (chr21) / 232 MB (tomato), of which 83 % / 89 % is already paid at the floor. **No
struct-size change in this scope can move peak RSS.**

**N5 — 2026-07-31's L2 is measured dead.** `io_open_contig_calls=4` per contig,
`io_fai_parses=0`: `make_reference` fires once per *chromosome*, not per region. The cursor's D1
already fixed it.

**N6 — the alignment index is not deep-copied per cursor.** `open_bam.rs:439`'s
`self.index.clone()` looks like a copy of a 7.3 MB parsed BAI; `AlignmentIndex` is all-`Arc`, so
it is three refcount bumps. DHAT agrees: one parse, 41,293 blocks, no second copy. **No finding,
and none under the fan-out either.**

**N7 — there is no peak-RSS allocation lever.** At-gmax heap is 27.9 MB, dominated by the parsed
index (26.2 %, shared), the SAM header (12.6 %, behind an `Arc`) and `tandem_repeat`'s scoring
stacks (26.7 %, the recorded 2026-07-31 non-win). The walk's own per-locus state is **724 KB,
2.6 %** — spec §7's depth-boundedness confirmed again from a new direction.

**N8 — allocation count is not the currency, and now there is an exchange rate.** Removing
751,447 allocations removed 255.6 M instructions = **340 instructions per allocation**. Against a
25.2 G-instruction chr21 walk making 5.65 M allocations, **1 M allocations ≈ 1.35 % of the walk**;
moving it 5 % needs **65 % of every allocation deleted**. That rate independently *predicts* the
2026-07-31 lazy-`Record` result (36.4 % of allocations → −0.5 %). **The hostile prior is no longer
a caution, it is arithmetic** — and it prices the top eight remaining allocation sites at under
1 % of the walk each, none of which removes a hash, a sort or a group.

**N9 — a QoS hint on the verification thread is not recommended.** Considered and measured; the
justifying measurement came back null, and it points the wrong way — it would lengthen the tail
in exactly the exposed case H3 is about.

**N10 — `available_parallelism()` returns 18 on this host and the usable worker count is 6.**
Measured directly. `src/var_calling/pipeline.rs:128` is the existing precedent in frozen
production; ng has no call site yet, which makes now the cheap moment to decide. State in ng's
fan-out spec that workers default to the high-performance core count.

**N11 — the 2026-07-31 non-wins N1–N6 stand and were not re-proposed.** Narrowing
`max_record_span`; caller-side region coalescing; the `folded_reads` map free-list (0 %); noodles
lazy `Record` (−0.5 %); `tandem_repeat` RT-stack reuse (≈0 %); `segment_criteria`'s per-candidate
`upper()`; mimalloc (+2.4 % wall, +45 % RSS); `ActiveRead` correctly AoS at 184 bytes;
`WitnessedLocusPositions`' two-run inline capacity at the SmallVec floor — the last one confirmed
again this round by the walk's own occupancy data.

---

## 6. Out-of-scope observations

- **⚠ Owner decision — H1 makes `raw_chrom_reader.rs` diverge from the production file it exists
  to mirror.** Its header says it is "ng's copy of `fasta::ManualEvictChromRefFetcher`, with
  exactly one behavioural difference" and that "this copy stays diffable against it on purpose"
  ([raw_chrom_reader.rs:3-41](../../../../src/ng/raw_chrom_reader.rs#L3-L41)). It is **not**
  compile-time locked by `copy_fidelity.rs` (which covers only `chain_id_allocator.rs`,
  `cigar_cursor.rs`, `decompose.rs`), and 2026-07-31's H1 already made a second change to it. H1
  here is a third and much larger one. Either accept the divergence explicitly and update the
  header, or port the fix to production first. **This needs an owner call before H1 lands.**
- **Production almost certainly has the same 357× read amplification.**
  `ManualEvictChromRefFetcher` has the identical unbuffered read-and-discard shape. That is a
  real opportunity in `src/pileup/`, which is frozen for ng and out of scope here — worth its own
  issue.
- **`md5::compress` is scalar** and would be a SIMD or parallel-hash candidate *if* per-run
  verification survives H3. Prefer H3's caching, which removes the work rather than speeding it.
- **The probe issues two verification calls**, which is why the profile shows a third thread
  parked in `__psynch_mutexwait` for its whole life. Not a defect — two entry points into the
  same setup sharing one cache, so the second is a parked thread rather than a second 3.1 GB read,
  and the probe's whole setup costs ≲9 ms. Recorded because it was misread once already.
- **`data_layout` routed a redundant `sort_by_key` over 9.65–23.6 M already-sorted elements**
  inside the hottest function. H2 cuts that population by 72 %, so re-measure after H2 rather than
  before.

## 7. What's already good

- **The alignment cursor did exactly what it was deferred to do, and the evidence is
  unambiguous** — region grain went from a 3.33× penalty to 1.08× and BGZF decode from 43 % of
  self time to 0.8 %, with `loci=256391` identical at every grain. Deferring H3/H4 to a designed
  redesign rather than patching the query was the right call.
- **`RegionRawAlignedReads::read_next` documents the precondition that makes L1's optimisation
  legal** — *"After `Ok(false)`, `buf` holds an unspecified record the caller must not read"* —
  written before anyone needed it, and the single caller honours it. The contract is why a
  measurable win was available without an invariant change.
- **The probe's counters are a model for how to make a fast run checkable.** They come from the
  dispatcher rather than the instrument's own bookkeeping, and the eight generator counters are
  destructured by name so a transposed pair cannot compile — which is what let every agent this
  round prove its variant did the same work rather than less of it.
- **Spec §7's "bounded by depth, not by region length" keeps surviving new attacks on it.** This
  round it was confirmed from a third direction: the walk's own per-locus state is 724 KB of a
  27.9 MB at-gmax heap, and the cursor's live-read set is 104 KB.

---

## Author responses

Owner decisions taken 2026-08-04, in conversation, and applied the same day.

| finding | response |
|---|---|
| **H1** — reference re-read once per locus (357× amplification) | **applied.** The `BufReader` + logical-offset fix is in the working tree. Re-measured after landing: chr21 walk **1.836 → 1.281 s, −30.2 %**, matching the sub-agent's −30.1 % independently. |
| **H1's owner decision** — divergence from the production copy | **resolved: diverge, and do not fix production.** *"Don't try to fix production, we are betting for the ng caller."* The file header no longer claims to be a diffable copy; it now names both divergences, states that production has the same waste and is deliberately not being fixed, and warns against reading a diff against `ManualEvictChromRefFetcher` as a bug list. |
| **H3** — the FASTA verification tail | **the proposed sidecar cache is rejected; a switch was built instead.** Owner: *"add a configuration parameter to the fasta reader, by default the checksum is done, but the user can choose to skip it."* This is the better design and the report's recommendation was worse — see below. |
| **H2** — scanner asked for intervals its consumer discards | **applied**, with the safety argument corrected (see below) and pinned by two new tests. |
| **H4** — two-entry `BTreeMap` | **applied.** |
| **L1** — per-region `RecordBuf` clone | **applied.** |
| **H5** — `folded_reads` as a sorted `Vec` | **built, applied, then REVERTED in `5d35490`.** The −3.7 % price was a calibration on one depth. On the GIAB depth sweep it is **−3.2 % at 30× and +16.4 % at 300×** — the same cliff, and the same mechanism (`_platform_memmove` second overall at 300×), that made production revert this shape in `perf_pileup_2026-05-12.md`. Not worth a 16 % cliff for a 3 % gain. The **open-record table stays a sorted `Vec`**: it holds at most two entries at any depth, so it has no depth term. |

### All five changes together, re-measured serially on a quiet host

Baseline is a binary built from `d19d4ab`; binaries alternated within each pair.

| fixture | baseline walk | applied | change |
|---|---:|---:|---:|
| chr21 (8 pairs) | 1.837 s | **1.071 s** | **−41.7 %** |
| chr1 (4 pairs) | 11.556 s | **6.599 s** | **−42.9 %** |
| tomato CRAM (4 pairs) | 6.819 s | **2.679 s** | **−60.7 %** |
| chr21 peak RSS | 21.39 MB | **18.68 MB** | **−12.7 %** |

**⚠ That table includes H5, which was then reverted.** The shipped state, re-measured across
the GIAB depth sweep on chr21, is a win at every depth with no crossover:

| depth | baseline | shipped | change |
|---|---:|---:|---:|
| 5× | 1.353 s | 0.765 s | **−43.5 %** |
| 30× | 1.831 s | 1.092 s | **−40.4 %** |
| 50× | 2.178 s | 1.373 s | **−37.0 %** |
| 300× | 7.267 s | 6.070 s | **−16.5 %** |

The benefit falls with depth because most of it is fixed cost (the reference reads and the
repeat scan), which deeper coverage amortises against more per-read work.

Raw values, chr21: base `1.846 1.833 1.833 1.844 1.834 1.828 1.840 1.841`;
variant `1.072 1.076 1.074 1.067 1.069 1.093 1.063 1.064`. **Ranges disjoint on every
fixture.** All four dumps `cmp`-identical; probe counters exact; `cargo test --lib` **2,874
passed**; clippy `-D warnings` clean; the cross-process determinism test
(`ng_emits_the_same_bytes_in_a_second_process`) run explicitly and green.

### H2's safety argument was overstated in the patch, and is now correct

The patch's doc comment claimed `prefilter` is the only reader of `ScannedWindow::detections`.
**That is false in general** — `RegionScanner` reads them raw, with no copy floor
([tandem_repeat.rs:1121](../../../../src/ng/tandem_repeat.rs#L1121)) — and true only on the
path the change touches. The comment now says so, warns against moving the raised floor into
`scan_window` or `TypedRegionConfig::scan` (which would silently change what `RegionScanner`
emits), and the three premises are each cited to the line that implements them.

Two tests pin it, and the second exists because the first could otherwise pass vacuously:
`raising_the_scan_floor_cannot_change_the_prefiltered_set` compares raised against un-raised
through `prefilter`, asserting first that the raise *removes something* from the raw set;
`raising_the_scan_floor_too_far_does_change_the_prefiltered_set` raises one copy past the
consumer floor and requires the comparison to notice. **The fixtures already in this module
would not have exercised it** — they are built from long clean tracts over a cycled 16-mer
filler, which contains no short intervals for the floor to drop, so a new
`contig_with_two_copy_micro_repeats` fixture plants exactly-two-copy repeats of every scanned
period.

### H5 was a re-implementation, not a rebase, and it had to avoid a known failure

The same swap was tried in **production's** walker and reverted
(`perf_pileup_2026-05-12.md` H1: +1.3 % mean, four of eight fixtures 3–12 % worse, worst
`multi_op/5000` at +11.8 %) for two named reasons — `Vec::remove`'s shift on re-fold, and byte
inflation from doubling past capacity 32. The constant's own doc comment recorded that history
and said to keep the hash map. Both causes are addressed rather than ignored:

- **The re-fold no longer shifts.** Production's attempt kept the map's `remove`-then-`insert`
  shape, which in a `Vec` shifts the tail down and then straight back up for a key that was
  already present. `FoldedReads::fold` replaces in place; `fold_read_into_record` now *reads*
  the previous state (both fields it needs are `Copy`) instead of removing it. A shift happens
  only when a genuinely new read arrives.
- **Depth is the axis that decides it**, and ng's is two orders of magnitude below the fixture
  that regressed: the shift is `O(reads per record)`, and ng's mean is **17.4**. That number is
  now in the type's doc comment, with an explicit note that **if ng's depth regime changes, the
  decision changes with it** — the shape is not universally better and should not be copied
  back to production on the strength of this result.

The ordering that the two deleted `sort_unstable` calls used to restore is now **structural**:
`FoldedReads` keeps entries in ascending `read_id`, which is exactly what both consumers need
(`q_sum` is an `f64` sum, so an order that varies run to run moves emitted bytes). The
cross-process determinism test therefore now guards against a future unordered container rather
than against deleting a sort, and its doc says so.

### Why the switch beats the cache the review proposed

The review measured a persisted digest keyed on `(len, mtime)` and flagged it as a correctness
weakening. The owner's alternative removes the weakening rather than trading it down, and it is
strictly better on three counts:

1. **No guess.** The cache asserted "size and timestamp suggest nothing changed", which can be
   wrong silently. The switch asserts "you told me not to check". When the check runs it is the
   full check; when it does not, someone asked for that.
2. **No second file and no stale-state rule.** The spec §3.11 conflict and the read-only /
   shared reference-directory problem both disappear, as does the invalidation logic — the part
   of the cache least likely to stay correct.
3. **It targets the case that actually hurts.** The 11 s is *fixed*, so it rounds to nothing on a
   real pipeline run and dominates only short ones — which are overwhelmingly measurement runs.
   The cache would have charged every real run a weaker guarantee to buy back time real runs do
   not need. Owner: *"If we create a SNP caller in which the time is defined by the time we need
   to read the FASTA file we have succeeded in a spectacular way."*

### What was built

- **[`ReferenceCheck`](../../../../src/ng/reference_info.rs)** — `VerifyAgainstIndex` (the
  `Default`) or `TrustIndexWithoutChecking`. An enum rather than a `bool`, because `verify: false`
  at a call site says nothing. Both variants document what is given up, in the reader's terms.
- **The parameter is on `read_reference_verifying_or_creating_fai`**, the batteries-included
  entry point, which already returned `Option<VerificationHandle>` — so `None` keeps meaning
  "nothing to await" and no caller shape changed. The primitive named
  `read_fai_verify_in_background` was left alone: its name is a promise.
- **Skipping is unavailable when the `.fai` is missing**, and that is documented as deliberate,
  not an oversight — writing the index requires reading the FASTA, so there is nothing to skip.
- **A CLI flag** on the experimental typed-regions tool, `--trust-reference-index`.
- **One shared rule for the dev tools**,
  [examples/shared/reference_check.rs](../../../../examples/shared/reference_check.rs), read from
  `PVC_TRUST_REFERENCE_INDEX`. **The default there is still to check**, because those tools
  produce the byte-identity dumps this review's gate depends on. A misspelled value is a usage
  error (exit 2), not a silent default — a typo that quietly meant "check" would show up only as
  an unexplained 11 seconds.
- **The probe now prints `reference_check=`**, because `seconds` is not interpretable without it.
  This also closes part of §4 B6.

### Measured after landing, chr21, 4 alternating pairs

```
1  CHECKED reference_check=verified_against_fai seconds=1.292 real 10.42
   SKIPPED reference_check=trusted_unverified   seconds=1.261 real 1.26
2  CHECKED …                                    seconds=1.287 real 10.40
   SKIPPED …                                    seconds=1.267 real 1.27
3  CHECKED …                                    seconds=1.281 real 10.47
   SKIPPED …                                    seconds=1.283 real 1.29
4  CHECKED …                                    seconds=1.281 real 10.39
   SKIPPED …                                    seconds=1.254 real 1.26
```

**Whole-run wall 10.42 s → 1.26 s (−88 %); the walk itself is untouched**, which is the point —
the switch removes a fixed tail, it does not make the generator faster.

**Gates.** All four dumps byte-identical in *both* modes (`cmp`); `cargo test --lib` 2,872 passed
(up 3 — see below); `cargo test --examples` 33 targets green; `cargo clippy --all-targets
--all-features -- -D warnings` clean; `cargo doc --no-deps` 12 unresolved links (baseline).

**Three tests, and the decisive one proves a negative.** Showing that a skip *skipped* needs a
fixture where reading the FASTA would fail: `pair_whose_fai_lies_about_the_wrapping` gives a FASTA
whose sibling `.fai` describes a different line wrapping. `TrustIndexWithoutChecking` returns `Ok`
on it — which it could only do without reading the bases — while `VerifyAgainstIndex` on the *same
fixture* raises `FastaFaiMismatch`. Without that paired control, a skip that silently did nothing
and a skip that correctly skipped would look identical.

## Author response convention

Address each finding by its identifier (H1, L1, …) with one of: `applied in <commit>` /
`experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` /
`won't fix because …`. The "no gain" path is expected — §5's Note section already carries eleven
of them, and two (N1, N2) were implemented and refuted inside this review.
