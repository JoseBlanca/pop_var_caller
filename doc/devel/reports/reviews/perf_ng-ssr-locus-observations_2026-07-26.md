# Performance Review: ng-ssr-locus-observations
**Date:** 2026-07-26
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** the ng STR **locus observation generation** — `SsrGenerator::next_locus` and its four steps (build the `SsrLocus`, fetch + depth-cap the reads, `classify_read` per read, tally), plus the default delimiter `SsrUnitRobustAligner`
**Verdict:** Apply the listed wins
**Hot-path evidence:** four sampling profiles + two three-point wall-clock sweeps + two counter experiments + two `cargo asm` listings, all collected for this review (full log: `tmp/perf_review_2026-07-26_ng-ssr-locus-observations/measurements.md`)

---

## 1. Scope and constraints

**What was reviewed.** Everything `next_locus` does between "here is a tract" and "here are the
observations": the tract ± flank reference fetch, the per-locus read query and depth cap, the
per-read classify pipeline and its pair-HMM tract delimiter, and the tally. This is the second
review of this code today — the first
([perf_ng-ssr-cohort-stutter_2026-07-26.md](perf_ng-ssr-cohort-stutter_2026-07-26.md)) reviewed the
same path from the cohort-walk end, landed H2 (the CRAM container cache, 129×) and left the
observation core mostly unexamined at 6.3% of a decode-bound profile. **H2 changed the shape of the
profile, so this review re-measures rather than re-reads.**

**Reviewed against:** branch `ng-aligner-bakeoff`. Measurements were taken at `5a9e71c`; the tree
moved to `e17be9b` mid-review (a parallel session committed
[examples/ng_ssr_cohort_stutter.rs](../../../../examples/ng_ssr_cohort_stutter.rs) input grouping
and two bench scripts). `git diff --stat 5a9e71c..HEAD -- src/` is **empty**, so every finding
below stands unchanged at HEAD.

**Throughput / latency targets, input sizes, hardware.** Two workloads with opposite shapes, and
the ranking depends on which one you mean:

| | reads / covered locus | measured | target |
|---|---|---|---|
| **shallow cohort** — 51 tomato samples × BED-restricted SL4.0ch01, CRAM | ~7 | 82.68 s cold / 73.78 s warm | cohort in minutes |
| **deep single sample** — HG002 300× BAM, chr16–20 Tier BED | ~412 | 51.71 s | the same code Stage-1 `ssr-pileup` runs per sample |

Hardware: macOS / Apple Silicon host (native release, `apple-m1`) and the Debian 12 aarch64 dev
container. **The PMU is not virtualised in either**, so `perf stat -e cache-misses`, `perf c2c` and
`perf sched` are unavailable; user-space sampling works on both and was not blocked.

**Hot-path evidence available.** Collected today, all quoted verbatim in the audit directory:

- three `sample` profiles — tomato 51-sample release (38,054 samples), the same workload under
  `[profile.profiling]` (37,957), HG002 300× 5-contig release (29,884 main-thread) — plus a fourth
  (HG002 chr20) that turned out startup-bound and is reported as such;
- a three-point `flank_bp` sweep on **both** workloads (a new experiment, temporary env knob,
  reverted);
- a delimiter call/retry counter experiment on both workloads (temporary, reverted);
- two independent `cargo asm --simplify "classify::delimit"` listings (1,781 and 1,681 lines).

**In-scope files.**

- [src/ng/locus_generation/ssr.rs](../../../../src/ng/locus_generation/ssr.rs) — the generator, `classify`, `read_region`, `tally`
- [src/ng/locus_generation/mod.rs](../../../../src/ng/locus_generation/mod.rs) — the shared output types it fills
- [src/ng/alignment/ssr_unit_robust.rs](../../../../src/ng/alignment/ssr_unit_robust.rs) — the default delimiter (algorithm 4u)
- [src/ng/alignment/emission.rs](../../../../src/ng/alignment/emission.rs), [src/ng/alignment/mod.rs](../../../../src/ng/alignment/mod.rs), [src/ng/alignment/stutter.rs](../../../../src/ng/alignment/stutter.rs) — the delimiter's inputs
- [src/ng/read/filtering.rs](../../../../src/ng/read/filtering.rs) — the per-read filter chain the fetch drives (filter #8 is the reference-touching one)
- [src/ng/read/input/region_query.rs](../../../../src/ng/read/input/region_query.rs), [src/ng/read/input/mod.rs](../../../../src/ng/read/input/mod.rs), [src/ng/read/input/open_bam.rs](../../../../src/ng/read/input/open_bam.rs) — **in scope only where the profile lands inside `next_locus`'s own call tree**
- [src/ng/ref_seq.rs](../../../../src/ng/ref_seq.rs), [src/ng/raw_chrom_reader.rs](../../../../src/ng/raw_chrom_reader.rs), [src/ng/reference_info.rs](../../../../src/ng/reference_info.rs)
- [examples/ng_ssr_cohort_stutter.rs](../../../../examples/ng_ssr_cohort_stutter.rs) — the driver

**Deliberately out of scope.** Region typing (`src/ng/region_typing/` — reviewed this morning;
`scan_set`'s whole-contig scan is an owner decision), production `src/ssr/` (frozen), patching
vendored noodles, and the BGZF/CRAM decode internals below the sites named above.

**Categories dispatched.** `methodology` (always), `hot_loops` (the DP is the top symbol in both
regimes), `allocations` (per-read `Box`/`Vec`/`HashMap` churn at 412 reads/locus), `data_layout`
(the DP ring and backpointer matrix), `io_and_syscalls` (the profile puts 11.7% of the run in
`open`/`read`/`lseek`), `concurrency` (single-threaded on an 8-core box; per-sample fan-out is a
live plan).

## 2. Verdict

**Apply the listed wins.** Three facts reorder everything the previous review left standing.

**1. The delimiter is now the code.** `classify::delimit` is the top self-time symbol in every
profile taken today: **36.0%** of the shallow cohort run (13,705 / 38,054), **37.9%** of the same
workload under `[profile.profiling]` where it appears as its own un-inlined symbol
(`SsrUnitRobustAligner::align`, 14,403 / 37,957), and **67.0%** of the deep single-sample run
(20,018 / 29,884). This morning's review measured it at 6.3% of a profile that was 71% zlib
inflate; H2 removed the inflate, and what was 6.3% of the per-locus stage is now a third to two
thirds of the whole process. **Its five DP findings (L3–L7) were filed with codegen evidence and no
hot-path evidence. They now have hot-path evidence, all five are still present in the code, and
none was made stale by H2.**

**2. The single largest lever is a config knob, not a code defect — and it is already a declared
open question.** `flank_bp` defaults to `DEFAULT_BUNDLE_THRESHOLD = 30`, so the DP's reference frame
is `tract + 60`. Measured reference-tract lengths are a median of **8 bp (tomato) and 14 bp
(HG002)**, so **~88% of the DP's columns are flank**. Narrowing the flank cut wall clock by
**31% (tomato) and 49% (HG002)** in a three-point sweep — while the reads behind *complete*
observations moved by **0.01–0.03%**. What shrinks is the partial-observation set. That is exactly
the trade `arch/locus_generation_ssr.md` §4 already lists as open ("Should `flank_bp` equal the
bundle threshold?" and "Do partial observations pay *for genotyping*?"). **This is a science
decision with a measured price tag, not a merge — the review's job is to put the number on it, and
the number is large.**

**3. H2's successor costs are the whole remaining fetch bill, and both are removable in place.** The
fetch branch is still 58.2% of the shallow-cohort run, but the decode is no longer in it. What is
in it: the container-cache **hit** path linear-scans every record of the held container and
re-walks each record's CIGAR (**≈15%**), and the per-query reference-reader factory re-opens the
`.fai` *and* the FASTA for every (locus, sample) — **14.0%, of which 1,306 samples are inside
`open(2)`**. Neither needs H1's grouping, its new generator entry point, or its coordinate-order
invariant. **H1's deferral is therefore strengthened, not reversed.**

Two candidate levers were measured and **closed**: the widen-and-retry second DP, which two
categories flagged as a possible 2× on the hottest function, fires on **0.10% / 0.17%** of reads;
and the "reads admitted on overlap pay a wasted DP" hypothesis is refuted — 99.7% / 98.5% of
fetched reads yield an observation. The prior review's H4 (margin-fetch cache) is confirmed dead at
**0.16%**, and the depth cap never binds (`reads_capped = 0` at 412 reads/locus against a cap of
1,000).

Ordering advice: **do H2/H3 (contained, bit-identical, ~29% of the shallow run between them) and H4
(one field, one file) before touching the DP's structure**, because H4 changes the register
allocation every other DP finding will be measured against. And build the delimiter bench first —
without it, nothing in H4–H6 is measurable (see M1).

## 3. Measurement plan

**1. The delimiter bench — the first deliverable, and it is unblocked today.** The prior review's
M1 ("`src/ng/` has no committed benchmark at all", re-verified: `grep -rn "pop_var_caller::ng"
benches/` returns nothing) was blocked on lifting `read/input/test_fixtures.rs` out of
`#[cfg(test)] pub(crate)`. A **delimiter-only** bench needs none of that — `SsrUnitRobustAligner`,
`delimit`, `align`, `UnitRobustScratch`, `PerQualityEmission`, `ReadBases`, `RepeatGeometry`,
`StutterModel::hipstr_shipped` and `Motif::new` are all `pub`, and the DP takes two byte slices.

```
benches/ng_ssr_delimiter_perf.rs   (harness = false; criterion 0.8.2 is already a dev-dep)
  groups: frames{period 1..=6 × tract 10/20/40} × depth{1, 7, 30, 100, 412}
  cargo bench --bench ng_ssr_delimiter_perf -- --save-baseline pre-<finding>
```

Three design rules, each tied to a real trap in this code: **`black_box` both ways** (fat LTO will
otherwise hoist exactly the per-locus setup H5 proposes to hoist, and report the win as already
present); **an in-bench verification assertion**, because `delimit` returns `None` immediately when
`reference_len == 0` and a mis-sliced fixture yields a plausible, meaninglessly fast number; and
**a multi-frame group**, because `UnitRobustScratch` is grow-and-keep and never re-zeroed, so a
single-input bench measures a permanently-warm L1 that real loci never see. The fixture builder
already exists at `examples/ng_ssr_synthetic_bakeoff.rs:41-335` in a binary `benches/` cannot
import — copy the four scenario constructors, or lift them into a shared module, and say which.

**Why it is load-bearing:** the DP findings are gated at ≥10% *on the delimiter*, which Amdahl
translates to 3.6% (tomato) / 6.7% (HG002) end to end — straddling the line where this project has
already recorded a **+6% thermal-drift false positive** (prior review, H4).

**2. The identity gate — free, it already exists.** Every change in section 5 marked
"bit-identical" must leave the driver's TSV byte-identical:

```
tmp/perf_review_2026-07-26_ng-ssr-locus-observations/cohort_rows.tsv    950,637 rows
    md5 = 9409dc94253d397155222b231be7afa3     (tomato, 51 samples, ch01 BED, CRAM)
tmp/perf_review_2026-07-26_ng-ssr-locus-observations/hg002_5c_rows.tsv
    md5 = 699d50d5379107127db0773b4088cd26     (HG002 300×, chr16–20 Tier BED)
```

plus the counted invariants for M7's asserts: tomato 1-sample `covered 5534 / obs 5207+13620 /
reads 19314+20907+116` (identical to the prior review's W2 row, two independent runs); HG002 chr20
`720 / 2976+27114 / 161862+130588+4098`; HG002 5-contig `4662 / 18886+174356 /
1033458+860551+29166`. Note `delimit_parity.rs` pins **algorithm 3, not 4u**, so it is not the
oracle for any DP change — the TSV md5 is.

**3. The cross-commit protocol.** Criterion `--save-baseline`/`--baseline` for the delimiter group;
**interleaved** end-to-end A/B (two binaries, alternating, ≥3 runs per side — never back-to-back
batches); and a **revert experiment** for any accepted win (revert, re-measure, confirm sign and
magnitude both flip). Write it into the bench file's module doc, not into a review document — that
is precisely what was lost when the prior review's probe was deleted.

**4. The `flank_bp` concordance study (gates H1).** The perf half is done (below). The missing half
is scientific: at `flank_bp ∈ {30, 20, 16, 12}`, compare against HipSTR / the GIAB HG002 truth set
what is actually lost — the partial observations, the 1.9% of (locus, sample) cells that stop being
covered, and whether the anchor rule still holds (`AnchorRule::window = 5`, `min_matches = 3`,
`min_support = 5`, and `JunctionGuard::flank_width = min(8, flank_len/2)` changes once
`flank_len < 16`). **Threshold: adopt the narrowest flank at which complete-observation
concordance is unchanged and the partial-fed likelihood does not degrade.** `flank_bp` is already a
config field, so the experiment is a flag flip plus a comparison.

**5. Re-profile the deep regime under `[profile.profiling]`.** Its second-largest main-thread symbol
is `run_cohort` at **5,456 / 29,884 = 18.3%** — a fat-LTO-fused bucket spanning in-scope work
(classify's non-DP part, the tally, the `outcomes` `Vec`) and out-of-scope work (TSV formatting).
The tomato pair proves the profiling build dissolves that bucket. Two minutes of machine time, no
code. **Run it before acting on any tally or allocation finding** (L7–L9): those were filed at Low
wall-clock confidence from a *shallow* profile, and the deep regime feeds them two orders of
magnitude more items.

**6. Stage timers in the driver (M4, still open).** Three `Instant`s and a moved `join` of the
verify handle. This review had to attribute the walk/fetch/classify split from profile trees rather
than read it, and one 40 s profile was discovered to be startup-bound only after it was taken.

## 4. Build / toolchain configuration

The build lever is still pulled (`lto = "fat"`, `codegen-units = 1`, `panic = "abort"`,
`opt-level = 3`, pinned toolchain, an existing `[profile.profiling]`). Two items are **already
filed** by the prior review and are not re-argued: **M9** (no `target-cpu` for
`aarch64-unknown-linux-gnu`; `.cargo/config.toml` covers only `x86_64-linux` and `aarch64-macos`)
and **M10** (container cache mounts). M9's stakes went up: it was correctly discounted on a
zlib-bound path, but the code that now matters is a scalar `f64` DP full of `%`, `min`/`max` and
`fcsel` chains, so host and container are two configurations for exactly the function under
review.

**M11 splits, and its own revisit condition is met.** It deferred PGO *and* an allocator A/B until
the profile was no longer inflate-dominated. It no longer is (`inflate_fast_help` is 115 of 38,054
on tomato, down from 934 of 2,162 pre-H2).

- **PGO: indicated.** The run is dominated by one function that is a nest of data-dependent
  branches (`best_of`'s comparison chains, `gap_open`/`guarded`'s `ccmp`/`cset` predicates, the
  terminal-route tests), whose outcomes are skewed by position within the frame and therefore
  highly predictable across runs, on a workload that is the definition of stable and repeatable.
  **Sequence it after the delimiter bench and after the cheap DP wins**, so PGO is not credited with
  a code win; threshold ≥5% on the tomato cohort *and* a non-regression on HG002, TSV byte-identical
  (if PGO changes results, something in the DP is order-sensitive and that is a much bigger
  finding). Cost: a three-stage build behind a documented script, never a default profile.
- **The allocator A/B: do not run it.** Allocator self-time is 1.9% (tomato) and **1.2% in the deep
  regime that allocates ~60× more per locus**, the path is single-threaded, and the function that
  replaced inflate at the top does zero allocations per call. The crate already carries an
  `alloc-mimalloc` feature if anyone wants the null result on the record.

## 5. Code-level findings

Full per-finding detail — diffs, verbatim `cargo asm`, the noodles source that confirms each
mechanism, and the layout probe output — is in the six per-category files under
`tmp/perf_review_2026-07-26_ng-ssr-locus-observations/`. Below is the synthesis.

### Hot-path

#### H1: src/ng/locus_generation/ssr.rs:115 (with src/ng/region_typing/segment_criteria.rs:252) — `flank_bp = 30` makes ~88% of the DP flank, and narrowing it is worth 31–49% of wall clock

- **Confidence:** High on the cost; the *decision* is not a performance decision.
- **Hot-path evidence:** my own three-point sweep, both workloads, `/usr/bin/time -p`, temporary
  `PVC_NG_FLANK_BP` override (reverted). HG002 300× chr16–20, 1 sample, covered loci identical
  (4,662) at all three widths:

  | `flank_bp` | real | reads → completes | reads → partials |
  |---|---|---|---|
  | 30 (default) | **52.74 s** | 1,033,458 | 860,551 |
  | 20 | **36.30 s** (−31%) | 1,033,298 | 732,988 |
  | 12 | **27.13 s** (−49%) | 1,033,180 | 591,652 |

  Tomato, 51 samples, aggregated: **73.78 → 59.71 → 50.58 s** (−19%, −31%); reads behind completes
  955,714 → 955,656 → 955,586 (**−0.01%**); reads behind partials −41%; covered (locus, sample)
  cells −1.9%. Supporting measurement: reference-tract lengths are **p50 8 bp (tomato), 14 bp
  (HG002)** against a frame of `tract + 60`.
- **Pattern matched:** the loop bound itself is the invariant — cell count is linear in frame width,
  and the query span (hence the number of reads that pay a DP at all) shrinks with it too.
- **Mechanism:** two effects at once, both real. The DP matrix is `(read slice + 1) × (tract + 2 ×
  flank + 1)`, so at the median tomato locus 60 of 67 columns are flank; and the read query spans
  the same window, so a wide flank admits marginally-overlapping reads that can only ever produce
  partials. What the delimiter's *rules* consult is far narrower: the anchor test looks at 5
  columns each side, the junction guard at `min(8, flank_len/2)`.
- **Measurement plan:** measurement item 4 — the concordance study. The perf side needs no further
  work.
- **Complexity cost:** **zero code** (it is a config field). The cost is scientific: fewer partial
  observations, 1.9% fewer covered cells, and a changed `flank_width` below `flank_len = 16`. It
  moves the measurement, so it must not be taken as a perf change.
- **Note:** this also caps every DP finding below. If `flank_bp` drops to 16, H4–H6 are optimising a
  matrix roughly 2.9× smaller.

#### H2: src/ng/read/input/region_query.rs:458-466 — the container-cache *hit* path rescans every record of the held container and re-walks each record's CIGAR, once per (locus, sample)

- **Confidence:** High
- **Hot-path evidence:** `[profile.profiling]` build, 37,957 samples, tree attribution verbatim:
  ```
  + ! : | 5945   <CramRegionSource as RecordSource>::read_next
  + ! : | + 3071   <Vec<T> as SpecFromIterNested<T,I>>::from_iter
  + ! : | + ! 1904   Cloned<I>::next
  + ! : | + ! : 1904   RecordBuf::alignment_end
  ```
  **3,559 (9.4%) in the overlap test + 2,120 (5.6%) in the filter-and-clone ≈ 15% of the run.** In
  the release profile the same site appears as `Cloned<I>::next` **7,458 = 19.6%**, second in the
  self-time ranking, with `alignment_end` inlined into it — treat 15% as measured and 19.6% as the
  upper end.
- **Pattern matched:** decode once into a reusable buffer — one level up from H2 of the prior
  review. The decode is cached; the *interpretation* of the cached buffer is not.
- **Mechanism:** on a cache hit, `refill` runs `held.records.iter().filter(…).cloned().collect()`
  over **every** record of the container. `overlaps` calls `RecordBuf::alignment_end`, which is a
  fresh CIGAR walk (`noodles-sam-0.85.0/src/alignment/record_buf/cigar.rs:36-41`). Measured from the
  actual index: the bench CRAM's `.crai` has **89 entries and 89 distinct offsets** (one slice per
  container), htslib's default is ~10k records per slice, and ch01's ~282k queries each walk that
  container to select ~7 reads. The records are coordinate-ordered and their footprints never change
  after decode, so both can be resolved **once per decode**.
- **Measurement plan:** step 1 (footprint side-table, no new invariant) first, on the 51-sample
  workload against 82.68 s cold / 73.78 s warm, gated on the `cohort_rows.tsv` md5. **Threshold:
  ≥10% off the wall**; below 5%, the scan was not the cost and the finding is refuted. Re-profile to
  confirm `Cloned<I>::next` and `alignment_end` have left the top of the ranking.
- **Complexity cost:** one `Vec<(u32, u32, u32)>` on `DecodedContainer` (~120 KB per pooled reader,
  on top of the +34 MB the container cache already costs) and a `debug_assert` that it stays the
  same length as `records`. Step 2 (a `partition_point` cursor) adds an invariant and should only
  follow if step 1 leaves the scan visible.
- **Silent-wrongness risk:** a mis-bounded cursor returns **fewer** reads — lower depth, no error.
  The trap is specific: records are ordered by *start*, but a record starting before the region can
  still overlap it (a long deletion), so the cursor needs a container-wide `max_span` look-back. The
  CRAM path's only structural oracle is
  `t8_a_cram_yields_the_same_ordered_reads_as_the_same_bam` (`open_bam.rs:1717`) — the BAM-only T5
  linear-scan oracle does not cover it — and its fixture needs a long-`D` record added.

#### H3: src/ng/locus_generation/ssr.rs:1615 (via src/ng/read/filtering.rs:291, src/ng/raw_chrom_reader.rs:111-141) — the `make_reference` factory re-opens the `.fai` *and* the FASTA for every (locus, sample)

- **Confidence:** High
- **Hot-path evidence:** release profile, 38,054 samples. Self-time `__open` **2,331**, `read`
  **1,856**, `__lseek` **257** = **11.7% of the run in syscalls**. Tree:
  ```
  + 21150  OrderVerified<I>::next
    + 5318   WindowedRefSeq::fetch_raw_into        (ref_seq.rs:648)   14.0% of total
    |  + 1674  RawChromReader::for_contig          (raw_chrom_reader.rs:235)
    |  |  + 1312  noodles_fasta::fai::fs::read     → 1,306 in open()
  ```
  Contrast: the generator's *own* margin fetch (`RefSeq::fetch_into`) is **61 samples = 0.16%**.
- **Pattern matched:** per-call syscall — here `File::open` per logical operation, and the buffering
  that exists (`fai::fs::read` wraps its own `BufReader`) is thrown away with the reader at the end
  of every query.
- **Mechanism:** the driver's factory is `move || WindowedRefSeq::new(fasta.clone(),
  contigs.clone())`, called once per file per query. `WindowedRefSeq::new` is lazy, so the cost
  lands on filter #8's first `fetch_raw_into`: `RawChromReader::for_contig` → `open_contig` does
  `fai::fs::read` (a `File::open` plus a parse of the whole index into an owned `Vec<Record>`), a
  linear `find`, then `File::open` on the FASTA. **Two `open(2)`s and a full `.fai` parse per
  (locus, sample) that has any read** — ≈282k queries, ≈564k opens on the measured chromosome. The
  window is discarded when the query's stream drops, so the next locus starts cold again.
  **This confirms the prior review's L11 and shows it was under-scoped**: L11 costed the `.fai`
  parse (which scales with contig count) and missed the two opens, which is where 1,306 of the 1,312
  samples actually land.
- **Measurement plan:** confirm the mechanism first with a counter in `for_contig` (it must print
  ≈ covered loci × samples, and drop to ≈ contigs × samples after the fix). Then **≥10% off the
  51-sample wall**, TSV md5 unchanged, `__open` gone from the self-time top-10.
- **Complexity cost:** the recommended shape is the prior review's own suggestion and closes the
  "Arc gap" the `SsrGenerator` type doc names: three blanket impls (`RefSeq`/`RawRefSeq`/
  `ContigTable` for `Arc<T>`) plus one driver line handing the factory `Arc<ResidentRefSeq>`. Keeps
  `FnMut() -> R` intact. Its real cost is **residency** — a whole contig (SL4.0ch01 ≈ 90.8 MB,
  GRCh38 chr1 ≈ 250 MB) instead of a 64 KiB-class window. On CRAM that contig is *already* resident
  in noodles' own `fasta::Repository`; on BAM it is a genuine new cost.
- **Silent-wrongness risk:** a different accessor returning different bytes changes filter #8's
  mismatch count, silently changing which reads are kept. The oracle already exists and is
  exhaustive — `windowed_raw_returns_verbatim_bytes_matching_resident` (`ref_seq.rs:1110`) and
  `windowed_canonical_matches_resident_across_all_windows` (`ref_seq.rs:1192`) are precisely the
  "this swap cannot change an answer" proof.

#### H4: src/ng/alignment/emission.rs:171 and :215 — `PER_QUALITY_LN` is a `LazyLock`, so every table lookup inside the DP is an acquire-load plus a live call site, and the seven cold call sites are what spill the loop's constants to the stack

- **Confidence:** High on the mechanism; Medium on the size of the win. **The cheapest DP change in
  this review — do it first.**
- **Hot-path evidence:** `cargo asm --example ng_ssr_cohort_stutter --simplify "classify::delimit"`
  (1,781 lines, host native release): **7 `ldapr` and 7 `bl <std::sys::sync::once::queue::Once>::call`
  inside `delimit`**, and 15 reloads of loop-invariant values from the stack inside the cell body —
  including `ldur w10, [x29, #-216]`, which is **the read base**, constant for the whole row.
- **Pattern matched:** hoist the invariant — here the invariant is the initialisation check.
- **Mechanism:** `PerQualityEmission` is a unit struct; `scores_for` dereferences the `LazyLock` on
  every call. LLVM must treat the loop's live values as call-clobbered at seven points inside
  `delimit`, which is why the guard constants and the read base are re-read from the stack per cell
  rather than held in registers. Giving `PerQualityEmission` one field
  (`&'static [BaseScores; 256]`, resolved in `new()`) removes the acquire-loads, the call sites,
  and — the real prize — the register pressure.
- **Measurement plan:** `cargo asm` must show 0 `ldapr`, 0 `Once::call`, and a drop in the cell
  body's 15 stack reloads; then the delimiter bench, **threshold ≥5%**; then the 51-sample wall with
  the TSV md5. All 84 references to `PerQualityEmission` go through `::new()`, so the blast radius
  is one file.
- **Complexity cost:** low. The type stops being a ZST (8 bytes, carried by value in every aligner,
  all `Copy`), and the module doc's "stateless: constructing one of these is free" must be updated
  rather than left to become a lie. **Bit-identical:** yes — same table, same values.

#### H5: src/ng/alignment/ssr_unit_robust.rs:368-396 (used at :482, :512, :517, :536, :588) — every cell re-derives quantities that are constant along a whole stretch of the column axis

- **Confidence:** High. This is the prior review's **L4 and L5 in a stronger form that deletes
  both.**
- **Hot-path evidence:** the delimiter's 36.0% / 37.9% / 67.0%, plus the cell body verbatim from the
  `cargo asm` listing — a `ccmp`/`cset`/`fcsel` chain in which four locus constants are reloaded
  from the stack (`ldur x11, [x29, #-232]`, `ldp x12, x11, [x29, #-248]`, `ldur x12, [x29, #-256]`)
  — and the L5 range test as an unrolled `+1/+2/+3/+4/+5/+6` compare chain. 18 of the cell body's 23
  out-of-line branches are `panic_bounds_check` targets.
- **Pattern matched:** hoist invariants; avoid branches in the innermost loop.
- **Mechanism:** `guarded`, `gap_open`, `column_in_tract` and the `all(column_in_tract)` test are
  functions of `column` and the locus geometry alone. Read together, the column axis is **piecewise
  constant in at most five stretches** (left-flank body, left junction window, tract interior, right
  junction window, right-flank body). L4 proposed *tabulating* them per column, trading ~20
  instructions for a load; **splitting the column loop at those boundaries is strictly better**:
  `open` becomes a register, `column_in_tract` becomes a compile-time fact of which sub-loop you are
  in, and L5's range test becomes a loop bound. Two consequences the prior review did not spot:
  **two of the five `Match` candidates are provably `−inf` outside the tract** and can be dropped
  bit-identically (with every candidate at `−inf`, `best_of` returns `candidates[0]` either way),
  and in the flank sub-loops the slip row is never touched, removing one of L6's three ring
  indirections for ~88% of columns.
- **Measurement plan:** the delimiter bench, **threshold ≥15% on the period-1/2, tract-10–20 cases**
  (the tomato shape). Gate: TSV md5.
- **Complexity cost:** real. The cell body must be written once as
  `#[inline(always)] fn cell<const IN_TRACT: bool>(…)` and instantiated twice, so the scoring
  arithmetic stays single-source. The five stretch boundaries are the new invariant and every one
  can degenerate (guards clamp to 0, the tract can be shorter than one unit, either flank can be
  empty). A wrong boundary is a silently wrong measurement, not a panic — the test must be a
  **boundary sweep** (flank 0/1/8/30 × tract 0/period/2·period/60) diffed against the current
  implementation. **Bit-identical:** yes, provably — no arithmetic changes, candidate and summation
  order preserved.

#### H6: src/ng/alignment/ssr_unit_robust.rs:557-583 (with :400) — the whole-unit slip emission is resolved per tract *column* although it depends only on (read row, motif phase) — L3 confirmed

- **Confidence:** High
- **Hot-path evidence:** the same 36.0/37.9/67.0% attribution, plus the unrolled `k` chain verbatim:
  six `ldapr` of `PER_QUALITY_LN` and six `udiv`+`msub` pairs, one per `k`. Function census: 7
  `ldapr`, 11 `udiv` — six of each are these.
- **Mechanism:** `unit_emit` sums `period` emission scores; the **only** dependence on `column` is
  the motif phase, which takes `period ≤ 6` distinct values. It is currently recomputed at every
  tract column — `tract_len / period` times more work than needed, each instance carrying a runtime
  `%` (the `udiv`+`msub`, because `period` is a runtime value), an acquire-load and two bounds
  checks.
- **Measurement plan:** the delimiter bench, **threshold ≥10% on the long-tract cases (tract 40,
  period 3–6)**; expect little at tract 10 / period 1, which is why this ranks below H5 for tomato
  and above it for HG002. Re-`cargo asm` must show `udiv` fall from 11 to ~5.
- **Complexity cost:** one `[f64; 6]` filled at the top of each read row. **Bit-identical: yes,
  provided the `k`-ascending summation order is preserved** — `f64` addition is not associative.

#### H7: examples/ng_ssr_cohort_stutter.rs:274-279 — 51 independent samples are still serialized through one `SsrGenerator`, and the fan-out unit the prior review recommended is now the wrong one

- **Confidence:** High. **Confirms the prior review's H3 and corrects its recommendation.**
- **Hot-path evidence:** 51 samples = `real 1:22.68` at **100% cpu** on an 8-core box; one sample =
  4.599 s; the prior review measured 3 samples at exactly 3× one sample. No lock, futex or park
  frame appears anywhere in any profile taken today — the serializer is the `&mut self` borrow, not
  contention.
- **Mechanism, and what is new:** `next_locus` is `&mut self` over *pure per-worker state*
  (`align_scratch`, `margin_buffer`, `qual_buffer`, `counts`, `current_region`, `produced`). H3
  recommended **N single-sample processes** as the zero-code shape. Post-H2 that is the wrong unit,
  because the per-sample work collapsed while the per-process fixed cost did not — the walk is
  ~3.3 s against ~1.55 s of per-sample work, plus a ~2.4 s reference digest per process. Arithmetic
  over the measured numbers (labelled as arithmetic):

  | shape | wall | CPU |
  |---|---|---|
  | today: 1 process, 51 samples | ≈82.7 s (matches the measurement) | ~83 s + 1 digest |
  | one process per sample, 8 at a time | ≈34 s | ≈247 s + 51 digests |
  | **8 workers, ~6–7 samples each** | **≈13 s** | ~83 s + 8 digests |

  Per-worker cost: ~140 KB of generator scratch, plus the **measured +34 MB** decoded container per
  pooled reader. The `fasta::Repository` is per sample either way — and a batched fan-out is a peak
  **RSS win**, because today all 51 CRAMs are opened up front (51 repositories, up to 51 live
  containers) where 8 workers would hold 8.
- **Measurement plan:** shape A (a shell fan-out over *batches*, zero library code) against the
  82.68 s baseline; **threshold: wall ≤20 s at 8 workers, peak RSS no worse.** The oracle is the
  **union of per-worker rows sorted** by `(sample, contig, start, end, coverage, observed)` — row
  order changes, row content must not. If a threaded shape is not ≥1.3× better than batched
  processes, stop at the shell loop.
- **Complexity cost:** shape A is none in the library. A threaded shape needs L12 (`Send`), L13
  (per-worker output), the segment set materialised, and per-worker counters summed.

### Likely

- **L1: [ssr_unit_robust.rs:239, :606-607](../../../../src/ng/alignment/ssr_unit_robust.rs#L606-L607) — the backpointer cell is written as five separate byte stores per cell** (prior L7, confirmed and
  sharpened). Confidence Medium. The asm shows `add x12, x0, x0, lsl #2` (the ×5 stride) then five
  `strb` — the 5-byte record is never naturally aligned, so LLVM cannot merge them. Measured
  traffic, from the layout probe at real dimensions: **26,645 B per read, 10.98 MB per locus,
  7.90 GB over the 12.3 s HG002 chr20 run** — note this corrects the prior review's 99 KB/read by
  ~4×. Five states × 3 bits fit a `u16`: one `strh`, 2 bytes instead of 5, and an aligned stride.
  Threshold ≥8% on the delimiter bench; the risk is that the loop is ALU-bound, which is exactly why
  it is worth measuring before L2. Complexity: a `pack`/`unpack_one` pair plus
  `const _: () = assert!(STATES <= 5);` beside it — a sixth state silently corrupts a 3-bit field.
  (Decline the `u8` that would also fit: a 1-bit-short field truncates silently and a
  `debug_assert` compiles out of the release build this repo runs.) Bit-identical.
- **L2: [ssr_unit_robust.rs:238](../../../../src/ng/alignment/ssr_unit_robust.rs#L238) — the `Vec<Vec<[f64; STATES]>>` ring costs six addressing instructions per access, three of them loads**
  (prior L6, confirmed by two categories independently; data_layout filed it Hot-path, hot_loops
  Likely — the site is hot, the *fix's* gain is what is uncertain, so it lands here). The asm proves
  the row `ptr`/`len` are **reloaded per access** because the cell store may alias the outer buffer;
  18 of the cell body's 23 out-of-line branches are `panic_bounds_check`. Two structural corrections
  to the prior review: the three live rows are **ring-adjacent** (`cur−1, cur, cur+1`), and the
  `period == 1` "aliasing trap" is **benign** (both are reads). Prefer variant (a) — bind the rows to
  local slices outside the column loop with `split_at_mut`/`get_disjoint_mut` — over flattening;
  flattening turns a `stride` mix-up into a silent wrong answer instead of a panic. **Do this last:**
  H5 already removes one of the three indirections for ~88% of columns, so L2's measured value will
  be smaller afterwards. Threshold ≥10%. Bit-identical.
- **L3: [ssr_unit_robust.rs:349, :406-439](../../../../src/ng/alignment/ssr_unit_robust.rs#L406-L439) — every read re-derives the locus's frame: four `bl _log` calls and the entire full-width row 0**
  (new). Confidence Medium. `SlipCosts::from_model` takes logarithms of *run* constants
  (`StutterModel::hipstr_shipped()` is built once in `SsrGenerator::new`) per `delimit` call; row 0
  reads **nothing** from the read, so its `reference_len` cells of `best_of` chains, range tests and
  5-byte backpointer writes are recomputed identically for all 7 (tomato) or **412 (HG002)** reads at
  the locus, into the same addresses. Little value alone (~1% of the fill); real value as the
  *vehicle* for H5 and H6 at depth. Complexity: the cheapest shape memoises inside
  `UnitRobustScratch`, and **the memo key is the risk** — a stale hit returns another locus's guard
  widths, a silently wrong measurement, so the key must include flank lengths, `reference_len`, the
  motif and the stutter-model identity, and the test must be a cache-*miss* test.
- **L4: [raw_chrom_reader.rs:184](../../../../src/ng/raw_chrom_reader.rs#L184) — every `read_raw_bases` call zero-initialises a 64 KiB stack array and then reads 64 KiB to keep ~150 bases**
  (prior L9, confirmed with numbers plus a half L9 never named). `read` 1,351 (3.5%) + `__bzero` 613
  (1.6%) = **5.1% of the run**. `let mut read_buf = [0u8; FILE_READ_CHUNK]` is a fresh 65,536-byte
  stack array per call which Rust must zero before `&mut` reaches `read`; ~99.8% of the bytes read
  are then discarded by the newline-stripping loop. Fix: a grow-and-keep `Vec<u8>` field sized from
  the `.fai` geometry already in hand. Two lines plus one arithmetic helper; threshold ≥4%.
  **H3 removes this call site from the STR path entirely** — so take L4 only if H3 is rejected, or
  for the region-typing walk, which is `RawChromReader`'s other consumer.
- **L5: [region_query.rs:495-501](../../../../src/ng/read/input/region_query.rs#L495-L501) — the CRAM slice-MD5 is 8.8% of the run and our-side caching is already at its floor.** Confidence
  High on the mechanism and arithmetic; the fix is an input-encoding experiment, not code. noodles
  MD5s the slice's whole reference span on every `records()` call
  (`noodles-cram-0.93.0/src/io/reader/container/slice.rs:352-363`), which post-H2 is once per
  container decode. Measured from the `.crai`: **89 entries / 89 distinct offsets (one slice per
  container), 10 containers tiling ch01's 90.79 Mb** — so the MD5 total is one whole-contig digest
  per sample (51 × 90.79 MB ≈ 4.63 GB ⇒ ~634 MB/s, which matches the 8.8%). Three consequences:
  H2's un-applied sub-fix (ii) ("decode only the slice `landmark()` names") **buys nothing on these
  files — do not spend the commit**; no larger cache or loop reordering can reduce 10 decodes per
  sample; and the check cannot be bypassed from outside noodles (`Slice::header`,
  `ReferenceSequence` and `Record`'s fields are all `pub(crate)`). Converting to BAM does **not**
  win either (post-H2: 1.35 s CRAM vs 1.38 s BAM). The one strict lever is re-encoding with
  `embed_ref=1`; test on one CRAM, threshold ≥7% off the single-sample wall, and record *why* in
  `bench.config.sh` or someone regenerates them with the default and the 8.8% silently returns.
- **L6: [reference_info.rs:435-469](../../../../src/ng/reference_info.rs#L435-L469) — the startup whole-FASTA digest sets a short deep-run's wall clock, and the main thread blocks joining it**
  (prior L17, now **measured** on the single-process path rather than predicted for a fan-out).
  HG002 chr20, `real 12.274s` at 166% cpu: the 40 s profile is 3,587 samples of which
  **`__ulock_wait` 3,160** (the main thread in `join`), `md5::compress` 2,098, `read_fasta` 1,434.
  `read_fasta_bases` feeds every byte of 3.1 GB through a per-byte state machine. Three cheap A/Bs
  in order: a skip flag (establishes the ceiling — if it is not ≥30%, stop), a buffer-size bump (run
  it to *close* the I/O question, not to win it), and a chunked inner loop (`memchr` + slice-wise
  uppercase, threshold ≥25% of the digest thread's own runtime). The skip flag is a judgement call
  for the owner and **must not** become the default for `ssr-pileup`, whose `.psp` contract depends
  on the reference digest.
- **L7: [ssr.rs:756](../../../../src/ng/locus_generation/ssr.rs#L756) and [:791](../../../../src/ng/locus_generation/ssr.rs#L791) — one `Box<[u8]>` per observed read, ~90% of which `entry()` immediately frees as a duplicate key.**
  Confidence High on the allocation count, **Low on wall clock** (see the Note on the allocator).
  Derived from the W-B tallies: 292,450 boxes minted, 30,090 survive as distinct buckets. The prior
  review filed this Speculative with the objection "stable Rust has no `raw_entry`"; the fix answers
  it — carry the tract as a `Range` into `read.seq` (which the tally already holds, zipped) and look
  the bucket up with a borrowed slice via a nested
  `HashMap<ReadCoverage, HashMap<Box<[u8]>, Support>>` (`Box<[u8]>: Borrow<[u8]>`). A counting-
  allocator model puts the tally at **412 allocations per deep locus today vs 42** for the proposed
  shape. Measure with DHAT (`--features dhat-heap`) for the count and the deep workload for the
  wall; expect the count, not the clock.
- **L8: [filtering.rs:819](../../../../src/ng/read/filtering.rs#L819) and [:397](../../../../src/ng/read/filtering.rs#L397) — four heap blocks per kept read from `MappedRead` decode, one of them (`qname`) never read on this path**
  (new). ~1,650 blocks per deep locus — four times the whole tally path, and the largest count lever
  in scope. Take the cheap half first (skip `qname`), **but** `merge.rs:158` *does* read `qname` on a
  tie, so a name-free decode must be confined to the single-file arm. Same confidence split as L7:
  count High, wall Low.
- **L9: [ssr.rs:1004](../../../../src/ng/locus_generation/ssr.rs#L1004) — the tally's per-locus `HashMap` re-grows its table on every deep locus** (prior L12, re-judged — **its premise was
  wrong**: `HashMap::new()` is lazy, so the 97%-empty shallow visits pay nothing; the real cost is
  five table-growth reallocations per *deep* locus). Fix: hoist onto the generator and `drain()` (not
  `into_iter()`). Ride the `FxHashMap` swap along on the same two lines — SipHash over ~30-byte keys
  is measured at only ~0.08%, so neither is worth a commit alone. Bit-identical: the map's iteration
  order changes but `observed_sequences` is sorted on a total order over unique keys.
- **L10: [open_bam.rs:591-601](../../../../src/ng/read/input/open_bam.rs#L591-L601) — neither pooled reader is wrapped in a `BufReader`** (prior L8, re-scoped). The **CRAM side is now
  negligible** post-H2 (~10 container decodes per sample) — said explicitly so it is not re-filed.
  The **BAM** side is the deep regime's only per-query syscall item: each BGZF frame is two
  `read_exact` calls. Count syscalls rather than guess; threshold ≥5% off the per-locus stage.
- **L11: [open_bam.rs:556](../../../../src/ng/read/input/open_bam.rs#L556) / [:576-578](../../../../src/ng/read/input/open_bam.rs#L576-L578) — the decoded container is keyed to the *handle*, but handles are pooled anonymously by LIFO `pop`**
  (new; mechanism High, magnitude unmeasured — no threaded caller exists yet). With T callers on one
  file a worker can be handed another worker's container, and H2's 129× is a **hit-rate** win. The
  decomposition must therefore use contiguous locus blocks (cheap) or a handle lease (an API
  change). Scoped explicitly to *not* affect the sample-level fan-out of H7.
- **L12: [locus_generation/mod.rs:330-332](../../../../src/ng/locus_generation/mod.rs#L330-L332) and [:487](../../../../src/ng/locus_generation/mod.rs#L487) — what blocks a threaded fan-out inside the library today** (new, compiler-verified).
  `SsrGenerator<WindowedRefSeq, …>`, `WindowedRefSeq` and `SampleReads` are `Send` (and `SampleReads`
  is `Sync`); what is not is `GeneratorSlot`'s `Box<dyn LocusGenerator<S>>`, which lacks `+ Send`,
  making `GeneratorSet` and `SampleLocusObservationsIterator` `!Send` — and the iterator owns
  `SampleReads` by value, so locus-level parallelism is inexpressible through the public surface.
  One bound and one borrow.
- **L13: [examples/ng_ssr_cohort_stutter.rs:199-200](../../../../examples/ng_ssr_cohort_stutter.rs#L199-L200) — the shared `BufWriter<StdoutLock>` becomes one lock per *row*, not per locus** (prior L16,
  sharpened): `write_locus` issues one `writeln!` per row, ~10⁶ rows. Fix before threading; per-worker
  output files also restore per-locus row grouping.

### Speculative

Filed so the shapes exist; **do not act without an experiment that contradicts "this won't
matter"**. Details in the per-category files.

- `ssr_unit_robust.rs:238` — the 40-byte score cell straddles a 64-byte line in exactly 32 of every
  64 consecutive cells. Asm-first; the working set is L1-resident either way (see the Note).
- `raw_chrom_reader.rs:340`/`:357`/`:378` — `self.chrom_name.clone()` is evaluated **eagerly** on
  every reference fetch to build an error that is almost never constructed. Three lines, no API
  change; new.
- `ssr.rs:1592` — `segment.clone()` copies the segment's owned contig name into every locus, per
  sample. A variant exists that avoids the `pub`-field change the prior review rejected.
- `ssr.rs:1624` — the `outcomes` `Vec` (prior L13); becomes a two-line freebie once L7 makes
  `Classified` `Copy`, and only then.
- `ssr.rs:88`, `:1650`, `:1656`, `:1657` — the four `Box<[u8]>` copies of the reference window per
  locus visit. Re-judged **won't fix**: refuted as a time cost by the 0.16% margin-fetch attribution.
- `region_query.rs:245` — the BAM source's `overlaps` re-walks the CIGAR of every over-returned
  record. Same mechanism as H2 on the other container, but the record was just decoded and the
  filter walks the CIGAR again anyway. Fold into any change that already touches footprints.
- `ssr.rs:1331-1344` (with `ref_seq.rs:399`) — closing the "Arc gap" with **one shared**
  `fasta::Repository` puts noodles' `RwLock` on the per-locus reference fetch: a read lock per fetch,
  and a **write** lock held across a whole-contig read on a miss. So H3's fix and the sharing
  decision are one decision, and the sharing half is a concurrency question.
- The per-locus read cap is 1,000 and the deep workload runs at ~412 reads/locus, so **it never
  binds** (`reads_capped = 0`). The DP cost is exactly linear in this number, so the 67% is 67%
  *because* nothing subsamples it. Output-changing (a statistical-precision decision): flagged, not
  proposed.
- `f64` → `f32` or fixed-point scores would halve the ring and open SIMD lanes, and is **not**
  bit-identical — `emission.rs` makes bit-equality with production's `EMISSION_LN` the module's only
  hard oracle. Not a perf decision.

### Note

- **Closed by measurement: the widen-and-retry is not a lever.** Two categories flagged the second
  DP at `ssr.rs:698-713` as a possible 2× on the hottest function, unmeasurable by any profile
  (both calls inline to one symbol). Counters: **42 retries in 40,337 calls (0.10%) on tomato, 513
  in 296,548 (0.17%) on HG002.** `delimit calls == reads_fetched` exactly on both, confirming one DP
  per kept read and no hidden third call. Below hot_loops' own "close the finding" threshold of 2%.
- **Refuted: reads admitted on overlap do not pay a wasted DP.** 40,221 of 40,337 fetched reads
  (99.7%) and 1,894,009 of 1,923,175 (98.5%) yield an observation; `no_border_anchored` is **1 and
  16 reads**. `extract_region` already shrinks the matrix height for a marginal overlap. Partials
  are more than half the observations and are deliberate output, not waste.
- **Banding will not help this frame — do not port it across.** `ssr_best_path_flat_gap.rs`
  (algorithm 3) *is* banded and proven byte-identical over a 200,000-case soak; algorithm 4, **4u
  (the current default) and every round-2 variant are not** (`grep -c -i band`: 69 vs **0**).
  Porting looks free and would exclude almost no cells: that band's own `left_flank + right_flank`
  term is 60 while the median frame is 68 wide. The lever on the same cost is H1, not the band.
  (Independently: a band around the CIGAR diagonal *narrow enough to matter* is not provably
  bit-identical here — the long-allele reads this delimiter exists to measure are exactly the ones
  whose optimal path leaves it — so that would be a new algorithm for the bake-off, not an
  optimisation of this one.)
- **There is no cache-miss finding to make, and the prior review's sizing was ~4× out.** At real
  dimensions (`flank_bp = 30`, measured tract distributions) the DP is ~68×68 and its **entire
  working set is 28–42 KB at the median, 66 KB at the p99 tract** — L1-resident on the host, and
  grow-and-keep, so after the first read at a locus there are no cold misses at all. L1/L2's
  mechanism is **addressing and store µops**, not misses. The PMU is unavailable to check either
  way, which is stated in every finding that would need it.
- **The DP is not autovectorizable and SIMD is not the lever.** The column loop carries a serial
  dependency inside the row (`del` at `(i, j)` reads `cur[j-1]`, written by the same iteration's
  predecessor). The only vector-friendly reformulation is an anti-diagonal sweep, which is a rewrite
  of the fill, the traceback indexing and the scratch layout for a 150-column axis. Recorded so
  nobody re-derives it.
- **Fill dominates traceback by ~40×** (arithmetic from the code, not a measurement — no profile
  splits them). The traceback's real cost is paid *in the fill*, as L1's store traffic.
- **The quality gate and the CIGAR mapping are not hot.** `complete_or_low_quality` is 27/38,054 and
  33/29,884; `read_footprint`/`ref_to_read` do not appear. Nothing to do in `read_region`.
- **The margin fetch is confirmed dead as a cost** — `RefSeq::fetch_into` is 61 samples = 0.16%,
  consistent with H4's measured revert this morning. Anyone reading H4 in the prior review should
  read its *Applied* section too.
- **No allocation finding in this scope can be sold as a speedup.** Allocator self-time is 1.9%
  (shallow) and **1.2% in the deep regime that allocates ~60× more per locus**, and two allocation
  findings were implemented, measured at zero and reverted this morning. L7–L9 are filed with
  allocation-count/RSS confidence stated separately from wall clock for that reason.
- **`fai::fs::read` is internally buffered** — the finding against it (H3) is *how often it is
  called*, not how it reads. Filed explicitly because "add a `BufReader`" is the wrong reflex here.

## 6. Out-of-scope observations

- **`reference_info.rs:546`/`:626` — `FastaPass::push_byte` processes 3.1 GB one byte at a time**
  through a state machine before the MD5 sees it. That is the larger half of L6's digest cost and a
  `hot_loops` change in a module this review did not otherwise cover. Its geometry checks
  (non-uniform line length, duplicate names) are exactly what a chunked loop drops by accident.
- **The prior review's M5 is still live in the driver, verbatim.**
  `examples/ng_ssr_cohort_stutter.rs:16-17` still asserts the `--regions` rationale is "~8 minutes
  per sample" of walk waste (measured: 2.44 s — a `cargo run` artefact), and `:206-208` still asserts
  restriction "cannot change what a covered locus is" (measured: 6,489 vs 5,534). Neither has an
  owner response yet, and both are the kind of claim a future reviewer takes as a baseline.
- **`src/ng/read/input/test_fixtures.rs` is still `#[cfg(test)] pub(crate)`** (`mod.rs:28-29`), so
  M1's `fetch`/`locus`/`walk` bench groups remain blocked. Only the delimiter group is unblocked —
  which is why measurement item 1 proposes it alone.
- **`examples/ng_ssr_cohort_stutter.rs:193` — `contigs.clone()` per query** clones the whole
  `ContigList` (an owned `String` per contig): 13 per query on tomato, **2,580 on GRCh38**. Removed
  as a side effect of H3.
- **The region-typing walk is 14.7% of the deep single-sample run** (4,395 / 29,884 on five contigs)
  — out of scope here and already covered by the prior review's L1 and its `scan_set` note, but it
  is no longer negligible once the fetch branch shrinks.
- **The driver's `reads_no_border` column merges all three no-observation reasons**, which is why
  the retry rate needed a temporary counter. Splitting that column (or printing
  `generator.counts()`, as the sibling `ng_ssr_loci_dump.rs:232` already does) would have made the
  W-E experiment unnecessary.

## 7. What's already good

- **The scratch discipline holds where it matters most.** `UnitRobustScratch::resize` is grow-and-keep
  and never re-zeroed, and the whole DP + traceback was read twice today (by two different
  categories) to confirm **zero allocations per `align`** — on the function that is 36–67% of the
  run. The three generator-owned buffers (`margin_buffer`, `qual_buffer`, and the aligner's scratch)
  are the same pattern applied consistently.
- **Static dispatch where it counts, with the reasoning written down.** The `RepeatDelimiter` trait
  alias keeps the aligner a type parameter, and both `cargo asm` listings confirm it: `align` is
  fully inlined into `classify::delimit`, with no `dyn` indirection anywhere in the per-read loop.
- **The tally is order-independent by construction and says so.** `entry().or_default()` (one hash,
  not two), integer moment accumulation, a single `sort_unstable_by` on a total order, and an
  explicit doc note naming `q_sum`'s `f64` fold as the one order-sensitive field — which is what let
  this review clear the `FxHashMap` swap as bit-identical in one reading instead of one experiment.

## Author response convention

Address each finding by its identifier (H1, L2, M-item, …) with one of: `applied in <commit>` /
`experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` /
`won't fix because …`. The "experiment shows no gain" path is expected and welcome — this morning's
review produced two of them, and this review's own two closed candidates (the widen-retry, the
wasted-DP hypothesis) are the same thing found before the code was written.

---

# Applied — same-day experiments (2026-07-26)

**Six changes, every one measured, every one kept: the shallow tomato cohort is 1.51× faster and
the deep HG002 sample 1.19×, with byte-identical output on both.** H1 (`flank_bp`) was left alone
at the owner's instruction — it is an empirical question about partial observations, not a perf
change.

Final interleaved A/B, three rounds per side, alternating binaries (never back-to-back batches),
`md5` of the emitted TSV checked on **all twelve runs**:

| workload | before | after | |
|---|---|---|---|
| tomato, 51 samples, SL4.0ch01 BED, CRAM | 70.18 / 70.73 / 66.91 s → **69.27 s** | 46.03 / 46.02 / 45.95 s → **46.00 s** | **−33.6% (1.51×)** |
| HG002 300×, chr16–20 Tier BED, BAM | 51.48 / 50.41 / 50.79 s → **50.89 s** | 42.83 / 42.79 / 42.66 s → **42.76 s** | **−16.0% (1.19×)** |

`md5 9409dc94253d397155222b231be7afa3` (tomato, 950,637 rows) and
`md5 699d50d5379107127db0773b4088cd26` (HG002) on every run, before and after. Peak RSS is
unchanged (5,457 MB → 5,415 MB on tomato) and system time collapses from **8.15 s to 1.04 s** —
that second number is the syscall half of H3 landing.

## M1 — the delimiter bench, built first

[benches/ng_ssr_delimiter_perf.rs](../../../../benches/ng_ssr_delimiter_perf.rs), registered in
`Cargo.toml`. Two groups (`frame` over four (period, tract) shapes × reference/expanded reads;
`depth` over N ∈ {1, 7, 30, 100, 412} reads sharing one scratch), with the three guards section 3
called for: `black_box` on both sides, a measured-length assertion inside the timed body, and more
than one frame per group. It calibrates against reality — 24.9 µs per read at the HG002 median
shape, against ~35 s of `delimit` in a 52 s run over 1.9 M reads.

**It also answered a question for free.** The `depth` group is *flat*: 24.9 µs/read at depth 1 and
25.3 µs/read at depth 412. Nothing is being amortized across the reads of a locus today, which
independently confirms L3 (the per-locus prepared frame) is worth little — and it is why L3 was not
implemented.

## The six changes, in the order they were applied and measured

Each row is its own interleaved A/B against the *same* base binary, so the cumulative column is
measured, not summed. Delimiter-bench figures are criterion `--baseline`, p < 0.05 throughout.

| # | change | delimiter bench | tomato | HG002 |
|---|---|---|---|---|
| 1 | **H4** `PerQualityEmission` holds the resolved table | −1.3% (p1) … −7.7% (p6) | −2.3% | −2.5% |
| 2 | **H6** slip emission per row, phase carried not divided | −2.2% … −25.1% | (−2.3%) | −4.1% |
| 3 | **H2** per-record footprints resolved once per decode | — | **−15.5%** | (unaffected) |
| 4 | **H3** one shared reference reader + a reposition rule | — | **−28.4%** | −8.1% |
| 5 | **H5-lite** the column axis resolved once per read | −6.5% … −24.3% | −30.4% | −10.6% |
| 6 | **L1** backpointers packed into a `u16` | **−13.0% … −32.9%** | **−33.6%** | **−16.0%** |

Codegen, `cargo asm --simplify "classify::delimit"`, start → end: **`ldapr` 7 → 0**,
**`Once::call` 7 → 0**, **`udiv` 11 → 6**, `panic_bounds_check` 32 → 31, and the whole function
1,781 → 1,320 lines.

### 1. H4 — the emission table is a field, not a `LazyLock` deref

[emission.rs](../../../../src/ng/alignment/emission.rs): `PerQualityEmission` gains one
`&'static [BaseScores; 256]` resolved in `new()`. The DP's inner loop had **seven acquire-loads and
seven live `Once::call` sites**, and because LLVM must treat the loop's live values as
call-clobbered at each of them, the per-row constants — including the read base — were re-read from
the stack every cell. Bit-identical (same table, same values;
`per_quality_table_is_bit_exact` still pins them).

### 2. H6 — the whole-unit slip emission, once per row

[ssr_unit_robust.rs](../../../../src/ng/alignment/ssr_unit_robust.rs): the sum over the unit's
`period` bases depends on (read row, motif phase) only, so it is resolved into a `[f64; 6]` per row
instead of per tract column. The `k`-ascending accumulation order is preserved — `f64` addition is
not associative, and that is what keeps it bit-identical.

Then the modulo went too: the first version indexed by `(column - left_flank_len) % period` and the
fill used `(phase + k) % period`, which LLVM unrolled into **24 `udiv`s**. Carrying the phase
forward (one add and one compare per column, one conditional subtract in the fill) took it back to
6 and turned a **+0.5% regression at period 1 into −2.2%** — worth noting because period 1 is 66% of
tomato's reads and the naive form would have made the common case slower.

### 3. H2 — a record's footprint is resolved once per decode, not once per query

[region_query.rs](../../../../src/ng/read/input/region_query.rs): `DecodedContainer` gains a
`Vec<Footprint>` built beside the records. The container cache (this morning's H2) removed the
*decode* from the per-query path but left the *interpretation*: every query re-walked every
record's CIGAR through `RecordBuf::alignment_end` to select ~7 reads out of ~10⁴. **−15.5% on the
CRAM workload; the BAM workload is untouched, as it should be.** No new invariant — the full linear
scan is kept, only the CIGAR walk is hoisted, so the review's step-2 cursor (and its `max_span`
trap) was not needed to get the win.

### 4. H3 — one reference reader for the whole walk, and a reposition rule that makes sharing pay

Three pieces, and the third is what makes the first two safe:

- [ref_seq.rs](../../../../src/ng/ref_seq.rs): blanket `RefSeq` / `RawRefSeq` / `ContigTable`
  impls for `Arc<T>`. All three traits are already `&self`-only, so this is pure forwarding — it
  just gives `FnMut() -> R` an *owned* handle to hand out, which is the "Arc gap" `SsrGenerator`'s
  own type doc names.
- [ng_ssr_cohort_stutter.rs](../../../../examples/ng_ssr_cohort_stutter.rs): the driver builds one
  `Arc<WindowedRefSeq>` and hands clones to both the margin fetch and the per-query factory,
  instead of constructing a fresh reader per query.
- [raw_chrom_reader.rs](../../../../src/ng/raw_chrom_reader.rs): `fetch` now **repositions** when
  the requested window lies more than one read chunk clear of the buffered one, instead of reading
  through the gap. Without this, sharing one reader across BED-restricted loci tens of kb apart
  would read and buffer every base in between — the review's L10, which it correctly called a
  *prerequisite* for sharing rather than an independent win.

**−28.4% on tomato, −8.1% on HG002**, and it is what takes system time from 8.15 s to 1.04 s: the
~564k `open(2)`s per chromosome are gone, and `__open` / `read` have left the profile's top ten
entirely.

The reposition branch is the one silent-wrongness risk in this set (a wrong re-seek returns
plausible bases from the wrong coordinate), so it has its own test —
`a_fetch_clear_of_the_window_repositions_instead_of_reading_through_the_gap` — over a 400 kb contig
whose base at every position is a function of that position, asserting **both** halves: the bases
are right at the new coordinate (forward *and* backward), and `buf_len` stays under the gap. A
missing branch fails the second assertion; wrong arithmetic fails the first.

### 5. H5-lite — the column axis resolved once per read

The review proposed splitting the column loop into five constant stretches. What landed is the
contained half: a `Vec<ColumnPlan>` in the scratch holding `gap_open`, `gap_open_terminal`,
`in_tract` and `unit_deletable` per column, built once per `delimit` call and read once per cell.
That deletes L4's per-cell `ccmp`/`cset`/`fcsel` chain (with its four stack reloads) and L5's
unrolled six-step range test — **without** the five stretch boundaries the fission version would
have made load-bearing, each of which can degenerate and none of which fails loudly. Same
mechanism, a fraction of the risk. **−6.5% at period 1**, which is where the tomato workload lives.

### 6. L1 — five byte stores per cell become one halfword

`[State; 5]` (5 bytes, never naturally aligned, five separate `strb` and a `×5` address multiply)
becomes `Backpointers(u16)` — five 3-bit fields. **This was the single largest DP win: −13% at
period 1, −33% at period 6**, on an array written ~10⁷ times per deep locus to serve ~230 reads.
The packing's invariant fails the *build*, not a test, via
`const _: () = assert!(STATES as u32 * STATE_BITS <= u16::BITS)`, and `State::from_code` panics on
an out-of-range field rather than answering `Match` and hiding a broken packing behind a plausible
traceback.

## What was deliberately not done

- **H1 (`flank_bp`)** — owner's call: the flank width is to be settled on empirical grounds
  (what partial observations are worth), not on the 31–49% it is worth in wall clock.
- **H5's full loop fission, and dropping the two provably-`−inf` `Match` candidates outside the
  tract.** The contained version (5) took most of the win; the rest needs a five-boundary invariant
  whose failure mode is a silently wrong measurement. Not worth it on top of a 1.5×.
- **L2 (ring flattening).** The review said do it last, if at all, and that H5 would shrink it
  first. It did.
- **L3 (per-locus prepared frame).** Refuted by the bench's own depth axis before it was written —
  per-read cost is flat from depth 1 to depth 412.
- **L7–L9 (the allocation findings).** The category's own preamble said no allocation finding here
  can be sold as a speedup, and the post-change profile agrees: allocator frames are ~2.8% of a run
  that is now 51% delimiter.
- **H7 (the per-sample fan-out).** A driver/threading change and a separate piece of work; the
  single-threaded cohort now finishes in 46 s per chromosome, which changes its urgency.

## Where the time is now (tomato, 51 samples, post-change profile, 32,438 samples)

```
16570  classify::delimit                                    51.1%
 4211  md5::compress                                        13.0%   (noodles CRAM slice MD5 — L5)
 1923  noodles_cram::…::Block::decode                        5.9%
 1698  ng::tandem_repeat::scan_window                        5.2%   (region typing — out of scope)
 1288  noodles_cram::…::Slice::records                       4.0%
```

`Cloned::next` (19.6%), `RecordBuf::alignment_end` (9.4%), `__open` (6.1%) and `read` (4.9%) are
all gone. The delimiter is now more than half the run *because* everything around it shrank, and
what remains beside it is mostly noodles' own CRAM decode and reference-MD5 — which L5 measured as
already at its floor from our side, with `embed_ref=1` re-encoding as the one remaining lever.

**Test status:** full suite **2,450 passed / 0 failed** (was 2,431 this morning), including the
delimiter parity oracle and the 12 alignment-file integration tests; `cargo fmt` clean, `cargo
clippy --lib` clean.
