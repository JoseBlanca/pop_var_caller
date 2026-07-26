# Performance Review: ng-ssr-cohort-stutter
**Date:** 2026-07-26
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** the ng STR per-locus read path (region-typing walk → `SsrGenerator::next_locus` → per-locus indexed read query → `SsrUnitRobustAligner`), as exercised by a cohort walk
**Verdict:** Apply the listed wins
**Hot-path evidence:** two macOS `sample` profiles + twelve wall-clock experiments + one `cargo asm` listing, all collected for this review (full log: `tmp/perf_review_2026-07-26_ng-ssr-cohort-stutter/measurements.md`)
**Status:** H2 **applied** (CRAM 174.2 s → 1.35 s, 129×, output byte-identical); H1 **refuted in the proposed form** and deferred with its ceiling measured; H4 and L2 **experiments show no gain, reverted**. See *Applied — same-day experiments* at the end.

---

## 1. Scope and constraints

**What was reviewed.** The module-level path a cohort stutter dump drives: the region-typing
walk, the per-locus STR generator, the per-locus indexed read query behind it, the reference
accessors both use, and the default STR delimiter.

**Reviewed against:** branch `ng-aligner-bakeoff` at `935ac2e`, plus the uncommitted
[examples/ng_ssr_cohort_stutter.rs](../../../../examples/ng_ssr_cohort_stutter.rs) as it stood
during the review. The driver changed under me twice while the review ran (the owner was wiring
`--regions` into it), so every number below was taken with a stable probe instead — see
*Measurement harness* at the end of section 3.

**Throughput / latency targets, input sizes, hardware.** 51 tomato samples × ~14,455 targeted SSR
regions (`benchmarks/ssr_tomato1`, inputs sliced to `ssr_regions.bed`, SL4.0 reference ≈ 780 Mb,
~800k reads per sample). The whole cohort walk should finish in minutes to a couple of hours — it
is an exploratory analysis tool that gets re-run as the question changes. The same `next_locus`
path is what Stage-1 `ssr-pileup` will use per sample in production, so wins here are not
throwaway. Hardware: macOS / Apple Silicon host and the Debian 12 aarch64 dev container (Apple
`container` VM — **no PMU**, so hardware counters, `perf c2c` and `perf sched` are unavailable;
user-space sampling works on both).

**Hot-path evidence available.** Two symbolicated `sample` profiles (BAM and CRAM inputs), a
stage-subtraction experiment series, a contig-length sweep of the walk, and a `cargo asm` listing
of the delimiter at the caller the profile names. Sampling was **not** blocked in this
environment.

**In-scope files.**

- [examples/ng_ssr_cohort_stutter.rs](../../../../examples/ng_ssr_cohort_stutter.rs) — the driver
- [src/ng/region_typing/mod.rs](../../../../src/ng/region_typing/mod.rs) — `TypedRegionIterator`
- [src/ng/locus_generation/ssr.rs](../../../../src/ng/locus_generation/ssr.rs) — `SsrGenerator::next_locus`, `fetch_capped_reads`
- [src/ng/read/input/mod.rs](../../../../src/ng/read/input/mod.rs) — `SampleReads::reads_in_region`
- [src/ng/read/input/region_query.rs](../../../../src/ng/read/input/region_query.rs) — `BamRegionSource` / `CramRegionSource` (**added to scope by the call graph**: it is where the profile lands)
- [src/ng/read/input/open_bam.rs](../../../../src/ng/read/input/open_bam.rs) — `AlignmentFile` (reader pool, index, `fasta::Repository`)
- [src/ng/ref_seq.rs](../../../../src/ng/ref_seq.rs), [src/ng/raw_chrom_reader.rs](../../../../src/ng/raw_chrom_reader.rs) — `WindowedRefSeq`
- [src/ng/alignment/ssr_unit_robust.rs](../../../../src/ng/alignment/ssr_unit_robust.rs) — `SsrUnitRobustAligner`
- [Cargo.toml](../../../../Cargo.toml), [.cargo/config.toml](../../../../.cargo/config.toml), [scripts/dev.sh](../../../../scripts/dev.sh) — build and harness configuration

**Deliberately out of scope.** `src/ssr/` (production, frozen — ng is the from-scratch caller);
the correctness of what the delimiter measures (a correctness review's business); patching
vendored noodles.

**Categories dispatched.** `methodology` (always), `io_and_syscalls` (the profile puts 87–99% of
the per-locus stage in the read query), `allocations` (per-locus `Box`/`Vec`/`HashMap`/clone
churn), `hot_loops` (the DP delimiter and the per-window scan loop), `data_layout` (the DP's
`Vec<Vec<[f64; 5]>>` ring, the per-locus output types), `concurrency` (single-threaded today, but
the shared state and per-sample parallelism are both live questions).

## 2. Verdict

**Apply the listed wins** — after one correction that reorders everything else.

**The 7m47s baseline is not this code path.** Re-measured, the identical work — same flags, same
sample, same chromosome, byte-identical output counts — takes **10.5 s on the host and 12.3 s in
the container**, including when run from the owner's own container-built binary. The ~58-hour
cohort extrapolation is off by roughly 40×. What the 7m47s measured was the
`./scripts/dev.sh cargo run` harness: ~11 s of container start + a `rustup` channel re-sync + a
crates.io index update on *every* invocation (`dev.sh` mounts no `CARGO_HOME`/`RUSTUP_HOME`), plus
whatever release rebuild `cargo run` decided to do inside the timed region (27 s after a
one-file `touch`; a cold `target-container` release build of this crate at `lto = "fat"` /
`codegen-units = 1` is minutes). I could not reproduce anything near 7m47s from the code path, and
I/O is not a factor (342 MB of cold CRAM read in 0.08 s in-container). **The single most valuable
change from this review is to stop timing `cargo run`** — and to correct the driver's doc comment,
which has already recorded the artefact as a design rationale (M5, M10).

**There is a real 126× problem, and it is CRAM.** The cohort inputs are the sliced CRAMs
(confirmed in `benchmarks/ssr_tomato1/bench.config.sh`), and on CRAM the same BED-restricted
chromosome takes **171.60 s instead of 1.36 s** for byte-identical results. 35% of that run is
`md5::compress` — noodles re-verifying the slice's reference MD5 on every per-locus query — and
most of the rest is re-decoding whole CRAM containers per locus. At cohort scale that is ~21 hours
against ~12 minutes for the same data as BAM. **The input format, not the code, currently decides
whether this tool meets its target.**

**The dominant code-level defect is one pattern, and it is the same pattern in both formats:** a
per-locus indexed query re-decompresses compressed data it has already decompressed. `bgzf::io::Reader::seek`
is unconditional — it never checks whether the block it holds is the target block — so with one
typed segment every ~426 bp and a BGZF block covering ~64 KiB of records, the same block is
inflated hundreds of times. A whole-file linear decode of the same reads (`samtools view -c`) is
**0.33 s against 8.04 s**.

Everything else is small on this workload and the report says so: the DP delimiter is 6.3% of the
per-locus stage, the reference margin fetch 6.2%, allocator frames ~4% (and most of those are
noodles' per-record tag decode, out of scope), and **the per-query reference-reader construction
the owner asked about is ~1% — that lead is refuted.** Three categories filed nothing at Hot-path
and said so in a preamble rather than dressing count findings up as time findings.

## 3. Measurement plan

Ordered so each entry unblocks the next. Entries 1–3 are done, with results; 4–8 are the
experiments the findings gate on.

**1. Re-baseline off `cargo run`. (Done — do this before believing any future timing.)**

```
./scripts/dev.sh cargo build --release --example ng_ssr_cohort_stutter     # build separately
./scripts/dev.sh bash -c 'time ./target-container/release/examples/ng_ssr_cohort_stutter \
    --contigs SL4.0ch01 $HOME/genomes/.../S_lycopersicum_chromosomes.4.00.fa <sample.bam> > /dev/null'
```

Result: `real 0m12.283s`, tallies identical to the 7m47s run. Threshold: any timing quoted for
this path must come from the binary, with stdout redirected.

**2. Stage subtraction. (Done.)** Probe modes that stop after successive stages, host release,
BAM, 1 sample, SL4.0ch01: walk 2.44 s → +query and depth cap 7.21 s → +delimit and tally 8.04 s.
The read query is ~90% of the per-locus stage.

**3. Sampling profiles. (Done.)** `sample <pid> 90 -file …` against a **native host** build
(`cargo build --release --example …`; the container binary is Linux and cannot be sampled from
macOS). Both listings are in the audit directory. BAM: 2162 samples, read query 87.4%, BGZF
inflate ~71% of total. CRAM: 50,894 samples, read query 99.3%, `md5::compress` 35% of total.

**4. Confirm H1's mechanism cheaply, before writing the fix.** Add a temporary counter beside the
`seek` in `BamRegionSource::read_next` recording `(queries, seeks, distinct
chunk.start().compressed() values)` over the ch01 walk. H1 predicts `seeks >> distinct offsets`
(order 10× or more); if they are comparable, H1 is refuted and the inflate is genuinely new data.

**5. The decode-reuse experiment (gates H1, H2).** BAM: group queries so one seek+inflate serves a
group of segments. CRAM: hoist the decoded container onto the file handle, keyed by `.crai`
offset, and decode only the slice the entry names. Metric: per-locus stage wall + the three
tallies. **Thresholds:** BAM 8.04 s → ≤ 2.7 s (floor 0.33 s); CRAM 171.6 s → ≤ 20 s. **Gate:**
`covered` / `observations` / `reads_fetched` and the TSV rows byte-identical — extend the existing
T5 oracle (`region_query.rs:828`, "the indexed query returns exactly what a linear scan returns")
to the batched source.

**6. The walk experiment (gates L1).** Profile the *walk* stage alone (it was never profiled — the
90 s sample covered the per-locus stage). If `emit_into` shows up, add the monotone cursor; if it
does not, close L1. Metric: `walk` wall on `ssr_regions.bed`, 26.66 s baseline.

**7. CRAM→BAM conversion as the zero-code answer (gates H2's priority).** `samtools view -b` the
51 cohort CRAMs once, re-run. Expected from the measured per-sample numbers: ~21 h → ~12 min.
Cost: ~2× the disk of the CRAM set.

**8. Per-sample process parallelism (gates H3).** Run N single-sample processes over the same BED
instead of one N-sample process, per the project's "benchmark each caller with its native
parallelism" convention. Metric: cohort wall and peak RSS. Threshold: near-linear to core count;
if it is, no threading work inside the generator is justified at all.

**Measurement harness — deleted at the owner's request; rebuild it from this description.** All
numbers came from a throwaway `examples/ng_perf_probe.rs`. Its shape is the reusable part, and M1's
bench should reproduce it:

- Six modes, each a strict superset of the previous, so **subtracting adjacent wall times isolates
  one stage**: `walk` (region typing only, no samples) → `collect` (+ collect every `SsrSegment`
  into memory, so the nesting variants below share one typing pass) → `fetch` (+ per (segment,
  sample) `fetch_capped_reads` only, no delimiting) → `locus` (+ the driver's real body,
  `begin_segment` + `next_locus`) → `locus_sample_outer` (`locus` with the loop nesting swapped) →
  `fetch_grouped` (one indexed query per 32 kb bucket of adjacent segments, the H1 prototype).
- `--contigs a,b` / `--regions r.bed` matching the driver, and `--samples N` truncating the sample
  list so cost-per-sample is measurable without re-typing.
- `--dump <path>` — **the identity gate**: one row per observation (`sample_index`, contig, start,
  end, reference bases, observed bases, coverage tag, `num_obs`, `num_fwd`, `q_sum` to 6 places,
  `mapq_sum`, `mapq_sum_sq`) plus a row carrying `reads_without_observation` and
  `reads_discarded_by_cap`. Compared by `md5`. This is what made every "no gain" verdict below
  trustworthy: a change is only measured after its dump matches.
- It reproduced the driver's tallies exactly (`covered=6489`, `observations=20639`,
  `reads_fetched=43991` on W1), which is what licensed using it as the driver's stand-in.

**M1 still stands**: `src/ng/` has no committed benchmark, so this path currently has no repeatable
measurement in the tree at all.

## 4. Build / toolchain configuration

The build-configuration lever is **already pulled** — `lto = "fat"`, `codegen-units = 1`,
`panic = "abort"`, `opt-level = 3` on release, a pinned toolchain, and a `[profile.profiling]`
that already exists. Do not spend this review's budget on build flags; 71% of the stage is zlib
inflate driven by query amplification, and no flag closes that. **In this category the
*measurement infrastructure*, not the build, is the whole problem.** All eleven items are Likely
except M11 (Note).

**Measurement infrastructure**

- **M1 `[bench]` `src/ng/` has no committed benchmark at all.** All 7 registered benches
  ([Cargo.toml:140-167](../../../../Cargo.toml#L140-L167)) cover production paths;
  `grep -rn "ng::" benches/` finds no `SsrGenerator`, `TypedRegionIterator`, `SampleReads` or
  aligner. The delimiter bake-off commits already on this branch were merged with no harness.
  Proposed `benches/ng_ssr_locus_perf.rs` with four groups (`walk`, `fetch`, `locus`, `delimit`);
  every entry point needed is already `pub`, and a synthetic indexed-BAM fixture builder exists
  twice in-tree (`examples/dhat_ng_merge.rs:306-322`). Blocker to fix first: `test_fixtures.rs` is
  `#[cfg(test)] pub(crate)` (see section 6). Cost: one bench file, mostly fixture construction.
- **M2 `[bench]` the review's instrument is untracked and self-marked for deletion** —
  `examples/ng_perf_probe.rs`. Every number in this report came from it; its stage-subtraction
  design is the thing worth keeping. Proposed promoting it (or folding its modes into M1's bench)
  rather than deleting. **Owner decision: deleted** once the experiments were done — so its design
  is recorded in section 3 instead, and M1 is the only route back to a repeatable measurement.
  *Closed: won't fix.*
- **M3 `[bench]` nothing has ever been measured at the target workload, and the only measurement
  of the real input format lands past the target.** Confirmed in-tree: `ssr_regions.bed` = 14,455
  spans, 51 `.bench.cram` files, and `benchmarks/ssr_tomato1/bench.config.sh:27-31` confirms the
  cohort inputs are **CRAM**. Extrapolating the measured per-sample numbers (labelled as
  extrapolation, not measurement): BAM ≈ 12 min, CRAM ≈ 21 h. Run a `--samples {1,3,8,51}` sweep in
  both formats before any code change, so H1/H2 have a real baseline. Cost: machine time only.
- **M4 `[bench]` the driver doubles as the timing harness but bundles TSV formatting into the
  timed region and has no stage timers**
  ([:275-292](../../../../examples/ng_ssr_cohort_stutter.rs#L275-L292)). It reports one bundled
  `real`, which is how a 40× harness artefact went unnoticed. Cost: ~20 lines — two `Instant`s and
  a stage line on stderr.
- **M5 `[bench]` the 7m47s artefact is already baked into the driver's doc comment as a design
  rationale** ([:13-17](../../../../examples/ng_ssr_cohort_stutter.rs#L13-L17)) — see section 6 for
  the correction. Root causes: M10's missing cache mounts, plus `lto = "fat"` /
  `codegen-units = 1` making a cold `target-container` release build expensive.
- **M6 `[bench]` the BED walk is not cheaper than the contig walk** (2.88 s vs 2.44 s) despite
  typing 37× fewer segments, and the covered-locus count changes (6,489 → 5,534) — so the two modes
  are not comparable workloads and the documented rationale at
  [:216-219](../../../../examples/ng_ssr_cohort_stutter.rs#L216-L219) is wrong on both counts.
  Nobody noticed because there is no bench and the driver reports one bundled `real`. Cost: a bench
  group plus a one-off segment-diff script.
- **M7 `[bench]` the workload has a counted invariant and nothing asserts it.** The tallies
  (`covered`, `observations`, `reads_fetched`) are exactly what a decode-reuse change must not
  move, so M1's bench bodies should `assert_eq!` them against committed constants. Cost: two
  asserts and two constants per bench body — and it is what makes H1/H2 safe to merge.

**Build configuration**

- **M8 `[profile]` iterate under `[profile.profiling]`, not release.** The quoted profiles came
  from `debug = "line-tables-only"` + fat LTO ([Cargo.toml:45-49](../../../../Cargo.toml#L45-L49)),
  so call-site attribution is sound (that is how the `ssr.rs:1610` / `:1626` / `:1591` split was
  read) but inlined frames are collapsed — which is why the delimiter appears only as
  `classify::delimit`. `[profile.profiling]` ([Cargo.toml:58-62](../../../../Cargo.toml#L58-L62))
  already exists and was not used; use it for code-level iteration on L3–L7, and cross-check
  against release because its codegen differs.
- **M9 `[build]` no `target-cpu` for `aarch64-unknown-linux-gnu`.**
  [.cargo/config.toml:15-23](../../../../.cargo/config.toml#L15-L23) sets `target-cpu` for
  x86_64-linux and aarch64-macos but not for the aarch64-linux container the owner measures in, so
  the container builds at the armv8-a baseline while the host gets `apple-m1`. Host and container
  numbers are therefore two configurations, not one. Small on a zlib-bound path; will not stay
  small for `delimit`. Four lines of config; the honest risk is portability — prefer
  `neoverse-n1` or an explicit feature list over `apple-m1` if the container image must run on
  other aarch64 hardware.
- **M10 `[bench]` mount the cargo/rustup caches in the container.**
  [scripts/dev.sh:118-124](../../../../scripts/dev.sh#L118-L124) mounts no `CARGO_HOME` /
  `RUSTUP_HOME`, so every `./scripts/dev.sh cargo …` re-syncs the toolchain channel (it downloads
  5 components) and re-updates the crates.io index. That is the ~11 s floor under every container
  cargo invocation and part of what produced the 7m47s. Two `-v` flags plus a `.gitignore` entry.
- **M11 `[build]` (Note) do not pursue PGO or the allocator A/B yet.** Filed explicitly so the
  option is closed rather than forgotten: revisit after H1/H2 land, when the profile is no longer
  inflate-dominated.

## 5. Code-level findings

Full per-finding detail — including diffs, verbatim `cargo asm` excerpts and the noodles source
that confirms each mechanism — is in the per-category files under
`tmp/perf_review_2026-07-26_ng-ssr-cohort-stutter/`. Entries below are the synthesis: Hot-path and
Likely findings in full, Speculative and Note tiers condensed with pointers.

### Hot-path

#### H1: src/ng/read/input/region_query.rs:190 — every per-locus query re-inflates a BGZF block it may already hold

- **Confidence:** High
- **Hot-path evidence:** `zlib_rs::inflate::inflate_fast_help` 934 of 2162 self samples, top of
  the ranking; `copy_match_runtime_dispatch` 315; `Crc32Fold::fold` 160. Tree attribution:
  `ssr.rs:1610` (`fetch_capped_reads`) = 1887 samples (87.4%), of which
  `noodles_bgzf::io::reader::frame::parse_block` holds **1081 on the first-block path — the one
  reached immediately after the per-query `seek`** — plus 416 under `read_record`. Stage
  subtraction: the query alone is 7.21 s of the 8.04 s per-locus stage. Floor: `samtools view -c`
  decodes all 834,941 records of the same file in **0.33 s**.
- **Pattern matched:** compressed I/O — decompress once into a reusable buffer; here a fresh
  inflate of an *already-inflated* block per query, plus a random-access `seek`+`read` per query.
- **Mechanism:** [region_query.rs:197](../../../../src/ng/read/input/region_query.rs#L197) seeks
  the pooled reader to `chunk.start()` on taking a chunk. `noodles_bgzf::io::Reader::seek`
  (`noodles-bgzf-0.47.0/src/io/reader.rs:175`) is unconditional: `self.inner.seek(...)` then
  `self.read_block()`, with no check that the held block *is* the target. So each query pays an
  `lseek` + a frame read + a full DEFLATE inflate of up to 64 KiB **before it can look at one
  record**, and the block dies with the stream. The workload issues one query per typed segment —
  213,344 on a 90.9 Mb contig, one every ~426 bp — against BGZF blocks covering ~64 KiB of
  records, so every locus inside a block re-inflates that block from scratch. The reads themselves
  are decoded exactly once either way; only the repetition is removable.
- **Measurement plan:** plan entries 4 and 5 above — counter first (predicts
  `seeks >> distinct compressed offsets`), then the fix, gated at ≥ 3× on the per-locus stage with
  byte-identical tallies, then re-checked at 3 samples against the 24.73 s baseline for linearity.
- **Complexity cost:** non-trivial and honest. **(a) Group the queries** (recommended first step):
  one indexed query per BED span or fixed 16 kb bucket, decode the group's reads once into a
  `Vec<MappedRead>`, serve each segment by overlap. Adds a grouping loop in the driver, a
  `reads: &[MappedRead]` entry point beside `fetch_capped_reads`, and one live `Vec` per group
  (~3k reads at 30× for a 16 kb bucket). **(b) A per-sample forward window** reaches the 0.33 s
  floor but adds a monotonicity invariant (segments must arrive in coordinate order) that the
  current random-access contract does not have — do not take (b) before (a) is measured.
- **Suggested fix:** see `io_and_syscalls.md` §Hot-path 1 for both variants with diffs.

#### H2: src/ng/read/input/region_query.rs:424 — CRAM decodes a whole container per query, every slice of it, and MD5s the reference each time

- **Confidence:** High
- **Hot-path evidence:** same BED-restricted ch01 work, byte-identical results
  (`covered=5534 observations=18827 reads_fetched=40337`): **CRAM 171.60 s vs BAM 1.36 s.** CRAM
  profile, 50,894 samples: `md5::compress` **17,842** top of the self-time ranking (35% of the
  run), `noodles_cram::…::Block::decode` 8,568, `Slice::records` 5,273; the read query is 99.3% of
  the process and `CramRegionSource::read_next`
  ([region_query.rs:477](../../../../src/ng/read/input/region_query.rs#L477)) alone is 26,456.
- **Pattern matched:** compressed I/O, decode-once — three separate multipliers on one query.
- **Mechanism:** three costs stack, **all of them ours to fix, none requiring a noodles patch.**
  (i) `last_decoded_offset` is reset by `CramRegionSource::new`
  ([region_query.rs:349](../../../../src/ng/read/input/region_query.rs#L349)), so the
  "decode a container once" guard works *within* a query and never *across* queries — and
  consecutive loci sit in the same container. (ii) `decode_container_at`
  ([region_query.rs:434](../../../../src/ng/read/input/region_query.rs#L434)) loops
  `for slice in container.slices()`, decoding every slice although the `.crai` entry names one via
  `landmark()`. (iii) noodles MD5-hashes the slice's whole reference span on **every**
  `slice.records()` call (`noodles-cram-0.93.0/src/io/reader/container/slice.rs:359`) — and
  because these CRAMs are sliced to a scattered target BED, a slice's reference span is wide, so
  each query hashes far more reference than it reads. The fix is to stop calling `records()` again,
  not to defeat the check.
- **Measurement plan:** plan entry 5 (CRAM threshold 171.6 s → ≤ 20 s, tallies byte-identical),
  and plan entry 7 as the zero-code comparison — convert the cohort to BAM and confirm ~12 min.
- **Complexity cost:** the cache is one `Option<(u64, Vec<RecordBuf>)>` moved from the per-query
  source onto the pooled handle, plus a memory bound (one container's records, ~10k
  `RecordBuf`s per worker — measure RSS, it is not free). Decoding only the named slice is a
  smaller, independent change and should be its own commit.
- **Suggested fix:** see `io_and_syscalls.md` §Hot-path 2.

#### H3: examples/ng_ssr_cohort_stutter.rs:283-288 — 51 independent samples are serialized through one shared `SsrGenerator`

- **Confidence:** High
- **Hot-path evidence:** 3 samples over whole ch01 = **24.73 s** against 8.04 s for one — exactly
  linear, so the shared generator buys no cross-sample reuse whatsoever. Whole-genome BED, 1
  sample: `walk 26.66 s` + `perlocus 13.66 s`, i.e. the per-sample cost is ~13.7 s and the walk is
  a one-time ~27 s. Nesting order is immaterial (sample-outer 25.73 s vs segment-outer 24.73 s).
- **Pattern matched:** an embarrassingly parallel loop held serial by a `&mut self` borrow, not by
  a data dependency.
- **Mechanism:** `next_locus` is `&mut self` over pure scratch and counters
  ([ssr.rs:1345-1373](../../../../src/ng/locus_generation/ssr.rs#L1345-L1373)); the borrow
  checker, not any shared state, is what forces 51 samples through one core. The enabling fact,
  verified in source: region typing is a pure function of (reference, spans, config), so N
  processes each walking the same BED produce the **identical** segment set — the shared walk is a
  cost optimisation, not load-bearing for the `(contig, start, end)` join the dump relies on. One
  generator per worker costs T× scratch and nothing else, and it also removes a wart the driver
  already documents (it recomputes tallies from rows "because the generator's own counters are
  shared by every sample it serves").
- **Measurement plan:** plan entry 8. **Order matters:** do H1/H2 or the CRAM→BAM conversion
  first — at BAM speed the cohort is already ~12 min single-threaded, inside target, and no core
  count rescues CRAM's ~21 h.
- **Complexity cost:** the recommended shape (N single-sample processes) is **zero code** — a
  shell loop plus per-sample output files — at the price of repeating the 27 s walk per process
  (~4.3 min of wall at 8-way). The threaded shape (collect the 48,426 segments once, fan out over
  samples with a per-worker generator) recovers that but introduces two hazards flagged
  pre-emptively as L16 and L17. No channel: every sample needs every segment, so producer/consumer
  is the wrong topology despite `TypedRegionIterator` being built for one.

#### H4: src/ng/locus_generation/ssr.rs:1591-1597 — the reference margin fetch is sample-invariant work redone per sample

- **Confidence:** High. **Four categories flagged this independently** (hot_loops as a finding;
  concurrency, data_layout and allocations as cross-category notes), which is why it is promoted
  here.
- **Hot-path evidence:** `ssr.rs:1591` (`SsrLocus::fetch`) = 134 of 2162 samples (**6.2%**) at one
  sample. In the 51-sample target, 50 of every 51 of those fetches recompute identical bytes.
- **Pattern matched:** loop-invariant work inside the loop.
- **Mechanism:** in the driver's `for segment { for sample { … } }` nesting, `next_locus` re-runs
  `SsrLocus::fetch` — the tract ± flank reference read, the margin copy into a `Box<[u8]>`, and
  `segment.clone()` — for every sample, though none of it depends on the sample. It also drags in
  the per-locus reference-reader construction the walk otherwise avoids.
- **Measurement plan:** add a one-entry cache keyed on `(contig, segment)` and re-run the 3-sample
  ch01 experiment against the 24.73 s baseline. **Threshold:** ≥ 4% off the per-locus stage at
  N = 3, scaling toward ~6% at N = 51; tallies byte-identical. This caps at the 6.2% the profile
  measured — do not oversell it, and do not do it before H1.
- **Complexity cost:** one `Option<(SsrSegment, SsrLocus)>` field on the generator plus a
  cache-validity rule. The rule is the risk: a stale hit returns another segment's reference
  bases, which is a silently wrong measurement, not a crash — so the key must include the contig
  and the test must be a real cache-miss test, not just a hit test.

### Likely

- **L1: [src/ng/region_typing/mod.rs:1288-1329](../../../../src/ng/region_typing/mod.rs#L1288-L1329) — `emit_into` linear-scans every requested BED span for every resolved region.** Confidence Medium.
  Evidence: the BED walk is *slower* than the whole-contig walk (2.88 s vs 2.44 s) while typing 37×
  fewer segments; ≈8.8×10⁸ overlap tests is the right order for the ~0.4 s delta. The walk stage
  was never profiled, so **profile it first and close this finding if `emit_into` is absent**
  (plan entry 6). Fix: a monotone `usize` cursor on `SpanWalk` beside `emitted_upto` — spans and
  regions are both in genomic order. Complexity: one cursor field, plus the invariant that regions
  arrive monotonically (they do, and `emitted_upto` already relies on it).
- **L2: [src/ng/locus_generation/ssr.rs:271](../../../../src/ng/locus_generation/ssr.rs#L271) — `Reservoir::new` pre-sizes 164 KiB per locus visit regardless of depth.** Confidence High on the
  allocation, Low on the wall-clock win. **Two categories found this independently** (allocations,
  data_layout), and both measured `size_of::<MappedRead>() = 168 B` with a `rustc -O` probe rather
  than guessing: `Vec::with_capacity(1000)` × 168 B ≈ 164 KiB, allocated and freed once per
  (locus × sample) — ~10.9 M times on the intended workload, ~97% of them for loci with no reads
  at all (mean depth on covered loci is 6.8). It is also the most plausible source of the profile's
  `__bzero` 26 / `_platform_memset` 15, since data_layout confirmed the DP scratch is *not*
  re-zeroed. Fix: `Vec::new()` — Algorithm R branches on `held.len() < capacity`, not on the
  allocation, so the kept set stays byte-identical. Measurement: DHAT via the repo's
  `--features dhat-heap` pattern (expect the large-block class to vanish) plus the ch01 wall.
  Complexity: one line for the `Vec::new()` variant; the generator-owned-buffer variant adds a
  field and a clear-on-entry rule.
- **L3: [src/ng/alignment/ssr_unit_robust.rs:557-583](../../../../src/ng/alignment/ssr_unit_robust.rs#L557-L583) (with `:400`) — the whole-unit slip emission is recomputed per tract *column*.** Confidence
  High; **codegen-backed**, and found independently by hot_loops and data_layout. The `cargo asm`
  listing at the caller the profile names (`cargo asm --example ng_ssr_cohort_stutter --simplify
  "classify::delimit"`, 1681 lines aarch64) shows **six `udiv`+`msub` pairs, six `ldar`
  load-acquires of the `PER_QUALITY_LN` `LazyLock`, and six bounds checks per tract cell**. The sum
  depends only on (read row, motif phase), so it can be resolved `period` times per row instead of
  per column, removing the runtime `%`, the acquire-load and the quality bounds check from the
  cell. Measurement: criterion A/B in a new `benches/ng_ssr_delimiter_perf.rs` (criterion 0.8.2 is
  already a dev-dep) plus a re-`cargo asm` to confirm the `udiv`s are gone; threshold ≥ 10% on the
  delimiter bench. Complexity: one `[f64; 6]` per row (period is capped at 6). **Parity:** must
  keep the k-ascending summation order or scores move in their low bits.
- **L4: [src/ng/alignment/ssr_unit_robust.rs:368-396](../../../../src/ng/alignment/ssr_unit_robust.rs#L368-L396) (used at `:512`, `:517`, `:536`) — the junction guard and tract-aware gap-open are recomputed
  per cell though they depend only on the column.** Confidence High, codegen-backed: ~15
  instructions of `ccmp`/`cset`/`fcsel` per cell with two locus constants reloaded from the stack.
  Fix: hoist to a per-column array in scratch; values bit-identical. Complexity: three
  grow-and-keep `Vec<f64>` fields on `UnitRobustScratch`.
- **L5: [src/ng/alignment/ssr_unit_robust.rs:422](../../../../src/ng/alignment/ssr_unit_robust.rs#L422) and `:588` — `(column - period + 1..=column).all(column_in_tract)` compiles to an unrolled
  six-step branch chain per cell** and is provably equivalent to two comparisons, since
  `column_in_tract`'s true-set is an interval. Confidence High, codegen-backed. Complexity: near
  zero — one named helper. The cheapest of the DP findings; do it first.
- **L6: [src/ng/alignment/ssr_unit_robust.rs:237-240](../../../../src/ng/alignment/ssr_unit_robust.rs#L237-L240) — the DP ring is `Vec<Vec<[f64; STATES]>>`, so every predecessor read pays an outer index,
  a ptr/len reload and a bounds check.** Confidence Medium; both hot_loops (31
  `panic_bounds_check` blocks in the listing) and data_layout (two dependent loads per access,
  ~14 accesses per cell, ring rows separately heap-allocated so the ring is not contiguous)
  filed it. Fix: flatten to one `Vec<[f64; STATES]>` indexed `slot * stride + column`, or build the
  current row in a separate buffer and swap. Complexity: real — index arithmetic replaces double
  indexing, so a `stride`/`ring_len` mix-up becomes a silent wrong answer instead of a panic
  (mitigate with an `#[inline] fn cell(slot, column)`), and the `period == 1` aliasing case is a
  trap. The most invasive and least certain of L3–L7; do it last, if at all.
- **L7: [src/ng/alignment/ssr_unit_robust.rs:239](../../../../src/ng/alignment/ssr_unit_robust.rs#L239) — the backpointer matrix streams ~99 KB of stores per read (151×131×5 B) to serve ~300 traceback
  reads.** Confidence Medium. Five states × 3 bits fit a `u16` → 2.5× less store traffic.
  Complexity: a `pack`/`unpack_one` pair plus a new invariant (the 3-bit field width depends on
  the state count, so a sixth state silently corrupts it — assert it).
- **L8: [src/ng/read/input/open_bam.rs:579](../../../../src/ng/read/input/open_bam.rs#L579) — the BAM is opened as `bgzf::io::Reader<File>` with no `BufReader`,** so every BGZF frame costs two
  `read` syscalls (an 18-byte header read, then the body). Confidence Medium. Bounded by `read` 160
  + `__lseek` 25 + `__open` 17 = 9.3% of samples — the ceiling for all syscall-count findings
  combined, so it cannot be sold as a headline win. Complexity: a type change threaded through
  three declarations. Largely subsumed by H1, which removes most of the frames.
- **L9: [src/ng/raw_chrom_reader.rs:182](../../../../src/ng/raw_chrom_reader.rs#L182) — `read_raw_bases` always reads a 64 KiB chunk and discards the tail,** so a ~150-base reference
  fetch costs an `lseek` plus a 64 KiB copy. Confidence Medium. Complexity: one arithmetic helper
  and a `min` in the read loop.
- **L10: [src/ng/raw_chrom_reader.rs:288](../../../../src/ng/raw_chrom_reader.rs#L288) — a forward jump past the buffered window reads and buffers every intervening base instead of
  repositioning.** Confidence Medium. **Gated:** on today's code the walk asks for every window of
  the contig anyway (see the Note on `scan_set`), so this cannot be the walk's cost — it becomes
  material only *after* scan spans are narrowed, and it is a prerequisite for reusing one reader
  across scattered loci. Complexity: one branch plus a threshold constant.
- **L11: [examples/ng_ssr_cohort_stutter.rs:202](../../../../examples/ng_ssr_cohort_stutter.rs#L202) + [src/ng/ref_seq.rs:524](../../../../src/ng/ref_seq.rs#L524) — the owner's lead, judged: the
  per-locus `ContigList` + `PathBuf` clone is *not* material; what is real is the first fetch.**
  Confidence Medium. `WindowedRefSeq::new` does not appear in the profile at all and is genuinely
  cheap (it stores a path and sets `current: None`). But the reader it produces re-reads and
  linear-searches the **entire `.fai`** and re-`File::open`s the FASTA on first fetch
  (`open_contig`, [raw_chrom_reader.rs:111-141](../../../../src/ng/raw_chrom_reader.rs#L111-L141)),
  once per locus that has reads — ~6,489 on ch01, ~660k per cohort walk. Bounded by `__open` 17 +
  `__lseek` 25 + the ~1% in `fetch_raw_into`, so **~1–2% at most on tomato**. The reason it is
  filed at all: the cost is proportional to contig count, and the sibling `ssr_hg002` reference has
  **2,580 contigs**, where the same path would allocate ~2,581 blocks and search 2,580 entries per
  locus per sample. Smallest fix that keeps `FnMut() -> R` intact and needs no library change: hand
  out `ResidentRefSeq` over one shared `Arc`-backed `fasta::Repository`. Measurement: a
  contig-count sweep, not another profile. Complexity: none in the library for the
  `ResidentRefSeq` variant; the `Arc<ContigList>` variant is the `Arc` seam the generator's own doc
  already names at `ssr.rs:1331`.
- **L12: [src/ng/locus_generation/ssr.rs:1004](../../../../src/ng/locus_generation/ssr.rs#L1004) — a fresh `HashMap` per locus for a handful of alleles.** Confidence High on the allocation, Low on
  wall clock (the tally is 1 of 2162 samples). Hoist onto the generator + `clear()`. Complexity:
  `tally` becomes a method or takes the map by `&mut`.
- **L13: [src/ng/locus_generation/ssr.rs:1624](../../../../src/ng/locus_generation/ssr.rs#L1624) — the `outcomes` `Vec` is collected only to be zipped once.** Confidence High on the allocation,
  Low on wall clock. Honest cost: fusing classify into the tally cuts against the project's
  "split data-shaping from math" preference — probably **won't fix**, recorded for completeness.
- **L14: [examples/ng_ssr_cohort_stutter.rs:97](../../../../examples/ng_ssr_cohort_stutter.rs#L97) — `motif.as_bytes().to_vec()` per emitted row, feeding a `from_utf8_lossy` that takes `&[u8]`.**
  Confidence High. Pure waste, two lines, no complexity cost. Take it while you are in the file.
- **L15: [src/ng/read/input/mod.rs:449](../../../../src/ng/read/input/mod.rs#L449) — the k=1 arm builds and frees a one-slot `Vec` of streams per query,** contradicting the design note 60 lines
  below it that keeps `RegionReads` unboxed for exactly this reason. Confidence High on the
  allocation, Low on wall clock. Complexity: ~8 duplicated lines.
- **L16 (would be introduced by H3): [examples/ng_ssr_cohort_stutter.rs:208-209](../../../../examples/ng_ssr_cohort_stutter.rs#L208-L209) — a shared `BufWriter<StdoutLock>` becomes one lock per output row** across millions of rows
  the moment samples run in parallel. Fix before threading: per-sample output files. Complexity: an
  output directory instead of stdout — a driver UX change, and the reason process fan-out wants
  per-sample files anyway.
- **L17 (would be introduced by H3): [examples/ng_ssr_cohort_stutter.rs:172-173](../../../../examples/ng_ssr_cohort_stutter.rs#L172-L173), `:294-296` — `read_reference_verifying_or_creating_fai` spawns a CPU-bound whole-genome MD5
  thread per process,** so an N-process fan-out runs N concurrent 795 MB digests against the
  workers — the degenerate two-pool oversubscription case. Measure before deciding; the fix is a
  flag or a documented decision to skip a correctness check in an exploratory tool, which is a
  judgement call, not a diff.

### Speculative

Filed so the shapes exist; **do not act without an experiment that contradicts "this won't
matter"**. Details in the per-category files.

- `ssr.rs:756`/`:791` — a `Box<[u8]>` per read that `entry()` usually drops as a duplicate key;
  stable Rust has no `raw_entry`, so this is the finding most likely to cost more clarity than it
  buys.
- `ssr.rs:88` — `SsrLocus` copies the margin out of the buffer it was just fetched into (a
  lifetime on a public type is the price).
- `segment_criteria.rs:143` / `ssr.rs:1591` — `SsrSegment`'s owned contig name is cloned per locus
  per sample; `Arc<str>` would change a `pub` field **outside** this review's scope.
- `ssr.rs:1650`, `:1656`, `:1657` — three `Box<[u8]>` of reference bases per locus visit, built
  before anyone knows the locus has reads; explicitly **refuted as a time cost** by the profile's
  1-sample attribution and recorded only in case that changes.
- `ssr_unit_robust.rs:268` — `best_of`'s AoS `(f64, State)` (16 B, 7 padding); asm-first, do not
  bench blind.
- `ssr_unit_robust.rs:238` — a state-plane (SoA) ring; honestly mixed, since one of four per-cell
  reads uses all five lanes.
- `ssr.rs:1004` — replacing the tally `HashMap` with a linear scan (O(rows²) at a pathological
  locus — keep the map above a threshold).
- `locus_generation/mod.rs:136` — `ObservedSequence` is 72 B (crosses a cache line) and 24 of them
  are a `chain_ids` the STR path never fills; the generic path *does* fill it, so this is a shared
  type change.
- `ssr.rs:769-772` — `ln_p_err_sum`'s scalar `f64` chain where an integer sum plus one multiply
  would do; spends the "q_sum is soft" licence.
- `input/mod.rs:449` — short-circuit a query for a segment with no reads (~97% of queries) before
  borrowing a pooled reader; largely subsumed by H1.
- `open_bam.rs:110-125`, `:484-489` — the CRAM `fasta::Repository` is per sample, which is *why*
  there is no cross-sample contention; sharing one to save RAM would convert the memory win into a
  serialization point, because noodles' `Repository::get` holds the **write** lock across the
  whole-contig disk read. Measure per-worker RSS; prefer bounded residency over sharing.

### Note

- **Do not "fix" these** — three sites were checked and found already right, and the reasoning is
  worth keeping: `best_of` is fully scalarised by SROA (15 `fcmp` / 17 `fcsel` in the whole
  function) so a rewrite buys nothing; the `Mutex` reader pool guards only a `Vec` pop/push with
  the guard deliberately dropped before `open_reader`, and no lock or futex frame appears anywhere
  in either profile (no `parking_lot`, no sharding, and `Relaxed` on `readers_opened` is correct);
  and the driver's output side (1 MiB `BufWriter` over a held `StdoutLock`, streamed not
  accumulated) is already correct.
- **The DP scratch is not re-zeroed per read** (`ssr_unit_robust.rs:248` is grow-only and every
  cell of the current window is written before it is read), so the profile's `__bzero` /
  `memset` samples do not come from it — which is what points at L2 instead.
- **`scan_set` types every contig end to end, by design.**
  [mod.rs:1213-1237](../../../../src/ng/region_typing/mod.rs#L1213-L1237) sets every `ScanSpan`'s
  scan range to `Position(1)..Position(contig_len)`, and `emit_into` filters the output. That is
  why `--regions` **cannot** make the walk cheaper: a single 1,012 bp target still costs 3.28 s on
  ch01, and measured walk time is 35 ms per Mb of *contig* length (1.64 s / 2.33 s / 3.28 s for
  47.3 / 66.7 / 90.9 Mb contigs; 26.66 s for the whole-genome BED against a 27.3 s prediction).
  The whole-region-comes-back-whole rule is an owner decision (2026-07-17) and a satellite can be
  megabases, so this review does not reopen it as a perf question — it quantifies the price
  (~27 s per genome pass, one-time, ~2% of a 51-sample BAM cohort) and files L1 as the part that
  *is* actionable.
- `delimit_parity.rs` pins **algorithm 3, not 4r**, so it is not the parity oracle for L3–L7; the
  real gate for those is a `diff` of the driver's TSV before and after.
- `read_footprint` / `ref_to_read` (`ssr.rs:388-435`) are fine as written; `PER_QUALITY_LN`'s
  `LazyLock` acquire-load is only a problem inside the cell (L3), not in principle.
- The one allocation site the profile actually names — noodles' per-record tag decode
  (`read_string` → `malloc_tiny` 13, `Data::clear` → `_xzm_free` 28) — is **out of scope**; the
  in-scope lever is to stop re-decoding the same records (H1), not to touch noodles. Relatedly,
  `read_record_buf` fully decodes aux data this path never reads (3.1% self time).

## 6. Out-of-scope observations

- **The driver's doc comment records the mis-measurement as a design rationale.**
  [examples/ng_ssr_cohort_stutter.rs:13-17](../../../../examples/ng_ssr_cohort_stutter.rs#L13-L17)
  states that walking untargeted tomato ch01 is "~8 minutes per sample" of waste and makes that the
  case for `--regions`. The walk measures **2.44 s**; the ~8 minutes was the `cargo run` harness.
  The conclusion (`--regions` is worth having) survives — it cuts the per-locus stage 8.04 s →
  1.36 s — but for a different reason: fewer *emitted* segments means fewer read queries, not
  cheaper typing. Fix the comment before it becomes lore.
- **The same comment's claim that restriction "cannot change what a covered locus is" is false.**
  [:216-219](../../../../examples/ng_ssr_cohort_stutter.rs#L216-L219). Measured: whole ch01 finds
  **6,489** covered loci, the ch01 BED finds **5,534** — 955 loci have reads but lie outside the
  BED. The two modes are therefore not comparable workloads, and a cohort run that switches to
  `--regions` will produce a different (smaller) locus set. That is a correctness-of-analysis
  point, not a perf one, and it deserves a decision rather than a silent change.
- **`--regions` was documented but unparsed.** `run_cohort` took a `regions_bed: Option<&Path>`
  while `main` never parsed the flag (the example did not compile at review start). I added the
  five-line parse to unblock measurement — see *What I touched* below.
- **`src/ng/read/input/test_fixtures.rs` is `#[cfg(test)] pub(crate)`** (`mod.rs:28-29`), so a
  `benches/` file cannot reach the synthetic indexed-BAM fixture builder it needs — and that
  builder already exists twice in-tree (`examples/dhat_ng_merge.rs:306-322`). Lifting it behind a
  feature or into a `pub` dev-support module is the prerequisite for M1.
- **`BinningIndex::query` allocates a `Vec<Chunk>` per region query** inside noodles
  (`region_query.rs:88`, `:134`, `:152`) — one per locus per sample. Nothing to do while H1 stands;
  H1 removes most of the queries.

## 7. What's already good

- **Generator-owned scratch, applied consistently.** `margin_buffer` / `qual_buffer` on
  `SsrGenerator` ([ssr.rs:1361-1363](../../../../src/ng/locus_generation/ssr.rs#L1361-L1363)) and
  `UnitRobustScratch`'s grow-and-keep `resize`
  ([ssr_unit_robust.rs:248](../../../../src/ng/alignment/ssr_unit_robust.rs#L248)) mean the DP
  performs **zero allocations per `align`** — the reviewer read the whole DP and traceback to
  confirm it. L2 is notable precisely because it is the one place the pattern was not applied.
- **Static dispatch where it counts, with the reasoning written down.** The `RepeatDelimiter`
  trait-alias keeps the aligner a type parameter so `align` is a direct call in the per-read loop
  ([ssr.rs:1304-1309](../../../../src/ng/locus_generation/ssr.rs#L1304-L1309)), and
  `SampleRegionReads` is an enum rather than `Box<dyn Iterator>` with the inlining argument spelled
  out ([mod.rs:494-514](../../../../src/ng/read/input/mod.rs#L494-L514)). The `cargo asm` listing
  confirms it: the aligner is fully inlined into `classify::delimit`.
- **The CRAM index work that *was* done is the right shape** — `.crai` entries grouped by contig
  once at open, a container-level early stop, and a documented rejection of production's
  rescan-from-entry-0 ([region_query.rs:274-281](../../../../src/ng/read/input/region_query.rs#L274-L281)).
  H2 is the same lesson one level down: the grouping was hoisted out of the query, the *decode* was
  not.

## What I touched

Two changes outside the report, both flagged rather than assumed:

1. **`examples/ng_ssr_cohort_stutter.rs`** — added `--regions <path>` parsing to `main` and passed
   it to `run_cohort` (five lines). The flag was already documented and already a parameter; the
   example did not compile without it. Nothing else in the driver was changed.
2. **`examples/ng_perf_probe.rs`** — the measurement harness described in section 3. Written for
   this review and **deleted at the owner's request once the experiments were done**, so M2 is
   closed as *won't fix*. Its mode/dump design is recorded in section 3 for M1's bench to pick up;
   until that bench exists, this path has no repeatable measurement in the tree.

## Author response convention

Address each finding by its identifier (H1, L2, M5, …) with one of: `applied in <commit>` /
`experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` /
`won't fix because …`. The "experiment shows no gain" path is expected and welcome — that is what
the measurement plan is for. H1, H2 and H3 in particular are gated on experiments that could
still refute them.

---

# Applied — same-day experiments (2026-07-26)

Four findings were taken to a measurement. **One landed a 129× win; two were refuted by their own
experiments and reverted; one was refuted in the form the review proposed.** Every change was
gated on a byte-identity dump of all evidence per locus (`--dump`: one row per observation with
bases, coverage tag, `num_obs`, `num_fwd`, `q_sum`, `mapq_sum`, `mapq_sum_sq`, plus the
no-observation tallies), diffed by md5 against a pre-change baseline.

Baselines (host, native release, warm cache), all with identical tallies
(`covered=6489 observations=20639 reads_fetched=43991` for W1;
`covered=5534 observations=18827 reads_fetched=40337` for W2/W3):

| workload | baseline |
|---|---|
| W1 — BAM, whole SL4.0ch01, 1 sample | 8.36 / 8.38 / 8.45 / 8.50 s |
| W2 — CRAM, ch01 target BED, 1 sample | 171.46 / 174.22 s |
| W3 — BAM, ch01 target BED, 1 sample | 1.36 / 1.38 / 1.44 s |
| W1 × 3 samples | 23.61 / 23.70 / 23.78 / 23.93 / 24.00 s |

## H2 — **applied.** CRAM: 174.2 s → 1.35 s (129×), output byte-identical

The decoded container now travels with the pooled reader
([region_query.rs](../../../../src/ng/read/input/region_query.rs), `DecodedContainer`;
[open_bam.rs](../../../../src/ng/read/input/open_bam.rs), `ReaderHandle::container`), keyed on the
`.crai` offset. That makes it **per worker** — each concurrent caller holds its own reader and so
its own container — and adds no lock to the query path, which is the shape `arch/alignment_file.md`
§7 anticipated. `last_decoded_offset` is kept and now has a documented, distinct job: it dedups
*within* a query (a multi-slice container appears once per slice), while the container cache stops
a re-*decode* *across* queries.

| | baseline | applied |
|---|---|---|
| W2 per-locus stage | 174.22 s | **1.35 s** |
| W2 peak RSS | 227.9 MB | 262.2 MB (+34 MB) |
| W3 (BAM, unaffected) | 1.36–1.44 s | 1.38 s |
| W1 (BAM, unaffected) | 8.36–8.50 s | 8.42 s |
| dump md5 | `5388550a…` | `5388550a…` |

CRAM is now the same speed as BAM for identical output, and the BAM path is untouched (the cache
is `None` for a BAM). Whole-genome target BED, one sample: per-locus **11.62 s** against a ~1485 s
projection — so the cohort projection moves from **~21 h to ~10 min** (23.9 s walk once +
11.6 s × 51). Verified on the production arch too (container: 2.34 s for W2).

**Safety.** The failure mode of a cache like this is silent — a stale hit reports another region's
reads as this region's depth — so `a_stale_container_is_not_reused` drives three queries through
one pooled reader (container A → container B → back to A) and requires each answer to equal what a
cold reader returns for that region alone. **Mutation-verified**: replacing the offset comparison
with `self.container.is_none()` makes it fail on the second query (and also trips an existing
multi-container test). `readers_opened()` widened to `pub(super)` so the test can assert all three
queries shared one reader, hence one cache. Full suite 2431 passed / 0 failed (was 2430; the 9
release-profile `should_panic`-on-`debug_assert` failures are the documented pre-existing set,
none in `read/input`). fmt and clippy clean.

## H1 — **refuted in the proposed form.** Grouping helps only where the review assumed it would help least

Prototyped as the review's variant (a) — one indexed query per 32 kb bucket of adjacent segments,
each segment served by overlap from the bucket's reads — in the probe rather than the library, so
the ceiling could be measured before committing to a generator API. Fetch stage only, tallies
identical (`reads_fetched` matches in every row):

| workload | per-segment query | 32 kb bucketed |
|---|---|---|
| BAM, whole ch01 (untargeted) | 7.02 s | **32.89 s — 4.7× slower** |
| BAM, ch01 target BED | 0.74 s | **0.33 s — 2.2× faster** |
| CRAM, ch01 target BED (post-H2) | 0.70 s | 0.58 s |

Two things fall out. **The 0.33 s is exactly the `samtools view -c` floor**, so the review's
headroom claim was right about the ceiling. But the naive grouping *regresses the untargeted walk
4.7×*, because a 32 kb bucket holds ~75 segments and ~6,400 reads and the prototype rescans the
window per segment — an O(segments × reads) inner loop that swamps the inflate it saves. A real
implementation needs a position cursor, not a rescan, and the untargeted path needs the bucket to
shrink where segments are dense.

What is left after that is a **~1.4× end-to-end** win on the intended (BED-restricted) workload —
the fetch stage is 0.74 s of a 1.38 s per-locus stage. That is not worth a new `LocusGenerator`
entry point, a new invariant (segments must arrive in coordinate order), and a bucket-sizing policy
that has to be right in both regimes. **Deferred**, with the ceiling now measured rather than
assumed: revisit if the cohort grows an order of magnitude, and implement with a cursor if so.

## H4 — **experiment shows no gain, reverted.** The sample-invariant fetch is not recoverable

Implemented as a one-entry `Option<(ContigId, SsrLocus)>` on the generator, keyed on contig +
segment, with the two discriminating tests the review asked for (miss-on-different-segment,
hit-equals-cold-fetch) — the miss test **mutation-verified** to fail when the segment is dropped
from the key. Interleaved A/B, two binaries, alternating runs:

| | no cache | locus cache |
|---|---|---|
| W1 × 3 samples | 23.93 / 23.78 / 24.00 s | 23.67 / 23.23 / 23.73 s |
| W1 × 8 samples | 64.04 / 62.57 s | 65.01 / 61.68 s |

At N=3 that is ~1.5%; at N=8, where the hit rate is 87.5% and the win should be near its maximum,
it is **zero**. So the 6.2% the profile attributed to `ssr.rs:1591` is not the reference read —
skipping it saves an allocation and a copy that the allocator was already serving from a hot free
list, and the cache's compare-and-restore costs about the same. Reverted: a cache whose stale hit
is a silently wrong measurement is not worth a correctness-critical invariant for no gain. (First
run showed +6%, which turned out to be thermal drift from back-to-back heavy runs — the interleaved
A/B is the trustworthy form, and it is why the earlier number is not quoted as the result.)

## L2 — **experiment shows no gain, reverted.** The 164 KiB reservoir pre-size costs nothing in time

`Vec::with_capacity(1000)` of `MappedRead` → `Vec::new()`, one line:

| | baseline | `Vec::new()` |
|---|---|---|
| W1 (host) | 8.36 / 8.38 / 8.45 s | 8.36 / 8.38 s |
| W3 (host) | 1.38 / 1.44 s | 1.36 / 1.39 s |
| W1 (container, production arch) | 9.80–10.50 s | 9.99 / 9.86 s |

Checked in the container as well as on the host, in case a 164 KiB block behaved differently under
glibc than under macOS's allocator. It does not: the allocation is transient per locus and both
allocators serve it from cache. Reverted — `with_capacity` reads as more intentional for a
reservoir than a `Vec::new()` that needs a paragraph explaining itself, and the finding's own
confidence on wall clock was Low.

## What this changes about the review's conclusions

- **H2 was the whole cohort problem, and it is fixed in code** rather than by converting the CRAMs.
  Plan entry 7 (CRAM→BAM conversion) is no longer needed; keep it only as a fallback.
- **H1's ranking was right and its fix was wrong.** The re-inflate is real and the floor is real,
  but the contained version of the fix does not reach it and the reaching version costs an API.
- **Two of the three cheap "Likely" wins were noise**, which is the expected hit rate for
  allocation-count findings on a path whose profile is 71% decompression — and is why the review
  filed them at Low confidence on wall clock rather than as wins.
- The measurement harness earned its keep three times over — it produced the baselines, the
  identity gate, and the H1 refutation — and was then **deleted at the owner's request** (M2 closed
  as won't fix). Its mode and `--dump` design is written down in section 3 precisely because the
  binary is gone; **M1's bench is now the only way to re-measure this path**, and until it exists
  any future timing here starts from scratch.
