# Performance Review: ng-generic-pileup
**Date:** 2026-07-31
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** ng's generic (SNP/indel) locus generator — the pileup walk, the read/reference path it drives, and the typed-region stream that feeds it
**Verdict:** Apply the listed wins
**Hot-path evidence:** CPU sampling profile, DHAT heap profile, a region-grain sweep, a depth sweep, and per-finding A/B measurements — all on real HG002 data

---

## 1. Scope and constraints

**What was reviewed.** The module `src/ng/locus_generation/pileup/` and everything the walk
pulls: the read ingestion path (`src/ng/read/input/`), the reference accessors
(`src/ng/ref_seq.rs`, `src/ng/raw_chrom_reader.rs`, `src/ng/reference_info.rs`), the
typed-region stream (`src/ng/region_typing/`, `src/ng/tandem_repeat.rs`), and the
generator's own bench (`benches/ng_generic_pileup_perf.rs`).

**Reviewed against.** Commit `d95ce8b` on branch `ng-generic-perf` (worktree
`/Users/jose/devel/pop_var_caller-ng-generic-perf`). Every sub-agent detached its own worktree
at that SHA and confirmed the branch-only files were present before measuring.

**Targets and hardware.** The owner named the reference workload: **human WGS, one sample,
HG002 30×**. Host is an Apple M5 Pro (18 cores, 64 GB, macOS 26.5.2); all builds native on the
host, `--release` (fat LTO, codegen-units = 1, panic = abort). The generator is
single-threaded. Both wall time and memory were in scope; §2 explains why memory stopped being
the interesting half.

**Baseline, measured (this review's own instrument).**

| run | loci | wall | peak RSS |
|---|---:|---:|---:|
| chr1 | 1,541,788 | 33.2 s | 23.9 MB |
| chr2 | 1,536,964 | 41.4 s | 22.4 MB |
| chr21 | 236,081 | 6.8 s | 20.0 MB |
| **whole genome** | **18,524,066** | **453 s walk / 552 s wall** | **30.1 MB** |

**⚠ The fixture is not the target workload, and this governs the whole report.**
`HG002_TR_v1.0.1_Tier_30x.bam` is **tandem-repeat-targeted**: 30× depth *inside the TR
benchmark regions only*. The probe's own counters prove it — 1,541,788 loci over 240,227,974
generic bp on chr1 is **0.64 % of positions covered**. No true WGS alignment exists on this
host (the GIAB per-sample BAMs are 8–12 MB region-selected subsets). Consequence: **absolute
per-region costs in this report are properties of this fixture's depth, not of the generator.**
Three independent measurements say so — the fixture gives 27 µs/region on chr1 and 40.6 µs on
chr21, a 300× run of the same regions gives 97.2 µs, and the synthetic 30 %-covered bench gives
563 µs. **Ratios transfer; microseconds do not.** Closing this gap needs a real 30× WGS BAM and
is the first item in §3.

**Hot-path evidence available.** A `sample` CPU profile of the full chr1 walk (25,446
main-thread samples), a DHAT heap profile of a 300,000-locus prefix built on `[profile.profiling]`,
a five-point region-grain sweep, a 30×-vs-300× depth sweep, an allocator A/B, and per-finding
A/B measurements by each sub-agent in its own worktree. Raw artefacts and the full baseline
write-up are the audit trail in `tmp/perf_review_2026-07-31_ng-generic-pileup/`.

**Deliberately out of scope.**

- `parity.rs`, `tests.rs`, `copy_fidelity.rs`, `mock_reference.rs` — test-only.
- `src/pileup/`, `src/var_calling/`, `src/ssr/` — frozen production; ng must not edit it.
- `cigar_cursor.rs`, `decompose.rs`, `chain_id_allocator.rs` — byte-identical copies of
  production, enforced at compile time by `copy_fidelity.rs`. Editing them breaks the build, so
  findings against them cap at **Likely** and name the owner decision required.
- Narrowing `max_record_span` (the halo) — measured and **refuted**, see §5 N1.
- Designing the parallel fan-out — not asked for; only the taxes it would pay were catalogued.

**Categories dispatched.** All six. `methodology` (always), `allocations` (46.6 allocations per
locus), `io_and_syscalls` (43 % of self time in BGZF inflate), `hot_loops` (the walk's
arithmetic and the never-reviewed typed-region scan), `data_layout` (the types a cohort fan-out
would hold N copies of), `concurrency` (a single-threaded walk on an 18-core box, plus a
background verification thread).

**Method note.** The first fan-out was cancelled after producing nothing; the second ran to
completion. Sub-agents measured concurrently on one host, so **each agent's absolute seconds
are comparable only within its own table.** Every number this report *recommends acting on* was
re-measured by the orchestrator on a quiet machine with alternating binaries — and that
re-measurement changed one of them (§5 H8).

---

## 2. Verdict

**Apply the listed wins.** A combined patch of five independent changes is measured, correct
and ready:

| | chr1 | chr2 |
|---|---:|---:|
| baseline | 34.106 s | 34.917 s |
| **combined patch** | **29.604 s** | **30.315 s** |
| | **−13.2 %** | **−13.2 %** |

Alternating binaries, three rounds, quiet host; `loci=1541788 observations=1647161` identical
on every run; peak RSS 26.5 → 25.9 MB. **`cargo test --lib ng::` on the debug profile: 1,471
passed, 0 failed.** Patches are preserved in
`tmp/perf_review_2026-07-31_ng-generic-pileup/patches/`.

Beyond that, the largest lever in the review is **not** in the patch because nobody has built
it: the per-region BAM query re-decodes the same records **30.3×** (§5 H3/H4). Its measured
ceiling is a further **~3.4×** on the region-dominated part of the run.

**Two results reframe what "memory work" means here, and both are negatives.**

*Memory is already solved.* Milestone E's unattributed +68 % peak-RSS growth was **not the
generator's**: the identical walk retaining nothing peaks at 23.9 MB against the 859 MB E
measured on `ng_generic_loci_dump`, whose peak is its own whole-run `Vec<ObservationRow>`. The
whole genome runs in **30.1 MB**, linear in loci. Spec §7's "bounded by depth, not by region
length" was then tested directly: two DHAT runs at 25× different region grain give a
**byte-identical live high-water of 14,250,994 bytes in 74,722 blocks** while cumulative
allocation differs 10.5 %. The 68.6 MB seen at 400 bp grain is allocator high-water over churn,
not retention.

*Allocation count is not the currency.* Five experiments, three agents:

| removed | measured |
|---|---:|
| all `malloc` cost (mimalloc) | **+2.4 %** wall, **+45 %** RSS |
| **36.4 % of every allocation** (noodles lazy `Record`) | **−0.5 %**, noise |
| 13 % of all **bytes** (`folded_reads` free-list) | **0 %** |
| 41 % of all **bytes** (`tandem_repeat` stack reuse) | **≈0 %**, variant lost several pairs |
| the same allocation **plus** its hashing and two sorts | **−3.9 %** |
| a map **plus** its hashing (mate overlap) | **−1.2 %** |

**Removing an allocation pays only when it also removes the hashing, sorting or grouping the
structure existed to support.** Taking the mimalloc ceiling first is what makes this readable:
without it, "36.4 % of all allocations are noodles aux-tag `BString`s" would have been this
review's headline recommendation, and it is worth −0.5 %.

---

## 3. Measurement plan

Ordered by what unblocks what.

1. **Get a real 30× WGS BAM for HG002 and re-run the baseline.** Everything in §1's table is
   0.64 %-covered data. Command is the probe as given below with the new BAM. **Threshold:**
   if the per-region constant lands near the bench's 563 µs rather than the fixture's 27 µs,
   H3/H4's ceiling grows and they become the only findings that matter; if it lands near
   27 µs, the §2 patch is proportionally more of the answer. Nothing else in this report needs
   re-deciding either way, because every recommendation is a ratio.
2. **Extend the bench's grain sweep to 400 bp and 100 bp** (two array entries at
   `benches/ng_generic_pileup_perf.rs:420`). Its own numbers show the curve still climbing
   where it stops: `T ≈ 296 ms + 563 µs × regions`, with the two coarsest points statistically
   indistinguishable. **Threshold:** the added points are the only ones with ±0.4 % CIs.
3. **Raise `sample_size` off 10** (`:395-396, :416-417`). A natural revert experiment found the
   three *unchanged* grain points swinging 12 % across two runs of identical code, with one
   `p = 0.00` "significant" change. At `sample_size(10)` the bench cannot gate any finding in
   §5.
4. **Add a bench for the typed-region stream.** It is 18 % of main-thread self time, 23 % of
   the reference run and 41 % of all bytes allocated, and has **no benchmark at all**.
   `TypedRegionIterator::over_contig` is `pub` and needs no BAM.
5. **Build H3** (the BGZF re-decode) and measure against the ceiling already established: at
   10 kb query windows, 4,671 regions instead of 116,775, records read 975,142 → 72,727.
   **Threshold:** merge if it recovers ≥ 50 % of the 3.4× at unchanged output.
6. **Commit the probe** (§4) so the next question does not start from a throwaway again.

**The reproducible commands.**

```
# the review's instrument (new, uncommitted — see §4)
/usr/bin/time -l ./target/release/examples/ng_generic_walk_probe \
  ~/genomes/h_sapiens/gca_grch38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
  /Users/jose/devel/pop_var_caller/benchmarks/ssr_hg002/bam/30x/HG002_TR_v1.0.1_Tier_30x.bam chr1

# knobs: PVC_PROBE_MAX_LOCI=n  PVC_PROBE_WHOLE_CONTIG=1  PVC_GENERIC_REGION_CHUNK_BP=n
#        PVC_PROBE_MAX_RECORD_SPAN=n
# chr21 is the iteration fixture: 6.8 s per run against chr1's 33 s.
```

---

## 4. Build / toolchain configuration

**The `[profile.*]` audit is a clean pass.** Fat LTO, codegen-units = 1, panic = abort,
`debug = "line-tables-only"`; `debug-assertions` correctly off in release with a dedicated
`[profile.soak]` arming them; `target-cpu` set per-target in `.cargo/config.toml`; toolchain
pinned; `[profile.profiling]` documented in-file at `Cargo.toml:55-62`. No change recommended.

Three things do need doing.

**B1 — The feature-varying-build hazard is documented nowhere, and it destroyed a measurement
set in this review.** `cargo build --release --features dhat-heap --example X` writes to the
same `target/release/examples/X` path as a plain release build, so a subsequent plain run
silently executes the **instrumented** binary at ~5–6× slower. This review lost its first
per-contig and whole-genome numbers to it (chr2 read 201.7 s instead of 41.4 s) and caught it
only because the ratio was implausible. All eight `examples/dhat_*.rs` document the invocation
without the warning. **Fix:** add `--target-dir target-dhat` to every documented dhat
invocation and one line saying why.

**B2 — `cargo test --release` is red on a clean tree.** Four tests fail in release at `d95ce8b`
with no patch applied — `ng::alignment::left_align_repeated::tests::an_offset_past_the_reference_panics`,
the same in `left_align_structured`, `ng::alignment::ssr_marginal_sequence::tests::the_epsilon_endpoints_degenerate_without_flooring`,
and `var_calling::dust_filter::tests::sdust_mask_debug_asserts_on_tiny_window`. All four assert
on a `debug_assert!` message that release compiles out, so they pass only in debug. This is a
trap for exactly the workflow a perf review uses: it cost this review a false correctness alarm
against a good patch. **Fix:** `#[cfg_attr(not(debug_assertions), ignore)]` on the four, or run
them under `[profile.soak]`.

**B3 — the probe's allocator guard did the opposite of what its comment claimed** — found by
the `methodology` agent, **already fixed** in the working tree. It read
`#[cfg(all(feature = "alloc-mimalloc", not(feature = "dhat-heap")))]` under a comment asserting
that enabling both would fail the build; the `not(...)` is precisely what let them coexist and
silently measure dhat. It is now a `compile_error!` that names itself.

---

## 5. Code-level findings

### Hot-path

**H1: [src/ng/raw_chrom_reader.rs:207](../../../../src/ng/raw_chrom_reader.rs#L207) — a 64 KiB
stack buffer is zero-initialised on every reference refill.** *Confidence: High.*
`let mut read_buf = [0u8; FILE_READ_CHUNK];` is a 65,536-byte local zeroed on each call, and
there are 384,739 refills in the profiled slice alone. This attributes **1,035 of the profile's
1,090 `__bzero` samples (95 %)**, which the baseline had left as an open question. Hoisting it
to a `Box<[u8]>` field on the reader: **chr21 typed −7.4 %, whole-contig control −30 %** —
the single largest item in the applied patch. *Complexity:* one field, one extra parameter on a
private helper. **In the §2 patch.**

**H2: [src/ng/tandem_repeat.rs:497](../../../../src/ng/tandem_repeat.rs#L497) — the score loop
canonicalises every base 12 times.** *Confidence: High.* Two operands × six periods, each
through a *dependent* load pair. Precomputing the canonical window once: **−5.0 % of the whole
run**, and null in the typing-free control — which is what proves the time is in the typing
scan rather than elsewhere. `cargo asm` confirms the variant's inner loop is
`ldrb/ldrb/cmp/ccmp/csel`: no bounds checks, no table indirection. *Complexity:* one `Vec<u8>`
per `find_tandem_repeats` call (~102 kB per window). **In the §2 patch.**

**H3: [src/ng/read/input/region_query.rs:218-221](../../../../src/ng/read/input/region_query.rs#L218-L221)
— every region query re-seeks a BGZF reader that already holds the block it wants, and noodles
re-inflates it.** *Confidence: High.* Measured with counters on chr21: **82 % of all BGZF seeks
target the block the reader already holds** (27,731 of 33,671), 98 % of them backwards;
`SEEK_SAME_BLOCK_NS` is **30.1 % of the walk** at 68.6 µs each; **the same 35,228 records are
decoded 1,067,729 times — 30.3× re-decode**. Root cause is in the dependency: noodles-bgzf
0.47.0 `src/io/reader.rs:175-186` `seek` unconditionally `lseek`s, `read`s and runs
`parse_block` (inflate + CRC-32), with no "already in this block" fast path and no public way
to reposition inside the held block. **This is the largest unrealised lever in the review.**
*Complexity:* ~150–250 lines with a real correctness surface — the BAM analogue of the
`DecodedContainer` already in this file for CRAM. **Not built. §3 item 5.**

**H4: [src/ng/locus_generation/pileup/generator.rs:834-887](../../../../src/ng/locus_generation/pileup/generator.rs#L834-L887)
— the read query is opened at region grain, so query cost is paid per region rather than per
base pair.** *Confidence: High.* `PileupWalker::new` peeks the first read
([genome_walk.rs:104-111](../../../../src/ng/locus_generation/pileup/genome_walk.rs#L104-L111)),
and that peek is **50.6 % of the walk**: `OPEN_WALK_NS − READS_IN_REGION_NS` (3.14 s) ≈
`SEEK_NS + READ_RECORD_NS` (2.96 s). It is *not* planning — the in-memory BAI query is
0.43 µs/region (0.7 %) and all of `reads_in_region` is 0.58 µs/region (0.9 %). Measured ceiling,
identical `loci=256391` at every grain:

| regions | seconds | seeks | records read | re-decode |
|---:|---:|---:|---:|---:|
| 116,775 (400 bp) | 4.602 | 30,532 | 975,142 | 27.7× |
| 4,671 (10 kb) | 1.914 | 1,592 | 72,727 | 2.1× |
| 1 (whole contig) | 1.382 | 1 | 35,228 | 1.0× |

**A 10 kb query window recovers 72 % of the win at 1/25th of the query count.** The fix is to
widen the *query* without widening the *region*, which leaves the interface unchanged — see
N2 for why caller-side coalescing is the wrong shape. *Complexity:* as H3. **Not built.**

**H5: [src/ng/locus_generation/pileup/open_record.rs:389](../../../../src/ng/locus_generation/pileup/open_record.rs#L389),
`:545`, `:640`, `:1728` — the per-record fold state is a hash map whose 96-byte entry forces a
6,216-byte allocation per record, and whose seeded iteration order is then repaired by two
sorts.** *Confidence: High.* `size_of::<(u32, FoldedReadState)>() == 96`, and
`RECORD_FOLDED_READS_INITIAL_CAPACITY = 32` makes hashbrown round to 64 buckets. Two program
points with that stack are **465,357 blocks / 2.89 GB = 20.3 % of everything the run
allocates** — the largest ng-owned byte site in the profile. `keyed_observations` then builds a
`Vec<u32>` and sorts it *only* to recover determinism from a per-process-seeded map, and
`refold_live_reads` sorts the same ids again; reads are admitted in ascending `read_id`, so a
sorted `Vec<(u32, FoldedReadState)>` gives that order structurally. **Measured −3.9 % on chr21,
6/6 paired runs, output identical, 316 tests green including the cross-process determinism
test.** *Complexity:* one `FoldedReads` newtype with nine methods, four annotations, two sorts
deleted; it retires an invariant rather than adding one. **Not in the §2 patch** — it touches
`open_record.rs` alongside H6/H7 and needs a rebase; see §3.

**H6: [src/ng/locus_generation/pileup/open_record.rs:239](../../../../src/ng/locus_generation/pileup/open_record.rs#L239)
— `witness_of` builds three buffers, a sort and a merge to answer `Complete`, which is
1,646,289 of 1,647,161 answers (99.95 %).** *Confidence: High.* An 8-line fast path:
**whole-contig −5.4 %, typed −1.4 %.** *Complexity:* one early return. **In the §2 patch.**

**H7: [src/ng/locus_generation/pileup/open_record.rs:779](../../../../src/ng/locus_generation/pileup/open_record.rs#L779)
— `finalise` resolves every folded read's witness twice**, once in `keyed_observations` and
again in a second loop over the same map purely to feed four counters. *Confidence: High.*
Moving the tally into the loop that already has the witness: **whole-contig −3.2 %, typed
−1.9 %.** **In the §2 patch.**

**H8: [src/ng/locus_generation/pileup/genome_walk.rs:845](../../../../src/ng/locus_generation/pileup/genome_walk.rs#L845)
— the mate-overlap resolver allocates a hash map plus one `Vec` per distinct chain id, for
lists almost always of length 1.** *Confidence: High.* Replaced with one reused
`Vec<(ChainId, usize)>`, sorted, groups read off as runs. Output identity holds because each
contributor is in exactly one group, `bq_updates` touch disjoint indices, `to_remove` is
sorted and deduped, and within a group indices stay ascending.

**⚠ The author measured −2.6 % (chr21) / −2.4 % (chr1); on a quiet host it is −1.2 % on chr1,
and −0.85 % marginal on top of the §2 patch.** The original number was taken while five other
agents were running on the same machine. It is still a win and it is in the patch, but it is
half its reported size — this is the one place the orchestrator's re-verification changed a
result, and it is why cross-agent numbers do not enter a report unre-run. **In the §2 patch.**

**H9: [src/ng/reference_info.rs:1018](../../../../src/ng/reference_info.rs#L1018) — the
background FASTA verification is a fixed ~11 s tail, not an overlap.** *Confidence: High.*
A/B with a skip knob, chr21, two runs each way, output identical on all four:

| | walk `seconds=` | `real` | `user` |
|---|---|---|---|
| verify on | 5.877 / 5.197 | 12.07 / 10.83 | 16.28 / 15.16 |
| verify off | 5.759 / 5.156 | 5.76 / 5.16 | 4.87 / 4.46 |

The md5 of the 3.1 GB reference costs ~11 s of CPU and **does not contend** with the walk
(`seconds=` moves +0.8–2.0 %); the finished walk then **blocks 5.9 s in `handle.join()`**. It
does not shrink with the work: ~0 % of the whole-genome run, **~50 % of any run under ~11 s** —
which is every worker in a future per-region or per-sample fan-out, and every one of this
review's own chr21 iterations. This also explains the whole-genome `552 s real` against `453 s`
of walk. *Fix:* persist the digest beside the `.fai` under the same `(path, size, mtime)` rule
`cache_key` already uses in memory. *Complexity:* one small on-disk artefact and its
invalidation rule.

### Likely

**L1: [src/ng/locus_generation/witness.rs:258-294](../../../../src/ng/locus_generation/witness.rs#L258-L294)
— `ReadWitness` costs 32 bytes in every observation so that 0.05 % of them can carry a set.**
Boxing the payload was measured for size, not for time: `ReadWitness` 32 → 8,
`SequenceObservation` 112 → 88, `ObservationKey` 64 → 40, `KeyedObservation` 120 → 96. Its real
argument is the fan-out — a 21 % cut in the type a per-sample cohort merge holds N copies of —
and the agent said so rather than claiming a wall-clock win. **Measurement plan:** paired
binaries on chr21; merge on the fan-out argument even at neutral wall.

**L2: [src/ng/raw_chrom_reader.rs:125-195](../../../../src/ng/raw_chrom_reader.rs#L125-L195) +
[generator.rs:554-582](../../../../src/ng/locus_generation/pileup/generator.rs#L554-L582) —
`make_reference` is called per region; its first fetch does an `open(2)` plus a linear scan of a
2,580-record `.fai`.** 9,820 FASTA opens on chr21, `OPEN_CONTIG_NS` 2.2 %. **The obvious fix the
source itself suggests is a trap:** `RawChromReader::fetch` extends its window whenever the gap
is under 64 KiB and only `evict_before` shrinks it, so a shared accessor walking a contig grows
monotonically to ~46 MB on chr21 and ~250 MB on chr1 against a 20 MB baseline — it would trade
this review's best result away. The safe fixes are a name→index map instead of the linear
`find`, and `pread` over a shared `Arc<File>`.

**L3: [src/ng/locus_generation/pileup/generator.rs:593](../../../../src/ng/locus_generation/pileup/generator.rs#L593)
+ chain_id_allocator — the run-lifetime chain-id allocator is a cross-region serial dependency,
and the constructor that would break it is `#[cfg(test)]` in a copy-fidelity-locked file.**
`next_id` is `u64` and survives `reset()`; the invariant is *disjointness*, not ordering, so
per-worker id blocks are free — but promoting `with_next_id_for_testing` needs the production
edit first. **This is the item that is cheap now and expensive later**, and it is the reason to
record it before a fan-out exists.

**L4: [src/ng/read/input/reference.rs:238-252](../../../../src/ng/read/input/reference.rs#L238-L252)
— `OpenReference`'s one-contig bound and noodles' repository write lock are where a fan-out
changes cost.** Per-region fan-out within a contig is safe; contig-parallel thrashes the
one-contig clear and must use `unbounded`. Zero cost today — a BAM run never opens the bases.

### Speculative

**S1: [src/ng/locus_generation/pileup/open_record.rs:640,643,670,672](../../../../src/ng/locus_generation/pileup/open_record.rs#L640-L672)
— `keyed_observations` builds four `Vec`s per record**, 11 % of blocks. Filed with an explicit
"the prior is that it won't pay", which §2's table supports. Largely subsumed by H5.

**S2: [src/ng/read/aligned_read.rs:71,102,103,109](../../../../src/ng/read/aligned_read.rs#L71-L109)
— `decode_record` builds four `Vec`s per read.** Same prior.

**S3: [src/ng/locus_generation/pileup/open_record.rs:2271-2320](../../../../src/ng/locus_generation/pileup/open_record.rs#L2271-L2320)
— `ReadContribution` is 104 bytes, 80 of them one inline `SmallVec`, and the vector of them is
rebuilt once per reference base.**

**S4: [generator.rs:422,553,585](../../../../src/ng/locus_generation/pileup/generator.rs#L422)
— `PileupGenerator` is `!Send` and so is `Arc<WindowedRefSeq>`**, so a fan-out must construct
the whole per-worker stack inside the worker; only the `Send + Sync` triple
`(PathBuf, Arc<ContigList>, Arc<fai::Index>)` can cross. The design is right; the signature does
not say so. Worth a constructor contract before a caller assumes otherwise.

**S5: `cigar_cursor.rs` per-position event rebuild** — capped here because the file is
copy-fidelity frozen and releasing it is an owner decision.

### Note

**N1 — the halo (`max_record_span`) is not a lever, measured and closed.** Narrowing it 12.5×
(5,000 → 400) left `records_outside_region` (883,083) and `reads_admitted` (374,437)
**bit-identical** and moved wall 7 %. The BAM index's 16 kb bin resolution is coarser than
either halo, so both queries touch the same BGZF blocks. It is also a correctness knob — it
bounds what a long deletion can reach. Recorded so it is not re-proposed.

**N2 — caller-side region coalescing is the wrong shape.** Generic regions are **not
adjacent**: 96.6 % of chr21 is generic across 102,938 regions, so ~10⁵ short STR regions sit
between every pair. Merging at the caller means merging across regions that belong to another
generator, which breaks spec §2 disjointness. The win lives in widening the *query*, not the
*region* (H3/H4).

**N3 — four measured non-wins, recorded so the byte ranking does not resurrect them.**
`open_record.rs:545`'s map free-list (13 % of all bytes) → **0 %**; `region_query.rs:236`'s
noodles lazy `Record` (36.4 % of all allocations) → **−0.5 %**; `tandem_repeat.rs:383`'s RT
stack reuse (41 % of all bytes) → **≈0 %**, variant losing several of ten interleaved pairs;
`segment_criteria.rs:889`'s per-candidate `upper()` → null, because `prefilter` runs first.
Patches for the first two are preserved in `patches/`.

**A direct contradiction between two agents, resolved by measurement.** `data_layout` closes
H5 by proposing a `folded_reads` free-list as "an independent second lever, three lines and
untested". `allocations` built exactly that and measured **zero**. The pair is the cleanest
illustration in the review, because it is the *same struct* under two interventions: pool the
allocation → 0 %; replace the structure so the hashing and both sorts go → −3.9 %.

**N4 — a first `BString` rewrite that removed 5.09 M allocations was *slower* than baseline.**
It differed from the successful version only in unpacking the sequence through the lazy
per-base iterator instead of a nibble-pair table. Worth knowing before anyone assumes the lazy
API is free.

**N5 — the reader-pool and counts mutexes are not a cost.** Taken ~3× per region query
(~1.8 M lock pairs on chr1) but `pthread_mutex_unlock` is **3 samples of 25,446**, with no
`__psynch_*` / `__ulock_wait` / `swtch_pri` frames on the main thread at all. Recorded so the
pool is not read as free if the per-region query count ever *grows*.

**N6 — `ReferenceInfoCache`'s single-flight is uncontended and doing its job.** The profile's
"third thread, 9,429 samples in `__psynch_mutexwait`" is not an idle thread: it is the
**second** verification call blocked on the single-flight slot and then hitting the cache —
the cache saving a duplicate 3.1 GB read. This corrects a misreading in the orchestrator's own
first pass.

---

## 6. Out-of-scope observations

- **`ng_generic_loci_dump` buffers every row of the whole run** before rendering, which is what
  made Milestone E's RSS number un-attributable. Not a defect in a dump tool, but it should
  carry a line saying it is not a memory instrument — the probe is.
- **`md5::compress` is scalar** and would be a SIMD or parallel-hash candidate *if* per-run
  verification survives H9. Prefer H9's caching, which removes the work rather than speeding it.
- **The whole-genome run shows `user + sys` = 447 s against 552 s `real`** — ~105 s with no CPU
  running at all, i.e. I/O-blocked, and only 0.26 s of it is startup on chr21. Unexplained;
  worth a look if whole-genome wall becomes the target.
- **`chrom_name.clone()` per refill** — 2.8 % of allocations, for an error value never built.

## 7. What's already good

- **`with_shared_index` is the right seam and it predates the need.** `WindowedRefSeq` holds
  `Arc<ContigList>` and `Arc<fai::Index>` so k accessors share one parsed `.fai` while holding
  k cursors — which is exactly what a per-worker fan-out needs, and the reason a fresh
  `WindowedRefSeq::new` per region (189 µs of `.fai` parsing) is not what the callers do.
- **`ActiveRead` is correctly AoS at 184 bytes** — the per-position loop reads essentially every
  field, so the SoA rule genuinely does not apply. Checked and rejected rather than skipped.
- **`WitnessedLocusPositions`' two-run inline capacity is already at the SmallVec floor**
  (24 bytes; four runs would cost 32), so "widen it, it's free" is wrong. The existing choice is
  the measured one.
- **`copy_fidelity.rs` makes the production-copy invariant a build failure**, which is why this
  review could bound its own edit surface precisely instead of guessing.

---

## Author response convention

Address each finding by its identifier (H1, L2, …) with one of: `applied in <commit>` /
`experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` /
`won't fix because …`. The "no gain" path is expected — §5 N3 already carries four of them.
