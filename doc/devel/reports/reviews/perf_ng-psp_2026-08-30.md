# Performance Review: ng-psp
**Date:** 2026-08-30
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** ng's psp store, `src/ng/psp/` — the whole module, reader and writer
**Verdict:** Apply the listed wins
**Hot-path evidence:** four sampling profiles, a new criterion bench, a new allocation oracle, and deterministic instruction counts

---

## 1. Scope and constraints

**What was reviewed:** the whole of `src/ng/psp/` — `block.rs`, `chain_ids.rs`, `footer.rs`,
`header.rs`, `index.rs`, `mod.rs`, `reader.rs`, `record.rs`, `trailer.rs`, `walk.rs`, `writer.rs`.
25,016 lines, of which the library is 4,936 and the rest is `#[cfg(test)]`.

**Reviewed against:** branch `ng-psp-perf`, commit `b22860ea` (the module as merged to `main` at
`504d172d`, plus this review's own benchmark).

**Targets.** A psp holds one sample's evidence at every reference position a run analysed — one
record a position, at three reads a position and at three hundred, for a cohort of one sample and of
several thousand. A caller opens one per sample and holds them all open for a whole run. The
design's headline claim ([spec §5.4](../../ng/spec/psp_file_format.md)) is that a 62-sample walk over
471,520,156 records takes 23.1 s against production's `.psp` 42.4 s — **1.8× faster on the same
records**, with 35 % fewer bytes and 7.7× less memory an open sample. The memory budget is 500 kB an
open sample (§1.1), met at 480 kB
([milestone H4](../implementations/ng_psp_h4_2026-08-30.md)). **Throughput was the target under
review; the memory budget was treated as settled and not tradeable.**

**Hardware:** Linux aarch64 in the dev container (the production build target) and an
Apple-silicon macOS host. **No hardware performance counters exist on either** — the PMU is not
virtualized — so there is no cache-miss, branch-miss, `perf c2c` or `perf sched` evidence anywhere in
this review, and every branch claim rests on callgrind's simulated counts or on a behavioural sweep.

**Hot-path evidence available.** All of it was created by this review; none existed beforehand.

- **`benches/ng_psp_perf.rs`** — a full walk, a head-only walk, a one-in-a-hundred walk and a write,
  at 10 reads a position and at 280. Synthesised records, so it needs no corpus under `tmp/`.
- **`examples/ng_psp_against_production.rs`** — both readers over the same records, interleaved.
- **`examples/dhat_ng_psp.rs`** — allocator calls a record, per workload. Identical run to run.
- **Four `sample(1)` profiles** and one un-collapsed `--profile profiling` profile.
- **Deterministic instruction and simulated-branch counts** from callgrind in the container.

**In-scope files:** the eleven above, plus `src/psp/varint.rs` and
`src/ng/locus_generation/mod.rs` where a psp walk calls into them.

**Deliberately out of scope:** `src/psp/` (production's store — frozen; read for comparison, nothing
proposed in it), `src/ssr/`, and everything else outside `src/ng/psp/`.

**The baseline the review established**, which the store did not have before it. Dev container,
Linux aarch64, mimalloc, synthesised records naming every read, **cheapest of three complete
passes**:

| workload | shallow — 10 reads a position | deep — 280 reads a position |
|---|---:|---:|
| full walk | 33.561 ms · **111.9 ns a record** | 36.571 ms · **609.5 ns** |
| head-only walk, no body built | 7.769 ms · **25.9 ns** | 19.935 ms · **332.3 ns** |
| one body in a hundred | 8.318 ms · 27.7 ns | 20.187 ms · 336.5 ns |
| write | 38.795 ms · 110.8 ns | 101.419 ms · 1,448.8 ns |

*Walks cover 300,000 and 60,000 records, writes 350,000 and 70,000; all four arms cut four blocks.*

**What the head is worth, on records that name every read: 4.32× at 10 reads a position and 1.83× at
280.** That is the full walk over the head-only walk, and it answers the question
[arch §7](../../ng/arch/psp_file_format.md) narrowed and left open — *how much of the head's speed-up
survives at 300 reads a position*. [Milestone H5](../implementations/ng_psp_h5_2026-08-30.md)
measured 2.51× and 2.40× on stores built from a production `.psp` and reported the erosion with depth
as about 5 %, while saying its figures were upper bounds and that closing the question needed a store
naming every read. **On records that name every read the erosion is 58 %, not 5 %** — the head grows
with depth while the body the skip avoids does not, which is exactly the mechanism H5 predicted and
could not size. ⚠ These records are synthesised, not ng's own output; the question closes properly
when ng can write a store.

**Categories dispatched:** `methodology` (always), `allocations` (the decode path builds a record per
position), `data_layout` (the live read-identifier set is a sorted array rewritten per record),
`hot_loops` (a row-shaped field parser runs ~16 times a record), `io_and_syscalls` (a 16 kB read
chunk with no `BufReader`). **`concurrency` was not dispatched:** `grep` finds no `Arc`, `Mutex`,
`RwLock`, atomic, `rayon`, `thread::spawn` or channel anywhere in `src/ng/psp/`.

**⚠ One constraint shapes what any of this can conclude.** Nothing in the run reads or writes a psp
yet — `grep` for the reader's entry points across `src/ng/` outside `src/ng/psp/` returns nothing. So
every call frequency below is *design intent* from the specs, not an observed rate, and the corpora
that exist are all written from a production `.psp`, which names about **3.4 %** of the reads ng will
name (the owner's ruling of 2026-08-17). Measured here: the live read-identifier set on those stores
averages **2.0 identifiers a record on hg002 and 0.1 on tomato**, against the ~280 and ~10 that ng
will carry. **They under-weight the chain-id work by far more than 3.4 % suggests**, which is why
this review's benchmark synthesises records that name every read.

---

## 2. Verdict

**Apply the listed wins.** Three sites carry measured fixes with contained complexity, and the
library test suite stays green on all of them, byte-exact record fixtures included: the live-set
apply (H1), the writer's choice of residual observation (H2), and the reader's residual list (H3).
Together they were measured at roughly a third off a write and a fifth to a half off a walk **at 280
reads a position**.

**And one thing they do not do, which is the more important half of this report.** The claim §5.4
sells the format on — 1.8× faster than production — **holds at 280 reads a position and inverts at
10**, where the shipped reader is 2.2× *slower*. §5.4's own walk is 62 tomato accessions, which is
the 10-reads corner. None of H1–H3 addresses that: on production-shaped tomato records the live set
averages 0.1 identifiers, so the chain-id work they fix is not what is slow there. What is slow there
is row-shaped field parsing, and the only finding aimed at it (L1) is the one candidate nobody got to
measure. **Section 5's closing note puts the decision that follows to the owner.**

---

## 3. Measurement plan

Everything here is built and committed on `ng-psp-perf`; this section is how to re-run it, and what
each number gates.

**The machine will lie to you if you let it.** Two independent agents measured the *same binary*
twice and got 76.5 ms against 46.7 ms, and 11.37 / 16.95 / 20.37 ms, on one benchmark. During this
review a full criterion pass reported `walk_one_in_100/deep` at 31.96 ms; re-run alone it was
20.19 ms, and the true figure — where it sits 0.28 ms above the head-only walk, exactly 600 bodies at
470 ns — is 20.19. **Criterion's own `change:` percentages compare across runs and must not be
quoted.** Three rules, all of which this review's numbers follow:

1. Gate on a **count** — allocator calls, syscalls, instructions, branches — because it is identical
   run to run. Use wall time only as a secondary check.
2. When wall time is the only instrument, **interleave the arms inside one process or one container
   session** and take the cheapest round of each.
3. Take **the cheapest of at least three complete passes**, never one.

| what to run | what it gates |
|---|---|
| `./scripts/dev.sh cargo bench --bench ng_psp_perf` ×3, cheapest of each arm | walk and write throughput at both depths |
| `./scripts/dev.sh cargo run --release --example dhat_ng_psp --no-default-features --features dhat-heap` | **allocator calls a record** — the gate for H2, H3 and every L-level allocation finding |
| `./scripts/dev.sh bash -c 'valgrind --tool=callgrind --branch-sim=yes <bench binary> walk_full/shallow'` (filter *without* `--bench` → criterion runs one iteration → deterministic) | **instructions and simulated branches** — the gate for H1 and L1 |
| `target/release/examples/ng_psp_against_production <prod.psp> <ng.ngpsp> --rounds 9` | the §5.4 claim, both arms interleaved |
| `target/release/examples/ng_psp_against_production … --only ng`, sampled with `sample <pid> 25` | where a walk's time goes |
| `cargo build --profile profiling …` before sampling | **required for a breakdown** — see section 4 |

**Two workloads the bench still lacks**, and three findings cannot be gated without them: a walk that
opens a store and walks it **more than once** (L6 prices per-walk setup at zero today), and a corpus
with **more than two observations a record** (L3 and S1 are invisible at one observation).

---

## 4. Build / toolchain configuration

**A `#[global_allocator]` is per binary, and `alloc-mimalloc` being a default feature does not put
one in a bench or an example.** A target without the declaration links the system allocator whatever
the feature says. Measured on `examples/ng_psp_skip_value.rs`, two builds of one source differing
only in that declaration, alternated over 7,687,686 tomato records:

| arm of that harness | system allocator | mimalloc |
|---|---:|---:|
| full walk | **2.594 s** | **1.837 s** |
| body-skipping walk | 0.649 s | 0.628 s |
| **the ratio it prints** | **3.99** | **2.93** |

The error is not neutral for a ratio: the full walk builds a record a position and the skipping walk
builds one in a hundred, so a slower allocator taxes the numerator and leaves the denominator alone —
41 % against 3 %.

**At the reviewed commit, 7 of 11 benches and 25 timing-or-memory examples had no declaration.**
Three consequences, and the third is the one that survives:

- **[Milestone H5's whole result](../implementations/ng_psp_h5_2026-08-30.md) was taken under the
  wrong allocator** — all six cells, from `examples/ng_psp_skip_value.rs`. Its tomato
  one-in-a-hundred cell reads 3.038× (3.025–3.040); the same quantity under mimalloc measured 2.93.
  The report's ±0.5 % spread is a fact about repeats of one misconfigured binary. **Re-take the table
  rather than arithmetically correcting it** — the original runs' load conditions are not
  recoverable.
- **Spec §5.4's own 1.8×** came from `examples/psp_row_stream_roundtrip.rs`, also without the
  declaration. Being a ratio does not protect it: the two readers allocate differently.
- **[Milestone H4's 480 kB an open sample survives](../implementations/ng_psp_h4_2026-08-30.md).**
  Checked rather than assumed: two builds differing only in the declaration, in fresh processes at
  1, 64 and 256 open samples, give a slope of **111.8 kB a sample either way** — the two rises agree
  to the kilobyte. The harness measures a slope across sample counts precisely so the allocator's
  fixed arena cancels, which its own doc predicted. Only the fixed part moves. **The budget stays
  settled.**

**Fix, and the gate that keeps it fixed.** Add to every bench and every timing example:

```rust
#[cfg(all(feature = "alloc-mimalloc", not(feature = "dhat-heap")))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
```

and to `ci.yml`, a check that is deterministic and takes no time:

```sh
test -z "$(grep -L global_allocator benches/*.rs)" \
  || { echo 'a bench with no #[global_allocator] times the wrong allocator'; exit 1; }
```

The `not(feature = "dhat-heap")` guard is required, not decoration: CI lints and tests with
`--all-features`, which turns both on, and two `#[global_allocator]`s do not compile.

**`[profile.release]` is right and is worth 1.48× — do not weaken it casually.** Six alternated
pairs of one example built under `release` (`lto = "fat"`, `codegen-units = 1`) and under
`profiling` (`lto = false`, `codegen-units = 16`), cheapest of each: **0.960 s against 1.421 s**,
release ahead in every pair.

**But a `--release` profile cannot say where a walk's time goes.** Fat LTO inlines
`read_record_head`, `decode_record_body` and `BlockStream::next_record_where` into
`RecordIter::next`, which then appears as one 35.4 % entry. Under `--profile profiling` the same walk
resolves:

| site | share of self time |
|---|---:|
| `FieldReader::read_varint` | **21.2 %** |
| `record::decode_record_body` | **16.2 %** |
| `ZSTD_decompressSequences` | 10.2 % |
| `record::decode_the_body_of` | 7.1 % |
| `block::BlockStream::next_record_where` | 5.3 % |
| `mi_free` | 3.7 % |
| `walk::RecordIter::next` *(itself)* | **3.0 %** |
| `record::read_record_head` | 2.9 % |
| `chain_ids::LiveSetReader::parse_changes` | 2.5 % |
| `psp::varint::decode_u64_leb128_cold` | 2.2 % |
| zstd, every symbol together | 13.4 % |

**Field decoding is 26.3 % of this walk against the whole zstd decompressor's 13.4 %** — twice as
much. Take both profiles from now on: `--release` for the shares that describe the shipped binary,
`--profile profiling` for the breakdown, and never compare a timing from one with a timing from the
other.

**One gap, filed but not measured.** `.cargo/config.toml` sets `target-cpu` for `x86_64-linux`
(`x86-64-v3`) and `aarch64-macos` (`apple-m1`) and **has no entry for `aarch64-linux`** — the dev
container, which is the production target. So the container's numbers are taken at the generic
aarch64 baseline and the host's at `apple-m1`. On aarch64 NEON is already in the baseline, so expect
little; the reason to record it is that the two tables this review reasons across were built at
different ISA floors and nothing said so.

---

## 5. Code-level findings

### Hot-path

#### H1: [chain_ids.rs:241](../../../../src/ng/psp/chain_ids.rs#L241) and [chain_ids.rs:222](../../../../src/ng/psp/chain_ids.rs#L222) — the live read-identifier set is rewritten end to end at every record that has any change

- **Confidence:** High
- **Hot-path evidence:** `sample(1)` on the head-only walk at 280 reads a position — a walk that
  builds no body at all — puts `chain_ids::apply_arrivals` at **44.2 % of self time**, with
  `apply_departures` folded into the 45.0 % walk-loop bucket. Callgrind on the shipped binary
  attributes 174,110,105 instructions (3.88 % of a run) to `apply_arrivals` plus 68,540,008 to
  `vec/mod.rs` beneath it.
- **Mechanism:** `apply_arrivals` merges into a scratch `Vec` and swaps it in, reading all 280 live
  identifiers and writing 280 more; `apply_departures` runs `Vec::retain`, evaluating a closure over
  all 280. At 280 reads a position about two identifiers arrive and two depart at every position, so
  **560 identifiers move to change four**, at nearly every record.
- **Two agents fixed this independently, with different designs, and neither has been measured
  against the other.** That comparison is the one piece of work this finding still needs.

  **(a) In-place, no scratch** — departures compacted from the first departure's index with one
  `copy_within` for the tail; arrivals appended outright when the lowest sorts above every live
  identifier, and otherwise merged backwards into the tail alone. Measured in one process with the
  arms rotating, cheapest of 7 rounds, 60,000 records at 280 reads:

  | | full walk | head-only walk | write |
  |---|---:|---:|---:|
  | as shipped | 1,266 ns a record | 661 ns | 4,224 ns |
  | in place | 816 ns (**0.645×**) | 152 ns (**0.231×**) | 3,141 ns (**0.744×**) |

  Decomposed, head-only: 415 ns → 217 (append path alone) → 96 (adding departure compaction) → 94
  (adding the tail merge).

  **(b) One merged pass** — departures skipped and arrivals interleaved in a single traversal,
  keeping the scratch buffer. Interleaved wall clock, head-only walk at 280 reads, cheapest of 3
  alternating rounds: **28.872 ms against 23.609 ms, 18.2 % faster**, variant ahead in every round.
  Deterministic: −34,379,123 conditional branches and −992,119 simulated mispredicts.

- **⚠ How much this is worth depends on a corpus property that neither the benchmark nor any real
  store had until this review.** The saving comes mostly from the append fast path, which applies
  only when no arriving identifier sorts below one already live. A chain id names a read *pair*, and
  spec [`psp_record_encoding.md`](../../ng/spec/psp_record_encoding.md) §6 measures **83 % of
  identifiers on the human sample and 91 % on tomato covering two stretches** — so identifiers go
  absent and come back, and the fast path does not apply to them. Across that axis:

  | corpus | full walk | head-only walk |
  |---|---:|---:|
  | 280 reads, **no** returning identifiers | 0.614× | 0.251× |
  | 280 reads, 0.69 returning a record *(about the real rate)* | 0.763× | 0.546× |
  | 280 reads, 2.76 returning a record *(three times it)* | 1.011× | 1.040× |
  | 10 reads a position | 0.984× | 1.025× |

  **At the realistic rate it is worth 0.55× to 0.76× at depth; at 10 reads a position it is a wash;
  it never costs more than 4 % anywhere measured.** `benches/ng_psp_perf.rs` has since been corrected
  to build paired reads and now reports its own rate — 0.933 returning identifiers a record at 280
  reads, 0.033 at 10.
- **Measurement plan:** build both designs and run them against each other **in one process**,
  rotating arms, on the corrected corpus at both depths. Gate on conditional-branch count from
  callgrind (deterministic); require the cheapest interleaved round not to regress at 10 reads a
  position. Merge whichever wins; they are not combinable.
- **Complexity cost:** design (a) is about 20 lines with a backwards merge whose invariant
  (`slot = unmoved + still_to_place - 1`) has to be written out, and it hand-rolls what `retain` did.
  Design (b) is one function replacing two, holding the same preconditions in one place. **On memory
  they differ and (a) wins: it deletes `merged_ids` from both `LiveSetReader` and `LiveSetWriter`,
  giving back 2,264 bytes an open sample at 280 reads** — 6.8 MB across three thousand samples —
  where (b) keeps it. The proptests already reach the interleaving case.

#### H2: [record.rs:701](../../../../src/ng/psp/record.rs#L701) — a one-observation record does the full general path to reach the only answer it could have reached

- **Confidence:** High
- **Hot-path evidence:** `sample(1)` on the deep write puts `record::encode_record_body` at 23.6 % of
  self time and `core::slice::sort::unstable::ipnsort::<u64>` at **9.7 %** — both sorts are this
  function's. Callgrind on the shipped binary: `ipnsort::<u64>` is **234,232,374 instructions,
  5.22 % of a run**, the fourth-largest entry, with `bcmp` a further 58,470,053.
- **Mechanism:** at one observation the residual is observation 0 or nothing, and which one depends
  only on whether that observation's identifier list is already a strictly ascending set. The general
  path reaches that answer by allocating `named_elsewhere` (empty here) and `every_read`, sorting and
  deduplicating both, building a `LiveSet`, deriving the residual and comparing it to the original.
  One pass of `windows(2)` decides it. **This is about ninety-nine positions in a hundred** at the
  measured corner.
- **Measured, in the same binary as H1's design (b) and not separable from it in wall clock:**
  cheapest of 3 interleaved rounds, `write/shallow` **6.377 ms → 4.346 ms (−31.8 %)`, `write/deep`
  **84.620 ms → 53.528 ms (−36.7 %)**, variant ahead in every round at both depths. Deterministic,
  and this one *is* separable: total instructions 4,489,248,708 → 3,447,947,891 (**−23.2 %**),
  simulated mispredicts −30.2 %, and `ipnsort::<u64>` **halved** — 234,232,374 → 118,429,059, which
  is the signature the mechanism predicts, two sorts a record becoming one.
- **Not a format change:** 5,330 library tests pass, byte-exact record fixtures included.
- **Measurement plan:** callgrind, gating on total instructions falling more than 5 % with the
  fixtures still byte-exact; `dhat_ng_psp`'s `write/*` allocator calls a record as the second gate.
- **Complexity cost:** one branch at the top of the function, and a second place stating the rule
  about what makes a list derivable — the two could drift. Any corpus contains both shapes, so the
  existing round-trip tests cover both.

#### H3: [chain_ids.rs:334](../../../../src/ng/psp/chain_ids.rs#L334) — the residual read list is built by pushing into an empty `Vec`, so it reallocates its way up to the depth at every record

- **Confidence:** High
- **Hot-path evidence:** `sample(1)` puts `chain_ids::residual_reads` at **20.1 %** of the deep full
  walk and 6.3 % of the shallow one, with `RawVecInner::finish_grow` at 3.4 % and `RawVec::grow_one`
  at 1.2 % beside it. `examples/dhat_ng_psp.rs` names the count: **11.04 allocator calls a record**
  at 280 reads, of which this is 7.
- **Mechanism:** `residual_reads` clears `into` and pushes one identifier at a time. `into` reaches
  the caller as the residual observation's `chain_ids`, so it is a fresh `Vec` at every record: from
  no capacity, 280 pushes cost eight allocations and copy 508 identifiers between buffers. The length
  is not a guess — every identifier in `named_elsewhere` is live, so the answer is exactly
  `live.len() - named_elsewhere.len()`. One `reserve` makes it one allocation, and it *shrinks* the
  record handed to the caller: 280 identifiers in a buffer of 280 rather than of 512.
- **Measured** (`examples/dhat_ng_psp.rs`, allocator calls a record — identical run to run):

  | workload | as shipped | with H3 | with H3 and H4 |
  |---|---:|---:|---:|
  | `walk_full/deep_280_reads` | 11.04 | **4.13** | 4.13 |
  | `walk_full/shallow_10_reads` | 6.02 | **4.05** | 4.05 |
  | `write/deep_280_reads` | 8.97 | 2.11 | **0.088** |
  | `write/shallow_10_reads` | 4.00 | 2.04 | **0.017** |

  Bytes asked for a record fell with them: `walk_full/deep` 7,933 → 2,393. Interleaved wall clock,
  paired rounds: the write was cheaper in **4 of 4** rounds at both depths; the deep full walk in 9
  of 12; the shallow full walk showed nothing above that machine's noise.
- **Measurement plan:** `dhat_ng_psp`; merge when `walk_full/deep` falls below 5 allocator calls a
  record. Wall time secondary, interleaved.
- **Complexity cost:** one line. The invariant it leans on — `named_elsewhere` ⊆ the live set — is
  what the caller's own guard `check_a_derived_read_list` already enforces, and `saturating_sub`
  makes a record that breaks it cost one growth rather than a panic. **Zero bytes added to the
  480 kB**: the buffer goes to the caller inside the record and is now smaller than it was.

#### H4: [record.rs:701](../../../../src/ng/psp/record.rs#L701) — `residual_observation_of` builds three fresh `Vec`s a record on the write path

- **Confidence:** High
- **Hot-path evidence:** the same profile entries as H2, plus `dhat_ng_psp`: **8.97 allocator calls
  and 10,082 bytes a record** to write at 280 reads, for a record whose own bytes already go into a
  reused buffer.
- **Mechanism:** every record collects `named_elsewhere`, collects `every_read`, moves the second
  into a `LiveSet` by value and pushes the derivation into a fourth `Vec::new()`. All four die at the
  end of the call. Held on the encoder and cleared per record they allocate once for the file —
  which is the shape `RecordEncoder` already uses for `body_scratch` and `LiveSetWriter` for
  `now_live`. **This is the module's own reuse discipline applied to the one place that skipped it.**
- **Measured:** the third column of H3's table — **write at 280 reads goes from nine allocator calls
  a record to one call every eleven records**, and bytes a record from 10,082 to 176. Interleaved
  wall clock: cheaper in 4 of 4 paired rounds, mean 184.0 → 118.0 ms deep and 15.7 → 9.7 ms shallow.
- **⚠ H2 and H4 overlap and should be decided together.** H2's fast path skips the general path
  entirely for one-observation records, which is 99 in 100; H4 makes the general path cheap for the
  other 1 in 100. Taking H2 first shrinks what H4 is worth, and the write measurements above cannot
  be added.
- **Complexity cost:** a `ResidualScratch` of three `Vec`s, a field on `RecordEncoder`, and a second
  entry point beside `encode_record_body` so the public API is untouched. **About 7.2 kB a *writer*
  at 300 reads a position, and nothing on the reader side** — a writer is one per sample being
  written, not one per sample held open, so the 480 kB budget is not touched.

### Likely

#### L1: [record.rs:1123](../../../../src/ng/psp/record.rs#L1123) — the error type is too wide to return in registers, so a six-instruction fast path carries a 96-byte frame

- **Confidence:** Medium — **this is the one candidate nobody measured**, and it is the only finding
  aimed at what actually makes ng slow at 10 reads a position.
- **Hot-path evidence:** `FieldReader::read_varint` is **21.2 %** of the un-collapsed real-corpus
  walk, `read_u32` 1.7 %, `read_locus_kind` 1.2 %. `objdump` on the shipped binary shows the whole
  out-of-line copy: ten memory operations of prologue and epilogue, six callee-saved registers
  spilled and a 96-byte frame, wrapped around six instructions that do the work.
- **Mechanism:** `x0` is a hidden return pointer, so `Result<u64, RecordDecodeError>` comes back
  through memory rather than in a register pair; and because `Malformed` owns a `String` the frame is
  96 bytes. `Cargo.toml` already carries `result_large_err = "allow"` naming exactly this trade for
  the cohort driver's errors — **this is the same trade on a path that runs about sixteen times a
  record**.
- **Measurement plan:** split the primitive — an inner `fn next_varint(&mut self) -> Result<u64,
  VarintError>` (a fieldless two-variant enum, so the `Result` is 16 bytes and returns in registers)
  and a `#[cold] #[inline(never)]` converter to `RecordDecodeError`. Gate on **total instructions**
  from callgrind falling more than 2 % with `read_varint`'s own attribution falling more than 20 %;
  interleaved wall clock as the secondary check. The risk the gate exists to catch is that LLVM
  already flattens this at the inlined call sites and only the out-of-line copy benefits.
- **Complexity cost:** a second, thinner fault type inside the reader and a two-line wrapper on every
  `read_*`. The public error class and its messages are unchanged, which is what Milestone D's
  restartable-read contract depends on.

#### L2: [record.rs:1181](../../../../src/ng/psp/record.rs#L1181) — a 64-bit division about five times a record, to compare against a constant

- **Confidence:** Medium — measured only bundled with H2, so its own share is not separated.
- **Mechanism:** `MOST_BYTES_A_BODY_CAN_DECLARE / least_bytes_each as u64` is evaluated on every
  call before anything is compared. `least_bytes_each` is a constant at every call site, so fat LTO
  may fold it — but `read_count` is `pub(super)` and reached from `chain_ids.rs` too, and where the
  fold does not happen an aarch64 `udiv` is tens of cycles. Multiplying instead is exact: for
  positive integers `declared > M / L` and `declared * L > M` are the same statement, because the
  truncation in the division is precisely the slack the multiplication keeps.
- **Measurement plan:** build it **alone**, not bundled, and gate on the callgrind total dropping at
  all — the change cannot cost anything.
- **Complexity cost:** one `saturating_mul` and a comment carrying the equivalence argument, which is
  not obvious. The refusal keeps its wording by computing the quotient inside the cold branch.

#### L3: [chain_ids.rs:308](../../../../src/ng/psp/chain_ids.rs#L308) — `decode_read_list` reads a count and then pushes that many identifiers without reserving

- **Confidence:** High that the allocations are there; Medium that they matter, because the record
  that exercises them has more than one observation, which is 1 in 100 in the only corpus that exists.
- **Mechanism:** the count is read before the loop, and the module already owns the safe way to turn
  a declared count into a reservation — `entries_to_reserve`, which bounds the declared number by
  what the remaining bytes could hold, exactly so a hostile body claiming a million reads in eleven
  bytes cannot drive the allocator. An observation naming 70 identifiers costs five allocations
  without it and one with.
- **Measured** on the corpus that exists: `walk_one_in_100/deep` 0.278 → 0.170 allocator calls a
  record. Real, and small.
- **Measurement plan:** **the corpus comes before the gate** — add a shape with three or four
  observations a record to `benches/ng_psp_perf.rs` and `examples/dhat_ng_psp.rs`, then gate on
  `walk_full` falling by about three allocator calls a record there.
- **Complexity cost:** one call, plus widening `FieldReader::bytes_left` and `entries_to_reserve` to
  `pub(super)`. Both already live in `record.rs` and are used by its own decoders.

#### L4: [record.rs:1239](../../../../src/ng/psp/record.rs#L1239) — the decoder proves the witness runs canonical, then hands them to a constructor that sorts and merges them again

- **Confidence:** Medium
- **Hot-path evidence:** `locus_generation::witness::canonicalise_runs` is **1.4 %** of the real
  tomato walk's self time — small, but entirely redundant.
- **Mechanism:** `read_witness` already refuses a run covering no position, one starting at or before
  the previous run's end, and one ending past `u16::MAX` — which is exactly `canonicalise_runs`'s
  postcondition. It then builds a `Vec` and passes it to `from_half_open_runs`, which collects into a
  `SmallVec`, checks emptiness again, sorts a list already known sorted, and runs a merge that can
  never merge anything.
- **Measurement plan:** **`benches/ng_psp_perf.rs` writes `ReadWitness::Complete` for every
  observation and cannot see this change at all** — measure on the real tomato corpus with
  `ng_psp_against_production --only ng`, gating on `canonicalise_runs` leaving the profile.
- **Complexity cost:** a second constructor on a type whose whole documented reason for existing is
  that only its constructors may build it, and whose invariant fails *silently* when broken — two
  spellings of one witness stop merging into one observation. **That is a real widening of a
  deliberately narrow interface for 1.4 %; take H1–H3 first and revisit only if a profile still names
  it.**

#### L5: [reader.rs:70](../../../../src/ng/psp/reader.rs#L70) and [index.rs:24](../../../../src/ng/psp/index.rs#L24) — the block index is 24 bytes an entry with 4 of them padding, and at a whole genome it is 336 kB of the 480 kB budget

- **Confidence:** High on the arithmetic, Medium that the complexity is worth spending.
- **This is a memory finding, not a speed one.** The only search over the index,
  `records_from`'s `partition_point`, runs once per walk.
- **Mechanism:** measured, `size_of::<BlockIndexEntry>() == 24` and `size_of::<GenomePosition>() ==
  16` — a `u32` beside a `u64`, so 4 bytes of padding that field order cannot remove. Spec §3.3 puts
  a whole genome at **roughly 14,000 blocks**: 14,000 × 24 = **336,000 bytes**, held per open sample
  for the whole run. A block never crosses a contig (§3.2) and entries are in genome order, so the
  contig column is a run per contig and can leave the per-entry record entirely, taking the entry to
  16 bytes.
- **⚠ And it says something about the budget that H4's measurement could not.** H4's 480 kB was taken
  on stores holding **one chromosome** — `tomato.ngpsp` has 160 blocks and `hg002.ngpsp` 281, so
  their whole indexes are under 7 kB. A whole-genome store adds ~336 kB on top. **The 500 kB budget
  is met on the corpora that exist and has not been tested on a store of the size the run will
  build.**
- **Measurement plan:** the `size_of` arithmetic is the first gate, exact and deterministic. End to
  end needs **a store with ~14,000 blocks built for the purpose** — the corpora under
  `tmp/perf_review_2026-08-30_ng-psp/corpora/` cannot show it — then `ng_psp_open_cost` at
  `--samples 1, 8, 64`, one fresh process each, and the slope.
- **Complexity cost:** `block_index()` hands out `&[BlockIndexEntry]` today and would become an
  accessor; `decode_index` / `encode_index` can keep their current types, confining the change to
  `PspReader`. **Gives 112,000 bytes back an open sample — 23 % of the budget, 336 MB across three
  thousand samples.**

#### L6: [walk.rs:107](../../../../src/ng/psp/walk.rs#L107) via [block.rs:1543](../../../../src/ng/psp/block.rs#L1543) — every walk builds a fresh `BlockStream`, so its buffers and its zstd context are per walk rather than per open sample

- **Confidence:** High that the cost is there; Medium that it matters, because the caller that would
  pay it repeatedly does not exist yet.
- **Measured:** starting a second walk on an **already-open** reader costs **34.74 allocator calls
  and 34,975 bytes every time**, identical at both depths because it is a function of the reader and
  not of the data. For scale, walking one record with its body built costs 4.05. **⚠ The zstd decode
  context is not in that 34,975** — libzstd allocates through C `malloc`, which dhat's Rust hook does
  not see; spec §5.3 puts it at about 190 kB, so the true figure is roughly 225 kB a walk.
- **Mechanism:** spec §6.2's own usage shape is `records_from(at)?.building_only_where(…)` — **one
  walk a region**. A cohort of several thousand samples reading region by region pays this per sample
  per region.
- **Measurement plan:** **`ng_psp_perf` walks each store exactly once and prices this at zero — add a
  many-short-regions workload before merging anything here.** Gate on the per-walk allocator count
  from `dhat_ng_psp`'s `walk_start_again` phase.
- **Complexity cost:** honest framing is that it **moves** the cost rather than removing it — the
  reader would hold 32 kB of buffers plus a ~190 kB zstd context whether or not a walk is running. It
  wins when a sample is walked more than once and loses a little when a sample is walked once and
  then held. **Establish the caller's shape first**: if the merge does one walk a sample, this is a
  Note.

#### L7: [record.rs:895](../../../../src/ng/psp/record.rs#L895) and [record.rs:921](../../../../src/ng/psp/record.rs#L921) — the reference bases and each observation's bases are a `Box<[u8]>`, which at SNP shape is a heap allocation for one byte

- **Confidence:** Medium
- **Hot-path evidence:** with H3 and H4 applied, a full walk costs **4.05 allocator calls a record at
  10 reads and 4.13 at 280** — and two of those four are these. The other two are the observations
  vector and the residual list, which the record genuinely owns.
- **Mechanism:** an inline small-bytes type would keep short sequences in the struct and spill only
  for long ones. **This codebase already has exactly one such type** — `Motif { buf: [u8;
  MAX_MOTIF_LEN], len: u8 }` — which is why an SSR record's motif costs nothing while its flanks cost
  an allocation each. **It leaves the record owned and lifetime-free**, which is the contract
  `SampleLocusObservations`'s own doc calls load-bearing, so this is *not* the
  borrow-from-the-rolling-buffer idea and does not pay that price.
- **Measurement plan:** count first — put `smallvec::SmallVec<[u8; 16]>` behind the two fields on a
  branch, rerun `dhat_ng_psp` expecting `walk_full/shallow` to fall from 4.05 toward 2.0, then the
  interleaved wall check. **And re-take the cohort merge's peak heap**, because every observation a
  merge holds pays the extra width whether it uses it or not. Merge only if the count halves *and*
  the merge's peak heap does not rise.
- **Complexity cost:** **large, and it lands outside the store.** `SequenceObservation` and
  `SampleLocusObservations` belong to `src/ng/locus_generation/`; every producer and consumer
  changes. The inline capacity trades struct width against spill frequency and has to be chosen
  against real allele-length distributions rather than guessed. The benefit is measurable inside the
  store; the cost is paid everywhere.

### Speculative

- **S1: [record.rs:918](../../../../src/ng/psp/record.rs#L918) — three scratch `Vec`s a record in
  `decode_record_body`.** They cost nothing at one observation, because `Vec::new` does not allocate:
  the measured 4.05 calls a record at 10 reads is 4.00 for single-observation records plus 0.05 for
  the one in a hundred with two. **Build the multi-observation corpus before touching this**; the fix
  needs a scratch parameter through two public signatures for a gain nothing has yet measured.
- **S2: [block.rs:587](../../../../src/ng/psp/block.rs#L587) — closing a block copies every record
  byte into a second buffer** to put the block head in front of them. Per block, not per record, and
  no profile names it. `push` documents a rollback property that a reserved prefix would have to keep
  holding.
- **S3: [chain_ids.rs:78](../../../../src/ng/psp/chain_ids.rs#L78) — the live set could hold `u32`
  offsets from a per-block base** and halve every byte the merge and `residual_reads` touch. Measured
  obstacle: on `hg002.ngpsp` the live set spans **10.8 identifier values per live identifier**, so it
  is not dense — which also rules out the bitset variant.
- **S4: [block.rs:1832](../../../../src/ng/psp/block.rs#L1832) — `mmap` instead of `read`. Do not do
  this**, and the measurement is why: reading the whole 55.85 MB store costs about **45 ms against a
  988 ms walk**, so the ceiling on removing the copy is 4.6 %, and touched pages of a mapping count
  against the process's resident set — so a cohort holding three thousand mapped stores would have an
  RSS that grows with how much of each it has walked, which is precisely the property §1.1 exists to
  forbid. Plus `unsafe`, a dependency, and SIGBUS on a truncated file turning a `PspReadError` into a
  process kill. Recorded so the next reader of `fill_from_source` finds the measurement rather than
  re-deriving the idea.

### Note

- **The varint lead was refuted by measurement, and the `#[cold]` split should stay.** I proposed
  that `src/psp/varint.rs`'s one-byte-only fast path was costing us, since several ng fields exceed
  127 and `decode_u64_leb128_cold` is 2.2–3.1 % of a walk. An ng-local codec inlining every width made
  all six walk benchmarks slower; one inlining the two-byte case as well still executed **0.67 % more
  instructions** and lost every interleaved round (cheapest 55.261 ms against 63.757 ms), while
  `read_varint`'s own attribution grew **137 %**. Removing the cold call does recover 769,899
  simulated mispredicts — 6.4 % of the run's whole mispredict budget — and pays 29,981,553 extra
  instructions for them. The branch in front of the cold call is decided by *which field* is being
  read, which is the same sequence at every record, so the predictor learns it. **No freeze exception
  for `src/psp/varint.rs` should be requested.**
- **The `&'static str` field name threaded through every `read_*` is free on the taken path.** The
  disassembly shows the name's pointer and length touched only inside the two error arms, with no
  `format!` and no allocator call. Keep the convention.
- **The 16 kB read chunk is the right size for the wrong reason.** Re-measured on the shipped reader,
  five buffer pairs alternating in one process over the 7.69 M-record tomato store: 16 kB, 64 kB and
  1 MB came out **within 1.6 % of each other, unordered**, against a 12.7 % run-to-run spread — while
  read syscalls fell from 3,413 to 58. `READ_CHUNK_BYTES`'s doc comment and spec §4.4 both carry the
  *prototype's* sweep ("64 kB is 13 % slower … 256 kB is 40 % slower") as though it described this
  reader. It does not. **The right justification is memory:** 64 kB would take an open sample from
  480 kB to 528 kB, past the budget, for no measurable time. Correct the comment; keep the constant.
- **ng reads each byte of a store exactly once** — 3,413 read syscalls for 55,851,989 bytes — and all
  3,413 together cost about 0.7 ms, 0.07 % of a walk. **None of ng's deficit is I/O.**
- **The `__fcntl` in the write profile is `PspWriter::finish`'s `sync_all`**
  ([writer.rs:485](../../../../src/ng/psp/writer.rs#L485)), one per finished store — the durability
  contract §6.3 asks for, not a defect. Measured on this host, one `F_FULLFSYNC` is **4.9 ms**. It is
  worth recording in the bench's doc, because at the write arm's old size it was about a third of the
  shallow arm's whole number.
- **Writing 400,000 records costs 11 `write(2)` calls** — one a block plus two — and bytes written
  equal the file size exactly. `BufWriter`'s 8 KiB default is not in the path: a block is larger than
  the buffer, so std writes it straight through. Nothing on the write path is per record.
- **Opening a psp costs 4 read syscalls, and they are the cheap part.** On a human reference it reads
  **269 kB and takes about 11.8 ms**, almost all of it parsing the header's TOML — twice, once as a
  bare table to read the version and once as the real type. Across three thousand samples that is
  about **35 seconds before the first record of the run is walked**, all of it the same contig list.
  Sharing it is spec §11's, ruled outside this module; the parse being done twice is not.
  `header.rs`'s comment justifying the second pass says "a header is about a kilobyte", which on a
  human reference it is not.
- **The rest of the module is already careful, and a review that lists only what it found misreports
  it.** `BlockBuilder` reuses three buffers and swaps rather than copies at a cut; `BlockCompressor`
  keeps its frame scratch across blocks; `LiveSetWriter`, `LiveSetReader` and `RecordEncoder` all
  keep and clear their scratch; `BlockStream` holds two fixed buffers and shrinks the rolling one
  back at every block. `header.rs`, `index.rs`, `footer.rs`, `trailer.rs` and `reader.rs` allocate at
  open and in error paths and nowhere else. **H3 and H4 are the two places the module's own reuse
  discipline was not applied — not a pattern running through it.**

### The finding that is not in any category, and the decision it puts to the owner

**§5.4's 1.8× holds at 280 reads a position and inverts at 10.** Both readers over the same records,
cheapest of nine interleaved rounds, macOS host, mimalloc, reproducible to 0.5 % across repeats:

| corpus | reads a record | production's `.psp` | ng's store | |
|---|---:|---:|---:|---|
| tomato SRR7279481 — 7,687,686 records | 10.3 | **0.369 s** (48.0 ns a record) | 0.821 s (106.7 ns) | ng **2.22× slower** |
| HG002 chr21 — 74,623 records | 280.0 | 0.0261 s (350 ns) | **0.0142 s** (190 ns) | ng **1.84× faster** |

**§5.4's walk is 62 tomato accessions** — 471,520,156 records is 7.6 M each, one accession's worth —
so the corner its claim was measured on is the top row.

**Where the difference is.** Production reads a **columnar** block — one tight loop a column, then
materialise — and a `--release` profile of its arm puts **43.6 %** of its time in zstd and about 24 %
in its own decode. ng reads a **row** — one field reader a record, walking about sixteen fields in
sequence — and the `--profile profiling` breakdown of its arm puts **13.4 %** in zstd and **26.3 %**
in field decoding alone (`read_varint` 21.2, the varint cold path 2.2, `read_u32` 1.7,
`read_locus_kind` 1.2). *Those two shares come from different builds — section 4 says why a
`--release` profile cannot break ng's walk down — so read them as two descriptions of where each
reader's work sits, not as one subtraction.* ng's file is *smaller* (55.9 MB against 62.9 MB), so it
decompresses less and still takes twice as long.

**Two things make the gap worse than it looks, not better.** Production's reader re-reads its own
file **27.7 times over** — 1,743,807,569 bytes read from a 62,946,874-byte file, because its
`BufReader` is discarded on every seek — where ng reads each byte once. And ng's advantage at 280
reads is measured against that same handicapped reader.

**None of H1–H4 closes the tomato gap.** They fix chain-id work, and on production-shaped stores the
live set averages **0.1 identifiers a record on tomato**. What is slow there is the row-shaped field
parser, and the only finding aimed at it — L1, the error width — is the one nobody measured.

**So the decision, and it is the owner's because it is a design question and not a coding one:** is
the row layout acceptable at shallow depth, given that it is where the cohort actually is? The
format's other wins are real and measured — 35 % fewer bytes, an index 11× smaller, 480 kB an open
sample against production's 2.6 MB, and the head-skip that a columnar layout cannot offer at all. My
recommendation is **keep the row layout and measure L1 first**: it is the cheapest experiment that
addresses the actual gap, it is contained inside `record.rs`, and 26.3 % of the walk is sitting in
the functions it touches. Revisit the layout only if L1 and H1–H4 together leave ng behind
production at 10 reads a position.

---

## 6. Out-of-scope observations

- **`src/psp/reader.rs:2275` — production's reader re-reads its file 27.7× on tomato and 217× on
  hg002.** `seek_to_offset` uses `SeekFrom::Current` intending to preserve the `BufReader`'s buffer,
  but `BufReader`'s `Seek::seek` discards the buffer for `Current` as well; only `seek_relative`
  preserves it. This is the 11.1 % `read` line in production's profile. **`src/psp/` is frozen, so
  this is not a proposal** — it is recorded because every ng-versus-production ratio in this review is
  taken against a reader carrying it.
- **`src/ng/psp/header.rs:608` and `:628` — the header body is parsed as TOML twice per open**, and
  the comment justifying the second pass says a header is about a kilobyte. On a human reference it
  is 269 kB and the two parses are ~11.8 ms a sample. Follow-up: separate PR, and restate the comment
  whatever is decided.
- **`src/ng/psp/block.rs:1732` — `pump` memmoves the rolling buffer's tail to offset zero on every
  refill.** `_platform_memmove` is 5.5 % of the shallow full walk and 7.4 % of the deep one, but the
  same bytes pass through `pump` in the head-only walk where it is 0.9 % — so most of that share is
  record building, not this. A head offset that compacts only when the tail would not otherwise fit
  would remove what is left. Low priority; measure before acting.
- **`src/ng/psp/block.rs:1437` — `vec![0u8; READ_CHUNK_BYTES]` zeroes 16 kB that `fill_from_source`
  overwrites before anything reads it.** A `memset` a walk, not an allocation. The fix needs `unsafe`
  or the `read_buf` API. Small either way.
- **`src/ng/psp/reader.rs:116` and `trailer.rs:230` — `try_clone` for the header is a `fcntl` and a
  `close` an open sample** that the file's own comment two lines above argues against ("a cohort opens
  thousands of these, and a second `open(2)` per sample buys nothing"). Taking `&mut File` removes
  both. **Below what this environment can resolve — merge for tidiness, not for speed.**

---

## 7. What's already good

- **The head really is a skip, and the measurement proves it rather than asserting it.** A walk
  declining every body costs **0.020 allocator calls a record at 10 reads and 0.101 at 280** — and
  those are one walk's fixed setup spread over its records, not a per-record cost: the same ~590 calls
  appear whether the walk covers 30,000 records or 6,000
  ([`examples/dhat_ng_psp.rs`](../../../../examples/dhat_ng_psp.rs)). The design's central claim, that
  a cohort's first pass can advance past a body without building it, is structurally true in the code.
- **Scratch buffers are reused everywhere the module reaches for one, and swapped rather than
  copied.** `BlockBuilder` keeps `records`, `next_records` and `closed_block_payload` across blocks
  and swaps at a cut; `BlockCompressor` keeps its frame scratch; `RecordEncoder` keeps `body_scratch`
  ([block.rs:431](../../../../src/ng/psp/block.rs#L431),
  [record.rs:1484](../../../../src/ng/psp/record.rs#L1484)). H3 and H4 are the exceptions that prove
  the rule.
- **Nothing a reader holds is a function of the data, and it is enforced by construction.**
  `BlockCursor`'s two constructors are exhaustive literals, so a per-block field added beside them is
  a compile error rather than a field that silently never resets; `begin_next_block` shrinks the
  rolling buffer back to its budget after an outsized record
  ([block.rs:1303](../../../../src/ng/psp/block.rs#L1303),
  [block.rs:1674](../../../../src/ng/psp/block.rs#L1674)). This is what makes the 500 kB-an-open-sample
  budget a property of the code and not of the corpus.

---

## Author response — applied 2026-08-30, same branch

| finding | outcome |
|---|---|
| H1 — the live set rewritten end to end | **applied**, `1fda32e2`, after a three-arm head-to-head |
| H2 — the writer's one-observation record | **applied**, `be9955b5` |
| H3 — the residual list grown from nothing | **applied**, `7715bb92` |
| H4 — the writer's three per-record lists | **applied**, `be9955b5`, together with H2 |
| L1 — the error type's width | **experiment shows no gain — closing**, see below |
| L2 — the division in `read_count` | **experiment shows no gain — closing**, `830344fa` |
| L3 — `decode_read_list` not reserving | **applied**, `7715bb92` |
| L4–L7, S1–S4 | **deferred** — each needs a workload or a corpus the bench still lacks |

**H1 was two designs and neither had been compared.** Both were built into one binary with the
shipped pair, all three arms rotating round by round inside one process, cheapest round of each,
three repeats, every arm asserted to read identical identifiers. At 280 reads a position the
in-place design took the head-only walk to **0.567×** and the full walk to 0.762×, against the
single-merged-pass design's 0.851× and 0.920×; at 10 reads both were within 1 % of the shipped
pair. The in-place design also deletes the scratch buffer, giving **2,264 bytes an open sample**
back. It is the one that shipped.

**What the three applied fixes are worth**, against this report's own baseline, cheapest of three
complete passes either side, on an otherwise idle machine:

| workload | before | after | |
|---|---:|---:|---:|
| head-only walk, 280 reads | 19.935 ms | 11.030 ms | **0.553×** |
| one body in a hundred, 280 reads | 20.187 ms | 11.338 ms | **0.562×** |
| write, 280 reads | 101.419 ms | 57.700 ms | **0.569×** |
| full walk, 280 reads | 36.571 ms | 26.356 ms | 0.721× |
| write, 10 reads | 38.795 ms | 26.320 ms | 0.678× |
| full walk, 10 reads | 33.561 ms | 30.501 ms | 0.909× |
| head-only walk, 10 reads | 7.769 ms | 7.936 ms | 1.021× |
| one body in a hundred, 10 reads | 8.318 ms | 8.349 ms | 1.004× |

Allocator calls a record fell from 11.04 to 4.13 on the deep full walk and from 8.97 to **0.09** on
the deep write — a writer now calls the allocator about once every eleven records.

**And the headline comparison did not move, exactly as section 2 predicted.** Re-run after all
three fixes: tomato **0.452×** against 0.450× before, HG002 1.96× against 1.84×. The corpora are
built from a production `.psp` whose live set averages 0.1 identifiers a record on tomato, so there
was never anything on that path for these fixes to save. **The 2.2× gap at 10 reads a position is
untouched and remains the open question.**

**L1 was the one candidate aimed at that gap, and it fails.** Splitting the varint fault into a
`#[cold] #[inline(never)]` constructor does exactly what was predicted structurally — `read_varint`
loses its out-of-line copy entirely and only the cold handler remains — but the trade is bad. Two
binaries differing only in that split, run alternately, three paired rounds each:

| arm | shipped | cold-fault split | |
|---|---:|---:|---|
| full walk, 10 reads | 30.247 ms | 29.557 ms | 0.977×, ahead in 3 of 3 |
| **head-only walk, 10 reads** | **8.351 ms** | **8.738 ms** | **1.046×, behind in 3 of 3** |
| full walk, 280 reads | 27.236 ms | 28.395 ms | no clean signal |

Inlining sixteen field reads into the body path bloats code that the head path pays for and does
not use — and the head-only walk is the cohort's first pass and the whole point of the format's
head design. **Regressing it 5 % to gain 2.3 % on the full walk is the wrong way round**, so the
split was reverted.

**One methodological result worth keeping from L1 and L2.** Instruction count was the wrong gate
for both. L1 executed **0.073 % more** instructions and ran **2.3 % faster** on the full walk — the
frame's loads and stores are what cost, and callgrind counts them without weighting them. L2 saved
721,200 instructions out of 6,762,760,740 — 0.011 % — and was not worth the arithmetic identity it
put inside a guard. **Gate on a count when the mechanism is a count (allocator calls, syscalls);
gate on interleaved wall clock when the mechanism is the shape of a frame.**

## Author response convention

Address each finding by its identifier (H1, L2, …) with one of: `applied in <commit>` /
`experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` /
`won't fix because …`. The "experiment shows no gain" path is expected and welcome — three of this
review's own findings closed that way.

---

*Per-category findings, the profiles, and every harness written during this review are left as an
audit trail under `tmp/perf_review_2026-08-30_ng-psp/`.*
