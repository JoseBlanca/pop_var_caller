# What one contributor visit costs the ng generic walk

Measured on `6fbbd09` (round 3's last commit), host-native release builds, `instructions
retired` from `/usr/bin/time -l`, floor-subtracted, min of 3 runs. Wall clock is never quoted:
three other agents were measuring on this host throughout.

Vocabulary used below, defined once:

- a **column** is one covered reference base the walk processes — one turn of
  `WalkerState::process_position`;
- a **contributor visit** is one (column × read) pair that survives the mate-overlap collapse
  and the per-column depth cap — the unit the fold loop iterates over;
- **depth** here means contributor visits per column, which is what the walk actually pays
  for, not the sample's nominal coverage.

Two fixtures, and one of them turned out not to be what its name suggests:

| fixture | what it is | columns | visits | depth |
|---|---|---:|---:|---:|
| tomato `DRR000741.p1.cram` `SL4.0ch01`, first 1 M loci | CRAM, ~130× nominal | 1,449,893 | 123,613,443 | 85.3 |
| HG002 `HG002_TR_v1.0.1_Tier_30x.bam` `chr1`, whole | BAM, 30× nominal | 2,499,261 | 43,084,509 | 17.2 |

**Two instructions in five on the 30× fixture are spent getting to the reads, not walking
them.** Its BAM is restricted to the GIAB tandem-repeat tiers, so walking `chr1` means
traversing 613,682 typed generic regions to reach 1,541,788 covered bases — **2.5 covered bases
per region**. Running the same BAM with `PVC_PROBE_WHOLE_CONTIG=1`, which replaces the region
stream with one region per contig, costs 32,325 instructions per column at depth 20.6 against
46,414 per column at depth 17.3 through the region stream. Correcting the first for the depth
difference leaves **18,480 instructions per column of pure region overhead — 40 % of that
fixture's instructions, about 68,600 per generic region traversed**. Every 30× number below
states which mode it came from.

---

## The per-visit price: the arithmetic estimate was right, and there is 3–6× of headroom

**A contributor visit costs about 1,300–1,800 retired instructions, and the number barely moves
with depth.** The estimate from source that put it near 1,800 is confirmed; the hand-count of
the work actually done — 300–600 instructions — is what the walk is not achieving.

Marginal measurements, which cancel start-up and reference loading by differencing two runs of
the same fixture at two loci counts:

| fixture and mode | Δ instructions | Δ visits | depth | **instructions per visit** |
|---|---:|---:|---:|---:|
| tomato, region stream, 1 M → 2 M loci | 210.465 G | 118,793,400 | 82.5 | **1,772** |
| tomato, whole contig, 1 M → 2 M loci | 180.123 G | 98,178,421 | 97.8 | **1,835** |
| HG002, whole contig, 770 k → 1.54 M loci | 25.078 G | 16,002,805 | 20.6 | **1,567** |
| HG002, region stream, 770 k → 1.54 M loci | 58.047 G | 21,611,380 | 17.3 | 2,686 |

The last row is the contaminated one — the extra 1,119 instructions per visit over the row
above it is region setup, not walking.

Splitting the price into a depth-independent per-column term and a per-visit term needs two
depths on the same reference, the same contig and the same file format. Walking the sparse
tomato bench CRAM (`SRR5079860.p1.bench.cram`, depth 4.2) and the big one (depth 97.8) over
`SL4.0ch01` in whole-contig mode gives:

```
sparse:  12,667 instructions per column at depth 4.16
big:    179,409 instructions per column at depth 97.79
        → per column  5,260 + per visit 1,781 × depth
```

Feeding the HG002 BAM's whole-contig marginal (32,325 instructions per column at depth 20.6)
through the same 5,260-instruction column term implies **1,312 instructions per visit on a
BAM**, against 1,781 on a CRAM. The gap is the input format: CRAM record decoding is 26.0 % of
the profiling build's main thread (`decode_container_at`, 6,951 of 26,766 samples), a BAM pays
BGZF inflation instead.

So the walk pays, per covered reference base:

- **≈ 5,300 instructions** that do not depend on depth — opening the record, reading the
  reference byte, closing the record, building the emitted locus;
- **≈ 1,300 (BAM) to 1,800 (CRAM) instructions for every read at that base.**

At the 30× target that is 5,300 + 17 × 1,300 ≈ 28,000 instructions to produce one locus from
seventeen reads.

**Which side of the estimate-versus-hand-count question this settles: the estimate.** The
per-visit price really is near 1,800, and the work a visit performs — pull at most two events
from a cursor, compare one base against one allele bucket, add six scalars into a struct,
insert one hash-map entry — is not 1,800 instructions. There is a factor of three to six here,
and it is not hiding in cache misses (see the working set below).

---

## The two sorts: a tenth of the walk at 130×, a twentieth at 30×

Both sort a depth-sized list once per covered base:
`resolve_mate_overlap_at_pos` sorts one `(chain_id, index)` tuple per contributor
(`genome_walk.rs:1176`), and `keyed_observations_counting` sorts one `(read_id, &state)` pair
per folded read (`open_record.rs:763`).

### How much they sort

Counted exactly, over the whole run:

| | tomato ~130× | HG002 30× |
|---|---:|---:|
| mate-overlap sort: lists sorted | 1,449,893 | 2,499,261 |
| mate-overlap sort: elements, total | 123,613,482 | 43,084,689 |
| mean / p50 / p90 / p99 / max | 85.3 / 95 / 121 / 139 / 192 | 17.2 / 16 / 33 / 45 / 67 |
| close-time sort: lists sorted | 1,442,402 | 2,424,871 |
| close-time sort: elements, total | 123,408,535 | 42,998,240 |
| mean / p50 / p90 / p99 / max | 85.6 / 95 / 121 / 139 / 192 | 17.7 / 17 / 33 / 45 / 67 |

The two sort essentially the same number of elements: every contributor at a column folds into
the one record opened there, and that record closes one step later holding what folded into it.
The distribution is left-skewed at 130× — the mean of 85.3 is below the median of 95 because
low-coverage stretches pull the mean down, and no column exceeded 192 elements in the 1 M-locus
prefix.

### What they cost, in release instructions

Priced by running each sort **twice** over the same unsorted input (a copy is taken before the
real sort and sorted afterwards, then discarded) and differencing against the default build.
The doubled close-time sort also clones its `Vec`, so its number carries one extra allocation
per record, subtracted at 340 instructions per allocation where noted.

| | tomato ~130× | HG002 30× |
|---|---:|---:|
| mate-overlap sort, per column | 7,613 | 917 |
| close-time sort, per record (clone subtracted) | 6,780 | 1,178 |
| per element sorted (mate / close) | 89 / 79 | 53 / 66 |
| **both sorts, share of the walk** | **9.6 %** | **4.7 %** |
| both sorts, per contributor visit | 172 of 1,748 | 139 of 2,549 |

Note the shape: the cost per element sorted grows with the list, so the sorts are a **depth**
term that gets worse as coverage rises — 9.6 % at 130×, 4.7 % at 30×, and the ratio would
widen further at 300×.

### Which profile lines they are

From a `[profile.profiling]` build with an `#[inline(never)]` wrapper on each sort, 3 M loci of
the big tomato CRAM, 26,766 main-thread samples (raw output:
`sample_sort_attribution.txt`). Inclusive counts:

| site | samples | share |
|---|---:|---:|
| the close-time sort (`sort_the_folded_reads`) | 1,868 | 6.98 % |
| the mate-overlap sort (`sort_the_mate_overlap_list`) | 1,815 | 6.78 % |
| the `collect()` that builds the list the close-time sort orders | 602 | 2.25 % |

Mapping those onto the profile in the brief:

- **`small_sort_general` 4.4 % is the close-time sort.** Its named instantiation
  (`hfeae6316d89cb7ee`, 1,201 self samples = 4.49 % here) is reached only through
  `sort_the_folded_reads`.
- **`Vec spec_from_iter` 2.1 % is the `folded` collect** in `keyed_observations_counting`
  (`SpecFromIterNested::from_iter::h7544d…`, 563 self samples = 2.10 % here) — not a CRAM
  buffer, as its neighbours in the ranking might suggest. The *large* `from_iter` in this
  binary (2,279 samples inclusive) is a different instantiation and is CRAM block decode.
- **`quicksort` 1.9 % is both sorts**, in two instantiations: `h779456c8…` (498 self) under the
  close-time sort, `h83b728ce…` (414 self) under the mate-overlap sort.
- **The mate-overlap sort's `small_sort_general` is inside `<deduplicated_symbol>`** — see the
  next section.

---

## How often the general machinery is needed: almost never

Per column, on both fixtures, as natural frequencies. Everything here is a **depth** or
**breadth** observation, not a format one; the two fixtures agree closely, which is the point.

| the column… | tomato ~130× | HG002 30× |
|---|---:|---:|
| touches more than one record | **0 in 10,000** | **0 in 10,000** |
| has a record that widens | 14 in 10,000 | 10 in 10,000 |
| has some contributor carrying an insertion or deletion | 24 in 10,000 | 16 in 10,000 |
| has an affected record wider than one reference base | 33 in 10,000 | 30 in 10,000 |
| has two contributors sharing a chain id (mates overlap) | 1,048 in 10,000 | 1,053 in 10,000 |
| has a contributor carrying a mate-overlap reconciliation mark | 1,048 in 10,000 | 1,052 in 10,000 |

Per fold, and per contributor visit:

| | tomato ~130× | HG002 30× |
|---|---:|---:|
| the fold reuses `events_at_pos` instead of re-querying the cursor | 9,876 in 10,000 | 9,511 in 10,000 |
| the fold is a re-fold into a record the read was already in | 17 in 10,000 | 20 in 10,000 |
| the visit carries a reconciliation mark | 27 in 10,000 | 137 in 10,000 |
| the visit carries an indel event | 2 in 10,000 | 5 in 10,000 |

Three of these are worth reading twice.

**No column in either run touched more than one open record.** 1,449,893 columns on tomato and
2,499,261 on HG002, and `affected` had length 0 or 1 every time: 0.997 and 0.972 affected
records per column. The `Vec<u32> affected` allocated fresh at `open_record.rs:2325`, the
`affected.contains(&key)` membership scan and the `affected.sort_unstable()` that follows are
all machinery for a case that did not occur once.

**Mate overlap is the one general case that is genuinely common** — a column in ten has two
contributors from the same fragment. But it marks only 27 contributor visits in 10,000 at 130×
(137 in 10,000 at 30×), because a shared chain id in a column of 85 reads still involves only
two of them.

**Everything else is rarer than 35 columns in 10,000.** Taking widening, indels and
wider-than-one-base records as if they never co-occurred — the pessimistic reading — a column
that needs none of them is **at least 9,929 in 10,000 at 130× and at least 9,945 in 10,000 at
30×**. Re-folds are rarer still, at 17 and 20 folds in 10,000.

---

## The unattributed 8.7 % line: three merged functions, none of them the fold table

`<deduplicated_symbol>` is what `sample` prints when the linker folds two or more identically
compiled functions onto one address, so the symbol is ambiguous. Resolving the addresses
against `nm` on the profiling binary (load slide 0x7b4000) splits the 2,676 self samples in this
build three ways:

| what it actually is | samples | share of the line | reached from |
|---|---:|---:|---|
| `core::slice::sort::shared::smallsort::small_sort_general`, two instantiations merged | 1,321 | 49 % | the **mate-overlap sort** |
| `hashbrown::map::HashMap::retain` and `::remove`, merged | 1,060 | 40 % | `ChainIdAllocator::allocate_for_read` |
| `core::hash::BuildHasher::hash_one`, two instantiations merged | 295 | 11 % | noodles CRAM `Byte::decode_take` |

**It is not `OpenPileupRecord::folded_reads`, and it is not `ActiveReads::by_read_id`.** Half
the line is sort machinery already counted in the sorts above. The hash-map half is the chain-id
allocator's `pending_mates`, and specifically `evict_stale_pending`
(`chain_id_allocator.rs:354`), which runs `retain` over the **whole** pending-mate map on
**every read admission**: 954 of the 1,060 samples sit there, 3.6 % of the walk. That is an
O(reads × pending-map size) term — it grows as depth², and at 130× the map holds thousands of
entries. It is in `chain_id_allocator.rs`, which `copy_fidelity.rs` still holds byte-identical
to production, so it is reported here rather than proposed as a change.

The last eleventh is inside noodles' CRAM codec — input format, not the walk.

---

## The working set: it fits in L1, so cache misses are not the explanation

Sizes taken from the build (`size_of`), buffer lengths measured over the runs: tomato reads are
96 bases with 1.06 CIGAR operations each, HG002 reads 148 bases with 1.20.

This host: performance cores have **128 kB of L1 data cache** and share **16 MB of L2**
(`hw.perflevel0.l1dcachesize` 131072, `hw.perflevel0.l2cachesize` 16777216); the efficiency
cores have 64 kB and 8 MB.

Bytes one column touches, at depth D:

| | 30× (D = 17.2) | 130× (D = 85.3) | deepest column seen (D = 193) |
|---|---:|---:|---:|
| active-set `ActiveRead` array (184 B each) | 3.1 kB | 15.3 kB | 34.7 kB |
| cursor offset tables (17 B payload, its own 64-B line each) | 1.1 kB | 5.3 kB | 12.1 kB |
| `seq` + `bq_baq` (one byte read from each, two lines per read) | 2.2 kB | 10.7 kB | 24.1 kB |
| `contributors_buf` (104 B each) | 1.8 kB | 8.7 kB | 19.6 kB |
| mate-overlap `(chain_id, index)` list (16 B each) | 0.3 kB | 1.3 kB | 3.0 kB |
| `folded_reads` map (4-B key + 88-B value + control byte) | 1.6 kB | 7.7 kB | 17.5 kB |
| the `(read_id, &state)` list `finalise` sorts (16 B each) | 0.3 kB | 1.3 kB | 3.0 kB |
| **total** | **10.2 kB** | **50.4 kB** | **114.0 kB** |

**At the 30× target the whole column fits in one eighth of L1.** At 130× it is 39 % of L1, and
even the deepest column in the 1 M-locus prefix stays inside it at 89 %. Nothing here reaches
L2, let alone memory.

That closes the question the brief raised: a working set that missed cache would explain a
per-visit price far above the instruction count of the work, and this one does not. It also
explains why the three struct-narrowing attempts in earlier rounds measured null, and it
predicts the same for a fourth. **The 1,300–1,800 instructions per visit are instructions being
executed, not stalls.**

---

## Verdict

**There is headroom of the kind a structural change reaches, it is roughly a factor of three,
and it is in the per-column fixed cost and the two depth-sized sorts — not in the fold itself
and not in memory behaviour.**

A contributor visit costs 1,300 instructions on a BAM and 1,800 on a CRAM, against 300–600 for
the work a visit performs. The gap is not cache: the entire column, at either depth, sits inside
a 128-kB L1 data cache, and at the 30× target it uses an eighth of it. So the gap is executed
instructions, and the census says where they are. The two sorts are 9.6 % of the walk at 130×
and 4.7 % at 30×, and they are a growing term — 89 instructions per element sorted at depth 85
against 53 at depth 17, sorting what is essentially the same list twice per base because the
mate-overlap pass and the close pass each build their own. The chain-id allocator sweeps its
whole pending-mate map on every read admission, a further 3.6 %, and that term grows as the
square of depth. And beneath all of it sits the shape the frequencies make plain: **the walk
runs its general machinery on every column, and the general case occurs on fewer than 71
columns in 10,000.** No column in four million touched two records; 9,876 folds in 10,000
already skip the cursor and reuse the events the walk fetched. A path specialised for "one base,
one record, no indel, no re-fold" would cover at least 9,929 columns in 10,000 at 130× and 9,945
in 10,000 at 30×, and it is the one structural change these numbers argue for, because the two
things it can remove — the per-column 5,300-instruction fixed term and the two sorts of a list
that is always the same list — are between them a third of the walk.

---

## What was left in the tree, and what it cost

Four cargo features, all off by default, patch at
`walk_census_instrumentation.patch` (452 added lines, no line removed):

- `walk-census` — the counters. Every call site is `#[cfg]`-gated, and the release binary built
  **without** the feature has a **byte-identical `__text` section** to the one built before any
  of this was written, so the measured build provably contains none of it. The census build
  itself costs 1.8 % (tomato 217.396 → 221.285 G; HG002 110.170 → 112.181 G) and its instruction
  count is never quoted.
- `walk-sort-attribution` — `#[inline(never)]` wrappers on the two sorts, for the profiler.
- `walk-sort-doubling-mate` / `walk-sort-doubling-close` — each sort run a second time over a
  saved copy, which is how the sorts were priced in release instructions.

Gates, all passing on the default build:

- four dumps `cmp`-identical to the stored baselines (generic and SSR × chr21 and `SL4.0ch01`);
- probe counters exact on chr21 — `loci=236081 observations=251786 reads_admitted=54709` — from
  the census build as well as the default one;
- `cargo test --lib` 2,882 passed, 0 failed; `cargo test --examples` 33 targets ok;
- `cargo clippy --all-targets --all-features -- -D warnings` clean;
- `cargo doc --no-deps` 12 unresolved links, the baseline.

Raw outputs beside this file: `census_tomato_130x_1Mloci.txt`,
`census_hg002_30x_chr1.txt`, `sample_sort_attribution.txt`.

### One caveat on a number carried in from the brief

The start-up floor quoted for HG002 (1.900 G) does not reproduce: `PVC_PROBE_MAX_LOCI=1` on
`chr1` measures **0.349 G** here, three runs agreeing to two parts in a thousand. The tomato
floor does reproduce (1.313 G against the quoted 1.306 G). Every HG002 figure above uses the
measured 0.349 G; using 1.900 G instead would move the per-visit price by 1.4 %, which changes
no conclusion. The absolute walk totals also sit below round 3's — 217.4 G against 233.3 G on
tomato at 1 M loci — on the same commit; all comparisons here are internal to this session and
use one baseline throughout.
