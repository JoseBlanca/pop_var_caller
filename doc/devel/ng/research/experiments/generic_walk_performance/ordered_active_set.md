# The walk's lost order: three stages, measured

**Worktree:** `/Users/jose/devel/pop_var_caller/.claude/worktrees/agent-ab705bfff8c2f2e1d`, detached at
`6fbbd093764662ed2496acde39424c8ee234ea1c`. Both required markers confirmed
(`spans_only_its_anchors` × 8, `fn finalise_recycling` × 1) before anything was measured.
The tree is left with all three stages applied.

---

## The answer

**Yes, and by a lot — but the win is entirely in stage 3, and stage 3 only works because
stage 2 exists.**

Walk instructions retired, start-up floor subtracted, minimum of three runs a side:

| | 30× HG002 chr1 | ~130× tomato `SL4.0ch01` | 300× HG002 chr1 |
|---|---:|---:|---:|
| baseline `6fbbd09` | 108.541 G | 215.982 G | 382.564 G |
| **all three stages** | **99.522 G** | **191.453 G** | **331.694 G** |
| | **−8.31 %** | **−11.36 %** | **−13.30 %** |

Ranges disjoint on all three fixtures. The gain grows with depth, which is the shape the
previous attempt at this change did not have.

**The one measurement that settles the hypothesis.** The same fold-table rewrite, applied
*without* the ordered active set — so contributors still arrive in the `swap_remove`
permutation — reproduces the 2026-08-04 revert almost exactly:

| fold table as an ordered `Vec` | 30× chr1 | ~130× tomato | 300× chr1 |
|---|---:|---:|---:|
| with contributors **permuted** (stage 1 + stage 3) | −3.18 % | **+12.92 %** | **+23.15 %** |
| with contributors **ascending** (stages 1+2+3) | −6.46 % | −8.75 % | −10.31 % |

The historical record for the reverted attempt was *−3.2 % at 30×, +16.4 % at 300×*. The
left-hand row here is −3.18 % at 30× and +23.15 % at 300×. **The regression was caused by
the arrival order and by nothing else**, and removing that cause turns a 23 % loss at 300×
into a 10 % win. (Percentages in this table are against stage 1, which is the common base
of both rows.)

**The cost is not zero and is not only in instructions.** Stage 2 reorders the contributor
list, and three things downstream read that order:

- the `record_widen_events` run counter moves by ±2 on the gate fixtures (423 → 425 on
  chr21, 622 → 621 on tomato) — **every emitted locus line is identical**;
- two `parity::` tests that compare ng against production's walker now fail, because the
  mate-overlap tie-break's last resort is contributor index;
- **at 300×, where the per-column depth cap fires, the emitted loci genuinely differ** —
  measured, not inferred: 88,351 of 341,094 rows change on a chr21 300× dump.

The 300× divergence is the one that needs an owner decision. Sections 4 and 6 give it in
full.

---

## 1. What was measured, and with what

`instructions retired` from `/usr/bin/time -l`, minimum of three runs a side, floors
subtracted. Wall clock is not used anywhere in this report: three other agents were
measuring on the same six performance cores.

Floors (`PVC_PROBE_MAX_LOCI=1`), measured here rather than carried over where new:

- tomato `SL4.0ch01` — 1.306 G (given in the brief, used as given)
- HG002 chr1 30× — 1.900 G (given in the brief, used as given)
- HG002 chr1 **300×** — 0.355 G (**measured here**, from a single run:
  `355413950 instructions retired`)

`PVC_TRUST_REFERENCE_INDEX=1` on every run, so no run includes the FASTA checksum. Probe
counters on chr21 are exact and unchanged at every stage:
`loci=236081 observations=251786 reads_admitted=54709`.

**The 300× fixture exists and was used**: `benchmarks/ssr_hg002/bam/300x/HG002_TR_v1.0.1_Tier_300x.bam`
(2.7 GB, with its `.bai`), walked over chr1 with `PVC_PROBE_MAX_LOCI=1000000` — 1 M loci,
2.28 M reads admitted, ~13.5 s a run. So the depth sweep is 30× / ~130× / 300× and not two
points.

---

## 2. Stage 1 — the contributor carries its read

`ReadContribution` gains `active_index: u32`, filled from `iter().enumerate()` in
`WalkerState::process_position`; the fold indexes the active set with it instead of hashing
`read_id` through the secondary `read_id → index` map.

| fixture | before (stage 0) | after | change |
|---|---:|---:|---:|
| 30× HG002 chr1 | 108.541 G | 106.397 G | **−1.98 %** |
| ~130× tomato | 215.982 G | 209.802 G | **−2.86 %** |
| 300× HG002 chr1 | 382.564 G | 369.840 G | **−3.33 %** |

Raw, in the order run (instructions retired, floor **not** yet subtracted):

```
tomato 130x  baseline: 217782241927 / 217287764830 / 217687000994
tomato 130x  stage 1 : 211107961349 / 211460433481 / 211952306501
chr1   30x   baseline: 110449948476 / 110640570368 / 110441216583
chr1   30x   stage 1 : 108571181151 / 108296951553 / 108320930475
chr1  300x   baseline: 383149934480 / 383124947055 / 382919107865
chr1  300x   stage 1 : 370350423535 / 370194687429 / 370407758680
```

Ranges disjoint on all three.

**Gates: fully clean.** All four dumps `cmp`-identical to the stored round-3 copies
(251,792 / 4,406 / 1,718,914 / 11,945 lines); probe counters exact; `cargo test --lib`
**2,882 passed, 0 failed** — the same count as the baseline tree.

**This stage stands alone.** It is byte-identical, it keeps every test green, and it is
worth −2 to −3 % on its own. If the owner declines stages 2 and 3, this one should still
land.

**The invariant it promotes.** A contributor's `active_index` is only valid while the
active set has not been touched since the list was built. It was already true that nothing
admits or expires between the two — both calls are in `WalkerState::process_position` and
the walk's `admit_read` / `expire_passed_reads` sit outside it — but it was previously a
fact nobody depended on, and it is now load-bearing. `ActiveReads::at` panics rather than
returning `None`, so a violation stops the walk instead of silently dropping a read.

Diff: `stage1.diff` (118 lines).

---

## 3. Stage 2 — the active set stops being a bag

The brief proposed a slab plus an ordered side list. **That shape was built first and it
regressed**, so it is reported and discarded:

| shape, on top of stage 1 | 30× chr1 | ~130× tomato |
|---|---:|---:|
| slab (`Vec<Option<ActiveRead>>` + free list + `Vec<LiveRead>`) | **+1.05 %** | **+1.59 %** |
| `VecDeque<ActiveRead>`, no expiry guard | +0.77 % | +0.68 % |
| `VecDeque<ActiveRead>` + `min_alignment_end` guard — **adopted** | +0.03 % | +0.35 % |

The slab pays a dependent load on every access to a read, and reads are accessed about 130
times per covered base at ~130× — the contributor loop plus one `at()` per (contributor ×
affected record), so roughly 260 M extra indirections per million loci. What it removes is
one hash insert per admission and one hash remove-plus-insert per expiry, about 1.5 M of
each over the same million loci. It is the wrong trade by two orders of magnitude in call
count, and no amount of tuning fixes that.

**What was adopted instead.** `ActiveReads::reads` becomes a `VecDeque<ActiveRead>`:

- `push_back` on admission keeps ascending `read_id` — the push *is* the ordering
  guarantee;
- expiry removes in place with `VecDeque::remove`, which shifts whichever side is shorter.
  Reads leave in very nearly the order they arrived (`read_id` ascends with alignment start,
  and a read's end tracks its start), so the removal point is at or near the front and
  `remove(0)` only steps the head pointer — the reads themselves are not moved;
- the `read_id → index` hash map is deleted; `get_by_read_id` is a binary search;
- a new `min_alignment_end` field short-circuits `expire_passed` entirely at positions
  where nothing can expire. That guard is independent of the ordering and is what takes the
  cost from +0.77 % to +0.03 % at 30× — at that depth roughly three positions in four
  expire nothing.

`resolve_mate_overlap_at_pos` also changes `contributors.swap_remove(idx)` to
`contributors.remove(idx)`, because `swap_remove` would undo the ordering for the whole
column. That path fires only on a mate overlap with an indel on one side.

The doc comment claiming *"Iteration order is admission order"* is fixed — it is now true.

| fixture | stage 1 | stage 2 | change |
|---|---:|---:|---:|
| 30× HG002 chr1 | 106.397 G | 106.431 G | **+0.03 %, ranges overlap — neutral** |
| ~130× tomato | 209.802 G | 210.542 G | **+0.35 %** |
| 300× HG002 chr1 | 369.840 G | 374.088 G | **+1.15 %** |

Raw:

```
tomato 130x  stage 2 : 211848261861 / 212294571700 / 212236440926
chr1   30x   stage 2 : 108394715944 / 108330727903 / 108349023916
chr1  300x   stage 2 : 375233783602 / 374443089653 / 374496271103
```

At 30× the stage-1 range (106.397–106.671 G) and the stage-2 range (106.431–106.495 G)
overlap, so that row is a null result and is reported as one.

**Stage 2 is a cost, not a win. Its whole value is what it makes possible in stage 3**, and
on its own it should not be applied.

### Gates: two dumps move by one line each

The whole of both differences:

```
5c5
< # record_widen_events=423 column_depth_truncations=0 ... loci_emitted=236081
> # record_widen_events=425 column_depth_truncations=0 ... loci_emitted=236081

5c5
< # record_widen_events=622 column_depth_truncations=0 ... loci_emitted=1711775
> # record_widen_events=621 column_depth_truncations=0 ... loci_emitted=1711775
```

**One header line each; all 251,792 and 1,718,914 locus lines are identical.** The two SSR
dumps are `cmp`-identical. `record_widen_events` counts records that grew, excluding fresh
opens, so whether a record is opened at full width or opened narrow and then widened
depends on which contributor's event is processed first. Reordering contributors moves it.
Nothing else in the walk reads that counter.

### Two tests turn red, and they are the honest cost

```
test ng::locus_generation::pileup::parity::ng_agrees_with_production_where_production_fabricated_nothing ... FAILED
test ng::locus_generation::pileup::parity::every_divergence_from_production_is_one_of_the_six_named_classes ... FAILED
test result: FAILED. 2880 passed; 2 failed; 5 ignored
```

The failure, quoted from the first (fields elided for width):

```
seed 0x5eed0001 case 7: locus 5 ...
  left:  SequenceObservation { bases: [65, 84], ..., num_obs: 1, q_sum: -4.605170185988092, mapq_sum: 32, ..., chain_ids: [9] }
  right: SequenceObservation { bases: [65, 84], ..., num_obs: 1, q_sum: -4.0,               mapq_sum: 23, ..., chain_ids: [10] }
```

Same allele, same observation count, a different read behind it. The cause is
`pick_agree_keeper` / `pick_overlap_loser` in `genome_walk.rs`: after BQ and first-of-pair
and `alignment_start` all tie, the last resort is *which contributor has the smaller index
in the list* — which is exactly what stage 2 changes. ng could reproduce production's
choice only by reproducing production's `swap_remove` permutation, which is what it was
doing.

On real data this does not fire: a genuine mate pair has opposite `mate_role`, so
`pick_agree_keeper` resolves at that step and never reaches the index comparison. The
synthetic parity fixtures build reads that tie all the way down. **This is a divergence
from production that is not one of the six named classes, and naming it (or making the
tie-break order-independent) is a spec decision, not a performance one.**

`parity::ng_emits_the_same_bytes_in_a_second_process` — the determinism test — **passes**,
run explicitly.

Diff: `stage2.diff` (cumulative, stage 1 + stage 2). The two rejected shapes are
`stage2_slab.diff` and `stage2_queue_noguard.diff`, both `active_read_set.rs`-only and both
applying over stage 1.

---

## 4. Stage 3 — the fold table in `read_id` order, and the sort deleted

`OpenPileupRecord::folded_reads` becomes `FoldedReads`, a newtype over
`Vec<(u32, FoldedReadState)>` kept ascending by `read_id`:

- `locate(read_id)` answers the append case — a read the record has not seen, with an id
  above every id it holds — with **one comparison against the last key**, and falls back to
  a binary search otherwise;
- the fold calls `locate` once and `store_at` once, so a first fold is a `push` and a
  re-fold is an assignment in place;
- `keyed_observations_counting` iterates the table **in place**: the per-record `Vec`
  allocation, the collect, and the `sort_unstable_by_key` are gone. That is one sort and one
  heap allocation removed *per covered base*;
- `refold_live_reads` loses its `ids.sort_unstable()` — the keys already come out ascending.

| fixture | stage 2 | stage 3 | change vs stage 2 | change vs baseline |
|---|---:|---:|---:|---:|
| 30× HG002 chr1 | 106.431 G | 99.522 G | **−6.49 %** | **−8.31 %** |
| ~130× tomato | 210.542 G | 191.453 G | **−9.07 %** | **−11.36 %** |
| 300× HG002 chr1 | 374.088 G | 331.694 G | **−11.33 %** | **−13.30 %** |

Raw:

```
tomato 130x  stage 3 : 192972549723 / 192868179902 / 192759340942
                       192933807767 / 193025707409 / 193087048796   (rebuild, confirmatory)
chr1   30x   stage 3 : 101520854337 / 101698768333 / 101428585717
                       101424137821 / 101421584522 / 101428461099   (rebuild, confirmatory)
chr1  300x   stage 3 : 332088910406 / 332049128965 / 332070393279
```

Ranges disjoint from stage 2 on all three fixtures. **The win grows monotonically with
depth**, which is what a change that removes a depth-sized sort per covered base should do.

### The depth sweep the previous attempt failed

The same `FoldedReads` code, on top of stage 1 only — active set still a `Vec` with
`swap_remove`, contributors still permuted:

```
tomato 130x  stage 1 + stage 3: 238215820440 / 238431907945 / 238393878268
chr1   30x   stage 1 + stage 3: 104911437573 / 104953132679 / 104910386089
chr1  300x   stage 1 + stage 3: 455890252573 / 455818020176 / 455839842066
```

| | 30× chr1 | ~130× tomato | 300× chr1 |
|---|---:|---:|---:|
| walk | 103.010 G | 236.910 G | 455.463 G |
| vs stage 1 | −3.18 % | **+12.92 %** | **+23.15 %** |

This is the reverted change, reproduced: −3.2 % at 30× and a cliff at depth. **Every entry
that is not a push is a `Vec::insert` shifting the tail**, and at 300× that tail is ~300
entries of `FoldedReadState`. With contributors ascending the same code is −6.5 % / −8.8 %
/ −10.3 %. The delta between the two rows at 300× is 33 percentage points.

Diff for this variant: `stage1_plus_3.diff`, kept so the comparison can be re-run.

### What is load-bearing, and what is not

Worth stating precisely, because it is the opposite of what one might assume:

- **Correctness does not depend on ascending arrival.** `locate`/`store_at` insert in sorted
  position whatever order reads arrive in, so `FoldedReads` is ascending by construction and
  `keyed_observations_counting`'s `f64` summation order is fixed no matter what. The
  determinism guarantee moved from the sort into the container.
- **Only the speed depends on it.** Ascending arrival is what makes the insert a push. That
  is why the +23 % row above is a performance result and not a wrong answer.
- The remaining out-of-order arrival is real but rare: a read silent at one position of a
  record's footprint (in a deletion, or on an `N`) that speaks at a later one, in a record
  that stays open across positions. It costs an `insert` and is handled.

### Gates

`cmp` against the round-3 stored copies:

```
DIFFERS   g_chr21  (  251792 lines)     one header line — record_widen_events 423 -> 425
IDENTICAL s_chr21  (    4406 lines)
DIFFERS   g_tom    ( 1718914 lines)     one header line — record_widen_events 622 -> 621
IDENTICAL s_tom    (   11945 lines)
```

`diff` on the two that differ reports exactly one changed line each — the same header
carried in from stage 2. Every locus line on all four dumps is identical.

- probe counters on chr21: `loci=236081 observations=251786 reads_admitted=54709` — exact.
- `cargo test --lib`: **2,880 passed, 2 failed** — the same two stage-2 parity tests, no new
  failures.
- `cargo test --examples`: 33 targets, all `ok`.
- `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- `cargo doc --no-deps`: 12 unresolved links, none in the three files touched — the
  documented baseline.

### Peak RSS

Measured alongside every instruction count. At 30× the reduction is real and the ranges are
disjoint:

```
chr1 30x baseline: 21053440 / 21135360 / 21004288 bytes  (21.00-21.14 MB)
chr1 30x stage 3 : 19283968 / 18808832 / 20135936 bytes  (18.81-20.14 MB)
```

−4 % to −11 %, from the fold table losing a hash map's load-factor slack and the active set
losing `by_read_id`. At ~130× peak RSS is unchanged (386–394 MB at every stage), because
there it is the CRAM decode holding the reference contig, not the generator.

Diff: `stage3.diff` (cumulative, all three stages — this is the worktree's current state).

---

## 5. Where the change could have moved bytes, checked rather than assumed

The brief named four places. Each was checked:

- **`q_sum` summation order.** Unchanged: `keyed_observations_counting` took reads in
  ascending `read_id` before (by sorting) and takes them in ascending `read_id` now (by
  iterating an ordered container). `ng_emits_the_same_bytes_in_a_second_process` was run
  explicitly and passes.
- **Allele bucket creation order.** It *does* change at stage 2, because it follows
  contributor order and `find_allele_index` assigns indices first-seen. It does not reach the
  output: `finalise_recycling` sorts the emitted observations by `(bases, witness, read_group)`,
  and the grouping is by bucket **index**, which is internally consistent either way. Confirmed
  by the 30× and ~130× dumps, where every locus line is identical.
- **Per-bucket `AlleleSupportStats` totals.** Their `q_sum` differs in the last bits, as
  predicted. They do not reach the emitted locus — the observations are re-derived from
  `folded_reads`, and the only reader of the bucket totals is `evict_unsupported_alleles`'
  `num_obs > 0` test, which is an integer. Confirmed by the same identical dumps.
- **`refold_live_reads`' bucket-creation order.** It sorted ids to defend against hash
  iteration order; the ids now arrive sorted, so the sort was removed and the order it
  produced is preserved exactly.

The debug assertions in `finalise_recycling` (every folded read resolves to exactly one
witness class; no unsupported bucket survives) were armed throughout: the whole `--lib`
suite is a debug build, and it runs the walk over ~257,000 loci.

---

## 6. The finding that needs an owner decision

- `src/ng/locus_generation/pileup/genome_walk.rs` (`contributors.truncate(cap)`, step 2b of
  `WalkerState::process_position`) — **[Hot-path]** Reordering contributors changes which
  reads the per-column depth cap keeps, and a capped read opens and widens nothing, so
  record footprints change with it
- **Confidence:** High
- **Hot-path evidence:** measured, not inferred. `ng_generic_loci_dump` over chr21 of
  `HG002_TR_v1.0.1_Tier_300x.bam`, baseline against stages 1+2+3:

  ```
  341094 lines (baseline)   341111 lines (all three stages)
  removed lines: 88356   added lines: 88373   diff hunks: 61917
  # generic_loci=241038 ... locus_sum_reads_discarded_by_cap=13205   (baseline)
  # generic_loci=241030 ... locus_sum_reads_discarded_by_cap=12214   (staged)
  # column_depth_truncations=909                                     (both)
  ```

  The *set* of loci barely moves — 131 loci present only in the baseline, 123 only in the
  staged run, out of 241,038, or one in a thousand. What moves is the content: 88,351 data
  rows of 341,094 differ. A representative hunk:

  ```
  < chr21	13427610	13427611	CT	345	complete	0	CT	343
  ---
  > chr21	13427610	13427610	C	250	complete	0	C	249
  > chr21	13427611	13427611	T	344	complete	0	T	344
  ```

  One 2-base record in the baseline, two 1-base records after — a read carrying a deletion
  was kept by the cap in one ordering and dropped in the other.
- **Pattern matched:** an order-dependent selection rule downstream of a container whose
  order was changed.
- **Mechanism:** `contributors.truncate(cap)` keeps the *first* `cap` contributors. Under
  the old permuted order that is a scrambled subset; under ascending `read_id` it is the
  `cap` **leftmost-starting** reads. The code's own comment already flags the corner — *"a
  truncated read carrying a deletion would have widened a record — so dropping it changes the
  footprint, and with it every other read's witness"* — but until now the subset was
  arbitrary rather than systematically biased toward early alignment starts, which is a
  statistical property (`placed_left`, witness extent), not only a determinism one.
- **Measurement plan:** already run; the numbers above are it. What is *not* run is whether
  the calls that come out of those loci differ — that needs a variant-calling comparison, not
  a walk measurement.
- **Complexity cost:** none to add; the question is whether the changed sampling is
  acceptable. Two ways out if it is not: make the cap order-independent (e.g. a strided or
  hash-ordered subsample rather than a prefix), or admit the change and re-baseline. Both are
  decisions above a perf review.
- **Fix:** not proposed. **At 30× and ~130× the cap never fires**
  (`column_depth_truncations=0` on every gate fixture and on both probe fixtures), so the
  target workload in the brief — human WGS, one sample, ~30× — is unaffected and every locus
  line is identical there.

---

## 7. Recommendation

1. **Take stage 1 unconditionally.** −2 to −3 % across the sweep, byte-identical, 2,882
   tests green. It is independent of the other two.
2. **Take stages 2 and 3 together or not at all.** Together they are −6.4 % at 30×, −8.8 %
   at ~130×, −10.3 % at 300× on top of stage 1. Separated, stage 2 is a small cost and stage 3
   is a 23 % regression at 300×.
3. **Two things must be decided before they land**, and neither is a performance question:
   the mate-overlap tie-break's divergence from production (two `parity::` tests), and the
   depth-cap sampling change at coverages high enough to hit `max_snp_column_depth`.

## Files

All copied into this directory, and present in the worktree's `tmp/`:

| file | what it is |
|---|---|
| `stage1.diff` | stage 1 alone, over `6fbbd09` |
| `stage2.diff` | stages 1 + 2 |
| `stage3.diff` | stages 1 + 2 + 3 — the worktree's current state |
| `stage1_plus_3.diff` | stages 1 + 3, no ordering — the +23 %-at-300× reproduction |
| `stage2_slab.diff` | the rejected slab shape (`active_read_set.rs` only, over stage 1) |
| `stage2_queue_noguard.diff` | the queue without the expiry guard (`active_read_set.rs` only, over stage 1) |

Measurement scripts in the worktree: `tmp/measure.sh` (30× + ~130×), `tmp/measure300.sh`
(300×), `tmp/gates.sh <tag>` (the four dumps plus `cmp` against the round-3 stored copies).
Dump outputs under `tmp/stage_measure/`, including the two 300× chr21 dumps
(`g300_base.txt`, `g300_s3.txt`) and their diff (`g300.diff`).
