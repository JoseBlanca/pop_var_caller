# ng — read filtering in stages, C2: the cursor takes over the loop

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `7e8cfce` (C1b)
**Plan:** [`read_filtering_stages.md`](../../ng/impl_plan/read_filtering_stages.md) step **C2**.
**Design authority:** [spec](../../ng/spec/read_filtering_stages.md) §2, §3, §5, §7 ·
[arch](../../ng/arch/read_filtering_stages.md) §3.4, §4, §6.

---

## 1. Plan

The step the whole milestone is built around, and the one the plan insists lands alone:

> `AlignmentCursor` gains the record buffer, the reference, the fetch scratch buffer, the config,
> the tally and a `failed` flag, and calls the two filters and the conversion itself. `ReadFilter`,
> `FilterState`, `restart_after_end_of_input`, `has_failed`, `source_mut`,
> `with_validated_contigs` and `ReadFilterBuffers` are all deleted.

Plus the owner's widening at C1 (2026-08-03): the `failed` flag covers **both** routes into a
stopped cursor — a failed read *and* a failed reposition — and `region` / `last_region_start` are
committed only after the jump has succeeded.

## 2. Changes made

### `read/input/cursor.rs` — gains the loop

`AlignmentCursor` now holds `reads`, `buffer`, `reference`, `ref_buf`, `config`,
`read_group_counts` and `failed` where it held one `filter`. The loop lives in
**`next_filtered_read`**, in the order spec §2 fixes: read into the one reused buffer, reject on
flag and mapping quality, convert only the survivors, reject on length / CIGAR / mismatch
fraction, charge every drop to its read group.

`fail` sets the flag and yields the error once. `tally_for_buffered_read` keys the tally on the
buffer, because the first filter drops a read before anything owns it. `read_group_counts()` folds
in the `other_sample` rider against the first entry.

**`move_to_region` reposition ordering (the owner's ruling).** The jump now happens *before* any
of the new region's state is committed, and a failure sets `failed`:

```rust
if !reuse && let Err(source) = self.reads.jump_to(region) {
    self.failed = true;
    return Err(CursorError::ReadRecord { … });
}
```

Nothing restarts anything on a move any more, and the absence is the point: a filter living apart
from the cursor reached "end of input" at the end of *every* region, so a cursor had to undo that
each time or the first region silenced it for the whole chromosome. The cursor causes region ends,
so it never has to ask (spec §5).

### `read/filtering.rs` — loses the loop, keeps the rules

`ReadFilter`, `FilterState`, `restart_after_end_of_input`, `has_failed`, `source_mut`,
`reference()`, `counts()`, `tally_for_current_record`, `fail`, the `Iterator` and `FusedIterator`
impls, `with_validated_contigs` and `ReadFilterBuffers` are all gone — 806 lines out.

`verdict_pre_decode` → **`verdict_on_raw_read`** and `verdict_post_decode` →
**`verdict_on_aligned_read`** (arch §6), both now `pub(crate)` so the cursor can call them.
`ReadFilterCounts::record_drop` likewise. `ReadFilterError` and `ReadGroupCounts` stay, with
`ReadFilterError`'s doc rewritten: it no longer describes an iterator that does not exist, and says
the cursor raises it.

## 3. A scope addition, on the owner's ruling — and the checkpoint condition

C2 first left four things in `filtering.rs` that are not keep-or-drop rules, on the argument that
moving them was scope C2 had not been given. **The review contradicted that and was right:** C2 is
what moved their home, because the loop that raised `ReadFilterError` and keyed `ReadGroupCounts`
used to live in that file, and arch §6 already lists `ReadGroupCounts` in the "move to the cursor"
row. It also left a **module cycle** open — `filtering.rs` importing `read::input::read_groups`
while `region_raw_aligned_reads.rs` imported back.

**Owner's ruling (2026-08-03): fold the re-homing into C2**, so the milestone does not reach its own
checkpoint with the condition unmet.

| moved | to | why there |
|---|---|---|
| `resolve_read_group`, `unreadable_record` | `read/input/read_groups.rs` | beside the `ReadGroupResolution` and `RecordOwner` they consult — and it is what closes the cycle |
| `ReadFilterError` | `read/input/cursor.rs` | the cursor raises it, now that the loop is there |
| `ReadGroupCounts` | `read/input/cursor.rs` | arch §6's "move to the cursor" row |

**The checkpoint condition is now met.** `filtering.rs` imports production's reused predicates, the
read it judges, the reference, `types` and `io` — nothing from `read::input`. 1,787 → **961 lines**.

Still recorded, and still C3's:

- **The cursor imports `RecordSource`**, because `read_next` and `other_sample_records` are trait
  methods until **C3** makes them inherent. The trait is the last thing in `filtering.rs` that is
  not a rule, and C3 deletes it.
- **Three dead test helpers were deleted** (`rewind`, `records_consumed`, `named_fake`) because
  clippy is `-D warnings` and they had no callers once the `ReadFilter` tests went. `FakeSource`
  and `FakeRecord` survive to C3.

## 4. Tests

**Ten tests died with `ReadFilter`, and C2's first account of where their properties went was
wrong on two rows.** The review corrected both by experiment; this is the version that survived it.

- **Four have no successor, by design** — they drove `source_mut`, `restart_after_end_of_input`
  and the three-way `FilterState`, which is exactly what spec §5 says collapses to one flag.
- **Two covered the fatal stop.** The cursor's fatal-path tests cover it for `Source` and
  `Reference` — but one of the two drove the **`Decode`** arm, which no test now reaches. C1
  established no input can, so it is unreachable rather than untested: a better disposition than
  the one first written, but not the same one.
- **One was the tally surviving a reposition, and it was *not* covered** by
  `the_step_one_tally_accumulates_across_regions` as first claimed — §5.
- **Three were re-homed:**

| re-homed as | property |
|---|---|
| `an_unmapped_read_that_clears_the_mapping_quality_filter_is_charged_to_unmapped` | the #5 counter the drop fixture omits — and it needs a *placed* unmapped read, because one with no footprint is dropped below the filters, uncounted |
| `a_cursor_over_an_empty_script_yields_nothing_and_counts_nothing` | an empty walk still answers with a tally entry to fold |
| `the_tally_is_readable_before_the_walk_is_finished` | the tally is running, not final |

**Two more were written because the plan's stated oracle turned out to be insufficient** — see §5.

| new | property |
|---|---|
| `two_read_groups_are_tallied_apart_rather_than_summed` | spec §7's central claim, which nothing pinned |
| `the_other_sample_count_rides_on_the_first_entry_and_is_not_a_drop` | the rider's placement, which needed a fixture nothing had |

And `a_reposition_that_fails_is_refused_rather_than_answered` grew the two assertions its C1 note
promised, becoming `a_reposition_that_fails_is_refused_and_stops_the_cursor`.

**Suite: 2,860 → 2,857 (−3)** — ten deleted, seven added. Fully accounted.

## 5. The plan's oracle is insufficient, and mutation is how we know

The plan names C2's oracle: `a_walk_charges_every_drop_reason_by_hand_count` green either side,
plus the `other_sample` rider on the first entry. **Two mutations passed that oracle and the entire
suite:**

| mutation | before | after |
|---|---|---|
| a drop met but never counted | killed by the oracle and the three re-homed tests | — |
| **every read group folded into one tally entry** | **survived all 2,853** | killed by `two_read_groups_are_tallied_apart_rather_than_summed`, alone |
| **the `other_sample` rider on the last entry, not the first** | **survived all 2,855** | killed by the rider test, alone |
| a failed reposition no longer stops the cursor | — | killed by the grown C1 test, alone |
| **the tally cleared on the jump path** | **survived all 2,855** (the review's Blocker) | killed by `the_tally_survives_a_reposition_that_drops_everything_and_jumps`, alone |
| **the commit ordering reverted** — state committed above the fallible jump again | **survived all 2,855** | killed by `a_failed_reposition_leaves_the_cursor_untouched`, alone |

The last two came from the review, and the first of them is why the accounting table above is now
two rows different. **`repositioning_the_source_does_not_reset_the_running_tally` was not covered
by `the_step_one_tally_accumulates_across_regions`:** that test's second region begins at or after
its first, so it takes the **reuse** path and never repositions — instrumented, it reports
`jumping=1 reusing=1 replayed=3`, and its own doc comment ("a backward region, so nothing is
replayed") is wrong on both clauses. Nothing in the tree drove the tally across a jump.

**The owner's widening also landed half-pinned.** The `failed` flag masks the ordering, so
reverting the ordering alone left the whole suite green — the docstring's claim that the ordering
is what prevents the reuse-path wrong answer only held for a version without the flag. Worse, the
claim that a failed jump "leaves the cursor exactly as it was" was measurably false: the eviction
and two counters still ran *before* it (`kept 4 → 0`, `reads_evicted 0 → 4`,
`regions_jumping 1 → 2`). Fixed in the **code** rather than the prose — the reposition now precedes
everything, so the promise is true — and `with_failing_seek` became **positional**
(`with_failing_seek_at`), because an all-or-nothing fault cannot express "a seek that fails after a
region has been served", which is the only state in which "left as it was" differs from "was never
anywhere".

The read-group one is the important one. Spec §7's whole argument is that a drop rate is a *read
group's* property and summing across them erases the signal — and merging every library into one
entry changed no output, no dump and no read. The existing multi-read-group test
(`a_cursor_keeps_every_read_group_of_its_sample_not_just_one`) collects `(qname, read_group)` off
the *reads* and never looks at the tally.

The rider needed a fixture nothing had: **two** of our read groups *and* a foreign record. With one
read group the first entry is the last, so moving the rider is undetectable.

Both gaps are **pre-existing** — the deleted `ReadFilter` had the same shape — but C2 is the step
that moves the tally, so they are C2's to close.

## 6. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo test --lib` | **2,857 passed**, 0 failed, 5 ignored |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 errors, 0 warnings |
| four acceptance dumps, `cmp` | **byte-identical** |
| `ng_generic_walk_probe` chr21 | `loci=236081 observations=251786 reads_admitted=54709` |

**The oracle, green either side:** `a_walk_charges_every_drop_reason_by_hand_count` passed before
the change and passes after.

`cargo doc --no-deps` holds at the **12**-link pre-existing baseline. It reached 13 twice during
this step — a bad `AlignmentCursor` path, then `read/mod.rs` still pointing at
`filtering::ReadFilterError` after the move — and the count is what caught both.

**Timing, six runs a side, on the same machine.** The first draft of this section compared one
C2 run against a single 1.782 s baseline reading and called the difference noise; the review
refused that, correctly, because B1 set the bar at six runs a side on this exact probe when it
called **1.4 %** consistent. Measured properly, by rebuilding `6e22718` in a throwaway worktree:

| | runs | mean |
|---|---|---|
| baseline `6e22718` | 1.820 1.821 1.825 1.829 1.830 1.844 | **1.828 s** |
| C2 | 1.824 1.824 1.829 1.830 1.834 1.834 | **1.829 s** |

**No measurable change — under 0.1 %, well inside the spread.** The single 1.782 s reading was the
outlier, not the baseline. That is the expected result on structural grounds too: the loop is
line-for-line the old one bar renames, and C2 *removes* work (a three-state compare becomes a
bool, and the per-region restart call is gone).

## 7. Follow-ups

- **`RecordSource` is the last thing in `filtering.rs` that is not a rule**, and **C3** deletes it.
  After that the file is exactly spec §6's list: `ReadFilterConfig`, `DropReason`, `FilterVerdict`,
  `ReadFilterCounts`, and the two verdicts.
- **`ReadFilterError::Decode` is now unreachable *and* unpinned.** C1 established no input can
  reach it; C2 deletes the two tests that raised it through `FakeSource`. Rewriting its arm as
  `Err(_) => continue` will survive the whole suite. Recorded at the code.
- **`ReadFilterError::Source` still names two unrelated fatal conditions** (a failed read and an
  unresolvable `RG` tag), so arch §4's table is wrong — Checkpoint C.
- **The step-1 order spec §2 calls "the whole design" is pinned by nothing.** Hoisting the
  conversion above the first filter leaves all 2,857 tests green. That is what plan **D2** is for,
  and D2 cannot be built as written — Checkpoint C.
- **The fusing no longer latches.** The old iterator stopped on a clean end of input as well as on
  a failure; the cursor's guard is `failed` alone. Behaviour-preserving today, but only because all
  three reader arms latch `done` and the narrowing re-holds the early-stop record — a guarantee the
  cursor used to make itself and now assumes of two layers below, unstated. Raised by the review as
  a Minor; recorded rather than fixed, because changing it is a contract decision.
