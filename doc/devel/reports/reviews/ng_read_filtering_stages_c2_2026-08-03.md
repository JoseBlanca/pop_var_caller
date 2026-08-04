# Code review — ng read filtering in stages, C2 (the cursor takes over the loop)

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `7e8cfce` (C1b)
**Impl report:** [`ng_read_filtering_stages_c2_2026-08-03.md`](../implementations/ng_read_filtering_stages_c2_2026-08-03.md).

---

## 1. Scope

The uncommitted working-tree diff for C2 — 2 files, ~559 insertions / ~839 deletions:
`src/ng/read/input/cursor.rs` and `src/ng/read/filtering.rs`.

**Categories dispatched**, one `general-purpose` agent each, in its own git worktree:

| category | why |
|---|---|
| `reliability` | the tally is the plan's named silent surface |
| `refactor_safety` | C2 is a *move*; a byte-identical dump cannot see semantic drift |
| `module_structure` + `naming` | the milestone's checkpoint condition is a module-boundary claim |

## 2. Verdict

**Request changes** → all applied. One Blocker, three Majors, eight Minors.

The moved loop was correct. Every finding was about what the change *claimed*, what it left
behind, or what nothing tested.

## 3. Execution status

All three agents reproduced the gate independently: `cargo fmt --check` clean, `cargo test --lib`
**2,855 passed**, clippy exit 0, `cargo doc` 12 pre-existing unresolved links.

**Mutation coverage: 15 + 11 + 6 experiments across the three agents.** Five mutations survived
the whole suite; four were previously unknown.

## 4. Findings

### Blocker

#### B1: `cursor.rs` — a reposition may discard the whole step-1 tally, and no test notices
**Category:** reliability · **Confidence:** High (mutation-verified)

C2's accounting claimed the deleted `repositioning_the_source_does_not_reset_the_running_tally`
was covered by `the_step_one_tally_accumulates_across_regions`. **It is not.** The deleted test had
two halves and said so explicitly — non-erasure across a reposition, *then* continued accumulation,
because "non-erasure alone would also be satisfied by a tally that stopped counting". The claimed
successor covers only the second half, and only on the **reuse** path: its two regions are
`region(1, 40)` then `region(1, 100)`, so `1 >= 1` takes reuse and never repositions. Instrumented:
`jumping=1 reusing=1 replayed=3` — and the test's own doc comment ("a backward region, so nothing
is replayed and every read is filtered again") is wrong on both clauses.

Adding `read_group_tally.clear()` to the jump branch left **all 2,855 tests green**.

This is exactly the surface the plan singles C2 out for, and **C4 is about to add `reset_counts`**
— the first legitimate caller of "clear the tally" — into a file where clearing it on the wrong
edge is invisible.

### Major

#### M1: the owner's widening landed half-pinned, and its justification was false
**Categories:** reliability, refactor_safety, module_structure (all three agents)
**Confidence:** High

Two independent halves; **only the flag was tested.** Restoring the pre-C2 early commit of
`region`/`last_region_start` — reverting the ordering exactly — survived all 2,855, because the
flag masks it. So the docstring's claim that the ordering is what prevents the reuse-path wrong
answer holds only for a version without the flag.

And on unmutated code the comment *"a failed jump leaves the cursor exactly as it was"* was
**measurably false**: the eviction and the `reads_evicted` / `regions_jumping` counters ran *before*
the fallible `jump_to`. Measured — a served region then a failing backward seek gives `kept 4 → 0`,
`reads_evicted 0 → 4`, `regions_jumping 1 → 2`.

Neither half was reachable by test, because `with_failing_seek` was **all-or-nothing**: a reader
whose *first* seek fails has served nothing, so "left exactly as it was" and "was never anywhere"
are the same observation. That knob was the prerequisite.

#### M2: the checkpoint condition is not met, and not nearly
**Category:** module_structure · **Confidence:** High

The plan's Checkpoint C requires *"`filtering.rs` holding only the keep-or-drop rules and their
thresholds"*. It still held `resolve_read_group` + `unreadable_record`, `ReadFilterError` and
`ReadGroupCounts` — and the first kept a **module cycle** open (`filtering.rs` →
`read::input::read_groups`; `region_raw_aligned_reads.rs` → back).

The agent **contradicted the author's "out of scope" call** on the grounds that C2 is what moved
their home — the loop that raised `ReadFilterError` and keyed `ReadGroupCounts` used to live in
that file, and arch §6 already lists `ReadGroupCounts` in the "move to the cursor" row — and priced
it by building the whole re-homing: three import lines, no body changed, 2,855 passing, clippy
clean, 1,054 → 956 lines.

#### M3: two things called "counts" on one type, returning different values
**Category:** naming · **Confidence:** High

`AlignmentCursor` had a `counts` field (`CursorCounts` — what the cursor *did*), a
`read_group_counts` field (the drop tally) **and** a `read_group_counts()` method whose value
differs from the field's, because the method stamps the `other_sample` rider onto the first entry.
Field and method spelled the same, returning different things, both reachable in one function.

### Minor

- **Row 5 of the accounting overstates.** `a_failed_filter_is_not_restarted_and_says_so` drove the
  **`Decode`** arm; the cursor's fatal-path tests cover `Source` and `Reference` only.
  Decode-not-fatal survives. The arm is genuinely *unreachable* (C1's finding, re-traced) — a
  better disposition than the one written, but not the same one.
- **Fusing no longer latches.** The old `FusedIterator` stopped on a clean end of input too; the
  new guard is `failed` alone. Behaviour-preserving today only because all three reader arms latch
  `done` and the narrowing re-holds the early-stop record — a guarantee the cursor used to make
  itself and now assumes of two layers below, unstated and untested.
- **The one-buffer-per-pass invariant is documented and untested** — a fresh
  `NoodlesRawAlignedRead` per read survives the whole suite.
- **`"the walk stays stopped"` cannot fail** in the reference-failure test: its fatal record is
  last, so the script is exhausted either way. Only the `Source` test kills the fuse.
- **Three visibilities are wider than needed.** `verdict_on_raw_read`, `verdict_on_aligned_read`
  and `record_drop` compile at `pub(in crate::ng::read)` — worth tightening *because* spec §3's
  rule is "these must not migrate above the cursor", which the narrower visibility makes the
  compiler enforce instead of the spec merely asserting.
- **`buffer` / "buffered read"** mint two new names for Milestone A's *raw aligned read*.
- **`verdict_on_raw_read` does not take a raw read.** Arch §7 left the `(flag, mapq)` versus
  `&impl RawAlignedRead` choice open **for Milestone C**, so it was C2's to close.
- **Stale prose** in three places: `filtering.rs`'s module doc still advertised the `RecordSource`
  seam and the `ReadFilter` iterator; `RecordSource::Record`'s doc still named the iterator;
  `over_records`' doc still described building a filter.

### Nits

`FakeRecord::decode_fails`'s doc claims it exercises a path C2 removed; the accounting comment says
"These three are the rest" above a block of five.

## 5. What the review confirmed

- **The moved loop is behaviour-identical.** A line-by-line diff of `ReadFilter::next` against
  `next_filtered_read`: same order, same single buffer, same three fatal arms and variants, same
  two `record_drop` sites, same `kept += 1`, same `continue`/`return` shape. `fail`,
  `tally_for_buffered_read` and `read_group_counts` are body-identical modulo renames, and
  `over_records` reproduces `with_validated_contigs`'s initialisation exactly. The only deltas are
  the missing `EndOfInput` latch and `move_to_region`'s reordering — both filed above.
- **Deleted-code accounting is clean.** No live reference anywhere in `src/`, `examples/`,
  `benches/` or `tests/` to any of the seven deleted items.
- **Both of C2's own mutation claims reproduce**, kill set for kill set: read groups folded into one
  entry → 2 failed; rider on the last entry → 1 failed; swallowing the jump failure → 1 failed.
- **"The cursor causes region ends so it never has to ask"** verified by probe: a region end is not
  sticky, the tally is unchanged on re-asking, and the next region is correct.
- **C3 is not made harder** — slightly easier, since no generic crossed into the cursor.

## 6. The timing challenge, upheld

The author called a one-run-each-side `seconds` difference noise. The agent pointed out that B1
called **1.4 %** "consistent rather than noise" using **six runs per side on this same probe**, so
the evidence did not meet the bar this branch set. Re-measured by rebuilding `6e22718` in a
throwaway worktree: baseline mean **1.828 s** (six runs), C2 mean **1.829 s** (six runs) — under
0.1 %. The original single 1.782 s reading was an outlier, not the baseline.

## 7. Out of scope observations

- The step-1 order spec §2 calls "the whole design" is pinned by nothing, and plan **D2**, which is
  meant to pin it, cannot as written — its fixture is dropped by the region narrowing before the
  filters. Now confirmed by four independent measurements across C1 and C2.
- `region_raw_aligned_reads.rs` propagates with a bare `?` and never names the region it was
  serving, though it has it in hand.
