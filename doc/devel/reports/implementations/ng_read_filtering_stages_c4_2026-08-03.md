# ng — read filtering in stages, C4: the caller chooses the tally window

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `5e3b22f` (C3)
**Plan:** [`read_filtering_stages.md`](../../ng/impl_plan/read_filtering_stages.md) step **C4** —
the last of Milestone C.
**Design authority:** [spec](../../ng/spec/read_filtering_stages.md) §7 ·
[arch](../../ng/arch/read_filtering_stages.md) §3.4.

---

## 1. Plan

> `AlignmentCursor::reset_counts`, with a test that a fresh window starts empty and that nothing
> else on the cursor moves. Small and additive; kept out of C2 so that step stays about the loop.

## 2. Changes made

`AlignmentCursor::reset_read_group_counts()` clears the step-1 tally and records the point the new
window opens from. One new private field, `other_sample_at_window_start`, and one private helper,
`other_sample_this_window`.

## 3. The design was silent on one thing, and the plan's own test requirement settles it

`read_group_counts()` folds in the narrowing's `other_sample_records()` — the count of records
skipped as another sample's. **That number is cumulative for the life of the cursor, and the layer
below has no notion of a window at all.** So a reset that cleared only the tally would open a window
with every foreign record the cursor had ever stepped over already in it.

Spec §7 does not address it. The plan does, implicitly and decisively: *"a test that a fresh window
starts empty"*. A window carrying the whole cursor's other-sample history does not start empty. So
the reset records a baseline and the reported count is the delta.

The subtraction saturates — **and this report's first justification for that was aimed at the wrong
direction.** It said an underflow would mean the narrowing's count went backwards. That is
unreachable: every contribution to the count is a `+=`, and the baseline is sampled from the same
number, so it can never exceed it. `saturating_sub` guards nothing that can happen; it is kept
because it is free and a plain `-` would put a panic in an accounting helper for no benefit.

**The hazard that is real runs the other way, and it is not the window's.**
`CramAlignedReadsReader` adds a container's foreign-record count each time it decodes one, so a
container decoded twice — which a backward reposition can cause — contributes twice. A window
therefore starts at zero honestly and may **over**-report afterwards on CRAM, exactly as the
unwindowed number always has: that field's own doc has always called it container-granular and
warned it can run ahead of where a walk has reached. Recorded at the code rather than fixed, because
it is the reader's accounting and it pre-dates the window by two milestones.

## 4. Two deviations, recorded

**The method is `reset_read_group_counts`, not arch §3.4's `reset_counts`.** The cursor keeps two
unrelated tallies: this one, reached by `read_group_counts()`, and `CursorCounts` — what the cursor
*did* — reached by `counts()`. C2's review filed that collision as a Major and it was fixed there by
renaming the *field*; a method called `reset_counts` that resets only one of the two would re-open
it at the API. Arch's preamble states that signatures are illustrative and the contract is the
deliverable, so this is a naming choice rather than a design change — but it is a departure from a
name the architecture spells out, so it is recorded rather than absorbed silently.

**The new field is state the arch sketch does not have.** §3.4 lists the cursor's fields and does
not include a window baseline, because it did not anticipate §3.

## 5. Tests

| test | what it pins |
|---|---|
| `resetting_the_tally_starts_a_fresh_window_and_moves_nothing_else` | the window is empty afterwards — *every* counter of every entry, not merely `kept` — and the held reads, the walk's own counters and the cursor's ability to carry on are untouched |
| `resetting_the_tally_also_starts_the_other_sample_count_from_zero` | §3's half, which a plain `clear()` would miss |

**The first test was wrong on its first writing, and failing it taught me the property.** It
asserted that re-walking the *same* region after a reset refills the tally. It does not: a replayed
read is not filtered again, so it is not tallied again either — which is exactly spec §7's "a read
is filtered once, when first read off the file — never again when replayed". The test now uses a
**backward** region, which drops what is held and genuinely re-reads, and the comment says why a
forward one cannot serve.

**Suite: 2,858 → 2,865 (+7)** — two as built, five added by the review (§5.1). Fully accounted.

### 5.1 Five more tests, because six mutations survived

The review ran eighteen mutations and **six survived all 2,860 tests** — the most of any step in
this milestone, on its smallest diff. Two Blockers, and both are "my test could not reach it":

- **The window was never applied on the arm every real cohort walk takes.**
  `read_group_counts` folds the rider onto the first entry when the walk met read groups of its own,
  and into a fabricated entry when it met none. This report's own fixture consumes the whole contig
  before the reset, so afterwards only the empty arm ever ran — and reverting the `first_mut` line,
  *the single line the new field exists to protect*, left the suite green.
- **"Nothing else moves" was unpinned for the region in flight**, because no test reset mid-region.
  Both existing tests reset at a boundary and then call `reads_of`, which begins with
  `move_to_region` and re-establishes anything a mutation disturbed. `region = None`,
  `examined = 0`, `last_emitted = None` and `last_region_start = None` all survived. `examined = 0`
  **serves a read twice**, and the order guard cannot catch it: it compares with `<`, so a replay at
  the same position passes.

| test added | closes |
|---|---|
| `resetting_the_tally_mid_region_leaves_the_region_being_served_untouched` | `region = None`, `examined = 0` |
| `resetting_the_tally_leaves_the_order_guard_armed_for_the_region_being_served` | `last_emitted = None` |
| `resetting_the_tally_does_not_make_the_next_forward_region_reposition` | `last_region_start = None` — that one number *is* the forget rule |
| `read_group_counts_scopes_the_other_sample_rider_when_the_new_window_met_its_own_reads` | the `first_mut` arm |
| `read_group_counts_stays_empty_when_a_window_is_served_only_from_replayed_reads` | the replay property, promoted out of a comment |

### Mutations run

| mutation | killed by |
|---|---|
| the baseline is not recorded on reset | `resetting_the_tally_also_starts_the_other_sample_count_from_zero`, alone |
| the reset also clears the kept reads | `resetting_the_tally_starts_a_fresh_window_and_moves_nothing_else`, alone |
| the `first_mut` arm ignores the window | the rider-scoping test, alone (**survived 2,860** before) |
| reset sets `region = None` / `examined = 0` / `last_region_start = None` | two tests each (**each survived 2,860** before) |
| reset sets `last_emitted = None` | the order-guard test, alone (**survived 2,860** before) |
| `saturating_sub` → `-` | still survives, unreachable by construction — §3 |

## 6. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo test --lib` | **2,865 passed**, 0 failed, 5 ignored |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 errors, 0 warnings |
| `cargo doc --no-deps` | **12** unresolved links — the pre-existing baseline |
| four acceptance dumps, `cmp` | **byte-identical** |
| `ng_generic_walk_probe` chr21 | `loci=236081 observations=251786 reads_admitted=54709` |

C4 adds a method nothing in the walk calls, so the dumps could not have moved; they were run
because the plan checks every step against them.

## 7. Follow-ups

- **Nothing in the tree calls `reset_read_group_counts` yet.** It is spec §7's stated capability —
  the caller chooses the window — and the first caller will be whatever wants per-window drop rates
  (a run report, or a per-chromosome summary). Recorded so the absence is a decision rather than an
  oversight.
- **Arch §3.4 should be reconciled** with the method's name and the new field — for the checkpoint.
