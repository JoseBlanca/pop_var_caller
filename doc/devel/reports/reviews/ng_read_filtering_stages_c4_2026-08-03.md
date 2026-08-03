# Code review — ng read filtering in stages, C4 (the caller chooses the tally window)

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `5e3b22f` (C3)
**Impl report:** [`ng_read_filtering_stages_c4_2026-08-03.md`](../implementations/ng_read_filtering_stages_c4_2026-08-03.md).

---

## 1. Scope

The uncommitted diff for C4 — one file, ~187 insertions / ~2 deletions. `reliability` plus
`naming` on the new API, one agent in its own worktree.

## 2. Verdict

**Request changes** → all applied. **Two Blockers, two Majors**, plus Minors and Nits.

**Eighteen mutations run. Six survived all 2,860 tests** — the most of any step in this milestone,
on its smallest diff. Five tests added; the suite is green at 2,865 with all six killed.

## 3. Findings

### Blocker

#### B1: the window is never applied under test on the arm every real walk takes
**Confidence:** High (mutation-verified)

`read_group_counts` folds the other-sample rider in two places: onto the first entry when the walk
met read groups of its own (`first_mut`), and into a fabricated entry when it met none. **Only the
second was reached.** The author's fixture consumes the whole contig before the reset, so
afterwards the tally is empty forever and the `None` arm is the only one that runs.

Reverting the `first_mut` arm to its pre-C4 body — i.e. the single line the new field exists to
protect — left **all 2,860 tests green**. That arm is the one every real cohort walk takes.

#### B2: "resets the tally and nothing else" is unpinned for the region in flight
**Confidence:** High (mutation-verified)

`region = None`, `examined = 0` and `last_emitted = None` each survived all 2,860. The root cause is
structural: **no test called the reset mid-region.** Both existing tests reset at a region boundary
and then call `reads_of`, which begins with `move_to_region` and re-establishes everything a
mutation had disturbed.

`examined = 0` is the worst of the three — it serves a read **twice** (`["r0","r0","r1",…]`), and
the order guard cannot catch it because it compares with `<`, so a replay at the same position
passes.

### Major

- **`last_region_start = None` also survives.** That one number *is* the forget rule, so clearing
  it makes every later region jump instead of reusing — a total loss of the feature, invisible in
  the reads and in byte-identical dumps. Only the cursor's own counters can see it.
- **The `saturating_sub` rationale points the wrong way.** Replacing it with plain `-` leaves 2,860
  green, because the branch is unreachable: every contribution to the underlying count is a `+=`
  and the baseline is sampled from that same number. The hazard that *can* occur is the opposite —
  `CramAlignedReadsReader` adds a container's foreign-record count on each decode, so a container
  decoded twice (which a backward reposition can cause) contributes twice. Filed Medium confidence
  with a named CRAM verification step.

### Minor / Nits

- **The `None`-arm comment is falsified by C4.** "Every record was another sample's, so the filters
  met none at all" was true before the reset existed; an empty tally now also means *the window was
  just reset*. And the fabricated `None` key is spelled the same as the genuine one
  `ReadGroupCounts` documents — a read whose reader never stamped a group — so a caller cannot tell
  them apart.
- **The property the author found by failing a test deserves its own test**, not a sentence inside
  another test's comment: a window served entirely from replayed reads tallies nothing.

## 4. On the two deviations the brief asked about

- **The rename to `reset_read_group_counts` is right** — `reset_counts` would name the tally by
  `counts()`'s word and re-open the Major C2 filed. **But the author executed half his own
  rationale:** `grep -rn "reset_counts" doc/devel/ng/` still finds **six live references** across
  the spec, the plan and the arch, and the arch sketch still calls the field `counts` rather than
  `read_group_tally`. Three design documents now name a method that does not exist.
- **The new field was the correct resolution of the design's silence.** The baseline mechanism
  itself is well covered (four mutations killed); it is the *arm that uses it* that was not.

## 5. What the review confirmed

- Both of the author's own mutations are honest kills.
- **The replay property is real**: re-walking the same region after a reset serves five reads and
  tallies zero, confirmed by running it — so the author's first test, which asserted the opposite
  and failed, was wrong in the way he recorded.
- `other_sample_at_window_start` is correctly initialised; that mutation is killed by five tests.

## 6. Out of scope observations

The six stale `reset_counts` references in the design documents are the checkpoint's business —
this skill does not edit spec, arch or plan.
