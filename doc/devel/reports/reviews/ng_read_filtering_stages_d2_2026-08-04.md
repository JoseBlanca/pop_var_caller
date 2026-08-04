# Review — ng read filtering in stages, D2: the conversion is asked for nothing when the first filter drops

**Date:** 2026-08-04 · **Branch:** `ng-generic-perf` · **Base:** `5320fe4` (D1)
**Scope:** the working-tree diff for plan step **D2** — as submitted, 165 added lines in
`src/ng/read/input/cursor.rs`.
**Design authority:** [spec](../../ng/spec/read_filtering_stages.md) §2, §3, §6, §8 ·
[arch](../../ng/arch/read_filtering_stages.md) §3.4, §5 ·
[plan](../../ng/impl_plan/read_filtering_stages.md) D2, Checkpoint D.
**Fixes:** [`fixes_applied_2026-08-04_v2.md`](fixes_applied_2026-08-04_v2.md).

Three `general-purpose` agents, each in its own git worktree, detached at `5320fe4` with the
staged diff applied by `git apply`. Raw findings under
`tmp/review_2026-08-04_d2-conversion-not-asked/`: `reliability.md`, `naming-structure.md`,
`intent-refactor.md`.

---

## 1. Verdict

**4 Major, 12 Minor, 12 nits. No Blocker.** The mechanism works — all three agents independently
confirmed both tests have **unique detection power of exactly 1**, which is the thing D1 did not
have. What they took apart was the *reach* of one test and the *evidence* offered for the design.

## 2. The Majors

**M1 — the test pinned three of the first filter's six reasons, and the gap was reachable.**
*(reliability)* The submitted test 1 scripted duplicate, supplementary and secondary. The reviewer
converted the `Unmapped` and `LowMapq` drops before rejecting them — a per-reason hoist, and the
`Drop` arm is one `match`, so it is an ordinary mistake — and **all 2,869 tests stayed green**.
Filter #2 (low MAPQ) is the highest-volume drop on real data and was the one not measured. This is
the milestone's recurring defect in miniature: a test that looks like it covers a property and
covers a sample of it.

**M2 — the stated failure mode of a call-site increment is false, and the truth is worse.**
*(naming, reliability)* The diff justified counting inside the callee by predicting that an
increment beside the call, left behind by a hoist, would report *"zero conversions while every read
was being converted"*. Two reviewers built it. It reports **2** — exactly the number the test
asserts — because the dropped reads `continue` before reaching the abandoned increment. **The test
passes and the instrument silently stops instrumenting.** The conclusion was right; the argument
for it was wrong, and in the direction that made it sound more dramatic than it is.

**M3 — the counter's shape was the weakest part, and the alternative was already in the file.**
*(naming, and reliability's Minor (a))* The thread-local's justification — *"a field would have to
be threaded through every constructor"* — is refutable by one grep: `AlignmentCursor` has exactly
one struct literal, and it already holds `counts: CursorCounts`. Both reviewers built
`CursorCounts::reads_converted` incremented inside `convert_buffered_read`, deleted the
thread-local and its two free functions, and measured identical detection power. It needs no reset
protocol, is folded across a sample's files for free by the exhaustive-destructure `AddAssign`, is
per-cursor rather than thread-global, and ships. The submitted version was also the crate's **first**
`#[cfg(test)]` static and its first `#[cfg(test)]` statement inside a production function body.

**M4 — a compile-time pin was possible and the diff claimed none was.** *(intent-refactor)* The
doc said the method was *"the only way this design's central ordering can be tested at all"*. A
reviewer built a private zero-sized witness minted by the first filter's `Keep` arm and required by
the conversion; the hoist then fails to **compile**. It was not adopted — spec §1 says the design
"does not … add a type", and decisively the witness pins only half the property, since no type can
forbid someone writing a second copy of the length rule — but the claim had to go.

## 3. What all three agreed on, and it matters

Each agent ran the mutations independently and got the same result: **the hoist fails test 1 and
nothing else; hoisting filter #7 fails test 2 and nothing else.** No other test in the crate
notices either. After D1 — where three agents measured the submitted test's unique detection power
at zero — that is the result this step needed.

They also each verified, from the code rather than from the diff's prose, that **the plan's own D2
fixture is genuinely unbuildable**: `read_next` proves contig, footprint and read group before
yielding, which is exactly what `decode` refuses.

## 4. The Minors worth naming

- **"2,866 passed" appears three times and does not reproduce** — the figure is from `f5630f8`,
  before D1 added its test; the tree's baseline is 2,867. Two of the three copies also state it as
  the result of running the suite *with the hoist applied*, which on the finished step is
  `2868 passed; 1 failed` — the opposite of what the sentence says.
- **The same justification is written out three times**, plus a fourth pre-existing copy on
  `next_filtered_read`. Three of the four carried the stale count.
- **The two tests describe the same movement with opposite verbs** — "risen above" and "sunk
  below" — and test 2's doc named its mutation backwards. Both mutations move a step *earlier*.
- **Test 2's rationale cited spec §2 for a claim §2 does not make about #7.** §2's "writing the
  rule a second time" argument is about the mismatch rule and the CIGAR scan. A raw record carries
  a sequence, so #7's length *can* be had without converting — which is why it needs a test at all.
- **A cross-reference that could not support its claim** —
  `an_unmapped_read_that_clears_the_mapping_quality_filter_is_charged_to_unmapped` scripts an
  unmapped read **with** a footprint, so it says nothing about the footprint-less case.
- **`CursorCounts::reads_decoded` does not count decodes** — it counts reads surviving both
  filters. Pre-existing, made visible for the first time by this diff.

## 5. Not reproduced as findings

The reliability agent closed every "passes for the wrong reason" trap it was asked to check: the
tally assertions sum to the record count in both tests, so a record lost to the narrowing breaks
them; `only_tally` asserts a single read group; there is no early stop; the reset protocol was
belt-and-braces (libtest gives each test its own thread even at `--test-threads=1`).
`convert_buffered_read` is behaviour-neutral, `&self` was correct for the submitted shape, and
`#[inline]` would be noise — there are zero `#[inline]`s across `cursor.rs`, `filtering.rs` and
`aligned_read.rs`, including on hotter functions.

## 6. Carried to Checkpoint D

- **The witness type (M4)** — rejected here on spec §1 and on covering only half the property, and
  recorded on the method. If the owner wants the hoist to be a build error rather than a test
  failure, it is a spec §1 amendment and a small change.
- **Design-doc amendments** the reviewers specified: plan D2 and Checkpoint D, arch §3.4 and §5,
  spec §8's "every read" wording. All landed in their own commit.
