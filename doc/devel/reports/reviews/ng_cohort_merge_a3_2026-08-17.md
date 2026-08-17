# Review — ng cohort merge, step A3 (the reach walk)

*2026-08-17, branch `ng-cohort-merge`, working tree at stash commit `26c93281`. Five
category checklists, three sub-agents, two of them mutation-testing in isolated
worktrees. Per-category audit trail: `tmp/review_2026-08-17_ng-cohort-merge-a3/`.*

## 1. Scope

- **Reviewed:** `LocusCloser`, `ClosedLocus`, `SampleMembers` and their tests in
  [`close.rs`](../../../../src/ng/run/cohort_merge/close.rs); the `pub mod close;` line.
- **Categories dispatched:** `reliability` and `refactor_safety` (mutation, in worktrees),
  `naming`, `smells`, `idiomatic`. **Not dispatched:** `defaults`, `errors`,
  `unsafe_concurrency`, `tooling`, `module_structure` — no parameter, no error path, no
  concurrency, no manifest or tree change.
- **Both worktree agents were interrupted by a transient API error and resumed**; both
  finished, and the reliability run is explicitly reported as **partial** (7 of the 11
  mutations it identified).

## 2. Verdict

**Approve-with-changes.** No Blocker. Five Major across the two mutation agents, every one
applied.

## 3. Execution status

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --all-features -- -D warnings` | clean |
| `cargo test --lib ng::run::cohort_merge` | 16 passed, 0 failed (at review time) |
| `cargo clippy --all-targets --all-features` | 49 pre-existing errors, none in scope |

**Mutation testing:**

| agent | run | survived | changed-no-behaviour |
|---|---|---|---|
| refactor_safety | 8 | 0 | 2 (both proved inert, not assumed) |
| reliability (partial) | 7 | 2 | 1 |

## 4. The strongest result: the walk is production's, measured

The refactor_safety agent **ported production's `CohortSpanFold` + `derive_is_kept`
grouping into a temporary test and ran it against `LocusCloser` over 5,000 random
cohorts** — fixed seed, 1 to 4 samples, 0 to 8 observations each, starts in 1..=30 and
spans in 1..=6, so shared starts and differing spans at one start are the common case
rather than the rare one.

- **Loci start and end: identical in all 5,000.**
- **Membership: identical in all 5,000.**
- **The boundary is genuinely exercised:** it instrumented the corpus and counted **5,025
  members joining at exactly the running reach** — so the agreement is not reached on
  inputs that dodge the `<=`.

And the boundary question answered directly: production takes a position when
`positions[j] <= group_end`; ng breaks when `start > reach`, and `!(start > reach)` is
`start <= reach`. Same predicate. Bracketed from both sides by a deletion at 10 reaching
14 with a SNP at exactly 14 (one locus in both) and at 15 (two in both).

**The one divergence is deliberate and now measured:** ng sums every sample's
non-reference reads where production sums the per-position maximum. Over the same corpus
ng's total is never below production's and **strictly above it in 1,824 of the 5,000**
cohorts.

## 5. Findings

### Major

**M1 — `ClosedLocus::span()` walked straight back into the coordinate-ceiling defect that
A2 documented.** *(reliability)*
`span()` was `region.len()`, and `GenomeRegion::len`'s own doc — written one step earlier
in this same milestone — records that it overflows its `+ 1` at `end == u64::MAX`: a panic
in debug, and **0 in the release profile**, where overflow checks are off. The agent
verified the panic by constructing a `ClosedLocus` directly. Zero is the dangerous answer:
`span()` is what A4's width verdict compares against, so the widest locus expressible would
be judged narrow enough to build.
**Fixed:** subtract first, add after, both saturating; pinned by
`the_span_is_the_base_count_at_the_coordinate_ceiling`.

**M2 — an unsorted sample produced a locus that does not contain its own member.**
*(refactor_safety, and reliability independently)*
The precondition was stated in three doc comments and checked nowhere. Verified: `[50–50,
10–10]` in one sample closes one locus spanning 50–50 holding an observation forty bases
outside it. Structural rather than accidental — `reach >= start` always, so any backward
jump is absorbed into whichever locus is open — and undetectable downstream: A4 would judge
the wrong span and §4.2's projection cannot pad a member outside the span.
**Fixed** with a **release-level** assertion inside the walk. One comparison per
observation, against per-observation work that already includes a base-by-base sequence
comparison. Release rather than `debug_assert!` because the release profile is the one this
repo runs — a trap it has recorded hitting twice. Pinned by
`a_sample_out_of_coordinate_order_is_refused`.

**M3 — the members' documented "ascending sample order" had no test that could fail.**
*(reliability)* Emitting them in consumption order instead passed all 16 tests; a probe
proved the outputs differ. Every existing fixture happened to make the two orders agree.
**Fixed** by `the_members_are_in_sample_order_not_consumption_order`, whose sample 1 opens
the locus and sample 0 joins inside it.

**M4 — the `u32` saturation on the total was unpinned.** *(reliability)* Replacing
`saturating_add` with `+` passed all 16; the largest total any fixture reached was 9.
**Fixed** by `the_non_reference_total_saturates_at_the_u32_ceiling`.

**M5 — every consumed observation cost a scan of all k samples *through* k scattered
pointers, and one whole scan per locus was computed and thrown away.** *(smells,
idiomatic)* The sibling merge settled this question in the other direction:
`MergedCursors` holds "keys beside the heads, not read through them". Here k is the cohort
size, and the walk runs over every observation of every sample on every run.
**Fixed:** `keys` beside the cursors, refreshed only for the sample just advanced, so the
argmin scans one contiguous array; the head now comes back *with* its sample, so no caller
looks it up again; `cursors_at_open` is a reused scratch field rather than a `Vec` cloned
per locus; and the member vector is sized exactly during the walk instead of growing
through a `Filter` whose `size_hint` lower bound is 0.

### Minor — all applied

- **Two doc comments named mechanisms the code does not use.** The tie-break was credited
  with the walk's determinism; inverting it gives byte-identical output, because member
  order comes from the sample-indexed collection and the two aggregations are `max` and
  `+`. And "production's group walk with the columns removed" omitted that production also
  *judges* in the same loop and sums a different quantity. Both rewritten. *(smells,
  reliability, convergent)*
- **The test named for the non-reference total did not discriminate ng's rule from
  production's** — its three observations sit at three different positions, where a
  per-position maximum also sums to 9. The assertion that did discriminate was buried in a
  test named for member ordering. **Fixed** by adding spec §15's own fixture as
  `one_non_reference_read_in_each_of_two_samples_sums_to_two`, and by saying on the old
  test what it does and does not separate.
- **The test helpers inverted the vocabulary this module fixes** — `locus()` built a
  `SampleLocusObservations`, which spec §1.3 and A2's own doc call an *observation*, while
  `observation()` built a `SequenceObservation`. Renamed to `observation` and `sequence`.
  *(naming)*
- **The `9` that a test's "answers 18" argument turns on was hidden in a helper.** Now a
  parameter at the call site. *(naming)*
- `per_sample` → `observations_per_sample`; `earliest_head` → `sample_with_earliest_head`
  (it returns a sample, and now returns the head with it); `opened_at` → `cursors_at_open`,
  which says which of the two number spaces it lives in. *(naming, idiomatic)*
- **The progress guarantee was undocumented and unasserted.** A mutation of the comparison
  did not fail an assertion — it killed the test binary with SIGKILL, the one failure class
  a suite cannot catch by failing. The real code cannot spin, and the argument is now in the
  doc plus a `debug_assert!`. *(reliability)*
- Missing coverage added: an empty sample among covered ones; contig-before-position
  ordering in the merge key; and an observation starting exactly at the reach.

### Not applied, with reasons

- **Fold the four parallel per-sample vectors into one vector of a struct.** *(disputed)*
  It would let the compiler check their equal lengths, at the cost of the thing the layout
  exists for: the argmin scans `keys` and nothing else, and that is contiguous only while
  it is its own array. The sibling merge makes the same trade. Recorded in the type's doc.
- **A heap or loser tree for the argmin.** The scan is still O(k) per observation. The
  reviewer's own advice was to take this only after a measurement, and there is no caller
  to measure yet. Recorded as a follow-up.
- **A property test over the walk** (partition-invariance, disjointness, one-locus-per-
  observation) and a **k = 3,000 fixture**. Both are right; partition-invariance is spec
  §15's regression anchor and the plan puts it at **E4**, over the whole builder
  arrangement, where it can catch what a unit-level version cannot. Recorded there rather
  than duplicated here.

## 6. What's good

- The differential harness is the strongest evidence produced anywhere in this milestone,
  and it was built rather than argued: production's grouping ported into a test and run
  against ng's on 5,000 cohorts, with the boundary case counted rather than assumed.
- Both agents proved every "no behaviour change" verdict on a constructed fixture instead
  of inferring it from a green suite — 3 such verdicts, none reported as a finding.
- The reliability agent labelled its own run partial and listed the four mutations it did
  not execute, so the three numbers mean what they say.
