# Fixes applied — ng cohort merge, step A3

*2026-08-17, branch `ng-cohort-merge`. Input:
[the A3 review](ng_cohort_merge_a3_2026-08-17.md) — 5 Major, 11 Minor, 8 Nits over five
category checklists. Every finding is accounted for below.*

## Findings table

| ID | Title | Severity | Decision | Status |
|---|---|---|---|---|
| M1 | `span()` reintroduces the coordinate-ceiling defect | Major | Apply | **Applied** |
| M2 | an unsorted sample yields a locus not containing its member | Major | Apply | **Applied with adaptation** |
| M3 | members' "ascending sample order" untested (mutant survived) | Major | Apply | **Applied** |
| M4 | the `u32` saturating sum untested (mutant survived) | Major | Apply | **Applied** |
| M5 | k-sample scan through scattered pointers; one scan discarded per locus | Major | Apply | **Applied** |
| Mi1 | the tie-break doc names the wrong mechanism | Minor | Apply | **Applied** |
| Mi2 | "production's group walk" omits two real differences | Minor | Apply | **Applied** |
| Mi3 | the non-reference-total test does not discriminate the rules | Minor | Apply | **Applied** |
| Mi4 | test helpers invert the module's own vocabulary | Minor | Apply | **Applied** |
| Mi5 | the `9` behind a test's argument hidden in a helper | Minor | Apply | **Applied** |
| Mi6 | `per_sample`, `earliest_head`, `opened_at` misname their contents | Minor | Apply | **Applied** |
| Mi7 | `members` as a `Vec` — cost unnamed, growth unreserved | Minor | Apply | **Applied** |
| Mi8 | two `.expect`s restating one invariant | Minor | Apply | **Applied** |
| Mi9 | the progress guarantee undocumented and unasserted | Minor | Apply | **Applied** |
| Mi10 | `quiet_locus` rebuilds a field its base already set | Minor | Apply | **Applied** |
| Mi11 | parallel per-sample vectors unchecked by any type | Minor | Dispute | **Disputed** |
| Nits | eight | Nit | Apply / Defer | **5 Applied, 3 Deferred** |
| — | a heap for the argmin; a property test; a k = 3,000 fixture | — | Defer | **Deferred, with homes** |

## The five Majors

**M1 — `span()` no longer calls `region.len()`.** It subtracts first and adds after, both
saturating. The defect was one step old: A2 documented that `GenomeRegion::len` overflows
its `+ 1` at the coordinate ceiling — a debug panic and **0 in release** — and named that
as the reason `SampleLocusObservations::reach` avoids it. `span()` called it anyway, and
`span()` is what A4's width verdict reads, so the widest locus expressible would have been
judged narrow enough to build. Pinned by
`the_span_is_the_base_count_at_the_coordinate_ceiling`.

**M2 — the coordinate-order precondition is now enforced, at release level.** Handed a
sample whose observations arrive out of order, the walk used to close a locus over ground
it does not cover: `[50–50, 10–10]` came out as one locus spanning 50–50 holding an
observation forty bases outside it. **Adaptation:** the reliability agent proposed a
`debug_assert!` in `over` scanning the whole input; the refactor_safety agent proposed a
release `assert!` inside the walk. Took the second, for its reasons — the release profile
is the one this repo runs, the check is O(1) per observation against work that already does
a per-observation sequence comparison, and it fires on the first observation out of order
rather than after a scan of input that may be large.

**M3 and M4 — the two surviving mutants are now killed.** A test where the walk consumes
sample 1 before sample 0, and a test at the `u32` ceiling.

**M5 — the walk's layout follows the sibling merge.** `keys` beside the cursors, refreshed
only for the sample just advanced, so picking the next observation scans one contiguous
array instead of jumping into k observations across the heap; the head returns with its
sample, so no caller looks it up twice; `cursors_at_open` is a reused field rather than a
per-locus `Vec` clone; and the member vector is sized exactly during the walk. The O(k)
comparison count is unchanged — see *Deferred*.

## Disputed

**Mi11 — fold the four parallel per-sample vectors into one vector of a struct.** The
finding is right that nothing but `over` ties their lengths. But the reason `keys` exists
at all is that the argmin scans it and nothing else, and that is contiguous only while it
is its own array; folding it into a struct-of-fields per sample would undo M5. The sibling
`MergedCursors` makes the same trade with the same three parallel arrays. Recorded in the
type's doc so the next reader meets the reasoning rather than the risk.

## Deferred, each with a home

- **A heap or loser tree for the argmin.** The scan is still O(k) per observation, which at
  the top of the committed range is 3,000 comparisons per observation. The reviewer's own
  sequencing was: keys first, a heap only if a measurement says the scan still dominates.
  There is no caller to measure until C2. **Home: the region-width sweep (spec §14 Q1),
  which is the first thing that runs this at scale.**
- **A property test over the walk** — partition-invariance, disjointness, exactly one locus
  per observation. Right, and partition-invariance is spec §15's regression anchor. **Home:
  E4**, where it runs over the whole builder arrangement and can catch what a unit-level
  version cannot.
- **A k = 3,000 fixture.** The agent ran it as a probe and it passes. **Home: E4**, with the
  property test.
- Three nits: `#[must_use]` on the iterator struct; `span()` having no caller until A4;
  `over` taking `impl IntoIterator`.

## Validation

In the container after the fixes: `cargo fmt --check` clean;
`cargo clippy --lib --all-features -- -D warnings` clean;
`cargo test --lib ng::run::cohort_merge` **25 passed, 0 failed** (16 before). Full-suite
figures are in the commit message.

**Mutants re-run after the fix:** the "match production" mutation
(`saturating_add` → `max`) now fails three tests, including the new §15 fixture; it
previously failed only an assertion buried in a test named for something else.
