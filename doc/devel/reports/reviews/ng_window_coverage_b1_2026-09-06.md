# Review — window coverage B1: which positions a record reports depth at

**Date:** 2026-09-06
**Reviewed:** commit `da0f6746` on `ng-window-coverage` (plan step B1 of
[window_coverage.md](../../ng/impl_plan/window_coverage.md))
**Impl report:** [ng_window_coverage_b1_2026-09-06.md](../implementations/ng_window_coverage_b1_2026-09-06.md)
**Per-category files (audit trail):** `tmp/review_2026-09-06_window-coverage-b1/`

## 1. Scope

One new file and one `mod` line.

- **In scope:** `src/ng/window_coverage/depth.rs`, and the added line in
  `src/ng/window_coverage/mod.rs`. For the claims check, the implementation report and the
  `#### Window coverage` block of `PROJECT_STATUS.md`.
- **Out of scope:** `src/sample_summary/`, `src/paralog/`, `src/var_calling/` (production,
  frozen); `src/ng/window_coverage/{accumulator.rs,production_parity.rs}` (landed at Milestone
  A); `src/ng/run/cohort_merge/` and `src/ng/locus_generation/`, read as the sources of the
  types this step uses; the four checks red on `main` before this branch.
- **Categories dispatched:** three agents, each in its own worktree detached at `da0f6746` —
  one on `reliability`, one on `naming`/`idiomatic`/`errors`/`defaults`/`smells`, one on
  `module_structure`/`refactor_safety` plus the check of every number the diff's own prose
  claims.

## 2. Verdict

**Approve with changes.** Two Blockers, both of the form the rubric files at that level: a code
path that gives wrong results without panicking, and no test that says so. Both are fixed by
tests, one of them written by a reviewer who first proved it fails on a mutant.

## 3. Top three

1. **No test covered a *generic* record that arrived kept — the psp-mode path for every
   deletion-widened record.** Both multi-base kept tests built a tract. Rewriting the kept arm to
   spread depth over the whole span left all ten tests green, while a kept generic record over
   `3:10-12` then reported `[(10,8),(11,3),(12,3)]` instead of `[(10,8)]`. That is spec §6 trap 2
   firing in psp mode alone: every deletion's interior counted once by the widened record and
   again by the interior's own records, so the two modes' windows part company at exactly the
   loci direct mode gets right. **Fixed** — the test the reviewer proved against the mutant.
2. **No test made `build` fail.** Replacing `build(body)?` with a form that returns `Ok(())` on
   error left the suite green. In psp mode `build` is the decode of a stored body, and a swallowed
   decode failure takes that record's positions out of the window's denominator and the sample's
   histogram while the run reports success. **Fixed**, and the test also pins that nothing is
   reported before the failure, so a caller's accumulator cannot be left holding half a record.
3. **The guard making the `summary` parameter safe was a `debug_assert!`, which this repo's
   release profile does not compile.** With it replaced by a no-op — what release builds — a
   record on `3:9000-9002` handed a summary claiming `3:10-12` reported 9000–9002 and the suite
   stayed green, against a doc comment promising every position is on the summary's contig. The
   repo has recorded this trap twice (`locus_generation/mod.rs:90-93`) and put a *release*
   `assert!` on the cache's ordering check one file away. **Fixed:** both checks are release
   `assert_eq!`s now, each with a `#[should_panic]` test.

## 4. What's good

- **The three reviewers converged on one design point independently**, from three different
  checklists: the function was handed a summary while `Drawn::Kept` carries its own, and the two
  could disagree. Fixed by using the kept record's own — which is free, and widens the check from
  the region to the whole summary, so the head count a one-base record reports is checked too.
- **Both `match`es are wildcard-free.** A fourth `LocusKind` variant or a third `Drawn` shape is
  a compile error here rather than a silent fall-through — checked, not assumed.
- **The mutation harness the step shipped was reused by two reviewers unchanged**, and one of
  them re-ran all seven mutations and reproduced the exact failing set each was claimed to
  produce.
- **The file is where it belongs.** The reviewer who owns module placement checked the
  alternative and recommended against moving it: the rule is this module's spec §3.1, its only
  consumer is this module's accumulator, and `src/ng/calling/` already imports from
  `run::cohort_merge` at eight sites. What it asked for — and got — is that the dependency and
  its becoming mutual at C2 be written down rather than discovered.

## 5. Findings, and what was done with each

### Blocker

- **B1 — `depth.rs`: no test for a generic record arriving kept.** Fixed by
  `a_kept_generic_record_spanning_more_than_one_base_reports_at_its_anchor_only`, which asserts
  both the anchor-only rule and that the two draw shapes agree on it. Mutation
  `the_kept_arm_ignores_the_records_kind` now fails it.
- **B2 — `depth.rs`: no test for a failing `build`.** Fixed by
  `a_build_that_fails_stops_the_walk_and_returns_its_error`. Mutation
  `a_failed_build_is_swallowed` now fails it.

### Major

- **M1 — the summary/record agreement was guarded only in debug.** Both guards are release
  `assert_eq!` now, with the reason in a comment beside them and the precedent named. Three
  `#[should_panic]` tests cover the three ways they fire; three mutations, one per check, each
  fail exactly one of them. **Categories:** reliability, refactor_safety, errors (convergent).
- **M2 — the only generic multi-base fixture could not tell "the anchor's depth" from "the
  span's largest".** Its depths were `[8, 3, 3]`, where the anchor is also the maximum. Fixed by
  `a_generic_records_anchor_depth_is_its_first_position_and_not_its_deepest`, over `[0, 5, 5]`.
  Mutation `anchor_takes_the_deepest_position` now fails it, alone.
- **M3 — the kept draw's own summary was discarded for the parameter.** Now used, with the
  parameter checked against it whole. **Categories:** errors, reliability, refactor_safety.

### Minor

- **Mi1 — the callback reported a `Position`, dropping the contig** the accumulator it feeds
  takes as its first argument. Now reports `GenomePosition`, which is also what
  `LocusSummary::start_position()` already returns — so the contig is no longer discarded and
  rebuilt.
- **Mi2 — "the first base of this record" was derived twice**, from the summary at the entry and
  from the record inside the helper. The helper now takes it as an argument: one derivation.
- **Mi3 — an inverted region reported nothing, silently.** Pinned by
  `a_region_that_ends_before_it_starts_reports_no_position`, with the reason in the function's
  doc: such a record describes no ground, `GenomeRegion::len` saturates to zero, and nothing
  mints one today.
- **Mi4 — the name `for_each_depth_the_record_reports` is a clause.** Renamed
  `for_each_reported_depth`; the call site says `depth::` and the first argument says which
  record.
- **Mi5 — no property test over what is a pure function on a small domain.** Added, over span,
  witness layout, kind and read count: the two draw shapes agree, positions strictly increase,
  every position lies inside the record's region, and how many there are is decided by the span
  and the kind.
- **Mi6 — three doc claims.** "The great majority of positions" asserted a size with no number
  and no source: the mechanism is now stated without the size. "Which every record's head carries
  in both modes" was wrong — direct mode has no head, it derives the count by walking the record.
  And spec §6 trap 1 was cited as the source of an argument it does not make; it is now cited as
  the analogy it is.
- **Mi7 — `Drawn::Kept { body, .. }` absorbed the variant's own summary.** Now `summary: _`, so a
  third field on that variant is a compile error at the one site that walks every drawn record.
- **Mi8 — `pub mod depth;` diverges from the sibling `mod accumulator;` plus a targeted
  re-export.** The report's stated reason was inaccurate. Kept `pub` — the two items have no
  caller, so a private module is a `dead_code` warning — with the end state written down on the
  declaration itself: `mod depth;` plus a `pub(crate) use` at C2.
- **Mi9 — the module's dependency direction was undocumented.** One paragraph in `mod.rs` now
  says that `depth` names two cohort-merge types, that C2 makes the dependency mutual, and what
  that costs.

### Nits, applied

`along_the_locus` and `at_the_anchor` are half names — now `depth_along_the_locus` and
`depth_at_the_anchor`; the test helper `reported` was shadowed by a local of the same name in six
tests — now `collect_reported`; two clear-writing Rule-6 forms deleted; a comment naming the
bound that makes the position arithmetic safe.

### Not fixed, with reasons

- **A `num_obs_along_locus_into(&mut Vec<u32>)` variant**, so the cover could keep one buffer for
  the run instead of allocating per wide record. The allocation is in `locus_generation`, which
  this step does not touch; carried to plan step C2, where the caller and the buffer would live.
- **Recycling the body built for a kept record** back into the cache's spare list. Same home: the
  function's shape is what a caller would have to change, and whether it is worth anything is a
  C2 measurement.

## 6. The claims check

Every number the diff's prose makes about its own reach was re-derived by a reviewer, including
running the suite at `da0f6746~1` for the "before" counts. **Four were wrong and are corrected in
the implementation report:**

1. "Nothing outside `src/ng/window_coverage/` is touched at all" — the commit touches four files;
   two are documentation. True of source, and now says so.
2. "The cache already derives **and holds** that summary at the draw" — the summary at
   `observation_cache.rs:921` is a local for the ordering check and is dropped; the one the cache
   holds is derived a second time at `observation_cache.rs:875`. So a built record's sequences are
   walked twice today, and passing the summary in avoids a **third** walk, not a second.
3. "Which every record's head carries in both modes" — direct mode has no head; `LocusSummary::of`
   derives the count by walking the record's observations. The summary is available in both modes;
   a head is not.
4. "One summation cheaper" / "one pass over a short vector" — the two phrasings disagreed with
   each other and both understated: `num_obs_along_locus()` allocates a vector the length of the
   span and writes at every position for every observation.

Two smaller corrections: "3 fail, all three one-base tests" is loose — a fourth test uses a
one-base record and survives, because its head and body both say 0; and an in-code comment said a
head count of 99 was "none of 8, 3 or 11", where 11 is not a quantity anything computes.

**Everything else checked out**: both test counts (50→68 in the module, 6,325→6,343 in the
library, both re-run by the reviewer at the parent commit), every row of the test table, all seven
original mutations with their exact failing sets, the fixture depths, the cited line
`observation_cache.rs:921`, and the reasoning that the calling oracle cannot have moved.

## 7. Validation after the fixes

In the container, on this worktree:

- `cargo test --lib --all-features` — **6,343 passed, 0 failed, 15 ignored** (6,325 before B1).
- `cargo test --lib --all-features ng::window_coverage` — **68 passed, 0 failed** (50 before B1).
- `cargo clippy --lib --all-features` — 3 warnings, all `needless_lifetimes` in
  `src/ng/run/cohort_merge/`, none in `src/ng/window_coverage/`; all three predate this branch.
- `rustfmt --check --edition 2024` clean on both changed files.
- **Thirteen mutations, each caught**, five of them by exactly one test — the seven the step
  shipped with plus six for the behaviour the review added.

## 8. Deferred, with a home

- **`mod depth;` plus a `pub(crate) use`** — plan step C2, when the caller arrives and the
  dead-code argument expires. Written on the declaration.
- **The kept arm's `assert_eq!` on the built body's ground becomes a `RunError`** — plan step C2,
  where the source doing the building is a psp reader and a mismatch is a fact about the file
  rather than a bug in this crate. `observation_cache.rs:930-934` already describes the shape.
- **A buffer for `num_obs_along_locus`, and recycling the built body** — plan step C2, priced
  there.
