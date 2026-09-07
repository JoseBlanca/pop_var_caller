# window coverage — B1: which positions a record reports depth at, and how much

**Date:** 2026-09-06
**Plan step:** [window_coverage.md](../../ng/impl_plan/window_coverage.md) Milestone B, step B1
**Spec:** [window_coverage.md](../../ng/spec/window_coverage.md) §3.1
**Review:** [ng_window_coverage_b1_2026-09-06.md](../reviews/ng_window_coverage_b1_2026-09-06.md)
**Branch:** `ng-window-coverage`
**Builds on:** [A1](ng_window_coverage_a1_2026-09-06.md), [A2](ng_window_coverage_a2_2026-09-06.md),
[A3](ng_window_coverage_a3_2026-09-06.md)

## The answer

**The rule that turns one drawn record into covered positions with a depth at each now exists as
one function, and its three branches are pinned by eighteen tests.** A record spanning one base
reports the count its summary already carries and decodes nothing; a generic record spanning more
reports at its first base only, because the positions inside it have records of their own; a
repeat tract reports at every position of its span, because nothing else covers that ground.

**Which branch applies is decided by the record's span, never by which shape the draw arrived
in.** That is what makes direct mode and psp mode compute the same number: a one-base record is
answered from its summary even in direct mode, where the evidence is already in hand and could
have been summed instead. Two tests do exactly that comparison — one over a fixed record, one
over a swept domain of spans and witness layouts — and both assert the two shapes agree.

Nothing calls this yet. The cache learns to at plan step C2.

## The defects the review found

**Two paths could have given wrong results with nothing to say so, and both were missing tests
rather than missing code.**

- **No test built a *generic* record that arrived kept**, which is the psp-mode path for every
  deletion-widened record. Both multi-base kept tests built a tract. A reviewer rewrote the kept
  arm to spread depth over the whole span: all ten tests stayed green, while a kept generic record
  over `3:10-12` reported `[(10,8),(11,3),(12,3)]` instead of `[(10,8)]`. That is spec §6 trap 2
  firing in psp mode alone — every deletion interior counted twice, in one mode and not the other.
- **No test made `build` fail.** Replacing `build(body)?` with a form that returned `Ok(())` on
  error also left the suite green. A swallowed psp decode failure removes that record's positions
  from the window's denominator and from the sample's histogram while the run reports success.

**And one guard was not compiled in the profile this repo runs.** The `summary` parameter's
precondition — it must be this record's — was a `debug_assert!`, and `[profile.release]` leaves
debug assertions off. With it replaced by a no-op, a record on `3:9000-9002` handed a summary
claiming `3:10-12` reported 9000–9002 with the suite green. Both checks are release `assert_eq!`s
now, matching what the cache's own ordering check does one file away.

## Plan

The plan's B1 asks for a pure function taking a `Drawn` and, for a kept record spanning more than
one base, the source's `build`, yielding the positions the record speaks for and the depth at
each, with fixture tests over each record shape in both draw shapes. Delivered as one new file,
[`src/ng/window_coverage/depth.rs`](../../../../src/ng/window_coverage/depth.rs), and one `mod`
line. All six of the shape × draw-shape cells now have a test; the sixth — a generic record
spanning more than one base, arriving kept — was the review's first Blocker.

## Assumptions and deviations, all minor and all recorded

1. **It reports through a callback rather than returning a collection.** The signature is
   `for_each_reported_depth(drawn, summary, build, report)`, where `report` is called once per
   position with a `GenomePosition` and a depth. A returned `Vec` would allocate at every
   single-base record, or would need a scratch buffer threaded through the cache to avoid it; the
   callback needs neither, and C2's caller can feed the accumulator directly from inside it. The
   tests push into a `Vec`, so nothing about the shape costs testability.
2. **The summary is passed in for a built record and taken from the draw for a kept one.** For a
   built record, deriving the summary walks every observation it carries — and the cache walks
   them twice already, once for the ordering check at
   [`observation_cache.rs:921`](../../../../src/ng/run/cohort_merge/observation_cache.rs) and once
   for the summary it holds at
   [`observation_cache.rs:875`](../../../../src/ng/run/cohort_merge/observation_cache.rs) — so
   passing it in is what stops a third walk. A kept draw carries its own summary, which is this
   record's by construction, so that one is used and the parameter is only checked against it.
   Both checks are release asserts, for the reason above.
3. **A repeat cluster (`LocusKind::SsrBundle`) reports at every position, with the tract.** The
   spec's rule names two kinds and the bundle is neither; it is out of this plan's scope because
   bundle regions emit no loci today ([`locus_generation.md`](../../ng/spec/locus_generation.md)
   §11), so no such record can reach the function. Grouped with the tract because the argument for
   the tract is the argument for a bundle — region typing partitions the reference, so a repeat
   segment's record is the only record over its ground. The alternative, the generic arm, would
   silently under-cover every window overlapping a bundle if such records ever start arriving. A
   test pins which arm it takes.
4. **The depth even at a generic record's anchor comes from `num_obs_along_locus()`**, whose other
   entries are then discarded. Taking it straight from the observations would save an allocation
   the length of the span plus a write at every position for every observation — and would leave
   two derivations of "depth" that agree until a witness kind neither author thought about, which
   is spec §6 trap 1's failure one level down.
5. **A record whose region ends before it starts reports nothing.** `GenomeRegion::len` saturates
   to zero for such a region, so the depth vector comes back empty and both arms report nothing,
   while the span test — which normalises the ends — has already decided the record spans several
   bases. Nothing mints one today; the behaviour is pinned by a test rather than left to be
   discovered.
6. **The module is `pub` with no caller yet**, unlike its sibling `accumulator`, which is private
   behind one targeted re-export. Both of `depth`'s items would be `dead_code` warnings otherwise.
   The end state is written on the declaration: `mod depth;` plus a `pub(crate) use` at C2.

## Changes made

- **[`src/ng/window_coverage/depth.rs`](../../../../src/ng/window_coverage/depth.rs)** — new. The
  module doc states the three rules and why the span decides them; `for_each_reported_depth` is
  the entry point; `report_body_depths` is the private half that reads a built record's kind.
  Eighteen tests.
- **[`src/ng/window_coverage/mod.rs`](../../../../src/ng/window_coverage/mod.rs)** — one
  `pub mod depth;` line with the deferral written on it, and a paragraph in the module doc naming
  which way this module depends: `depth` uses two cohort-merge types, and at C2 the merge calls
  back into this module, so from C2 the two import each other.

Two of the four files this commit touches are documentation. **Nothing under
`src/sample_summary/`, `src/paralog/` or `src/var_calling/` is touched, and no source file outside
`src/ng/window_coverage/` is.**

## Tests added

Eighteen, all in the new file — ten written with the step, eight the review required:

| test | what it pins |
|---|---|
| `a_record_spanning_one_base_reports_the_head_count_and_not_its_evidence` | the summary says 7, the record's own observations sum to 9, and 7 is reported — so a branch reading the body here fails rather than agreeing by coincidence |
| `a_kept_record_spanning_one_base_builds_no_body` | the `build` handed in panics; the one-base rule never calls it |
| `a_generic_record_spanning_more_than_one_base_reports_at_its_anchor_only` | one position reported out of a three-base span, at the depth `num_obs_along_locus` gives there (8, against 3 at the two interior positions) |
| `a_generic_records_anchor_depth_is_its_first_position_and_not_its_deepest` | depths `[0, 5, 5]`: the anchor is reported at 0, so "the span's largest" fails where the fixture above cannot tell the two apart |
| `a_repeat_tract_reports_every_position_of_its_span` | three positions at depths 3, 7, 7 — asserted in order, so a flat spread, a reversal and an anchor-only rule all fail |
| `a_kept_record_spanning_more_than_one_base_is_built_through_the_source` | the body range the draw named (`4096..4224`) is the one `build` is asked for |
| `a_kept_generic_record_spanning_more_than_one_base_reports_at_its_anchor_only` | the psp-mode path for a deletion-widened record — spec §6 trap 2, in the one mode it can fire in |
| `both_draw_shapes_report_the_same_depths_for_the_same_record` | mode equivalence on one record: `Built` and `Kept` + `build` give identical output |
| `a_build_that_fails_stops_the_walk_and_returns_its_error` | the error comes back untouched, and nothing was reported before it |
| `a_record_with_no_reads_is_reported_at_depth_zero` | a covered position with no reads is still reported, at both spans — dropping it would shrink the window's denominator to the positions that happened to have depth |
| `a_partial_witness_raises_only_the_positions_it_witnessed` | a witness stopping inside a tract raises its last two positions and not the first four |
| `a_repeat_bundle_reports_every_position_of_its_span` | assumption 3 above |
| `the_span_decides_the_rule_and_not_the_records_kind` | a one-base tract is answered from its summary, so the two modes cannot answer differently by having a body in hand |
| `a_summary_that_is_not_the_records_own_is_refused` | the release check on a built record |
| `a_kept_draw_whose_own_summary_disagrees_with_the_parameter_is_refused` | the release check on a kept draw |
| `a_built_body_that_covers_other_ground_than_its_head_claimed_is_refused` | the release check on a stored head against the body behind it |
| `a_region_that_ends_before_it_starts_reports_no_position` | assumption 5 above |
| `every_reported_position_is_inside_the_record_and_the_two_draw_shapes_agree` | a property test over span 1–8, witness layout, kind and read count: the two shapes agree, positions strictly increase, each lies inside the region, and the count is decided by the span and the kind |

**Thirteen mutations, each caught, run in the container on this tree** — the harness is
`tmp/b1_mutations/` (a Python patch file plus a shell runner, because `scripts/dev.sh` forwards
only `CARGO_*` and `RUST*` variables, so a harness passing patterns through the environment
applies nothing and every mutation appears to survive):

| mutation | result |
|---|---|
| the one-base rule answers from the body instead of the summary | 3 fail — the three tests that pin the one-base rule |
| a generic record spreads its depth over its whole span (spec §6 trap 2) | 4 fail, the property test among them |
| a tract reports at its anchor alone | 7 fail |
| a repeat cluster joins the generic arm | `a_repeat_bundle_reports_every_position_of_its_span`, alone |
| the depth along a locus replaced by one pooled count spread flat | 4 fail |
| every reported position shifted one base right | 7 fail |
| the anchor's depth taken from the last position of the span | 3 fail |
| the anchor's depth taken as the span's largest | `a_generic_records_anchor_depth_is_its_first_position_and_not_its_deepest`, alone |
| a failing `build` swallowed, the record reporting nothing | `a_build_that_fails_stops_the_walk_and_returns_its_error`, alone |
| the kept arm ignores the record's kind and spreads over the span | 2 fail — the kept-generic test and the property test |
| the release check on a built record's summary removed | `a_summary_that_is_not_the_records_own_is_refused`, alone |
| the release check on a kept draw's summary removed | `a_kept_draw_whose_own_summary_disagrees_with_the_parameter_is_refused`, alone |
| the release check on a built body's ground removed | `a_built_body_that_covers_other_ground_than_its_head_claimed_is_refused`, alone |

## Validation results

In the container, on this worktree
(`/Users/jose/devel/pop_var_caller-window-coverage/scripts/dev.sh`):

- `cargo test --lib --all-features` — **6,343 passed, 0 failed, 15 ignored**, against 6,325 before
  this step. The eighteen above are the whole increase.
- `cargo test --lib --all-features ng::window_coverage` — **68 passed, 0 failed** (50 before).
- `cargo clippy --lib --all-features` — 3 warnings, all `needless_lifetimes` in
  `src/ng/run/cohort_merge/`; none in `src/ng/window_coverage/`, and all three predate this
  branch.
- `rustfmt --check --edition 2024` clean on both changed files.

**The standing calling oracle was not re-run for this step and cannot have moved**: nothing calls
the new function, so no record and no VCF byte can change. The plan runs it at every step of
Milestone C, which is where a byte can start to move.

**Four checks are red on `main` before this branch and are left alone**, three of them in
`src/ng/run/cohort_merge/` where another branch is working: nine files fail `cargo fmt --check`;
the three clippy lints above; the integration test
`a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele` fails; and the example
`ng_candidate_selection_probe` does not compile.

## Tradeoffs and follow-ups

- **A built body is dropped, not recycled.** In psp mode a record spanning more than one base
  builds a `SampleLocusObservations` that this function drops when it returns. The cache has a
  recycling channel for exactly such records — `SampleWindow::spare`, offered back to the source
  on the next draw — and the function's shape gives the caller no way to reach the built record.
  Whether that costs anything is a C2 question, and it is one allocation per wide record, not per
  position.
- **`num_obs_along_locus()` allocates a `Vec` per wide record**, and writes at every position of
  the span for every observation. It is production's API, reused as spec §7 asks; a
  `..._into(&mut Vec<u32>)` variant would let the cover keep one buffer for the run. Carried to
  C2, where the buffer would live.
- **The kept arm's ground check should become a `RunError` at C2.** A stored head claiming other
  ground than the body behind it is a fact about the file, not a bug in this crate, and
  `observation_cache.rs:930-934` already describes the shape such a check takes once a psp source
  is the one building.
- **The head-equals-body claim is still a claim.** The one-base branch rests on it, and B2
  measures it on a real store. A counter-example changes spec §3.1 to "build every body" before
  Milestone C begins.
