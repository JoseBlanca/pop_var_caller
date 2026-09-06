# Review — window coverage C1: the reference into the merge's cache

**Date:** 2026-09-06
**Reviewed:** commit `c8d7f82f` on `ng-window-coverage` (plan step C1 of
[window_coverage.md](../../ng/impl_plan/window_coverage.md))
**Impl report:** [ng_window_coverage_c1_2026-09-06.md](../implementations/ng_window_coverage_c1_2026-09-06.md)
**Per-category files (audit trail):** `tmp/review_2026-09-06_window-coverage-c1/`

## 1. Scope

- **In scope:** `src/ng/run/cohort_merge/observation_cache.rs`, `src/ng/run/mod.rs`,
  `src/ng/run/cohort_merge/{serial.rs,parallel.rs,mod.rs}`, `src/ng/run/{callers.rs,psp_caller.rs}`,
  and every number the step's prose claims.
- **Out of scope:** production (`src/sample_summary/`, `src/paralog/`, `src/var_calling/`);
  `src/ng/window_coverage/`, landed at Milestones A and B; the four checks red on `main`.
- **Categories dispatched:** two agents, each in its own worktree detached at `c8d7f82f` — one on
  `reliability` and `errors` plus the claims check, one on `naming`, `idiomatic`, `smells`,
  `refactor_safety` and `module_structure`.

## 2. Verdict

**Request changes, and they are made.** The step's arithmetic is right and the oracle green, but
the ground it fetched was narrower than the records the cache hands out — which the next step is
built on. Both reviewers found it independently, from different checklists, each with tests that
pass on the reviewed commit.

## 3. Top three

1. **The ground fetched was not the ground the cache holds, and it fails at both ends on ordinary
   input.** Past the region: every sample is drawn one record *beyond* the reach — the only way to
   know a record is past it is to draw it — and the ground stopped at the reach, so the last
   record every sample draws in every cover had no base. Before the region: a record drawn by an
   earlier cover is still held and handed to this cover's builders, and the buffer is refilled
   rather than appended to, so the bases the earlier cover fetched are gone. Three tests written
   against the reviewed commit demonstrated both. **Fixed:** the ground is now every position of
   every record any sample holds on the cover's contig, widened to the region either way, and a
   test walks the held records and asserts a base at each.
2. **The merge's accessor never released what it had walked past, and the type chosen made
   releasing impossible.** `Box<dyn RefSeq>` has no `evict_before`; the run's padding accessor
   calls it after every fetch for exactly this reason. A windowed reader extends its buffer while
   each request lands near the last — the merge's pattern precisely — so a merge that never
   released would end a contig holding every base it passed: about 250 MB on human chromosome 1
   against a walk that otherwise peaks near 25 MB, measured in `ref_seq.rs`'s own test. **Fixed:**
   the cache asks for `MergeReference: RefSeq + EvictableRefSeq` and releases before the ground
   each cover goes on to read, with a fixture reference that counts what it was told to release.
3. **Three of the five tests could not fail on what their names claim.** The failure's region was
   erased by the fixtures' conversion, so a mutation making a failure name a one-base region
   instead of the ground survived; the contig-guard fixture put its second cover at higher
   coordinates, where the position guard answers `None` on its own, so deleting the contig guard
   survived; and nothing drove `cover_in_parallel` at all, so deleting its fetch survived — in the
   path a real run takes on any multi-core machine. **All three fixed**, each with a fixture that
   separates what it is about from what sits beside it.

## 4. What's good

- **Both reviewers converged on the blocker from different checklists**, and both wrote the tests
  rather than arguing from the code.
- **The fixpoint mechanism the step rests on is correct**, checked against `sweep` and `draw_to`:
  the ground a cover reaches genuinely is not known until the drawing stops, so fetching after it
  is right. What was wrong was the inference drawn from it, not the mechanism.
- **`over` spells every field of both structs**, so a field added to either stops the build at
  construction. That half of the new-field discipline was already right; the evictors' `..` was
  the half that was not.
- **The `-> Result<GenomePosition, E>` annotation in `cover_in_parallel` was challenged and
  cleared**: removing it gives `E0277`, so the annotation and its comment stand.

## 5. Findings, and what was done with each

### Blocker

- **B1 — the fetched ground missed the overshoot record and the held prefix.** Fixed by
  `ground_this_cover_holds`, which widens the region to every record any sample holds on the
  cover's contig. `every_position_of_every_held_record_has_a_reference_base` walks a fixture with
  both shapes and asserts a base at each. `reference_base_at`'s doc no longer claims what it
  cannot: `None` now means a contig the merge has left, or ground no cover has reached.

### Major

- **M1 — the accessor never released.** Fixed with the `MergeReference` bound and one
  `evict_before` after each successful fetch; `a_cover_releases_the_reference_bases_the_merge_has_walked_past`
  pins that it happens and from where.
- **M2 — the failure's region was untestable.** The fixtures gained a failure type that keeps it,
  and `a_failed_fetch_names_the_ground_the_cover_held` asserts the whole region.
- **M3 — the contig-guard fixture could not separate the two guards.** Its covers now overlap in
  coordinate, which is the ordinary shape since positions restart at 1 on every contig.
- **M4 — nothing drove `cover_in_parallel`.** It now has its own ground test.
- **M5 — a failed fetch left the cache answering as though the cover succeeded**: `covered_to`
  advanced and the previous cover's buffer still readable. Fixed twice over — the buffer is
  invalidated before the fetch is attempted, and `covered_to` moves only after it succeeds — with
  a test that a failed cover leaves nothing readable and a retry restores it.
- **M6 — both evictors destructured `SampleWindow` with `..`**, which would silently absorb the
  two fields C2 adds and must drain, in one evictor or both. Every field is named at both sites.

### Minor

- **Mi1 — the shared tail of the two covers was written twice** and grows at every step of this
  plan. Extracted as `close_the_cover`, which is also where C2's observation pass lands.
- **Mi2 — `ReferenceUnreadable` was a public error type that was not an error.** Now
  `thiserror`-typed with a `#[source]` chain and `#[non_exhaustive]`.
- **Mi3 — `reference_base_at` was `pub` over state no public API can set.** Now `pub(super)`,
  matching `cover`, with a `not(test)` expectation that the compiler will ask back the moment C2
  supplies its caller.
- **Mi4 — the inverted-region `min` had no test.** The existing inverted-region test now asserts
  the ground starts at the region's lower bound.
- **Mi5 — a region whose first base is 0 fails the whole cover**, a new precondition on a
  `pub(super)` method. Said in the method's `# Errors`; unreachable from the shipped segmentation,
  which is 1-based.
- **Mi6 — "Both callers" was four callers**, two of which mint no padding accessor. The doc now
  says which mint one and why the two are not one accessor.
- **Mi7 — the 14.6 GB figure priced a different failure** — rebuilding a reader per region, not
  sharing one between callers. Replaced with what sharing actually costs: two callers at
  different offsets in one contig reposition each other, an `lseek` and a re-read per fetch.
- **Mi8 — "no test reaches this" was reached by a test in the same commit**, and the message it
  handed back asserted its own impossibility. Both corrected.

### Nits, applied

The subtraction in `reference_base_at` now names the guard that makes it total; the fixture
reference builds its four contigs with `vec![bases; 4]` rather than three clones.

## 6. The claims check

**Correct:** the module test count, 6,343 → 6,348 and "the five are the whole increase" (the
`#[test]` delta is exactly five), the four production construction sites, that nothing outside the
named files is touched, that `callers.rs` and `psp_caller.rs` were rustfmt-clean before and after,
and that the five reverted formatting hunks left every `cohort_merge` file's hunk count unchanged
(5→5, 4→4, 1→1, 0→0).

**Three were wrong:**

1. **"At the shipped 20 kb building region [the buffer] is about 20 kB per merge."** The shipped
   building region is **500 bases** (`DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN`, with no CLI
   override), so the buffer is about **500 bytes a cover** in the serial driver and `500 × threads`
   in the parallel one — about **4 kB at eight threads**. The figure also omitted the larger cost
   the step actually added, which is finding M1.
2. **"Thirty construction sites … two about the reference, twenty-eight that are not."** There are
   **54** fixture construction sites, and **five** tests about the reference. The argument for a
   test-only constructor is stronger at the right numbers, not weaker.
3. **"The three driver entry points carry the conversion bound."** **Four** public functions gained
   it: three in `serial.rs` and one in `parallel.rs`.

## 7. Validation after the fixes

In the container, on this worktree:

- `cargo test --lib --all-features` — **6,353 passed, 0 failed, 15 ignored** (6,343 before C1).
- `cargo test --lib --all-features cohort_merge::observation_cache` — 42 passed, 0 failed.
- `cargo clippy --lib --all-features` — 3 warnings, all `needless_lifetimes` in
  `src/ng/run/cohort_merge/`, all predating this branch; none in this step's code, and no
  dead-code or unfulfilled-expectation warning in either profile.
- `rustfmt --check` — every `cohort_merge` file is back at exactly the hunk count it had before
  this step (`observation_cache.rs` 5, `serial.rs` 4, `parallel.rs` 1, `mod.rs` 0), and
  `callers.rs` and `psp_caller.rs` are clean. **A whole-file `rustfmt` run recurses into a
  module's children**, which reformatted four files this step does not touch; that was reverted,
  twice, and the counts above are how it is checked.
- **The standing oracle, after the fixes**: 2,311 records, sha256 `84ad19c2…0590d` on both routes
  — unmoved, as it was before C1 and after C1's first cut.

## 8. Deferred, with a home

- **Observing each held record exactly once** — plan step C2. The ground now covers every held
  record, which is what C2 needs; what C2 still owes is that a record held across two covers is
  observed by one of them and not both.
- **A second FASTA handle per run** — three of the four merge entry points now open one. Harmless
  in count; worth a line in whatever prices file descriptors at three thousand samples.
