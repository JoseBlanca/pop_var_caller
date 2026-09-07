# window coverage — C1: the reference, into the merge's cache

**Date:** 2026-09-06
**Plan step:** [window_coverage.md](../../ng/impl_plan/window_coverage.md) Milestone C, step C1
**Spec:** [window_coverage.md](../../ng/spec/window_coverage.md) §3.2, §5
**Review:** [ng_window_coverage_c1_2026-09-06.md](../reviews/ng_window_coverage_c1_2026-09-06.md)
**Branch:** `ng-window-coverage`
**Builds on:** [B1](ng_window_coverage_b1_2026-09-06.md), [B2](ng_window_coverage_b2_2026-09-06.md)

## The answer

**The merge's cache now holds a reference accessor of its own and reads, once per cover, the
ground that cover's records lie on — into a buffer every sample will read by offset.** Nothing
reads that buffer yet except a test and the accessor beside it; the accumulator starts drawing on
it at C2.

**No VCF byte moved.** The standing oracle — six tomato accessions over the first two 100 kb
intervals, called from the CRAMs and from stored files — gives 2,311 records and sha256
`84ad19c22dd14de583cd85805dcd2e5169e799d7a63691c979b7fa43d400590d` on both routes, identical to
the baseline taken before A1 and re-taken immediately before this step.

## Assumptions and deviations, all minor and all recorded

1. **The fetch happens after the cover's fixpoint, not before it.** What the cache holds is not
   known until the drawing has stopped: an observation chaining past the region widens the reach,
   and each sample is then drawn one record past *that*. Fetching first would mean either a
   guessed margin or a re-fetch every time a sweep grew the reach. What this costs is that C2
   cannot observe a record at the moment it is drawn, which is how the plan's C2 text describes
   it; it observes after the fetch instead, in the same cover, which preserves the intent — every
   record observed exactly once, in coordinate order, against bases this cover read.
2. **The accessor is a `Box<dyn RefSeq + Send + Sync>` rather than a type parameter.** It is
   read once per cover and not once per position, so the indirection costs nothing measurable,
   and a second type parameter would have spread through every signature that names the cache —
   the two drivers, the parallel merge, both callers — as well as through thirty fixtures that
   want a hand-built reference rather than a FASTA on disk.
3. **A failed fetch is the cache's own failure, so it has its own type.** The cache is generic
   over what its readers refuse, and a reference fetch is not one of those. `cover` hands back
   `ReferenceUnreadable`, and the caller converts: a run converts it into
   `RunError::WindowCoverageGroundUnreadable`, which names the ground. The bound that asks for
   the conversion sits on `cover` rather than on the cache, so a caller that never covers never
   meets it.
4. **The fixtures reach the cache through `over_fixture`**, a test-only constructor that supplies
   a hand-built reference of four contigs of repeating `ACGT`. Passing a reference at each of the
   54 fixture construction sites would bury the ten tests that are about the reference among
   forty-four that are not; the alternative — letting a cache hold no reference — is the one shape
   that could measure no coverage without saying so.
5. **Each of the run's four merge entry points mints an accessor for the cache.** Two of them
   also write records and so mint it beside their padding accessor; the other two have no padding
   accessor at all. The two are not one accessor because two callers reading at different offsets
   in one contig would reposition each other's window — an `lseek` and a re-read per fetch. A
   fresh accessor shares the run's index and contig table and costs about 18 µs, the `open(2)`.
6. **One test's source error type changed from `Infallible` to `RunError`.** The cover can now
   fail on its own account, so a source that cannot fail no longer makes the merge infallible.

## Changes made

- **[`observation_cache.rs`](../../../../src/ng/run/cohort_merge/observation_cache.rs)** — the
  cache gains the accessor, the buffer and the position the buffer starts at;
  `fetch_the_ground_this_cover_reached` runs at the end of both covers; `reference_base_at` is
  how anything reads it; `ReferenceUnreadable` is what a failed fetch is; `over` takes the
  accessor and `over_fixture` is the fixtures' way in.
- **[`run/mod.rs`](../../../../src/ng/run/mod.rs)** — `RunError::WindowCoverageGroundUnreadable`,
  and the conversion from `ReferenceUnreadable` that a run makes.
- **[`serial.rs`](../../../../src/ng/run/cohort_merge/serial.rs)**,
  **[`parallel.rs`](../../../../src/ng/run/cohort_merge/parallel.rs)** — four driver entry points,
  three in the first and one in the second, carry the conversion bound through to their callers.
- **[`cohort_merge/mod.rs`](../../../../src/ng/run/cohort_merge/mod.rs)** — the fixtures' own
  reference, and the conversion their failure type needs in order to call `cover` at all.
- **[`callers.rs`](../../../../src/ng/run/callers.rs)**,
  **[`psp_caller.rs`](../../../../src/ng/run/psp_caller.rs)** — each of the four production
  construction sites mints the merge's accessor beside the padding one.

**Nothing under `src/sample_summary/`, `src/paralog/` or `src/var_calling/` is touched.**
Formatting-only changes that `rustfmt` made in code this step does not touch were reverted —
twice, because **a whole-file `rustfmt` run recurses into a module's children**, so formatting
`cohort_merge/mod.rs` also reformatted `build.rs`, `close.rs`, `serial.rs`, `parallel.rs` and
`observation_cache.rs`. What checks it is the per-file `rustfmt --check` hunk count against the
step's parent: `observation_cache.rs` 5, `serial.rs` 4, `parallel.rs` 1, `mod.rs` 0, all equal
before and after. `src/ng/run/cohort_merge/` is where three of the four checks red on `main` live
and where another branch is working, so tidying it unasked would widen the conflict surface for no
gain.

## What the review changed, and it was the ground itself

**The ground the first cut fetched was narrower than the records the cache hands out**, and it
failed at both ends on ordinary input. Each sample is drawn one record *past* the reach — the only
way to know a record is beyond the reach is to draw it — and the ground stopped at the reach, so
the last record every sample draws in every cover had no base. And a record drawn by an earlier
cover is still held and handed to this cover's builders, while the buffer is refilled rather than
appended to, so the bases the earlier cover read are gone. Both reviewers found it independently,
each with tests that pass on the first cut. **The ground is now every position of every record any
sample holds on the cover's contig**, widened to the region either way, and a test walks the held
records and asserts a base at each.

**And the accessor never released what the merge had walked past.** `Box<dyn RefSeq>` has no
`evict_before` at all, so the type made releasing unavailable rather than merely unused — while
the run's padding accessor calls it after every fetch for exactly this reason. A windowed reader
extends its buffer while each request lands near the last, which is the merge's pattern precisely,
so a merge that never released would end a contig holding every base it passed: about 250 MB on
human chromosome 1 against a walk that otherwise peaks near 25 MB. The cache now asks for
`MergeReference: RefSeq + EvictableRefSeq` and releases before the ground each cover goes on to
read.

Three of the first cut's five tests could not fail on what their names claimed, and all three are
repaired: the failure's region was thrown away by the fixtures' conversion, the contig-guard
fixture could not tell the contig guard from the position guard beside it, and nothing drove the
parallel cover — the path a real run takes on any multi-core machine. Two further fixes: a failed
fetch no longer leaves the previous cover's bases readable or `covered_to` advanced, and both
evictors name every field of `SampleWindow` rather than absorbing new ones with a `..`, which is
what C2's two new fields need.

## Tests added

Ten, all in the cache's own module — five with the step and five the review required:

| test | what it pins |
|---|---|
| `a_cover_holds_the_reference_bases_over_the_ground_it_drew` | every position of the region reads back the reference's own base — the fixture repeats `ACGT`, so a wrong offset lands on a wrong letter rather than on a plausible one |
| `the_ground_read_reaches_as_far_as_the_cover_drew_and_no_further` | an observation opening inside the region and reaching ten bases past its end: the base at the reach is there, the one past it is not, and neither is the one before the region's start |
| `a_position_on_the_contig_left_behind_is_not_read_from_the_new_ones_buffer` | two covers **overlapping in coordinate** on different contigs, which is the ordinary shape since positions restart at 1: the old contig's position comes back absent, and it is the contig guard rather than the position guard that says so |
| `a_cache_that_has_not_covered_yet_has_no_reference_base_anywhere` | before the first cover there is no ground, and asking gives `None` rather than an empty buffer's byte |
| `a_reference_that_cannot_serve_the_cover_ground_ends_the_cover_naming_it` | a cover on a contig the reference does not have fails through the conversion rather than being absorbed |
| `every_position_of_every_held_record_has_a_reference_base` | the property C2 reads: a fixture with both a record held from an earlier cover and one drawn past this cover's reach, and a base at every position of both |
| `the_parallel_cover_holds_the_reference_bases_over_the_ground_it_drew` | the production path on any multi-core machine reads its ground too |
| `a_cover_releases_the_reference_bases_the_merge_has_walked_past` | the release happens, and from the ground each cover went on to read — pinned with a fixture reference that records what it was told to release |
| `a_failed_fetch_names_the_ground_the_cover_held` | the whole region, through a fixture failure type that keeps it: the string-only one throws it away, so any change to what a failure names would pass unnoticed |
| `a_failed_fetch_leaves_no_ground_readable_and_a_retry_restores_it` | a cover that could not read its ground answers nothing rather than the previous cover's bases, and a cover made again restores it |

**The other thirty-two tests in the module, and every cache test in the two drivers, cover
against the fixture reference and all pass** — which is what says the fetch runs on every cover
without disturbing anything the cache already did.

## Validation results

In the container, on this worktree:

- `cargo test --lib --all-features` — **6,353 passed, 0 failed, 15 ignored** (6,343 before this
  step). The ten above are the whole increase.
- `cargo test --lib --all-features cohort_merge::observation_cache` — 42 passed, 0 failed.
- `cargo clippy --lib --all-features` — 3 warnings, all `needless_lifetimes` in
  `src/ng/run/cohort_merge/`, all predating this branch; none in the code this step adds.
- `rustfmt --check --edition 2024` — `callers.rs` and `psp_caller.rs`, which were clean before
  this step, are clean after it. `observation_cache.rs`, `serial.rs`, `parallel.rs` and
  `run/mod.rs`'s submodules were already failing before it and are left as they were.
- **The standing oracle, run before this step, after its first cut and after the review's
  fixes**: 2,311 records, sha256 `84ad19c2…0590d` on both routes all three times.

## Tradeoffs and follow-ups

- **The buffer is one cover's ground, not the genome's.** The shipped building region is 500
  bases, so it is about **500 bytes a cover** in the serial driver and `500 × threads` in the
  parallel one — about 4 kB at eight threads — plus whatever the held records widen it to. One
  buffer for the whole cohort, not one per sample, because the bases are the same for every
  sample and only which of them a sample's records fall on differs. The larger memory this step
  adds is the accessor's own window, which is why it now releases what the merge has walked past;
  plan step D3 measures the milestone.
- **A cover whose reach runs past the end of a contig would fail the fetch.** It cannot happen
  today: the reach comes from records, which are inside their contig. Plan step C3 pushes the
  reach half a window past the region's end deliberately, so the clamp against the contig's
  length belongs to that step, which is the one that creates the case.
- **What C2 still owes is that each held record is observed exactly once.** The ground now covers
  every record every sample holds on the cover's contig, so C2 will find a base at each of their
  positions — but a record held across two covers is in both covers' ground, and observing it
  twice would count its positions twice.
