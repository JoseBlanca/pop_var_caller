# Code-quality map: the cohort merge, the hidden-duplication filter's wiring, and its statistics

Area: `src/run/cohort_merge/` (9 files plus `observation_cache/round_view.rs`), `src/run/paralog_filter/`
(9 files) and `src/paralog/` (9 files). Read from code only; nothing was built or run.

## Summary

1. **`call-from-alignments` cuts its building regions narrower as the thread count rises, and then
   builds them one at a time.** `round_shape_for` sizes the regions for a driver that builds
   several at once, but the alignments path always takes the one-region-at-a-time driver and
   ignores the "regions in flight" count. At 48 threads a run pays for about 32 covers where 4
   threads pay for 4, with no parallel building in return. The width also varies with
   `--threads`, which is the condition under which the last thread-count bug in this area
   appeared. Unmeasured, but it is the finding with the most direct effect on runs.
2. **The calling layer imports its input types from `run::cohort_merge::build`.** Six data
   types (`CohortObservation`, `SampleSupport`, `SupportedAllele`, `AlleleSupport`,
   `PartialObservation`, `LocusBorder`) plus the `MinAltReads` threshold are the only reason
   `src/calling/` depends on `crate::run`: 12 files and about 25 import lines. Moving those
   definitions to a module below both layers removes that dependency. The constructors can stay
   in `build.rs`.
3. **About 900 lines of merge code exist only for tests, and the production frontier rule has no
   direct tests.** `organise.rs`'s `Organiser` is test-only, and its 36 tests check the
   locus-ownership rule, so that rule is tested on code no run executes. Production applies the
   same rule through `rounds::Frontier`, which no test calls directly. `rounds.rs`'s
   `merge_and_call_in_rounds` is not test-gated but is used only as a reference in tests.
4. **Comments from before the psp path still describe it as future work, and five places cite a
   deleted file.** Five comments say a corrupt `.psp` "must become a `RunError`" once
   observations are decoded from psp files. That decoding exists, and `psp_source` already
   refuses that kind of damage. Five other places refer to `parallel.rs` or
   `merge_cohort_in_parallel`, which no longer exist.
5. **Some code is duplicated:**
   - the window lookup appears twice: on the cache and on the sealed round;
   - the per-region build closure appears twice: in `rounds` and in `stream`;
   - the spill read-ahead thread appears twice: in pass two and pass three;
   - the hidden-duplication calibration is split across `calibrate.rs` and `calibration.rs`.
6. **The resume file format's histogram codec is in the wrong file.** `paralog_filter::resume`,
   a format a later run must read, borrows it from `cohort_merge::recorded_windows`, which is a
   file written only for a debug probe.

The two `crate::cli` references in `stream.rs` and `rounds.rs` occur only in doc comments; there
is no code dependency. `src/paralog/` is clean: it imports nothing from `run` or `cli`.

---

## Bugs

### B1. The alignments path narrows its regions by thread count but never builds them in parallel

- **Where**
  - `src/cli/calling_run.rs:485-519` (`round_shape_for`, the `Alignments` arm)
  - `src/cli/calling_run.rs:1505-1508`
  - `src/cli/call_from_alignments.rs:926-927`
  - `src/run/callers.rs:867` (always `call_cohort_from_sources_handing_each_record_over`)
  - Compare `src/run/psp_caller.rs:856-875`, which branches on `regions_in_flight > 1`.
- **What.** For alignments, `round_shape_for` returns
  `in_flight = clamp(ground / 500, 1, workers)` and `width = ground / in_flight`. Only the psp
  path reads `in_flight`. The alignments path runs the serial driver
  (`serial::merge_cohort_handing_each_locus_over_covering_samples_in_parallel`). That driver
  covers, evicts and builds one building region at a time.
  - The narrowing therefore buys no parallelism. It multiplies the number of covers by up to
    `workers`. Each cover sweeps every sample, fetches the reference, sorts contigs and folds
    window coverage.
  - The doc at `calling_run.rs:466-472` ("one region a thread … up to 5.5× faster than one region
    at a time") describes the round driver, which `rounds.rs:4-6` says no calling run uses any
    more.
  - The width also depends on `--threads`. Output is meant to be width-invariant, and tests
    pin that. The HG002 window-coverage defect recorded at `observation_cache.rs:659-667`
    (nine records unscored at two threads and scored at one) was a failure of exactly this
    invariance.
- **Scenario.** At 1 sample, `round_width_for` gives 16,000 bases.
  - At 48 threads that is 32 regions of 500 bases: 32 covers per 16 kb.
  - At 4 threads it is 4 regions of 4,000 bases.
  - The VCF should be identical; the 48-thread run does about 8 times the cover work.
- **Proposal.** Choose one of two fixes:
  - Make the `Alignments` arm return `(ground, ONE)`, so the region width does not depend on
    the thread count.
  - Route alignments through `call_cohort_from_sources_in_rounds` as psp does, if building
    regions in parallel pays there. The doc says the measurement on a 48-accession run went
    against it.

  Either way, fix the stale doc.
- **Downside.** Wider covers hold more ground at once. At high thread counts, the first fix
  raises peak memory back to what one region of `ground` costs, which is the pre-rounds
  baseline. The current narrow regions may be lowering memory by accident.
- **Effect.** No output change, if width-invariance holds. Speed probably improves at high
  thread counts. **Memory rises** at high thread counts with the first fix, so it needs an
  owner decision and a measurement.
- **Cost.** S for the code, plus one measurement (M).
- **Confidence.** High that `in_flight` is unused on this path. Medium on the size of the cost.

### B2. The tests of the frontier rule exercise a test-only copy of it

- **Where**
  - `src/run/cohort_merge/organise.rs:115-492`: `Organiser`, every item `#[cfg(test)]`; 36
    tests at `:511+`.
  - `src/run/cohort_merge/rounds.rs:269-324`: `Frontier`, the production rule used by
    `stream.rs:169`.
  - `src/run/cohort_merge/serial.rs:1095-1110`: `refuse_displaced_loci` checks the serial driver
    with `Organiser`.
- **What.** Spec §6.1's rule says a locus that begins on ground an earlier locus owns is dropped
  and counted. It is implemented twice.
  - **`Organiser::resolve_and_release`** is test-only and has 36 tests.
  - **`Frontier::hand_over`/`claim`** is what the streaming run uses. Nothing calls it directly
    in a test. The stream's tests (`stream.rs:609-1000`) and the round driver's tests reach it
    only through fixtures in which nothing is displaced. So the displacement branch of `claim`
    (`rounds.rs:313-322`) is never executed.
  - Production asserts that the count is zero (`callers.rs:7779-7785`).
    That assertion is the only thing standing behind the branch.
- **Scenario.** A regression that inverted the `<=` in `Frontier::claim` would drop every second
  locus at each region edge where loci touch. `callers.rs` would then panic on the non-zero
  count. The same regression in `Organiser` would be caught by its own tests. One regression is
  caught in tests, the other only in a run.
- **Proposal.** Port the displacement tests to `Frontier`; most of the fixtures carry over. Make
  `refuse_displaced_loci` use `Frontier`. Delete `Organiser`, `RegionIndex`, `MergeTally` and
  `RunEndedShort` (about 380 lines plus 36 tests, some of which test reorder-buffer behaviour no
  longer needed). Move `first_base_of` and `last_base_of` next to `Frontier`.
- **Downside.** The `Organiser` tests also document the reorder buffer and the ordered release,
  which the stream now does with channels and a `VecDeque`. Deleting them removes that
  documentation. The stream's own ordering tests cover the behaviour, so the loss is in
  explanation, not in coverage.
- **Effect.** None.
- **Cost.** M.
- **Confidence.** High.

### B3 (suspected, latent). The window search assumes each sample's reaches never decrease

- **Where**
  - `src/run/cohort_merge/observation_cache.rs:2099-2101` (`first_reaching_summary`, a
    `partition_point` on `reach_position() < position`)
  - Used by `with_observations` (`:1219`), `evict_before` (`:721`), `seal_round`
    (`round_view.rs:63`) and `window_starts_into`.
- **What.** `partition_point` gives a meaningful answer only if the predicate flips once. That
  holds when one sample's records are disjoint and ascending, which `draw_next` checks only for
  starts (`:2005-2011`).
  - On the psp path, `psp_source.rs:518-540` refuses overlapping records, so the path is safe.
  - On the direct path, the only overlap check is in `build.rs:156-167`, and it runs only on
    loci that get built.
  - A generator defect that produced `[100–130], [105–105]` in one sample would make the window
    start too late. A record would silently leave the window, and eviction would drop it early,
    before any assertion fired.
- **Proposal.** Either check `previous.reach < start` in `draw_next`, beside the existing start
  check (one comparison per record), or document the precondition at `first_reaching_summary`.
  The check moves the failure from silent to loud for direct mode. The comment at
  `observation_cache.rs:1999-2003` argues against it, but the test it cites to justify that lives
  in the deleted `parallel.rs`.
- **Downside.** A direct-mode generator bug would abort a run rather than produce slightly wrong
  windows. That is the intended trade everywhere else in this module.
- **Effect.** None on correct input.
- **Cost.** S.
- **Confidence.** Medium. No producer is known to violate the precondition.

### B4 (suspected, latent). A cohort mixing whole records and kept bodies would be indexed wrongly

- **Where**
  - `observation_cache.rs:1605` (`keeps_evidence` is cohort-wide)
  - `:1226` (records offered only if no sample keeps evidence)
  - `:1392` (`build_at` always indexes `held_bodies`)
  - `:1947` (`held_observations.get(index)`, else `held_bodies[index]`)
- **What.** `ObservationSource` allows each draw to return `Drawn::Built` or `Drawn::Kept`. But
  the cache stores them in two vectors that are indexed in step with `held_summaries`. A source
  that mixed the two shapes would misindex `measure_coverage_on`. A cohort with one source of
  each kind would make `build_at` index an empty `held_bodies` and panic. No source mixes today:
  walkers build, and psp readers keep.
- **Proposal.** Either assert in `draw_to` that a sample never switches shape, or say in the
  trait doc that a source must pick one shape for its whole life.
- **Effect.** None.
- **Cost.** S.
- **Confidence.** Low that this is reachable today; high that the contract is unwritten.

No defect was found that changes a genotype, QUAL or record today. The concurrency code
(`stream.rs`, both paralog passes) was read for shutdown order, error propagation and ordering
(see "Looked at, left alone").

---

## Architecture

### A1. Move the merge's output types below both `run` and `calling`

- **Where**
  - `src/run/cohort_merge/build.rs:465` (`REFERENCE_ALLELE`)
  - `build.rs:1289-1810`: `CohortObservation`, `SampleSupport`, `SupportedAllele`,
    `PartialObservation`, `LocusBorder`, `AlleleSupport`, and the one accessor
    `flank_deletion_of`.
  - `src/run/cohort_merge/mod.rs:117-383`: `MinAltObs`, `MinAltReadShare`, `MinAltReads` and
    their defaults.
  - Importers: 12 files in `src/calling/`, plus `tests/ng_calling_loop_calls_genotypes.rs` and
    `tests/ng_candidate_selection_truth_recall.rs`.
- **What.** These are the interface between two steps: the merge produces them and calling
  consumes them. They are data with no behaviour that depends on the merge.
  - The constructors (`CohortObservation::over`, `over_from_cursors`) need `ClosedLocus` and
    `WindowedCohort`. Rust lets those `impl` blocks stay in `build.rs` even if the types move,
    because both are in the same crate.
  - `MinAltReads` is a threshold that both the merge's keep rule and calling's candidate
    admission apply.
  - These imports are calling's only dependency on `crate::run`: `grep 'crate::run::' src/calling`
    finds nothing outside `cohort_merge`.
- **Proposal.** Put the definitions in a new top-level module, for example
  `src/cohort_observation.rs`, with `MinAltReads` beside them or in `types`. Leave the
  constructors and the `AlleleTable` machinery in `build.rs`. The question in the brief was
  whether these are calling-layer vocabulary: they are the shared vocabulary of two layers, so
  they belong below both, not inside `calling` either.
- **Downside.** Each type's doc then sits away from the code that fills it. Several docs explain
  the field by how `build.rs` computes it (for example `SampleSupport::supported`'s ordering). A
  reader changing the derivation must remember to update a doc in another file. Moving the whole
  of `cohort_merge` to the crate root instead would also work, but it imports
  `run::explain::{DroppedLoci, MergeDrop, ExplainRegions}` and `run::side_pool::SidePool`, so
  those would have to move too. That is a larger change.
- **Effect.** None.
- **Cost.** M (about 25 import lines, plus doc links).
- **Confidence.** High.

### A2. The resume file's histogram codec lives in a probe file

- **Where**
  - `src/run/cohort_merge/recorded_windows.rs:271` (`write_the_histogram`) and `:339`
    (`read_the_histogram`)
  - Used by `src/run/paralog_filter/resume.rs:37, 75, 177`
- **What.** `recorded_windows.rs` exists for one probe: an environment variable makes a run
  write its windows for `examples/ng_window_coverage_probe.rs`. The histogram codec is
  production state: a resumed run reads `coverage-histograms.tsv` through it. So the filter's
  resume path depends on a cohort-merge debugging file. The codec serialises
  `window_coverage::SampleHistogram` and knows nothing of the merge.
- **Proposal.** Move the two functions next to `SampleHistogram` in `src/window_coverage/`, or
  into `paralog_filter::resume`. Have `recorded_windows` import them from there.
- **Downside.** None of substance. The two formats must stay identical, which one shared
  function guarantees wherever it lives.
- **Effect.** None.
- **Cost.** S.
- **Confidence.** High.

### A3. The round driver is compiled into the release build but used only as a test reference

- **Where:** `src/run/cohort_merge/rounds.rs:94-262` (`merge_and_call_in_rounds`, `Round`,
  `Round::work`); `rounds.rs:4-8` says so itself.
- **What.** About 170 lines of production-compiled code with no production caller. Being `pub`
  in a `pub mod`, it raises no dead-code warning. Its only users are `stream.rs`'s tests
  (`:734`, `:872`), which compare the stream against it. `CalledRegion` and `Frontier` in the
  same file are live.
- **Proposal.** Gate it with `#[cfg(test)]`. Alternatively, delete it and compare the stream
  against `serial::merge_cohort_handing_each_locus_over`, which those tests already import
  (`stream.rs:530`) and which is the documented oracle. Move `Frontier` and `CalledRegion` into
  `stream.rs` or a small `frontier.rs`.
- **Downside.**
  - The round driver is simpler to read than the stream, and the module doc argues for keeping
    it for that reason.
  - Deleting it loses one comparison: whether the stream reports the same *first error*
    (`stream.rs:711`, `:850`). The serial driver may stop at a different point.

  `#[cfg(test)]` keeps both properties.
- **Effect.** None. The binary gets smaller.
- **Cost.** S.
- **Confidence.** High.

### A4. `crate::cli` in `stream.rs` and `rounds.rs` is a doc reference, not a dependency

- **Where:** `stream.rs:30`, `rounds.rs:27`.
- **What.** Both are backquoted paths in `//!` text explaining who sizes the rounds. No `use`.
  The real coupling runs the other way: `cli/calling_run.rs:496` reads
  `cohort_merge::stream::ROUNDS_HELD`, and that direction is correct.
- **Proposal.** Leave the code. When B1 is settled, check that the sentence is still true.
- **Effect.** None. **Cost.** None. **Confidence.** High.

### A5. Module visibility says "pub until the caller objects land", and they have landed

- **Where:** `src/run/cohort_merge/mod.rs:31-36`.
- **What.** The comment justifies blanket `pub` because "the two caller objects that will own it
  do not exist yet". They exist (`run::callers`, `run::psp_caller`). Of the outside users, the
  integration tests need the A1 types plus `close::{ClosedLocus, SampleMembers, Verdict}` and
  `observation_cache::WindowedCohort`. Nothing else outside `run` needs to see `organise`,
  `serial`, `rounds` or `stream` internals.
- **Proposal.** After A1, narrow `cohort_merge`'s submodules to `pub(crate)` where no
  integration test or example needs them, and rewrite the comment.
- **Downside.** Integration tests that build loci by hand (`tests/ng_calling_loop_calls_genotypes.rs`)
  keep `close` and `observation_cache` partly public.
- **Effect.** None. **Cost.** S. **Confidence.** Medium (examples were checked by grep only).

---

## Duplication

### D1. Window lookup written twice: on the cache and on the sealed round

- **Where**
  - `observation_cache.rs:1172-1258` (`with_observations`, `window_starts_into`)
  - `observation_cache/round_view.rs:112-163` (the same two methods)
- **What.** The two copies match line for line in shape: compute `left_edge`/`last_base` with
  min/max, assert against `covered_to`, call `first_reaching_summary` per sample, and offer
  records only when no evidence is kept. They differ only in the per-sample source fields. That
  is about 45 lines each. A rule change made to one copy (the inverted-region handling, the
  coverage assertion) has to be repeated in the other.
- **Proposal.** Write one free function that takes the per-sample slices (summaries,
  observations, windows), `covered_to` and `keeps_evidence`, and builds the `WindowedCohort`.
  Both methods would call it.
- **Downside.** The two owners keep their data in different structs, so the helper needs closures
  or three parallel slices. Somewhat less obvious to read than two plain loops.
- **Effect.** None. **Cost.** S. **Confidence.** High.

### D2. The per-region build-and-call closure, written twice

- **Where**
  - `rounds.rs:219-250` (inside `Round::work`)
  - `stream.rs:488-522` (`build_one_region`)
- **What.** Both make `calls`/`failed`, fill `starts` with `window_starts_into`, call
  `with_observations` and `build_region_handing_over_windowed` with the same `keep` closure and
  `build_at(sample, starts[sample] + index)`, and return `CalledRegion`. The `+ starts[sample]`
  offset is the subtle bug fixed in both: the comment at `mod.rs:140-154` (the `KeepingSource` fixture) records that a defect
  there passed every test. It now has to stay right in two places.
- **Proposal.** If A3 gates the round driver, `build_one_region` becomes the only live copy. If
  both stay, make `build_one_region` generic over a small "round reader" trait implemented by the
  cache and by `RoundView`, so the round driver calls it too.
- **Downside.** The trait is a new abstraction with two implementations. A3 makes it
  unnecessary.
- **Effect.** None. **Cost.** S (after A3). **Confidence.** High.

### D3. Two read-ahead threads over the spill

- **Where:** `src/run/paralog_filter/pass_two.rs:286-325` and `pass_three.rs:267-312`.
- **What.** The two are about 40 lines each. They share the same loop: take a refill buffer,
  resize it, fill it until full or the stream ends, send it, stop on an error or a closed
  channel. They differ in the item type (`CompressedSpillBody` against `CompressedSpillEntry`),
  the batch size (2,048 against 256), and pass three's running `first` index.
- **Proposal.** One generic `read_ahead<T: Default, R>(next: impl FnMut(&mut T) -> Option<Result<(), SpillError>>, batch: usize, …)`
  that returns batches carrying `first`. Pass two would ignore `first`.
- **Downside.** Pass two and pass three may need to tune batching differently, and they already
  use different sizes. The generic version keeps the size a parameter, so that freedom stays.
- **Effect.** None. **Cost.** S. **Confidence.** Medium (the value is modest).

### D4. Two width-verdict functions for one rule

- **Where:** `src/run/cohort_merge/close.rs:294-308` (`judge_here`) and `:310-338` (`judge`).
- **What.** `judge_here` repeats `judge`'s first test (quiet first), then adds the
  "deferred → Build" case and calls `judge` with `true`. `judge` is reached in production only
  through `judge_here`; its other callers are tests (`:1918-1930`).
- **Proposal.** Fold them into one function with a `deferred` argument. Keep the long comment,
  which belongs to the one rule.
- **Downside.** None.
- **Effect.** None. **Cost.** S. **Confidence.** High.

### D5. `paralog/calibrate.rs` and `paralog/calibration.rs`

- **Where:** `src/paralog/calibrate.rs` (one 20-line function) and `src/paralog/calibration.rs`
  (the types it fills).
- **What.** Two files whose names differ by four letters. `calibrate.rs:12-18` explains the split:
  `calibration.rs` was a byte-for-byte copy of production, checked by a textual guard. That guard
  was deleted at promotion step C3 (`paralog/mod.rs:24-25`), so the split no longer protects
  anything.
- **Proposal.** Move `calibrate_from_the_ratio_histogram` into `calibration.rs`.
- **Downside.** It loses the visible line between "production's code" and "ng's code". The
  numerical parity test (`production_parity.rs`) still pins the numbers, which is what matters.
- **Effect.** None. **Cost.** S. **Confidence.** High.

---

## Readability

### R1. Comments describe a psp path that has since been built

- **Where:** each of these says a check "must become a `RunError` when observations are decoded
  from a psp file":
  - `build.rs:103-106`
  - `build.rs:152-155`
  - `build.rs:1987-1990`
  - `build.rs:2316-2318`
  - `observation_cache.rs:1985-1990`

  Also `build.rs:2487-2493` ("this module has no `RunError` to raise yet").
- **What.** Observations are decoded from psp files, and `psp_source` refuses each of these
  shapes as damage before the merge sees it:
  - overlapping or out-of-order records: `psp_source.rs:523-540`, "found by fuzzing,
    2026-10-31";
  - reference-width mismatch: `psp_source.rs:1101-1109`;
  - empty read lists with a non-zero read count: `psp/record.rs:930-934`.

  The comments now tell a reader that a crash path is open when it is not. They also name
  `ObservationExceedsReachCeiling`, which does not exist anywhere in the tree.
- **Proposal.** Rewrite each in one sentence: the producer guarantees this; the psp reader refuses
  it as damage (`psp_source.rs:…`); this assertion catches a defect in the generator.
- **Effect.** None. **Cost.** S. **Confidence.** High.

### R2. References to deleted code and finished plans

- **Where**
  - `parallel` / `merge_cohort_in_parallel` / `parallel.rs`: `mod.rs:24`, `serial.rs:18`,
    `organise.rs:17`, `organise.rs:155`, `rounds.rs:10`, `observation_cache.rs:2002` (cites a
    test in that file).
  - `mod.rs:14-29`: the "what has landed / still to come" list, which says the caller objects do
    not exist.
  - `organise.rs:1-19`: describes the organiser as the component that releases loci.
  - `paralog/mod.rs:15-16`: says the run's wiring "has not landed".
  - `paralog/mod.rs:35-43`: the "Landed so far, Milestone A" list.
  - `doc/devel/ng/arch/cohort_merge.md` §4, and the layout at `:30`: still give `organise.rs`
    the production role that `stream.rs` and `Frontier` now have.
- **Proposal.** Rewrite the module headers to say what each file does today. Move the history to
  the reports, where it already exists.
- **Downside.** Some history explains a decision. Keep the reason; drop the dates and plan codes.
- **Effect.** None. **Cost.** S. **Confidence.** High.

### R3. Duplicated and interleaved comment blocks in the closing walk

- **Where**
  - `close.rs:653-664`: the `summaries_per_sample` doc has two openings, the second restating
    the first.
  - `close.rs:911-965`: four blocks overlap. One explains the ordering check, one the summary
    reading. Two others (`:940-955` and `:956-965`) both explain the never-mixed-kinds check.
    The ordering check's block sits above the summary block, not above the assertion it
    explains.
- **Proposal.** Keep one paragraph per check, each directly above its assertion.
- **Effect.** None. **Cost.** S. **Confidence.** High.

### R4. Comments outnumber code in the three largest merge files

- **Where.** Comment lines out of non-test lines:
  - `build.rs`: 1,374 of 2,548
  - `close.rs`: 547 of 1,083
  - `observation_cache.rs`: 1,077 of 2,103
- **What.** Much of the text is history: dates, review rounds, plan-step codes (`B1`, `B3`,
  `C1`, `D2`, `E3`, `A4`), and earlier drafts ("an earlier draft of this comment said…",
  `build.rs:347`). The clear-writing skill's Rule 5 asks for none of these in reader-facing
  text. The reasons a reader needs are present, but under the history.
- **Proposal.** Do not run a blanket trimming pass. That would risk deleting load-bearing
  reasons, and the owner values them. Instead, whenever a fix touches a function, rewrite its
  comment to the current reason, and drop plan codes and "until <date>" narratives. Start with
  R1 and R2.
- **Downside.** Slow, and depends on discipline.
- **Effect.** None. **Cost.** L if done wholesale; S each when done in passing.
- **Confidence.** Medium (taste).

### R5. Test-only methods inside production impls

- **Where**
  - `build.rs:491, 497`: `AlleleTable` fields marked `expect(dead_code)` outside tests.
  - `build.rs:604-629`, `864-911`: `#[cfg(test)]` methods.
  - `build.rs:914-1052`: `RegionOutcome`, `build_region`, `build_region_windowed`.
  - `close.rs:84-89`, `393-404`, `717-751`: test-only `span` and `over`.
  - `serial.rs:36-180`: two test-only drivers.
- **What.** A reader of the production path has to skip a `#[cfg(test)]` item every few hundred
  lines to find the code that runs. `AlleleTable` keeps a whole `LocusReferenceBases` and a
  counter alive in production that only tests read.
- **Proposal.** Group the test-only drivers and accessors in a `#[cfg(test)] mod oracle` per file,
  or in a `test_support.rs`. Where `AlleleTable` holds fields only tests read, return them from
  `assemble` to the test entry point instead.
- **Downside.** The oracle drivers are a deliberate part of the design: the serial driver is the
  oracle the module claims byte identity against (spec §15). Keeping each next to its production
  twin has some value.
- **Effect.** None. Removing the unused `reference` from `AlleleTable` saves one move per built
  locus. **Cost.** M. **Confidence.** Medium.

---

## Looked at, left alone

- `stream.rs` shutdown. Ending the dispatcher drops `ask_to_read`, which ends the reader. Every
  region task sends to an unbounded channel, so none blocks during unwinding. A failed or panicked
  `take` drops `arrived`, which turns the next `send` into `stopping`. No deadlock was found. The
  `catch_unwind` paths do nothing under `panic = "abort"`, but they cost nothing.
- `stream.rs` ordering. Regions are handed over in round order, and within a round in region
  order, through one `Frontier`, so the output does not depend on which thread finishes first.
- `cover` against `cover_in_parallel`. The serial sweep and the snapshot-and-max loop reach the
  same least fixpoint, because each draws only what starts at or before the current reach.
- Window truncation at `left_edge` cannot split a locus that an earlier region owns. Any record
  that a truncated record chained to has a reach at or past the region start, so it stays in the
  window (given B3's precondition).
- `seal_round` hands the round its finalised windows and copies back those at or past
  `next_from`. Loci that begin in the round have their windows in the view.
- Pass two and pass three shutdown. An early error drops the receiver, which ends the reader and
  the preparer. Ratios and records pair by index because the reader refuses a count mismatch.
- The window-based prior (`paralog/windowed_prior.rs`) uses `par_iter` per bin but sums
  serially, so it does not depend on the thread count.
- `ParalogScorePrecompute::new(…, samples_with_model)` builds its tables at the full cohort size
  and sets the frequency floor from samples that have a model. That is consistent with the
  scorer's cohort-size check.
- `ParalogScoringContext::observation_of` rejects infinite windows as well as `NaN` ones. This is
  deliberate and correct.
- `rewrite_filter_and_info` refuses a line with fewer than 9 columns. A sites-only line cannot
  occur, because a zero-sample run is refused earlier.
- `resume.rs` saves every field of `VcfHeaderMetadata` except the hidden-duplication provenance,
  which pass two recomputes.
- `refuse_malformed_analysed_regions` cannot fire on user input, because `regions.rs` sorts and
  merges `--regions` first.
- The two reference panics in `LocusReferenceBases::over` (an uncovered position, or members
  disagreeing on the reference) are covered at the run's start by the per-contig checksum checks
  (`psp_caller.rs:509, 532`).
- `recorded_windows::record_the_windows_at` panics on a write error. It runs only when
  `NG_WINDOW_COVERAGE_FILE` is set, as a probe.
- `src/paralog/` imports only `float`, `genetics` and `window_coverage`, which is the right
  direction. `production_parity.rs` (1,357 lines) is `#[cfg(test)]`.
- The scan of the `AlleleGroupTally` list in `assemble_capturing` (`build.rs:720-733`) is
  quadratic only in the distinct pairs per sample. The comment at `:715-719` already prices the
  worst case.
