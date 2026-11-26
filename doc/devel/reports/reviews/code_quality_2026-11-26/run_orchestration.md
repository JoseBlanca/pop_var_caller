# Code-quality map: `src/run/`, the pipeline driver (2026-11-26)

Scope: `src/run/` except `cohort_merge/` and `paralog_filter/`. The main files are `callers.rs`,
`psp_caller.rs`, `psp_source.rs`, the three `census_*.rs` files, `gatherer.rs`, `walker.rs`,
`records.rs`, `report.rs`, `explain.rs`, `edge_deletion_gathering.rs`, `fragment_keeper.rs`,
`calling_progress.rs` and `mod.rs`. Nothing was built or run; every finding is reasoned from the
code at `f84129cf`.

## Summary

1. **Two copies of "call one locus and build its record".** The calling loop has two drivers: one
   that works one building region at a time, and one that works a round of regions on the thread
   pool. Each carries its own 150-line copy of the per-locus work: re-read the reference, call,
   pad, correct the quality, apply `--min-site-quality`, gather read alleles, assemble the record.
   `call-from-alignments` always uses the first driver. `call-from-psps` uses the second whenever
   a round holds more than one region, which is the usual case. So the two commands write the same
   VCF only because two hand-kept copies agree. A third copy sits in the test-only oracle
   `call_cohort`. Merge them into one function. Output does not change.
2. **`call-from-alignments` computes a round shape it never uses.** `merge_parameters_for` divides
   the run's ground into as many regions as there are threads, but the direct-mode driver ignores
   the region count and works the narrower regions one at a time. So the region width depends on
   `--threads` for no gain, and the doc of `round_shape_for` claims a speed-up for
   `call-from-alignments` that its driver cannot deliver. This needs an owner decision and a
   measurement; it may change speed, not output.
3. **The record hand-over chain is written twice.** After the driver, each record goes through the
   edge-deletion gathering, then the fragment keeper, then the output. Direct mode and psp mode
   each build this chain, about 45 and 95 lines, with different batching.
4. **Five run settings are stored three times and set twice.** `explain`, `min_site_quality`,
   `low_base_quality_weight`, `gathers_fragments` and `depth_ceiling` are fields and builder methods
   on both callers and fields on `CohortCallingInputs`. Both subcommands set them with the same
   chain of builder calls.
5. **`run` imports `crate::cli` in only one place that matters**: `calling_progress.rs` takes three
   text-formatting helpers from `cli::run_output`. Every other import is in a test module or a doc
   link. Move the helpers down; `parameter_estimation/progress.rs` has the same dependency.
6. **The three `census_*` modules are `estimate-parameters`' pipeline, not the run's.** Only
   `cli/estimate_parameters.rs` uses them. One function, `fitted_inbreeding_of`, is called from
   `parameter_estimation`'s tests, which makes a test-only dependency cycle.
7. **No bug that changes a call was found.** Two small ones were: a refused record can be reported
   at the wrong locus (a wrong error message), and `--max-locus-span` above 65,535 can panic on a
   very long generic locus (suspected).
8. **The documentation has drifted.** Module and function docs in `callers.rs`, `mod.rs`,
   `psp_source.rs` and `walker.rs` describe a state several months old. One doc cites a test that
   has never existed.

---

## Bugs

### B1. A record refused by the output is reported at the wrong locus

- **Where** — `src/run/callers.rs:872-880` (direct mode), `src/run/psp_caller.rs:843-855` (psp
  mode). The error is wrapped at `src/run/callers.rs:1256-1265` and `:7744-7751`.
- **What** — The driver hands each finished record to a closure that runs the gathering and the
  fragment keeper and then calls the output's `hand_over`. Both stages hold records back: the
  gathering holds them for up to 256 bases, and the keeper until reads move more than
  `READ_REACH` past them. So a `hand_over` refusal returned from inside the closure belongs to an
  *earlier* record, but the driver wraps it as `RunError::RecordNotWritten { locus: region }`,
  where `region` is the locus just called. psp mode says so in a comment
  (`psp_caller.rs:832-834`); direct mode does not. The final flush after the driver
  (`callers.rs:890-897`, `psp_caller.rs:893-902`) names the right locus.
- **Scenario** — The disk fills while the keeper releases a record at chr1:1,000. The record just
  called is at chr1:1,180. The error says chr1:1,180.
- **Proposal** — Have the hand-over closure wrap its own error with the released record's region
  and return a `RunError`, rather than returning the bare `E` for the driver to wrap. This falls out
  of D2 below, where the chain becomes one type that knows which record it was writing.
- **Downside** — The driver's closure type changes from `Result<(), E>` to `Result<(), RunError>`
  for this path. The round driver's `prepare` keeps `E`.
- **Effect** — Output: none. Error message: corrected.
- **Cost** — S, or free inside D2.
- **Confidence** — high (likely: traced end to end).

### B2. `--max-locus-span` above 65,535 can panic on a long generic locus

- **Where** — `src/run/records.rs:358-361`; the flag is `u32` with only a non-zero check
  (`src/cli/call_from_alignments.rs:243-244`, `src/cli/calling_run.rs:1481-1490`).
- **What** — `reread_spellings_cut_short` converts the reference allele's length to `u16` with
  `expect("an ordinary locus is bounded by --max-locus-span, far inside u16")`. The default bound
  is 100, so this never fires at the default. Nothing stops a user from typing a bound above
  65,535, and a generic locus with three or more alleles longer than that then panics.
- **Scenario** — `--max-locus-span 200000` on long-read data with a 70 kb deletion that the
  merge joins into one generic locus with two other spellings. The run panics mid-genome.
- **Proposal** — Refuse `--max-locus-span` above `u16::MAX` where it is parsed, with a message.
  Check first whether other code also assumes a generic locus fits in `u16`, in which case the
  refusal fixes them too.
- **Downside** — A command line that used to be accepted is now refused. That is a CLI change and
  needs the owner's approval as a bug fix.
- **Effect** — Output: none on any run that did not panic.
- **Cost** — S.
- **Confidence** — medium (suspected: I did not confirm that the merge can build such a locus).

### B3. Nothing pins the round driver to the streaming driver, though a comment says something does

- **Where** — `src/run/psp_caller.rs:805-811`.
- **What** — The comment says the two drivers write the same VCF, "which is what
  `the_round_driver_writes_what_the_streaming_driver_writes` holds them to". No such test exists,
  and `git log -S` finds none in the history. What actually holds them together is
  `cli/mode_equivalence.rs`, which compares direct mode (always the streaming driver) with psp mode
  (the round driver whenever the test pool has enough threads). On a one-core CI machine psp mode
  takes the streaming driver too, and nothing compares the two drivers at all. The direct-mode
  thread-count test, `the_record_path_is_byte_identical_at_every_thread_count`, uses
  `MergeParameters::DEFAULT`, whose region count is one, so it never reaches the round driver.
- **Proposal** — Write the test the comment names: one psp cohort called at a region count of one
  and of several, VCF bytes compared. If A1 below is done first, the test mostly pins the
  scheduling, which is still worth pinning.
- **Downside** — None beyond test run time.
- **Effect** — none.
- **Cost** — S.
- **Confidence** — high.

---

## Architecture

### A1. `call-from-alignments` computes a region count its driver ignores

- **Where** — `src/cli/call_from_alignments.rs:696-701` asks
  `calling_run::merge_parameters_for(..., ObservationsReadFrom::Alignments)`.
  `src/cli/calling_run.rs:485-516` (`round_shape_for`) then sets
  `in_flight = min(ground / 500, threads)` and `width = ground / in_flight`. Direct mode always
  calls the streaming driver (`src/run/callers.rs:867`), which takes the width
  (`callers.rs:1105-1110`) and never reads `cohort_locus_builder_regions_in_flight`. Only psp mode
  branches on it (`src/run/psp_caller.rs:856-876`).
- **What** — Three consequences.
  1. The region width in direct mode depends on `--threads`. For example, with 63 samples and the
     default ground of 7,936 bases, one thread gives one 7,936-base region and 48 threads give
     fifteen 529-base regions, worked one after another. The VCF should not change, because the
     merge's answer does not depend on where regions end
     (`the_answer_does_not_move_when_the_ground_is_divided_ten_ways`). The number of covers, and
     so the number of times the cover waits for every sample's thread, grows with the thread count.
  2. The doc at `calling_run.rs:468-475` says `call-from-alignments` gets "the calling pass up to
     5.5× faster than one region at a time". That figure comes from commit `2137b8b3`, which added
     the round driver for psp mode. `git log -S` finds no version of `callers.rs` that called the
     round driver.
  3. The 76 s against 65 s measurement in the same doc, 125-base against 500-base regions at 48
     accessions, was therefore made on the streaming driver. It measures more covers, not a
     parallel gain.
- **Proposal** — Owner decision between two options.
  - **(a)** Return a region count of one and the full ground for `Alignments`. This makes the width
    independent of `--threads` and matches what the driver does. Fix the doc.
  - **(b)** Route direct mode through the round driver as psp mode does. That needs the walker to
    be `Sync`, or the stream to hold the walkers differently. It is the bigger change and the one
    the doc already claims.
- **Downside** — (a) may be slower or faster than today, depending on how cover cost scales with
  width, so it needs a timing on the tomato benchmark first. (b) is a real piece of work in
  `cohort_merge`, which is another agent's area.
- **Effect** — Output: none, as far as the merge's width invariance holds. Speed: **possible
  change either way**. Memory: (a) holds one region of the full ground width, which is what the
  round already holds today, so it should not change.
- **Cost** — (a) S plus a measurement; (b) L.
- **Confidence** — high that the region count is ignored; medium on the speed effect.

### A2. `run` imports `crate::cli` for three formatting helpers

- **Where** — `src/run/calling_progress.rs:23` imports
  `cli::run_output::{PROGRESS, elapsed_text, with_thousands}`. The same wrong-way import is at
  `src/parameter_estimation/progress.rs:25` (outside this area). `with_thousands` has a private
  duplicate at `src/neighbour_evidence/peak.rs:408`.
- **What** — Of the seven files the brief lists, this is the only production import. The rest are:
  - doc links only: `mod.rs:27`, `callers.rs:180`, `census_freshness.rs:47`,
    `calling_progress.rs:3`;
  - inside `#[cfg(test)]` modules: `census_cohort.rs:280-285`, `census_fit.rs:972-986`,
    `gatherer.rs:1614`, `test_fixtures.rs:137`. These build fixtures by running the real
    `generate-psps` command (`cli::generate_psps::run_generate_psps`) and the CLI's cohort
    fixtures.

  `PROGRESS` is a tracing target string (`"pvc::progress"`, `cli/run_output.rs:44`).
  `with_thousands` (`:480`) and `elapsed_text` (`:495`) are pure formatting functions. None of
  them needs anything from `cli`.
- **Proposal** — Move the three items to a small library module (for example `src/progress.rs`, or
  next to `error_render`) and re-export them from `cli::run_output` so the CLI's own uses do not
  change. Point `neighbour_evidence/peak.rs` at the shared `with_thousands`. Leave the test-module
  imports alone: building a fixture through the real command is a stronger test than a hand-built
  file, and the cycle exists only in test builds.
- **Downside** — One more small module. `run_output`'s subscriber reads `PROGRESS` to draw lines in
  place, so the const has to live where both can see it, which the new module allows.
- **Effect** — none.
- **Cost** — S.
- **Confidence** — high.

### A3. The `census_*` modules belong to `estimate-parameters`, not to the calling run

- **Where** — `src/run/census_fit.rs` (about 870 production lines), `census_cohort.rs` (about 265),
  `census_freshness.rs` (about 570). Their only production caller is
  `src/cli/estimate_parameters.rs:72-77`. `fitted_inbreeding_of` (`census_fit.rs:620`) is called
  from `parameter_estimation/joint/fit.rs:5881` and `:8567`, both inside test modules.
- **What** — These modules fit a cohort of censuses and turn the fit into a parameters file. They
  touch `run` only through `OpenPspCohort` (`census_cohort.rs:54`, `census_freshness.rs:29`) and
  one string constant (`census_fit.rs:68`). Living in `run` makes `run` the home of both pipelines
  and gives `parameter_estimation` a test-only dependency on `run`.
- **Proposal** — Two steps, the first cheap:
  1. Move `fitted_inbreeding_of` into `parameter_estimation` (next to `HomozygoteExcess`), and
     re-export it from `run` so the CLI does not change. This removes the cycle.
  2. Optionally, move the three modules into a sibling of `run`, such as `src/cohort_fit/`, which
     depends on `run::OpenPspCohort`.
- **Downside** — Step 2 is churn across about 1,700 production lines and their tests for a gain in
  findability only. Moving `census_fit` *into* `parameter_estimation` instead would add a new
  dependency from `parameter_estimation` on `calling::parameters_file`, which it does not have
  today; that is worse. I would do step 1 and leave step 2 unless the owner wants `run` to mean
  "calling run" strictly.
- **Effect** — none.
- **Cost** — step 1 S; step 2 M.
- **Confidence** — high on the facts; medium on whether step 2 is worth it.

### A4. Test-only accessors and two fields kept only for tests on `AlignedFilesVariantCaller`

- **Where** — `src/run/callers.rs:224-235` (`read_groups` and `read_filters`, both
  `expect(dead_code, reason = "only tests read it")`). Eleven `#[cfg(test)]` accessors are at
  `callers.rs:430-521`.
- **What** — The doc on `read_groups` (`:227-231`) gives a production reason for keeping the field,
  "writing the parameters file beside the output and reporting what the run used both need this
  table back". The CLI does both from its own copy (`call_from_alignments.rs:808-812`). So the
  field clones a cohort-sized table into every run only for tests. Most of the accessors exist for
  `every_setting_comes_back_as_it_was_given` (`callers.rs:2768`) and its neighbours, which check
  that a setting given to `open` can be read back. They do not check that the run uses it; the
  tests in `calling_joined_to_the_merge` (`:5109-5309`) already check that.
- **Proposal** — Drop the two fields and the accessors that only round-trip a value. Keep
  `sample_count`, `sample_reads` and `assembly_check` if tests use them for behaviour. Delete or
  rewrite the round-trip tests.
- **Downside** — A few tests are lost. They catch only a field mis-assignment in `open`, which the
  behaviour tests would also catch.
- **Effect** — none; memory falls by one `ReadGroups` clone a run, which is negligible.
- **Cost** — S.
- **Confidence** — medium.

---

## Duplication

### D1. The per-locus record build is written out three times

- **Where**
  - The streaming driver, `call_cohort_from_sources_handing_each_record_over`: per-locus body at
    `src/run/callers.rs:1111-1269`, record-building closure at `:1154-1222`.
  - The round driver, `call_cohort_from_sources_in_rounds`: per-locus body at `:7550-7706`,
    closure at `:7606-7667`.
  - The test-only oracle, `AlignedFilesVariantCaller::call_cohort`: `:698-755`, which builds no
    record.

  The two drivers' set-up is copied as well (`:1044-1077` against `:7503-7524`: destructuring,
  `explaining_drops`, the tract selection config), and so is their tail (`:1278-1305` against
  `:7786-7811`: dropped-loci rows, too-wide rows, histograms, the outcome struct).
- **What** — Both production copies do the same steps in the same order. Re-read cut-short
  spellings, read the reference beside the locus, call `call_one_cohort_locus`, and in its closure:
  skip loci with no alternative, read the base before the left flank and the padding base, build
  the evidence, apply the site-quality gate, gather read alleles, fill the window-coverage buffer,
  assemble the record. Then both complete the explanation rows. I compared the closures line by
  line and they are identical apart from how they borrow their scratch. **This is where the two
  commands could start writing different VCFs:** direct mode always uses the first driver and psp
  mode usually uses the second (B3, A1).

  The comment above the round driver calls it a prototype "reached only when `NG_ROUND_CALL` is set
  in the environment" (`:7424-7429`). That is no longer true: nothing reads `NG_ROUND_CALL`, and
  the round driver is psp mode's production path.
- **Proposal** — One per-worker scratch struct and one function:

  ```rust
  struct LocusWorker<S> {            // today's RoundCallScratch minus `totals`
      shaping, tract_shaping, scratch, padding_scratch, beside_scratch, window_coverage, tracts,
  }
  struct LocusCallingContext<'a, G> { // the read-only values both drivers destructure
      genotyper, explain, frozen, candidate_selection, tract_selection, calling_loop_config,
      run_sample_count, padding_reference, min_site_quality, low_base_quality_weight,
      gathers_fragments,
  }
  fn call_and_assemble<S, G>(
      context: &LocusCallingContext<'_, G>,
      worker: &mut LocusWorker<S>,
      observation: CohortObservation,
      reads: &LocusReads<'_, '_>,
  ) -> (LocusOutcome<Result<CalledLocus, RunError>>, Option<Vec<ExplainRow>>)
  ```

  The streaming driver owns one `LocusWorker` and calls this from its sink. The round driver keeps
  `RoundCallScratch { worker: LocusWorker<S>, totals }` for its `Drop` flush and calls the same
  function. Pull the set-up and the tail into two small helpers too. `call_cohort` can call
  `call_one_cohort_locus` through the same context, or be left as it is; it is a test oracle and
  builds no record.
- **Downside** — One more struct to read through. The streaming driver's closure today reads its
  surroundings directly, which some readers find easier than a context argument. The function sits
  on the hottest path of calling, so it must inline or cost no more than the closure. A call per
  locus through a struct reference costs nothing measurable beside a genotype fit, but confirm
  with the benchmark.
- **Effect** — Output: none, and an output-identity check would confirm it. Speed: none expected.
- **Cost** — M.
- **Confidence** — high.

### D2. The record hand-over chain is built twice, with different batching

- **Where** — Direct mode at `src/run/callers.rs:857-901`; psp mode at
  `src/run/psp_caller.rs:812-907`.
- **What** — Both modes run each finished record through `EdgeDeletionGathering`, then
  `FragmentKeeper::for_a_run`. At the end they drain the gathering into the keeper, finish the
  keeper, log the keeper's counts if it keeps fragments, hand everything released to the output,
  and call `after_gathering` on the tallies. psp mode adds batching (`KEPT_A_BATCH`,
  `PREPARED_A_BATCH`) and a side pool for preparing records. Direct mode does the same work one
  record at a time on the calling thread. The comments in each say "as in psp mode" and "as in
  direct mode": the sequence is meant to be one.
- **Proposal** — One `RecordRelease<P, E>` type in `run` that owns the gathering, the keeper and
  the two waiting lists. It has `accept(record, windows)` and `finish() -> EdgeDeletionGatheringCounts`,
  and takes a `prepare` and a `hand_over`. psp mode builds it with the side pool and its batch
  sizes. Direct mode builds it with a batch size of one and no pool, which is today's behaviour.
  This is also where B1's correct locus naming lives.
- **Downside** — Direct mode would carry the batching code as configuration it does not use. If
  direct mode should batch too (its calling thread does the same spill encoding with the duplication
  filter on), that is a speed change to measure separately, not part of this merge.
- **Effect** — Output: none. `FragmentKeeper::accept` is documented as `accept_batch` with a batch
  of one (`fragment_keeper.rs:207-216`). Speed: none if direct mode keeps a batch of one.
- **Cost** — M.
- **Confidence** — high.

### D3. Five run settings live in three structs and are set by two identical builder chains

- **Where**
  - Fields and defaults: `src/run/callers.rs:263-276` and `:371-377`; `src/run/psp_caller.rs:450-461`
    and `:554-560`.
  - Builder methods: `callers.rs:381-428`; `psp_caller.rs:564-606`.
  - Passed again as fields of `CohortCallingInputs`: `callers.rs:935-948`.
  - Set by identical chains in `src/cli/call_from_alignments.rs:835-844, 879-880` and
    `src/cli/call_from_psps.rs:728-737, 773-774`.
- **What** — `depth_ceiling`, `explain`, `min_site_quality`, `low_base_quality_weight` and
  `gathers_fragments`. Each has the same default and the same doc in both callers. The two
  `reference_path` computations (`callers.rs:337-342`, `psp_caller.rs:519-523`) are a smaller copy
  of the same kind.
- **Proposal** — A `CallingSettings` struct with a `Default`, held once by each caller and embedded
  in `CohortCallingInputs`. Build it once in `calling_run` from the shared arguments. Keep
  `with_fragments_gathered` as the one builder, because the CLI learns it only after building the
  header. Move the `reference_path` lines into a function beside
  `refuse_a_catalog_built_on_another_reference`.
- **Downside** — `depth_ceiling` is applied in different places in the two modes (by the walker in
  direct mode, by the psp source in psp mode), so the struct groups settings that are consumed at
  different depths. That is acceptable: it is what the run was configured with. It does touch both
  subcommands.
- **Effect** — none.
- **Cost** — S to M.
- **Confidence** — high.

### D4. The psp walk re-implements the walker's progress tracking and failure wrapping

- **Where** — `src/run/gatherer.rs:752-776` against `src/run/walker.rs:249-270`.
- **What** — Both wrap `SampleLocusObservationsIterator`, update `reached` from `reach_position()`
  and wrap a failure as `RunError::SourceFailed { sample, reached, source }`: about 15 lines, twice.
- **Proposal** — **Leave it.** The walker is deliberately not an `Iterator`
  (`walker.rs:40-45`) and the gatherer deliberately is one, and the census is fed at the gatherer's
  yield point. Sharing the 15 lines would need a third type or a trait-coherence workaround, which
  is more to read than the copy.
- **Effect** — n/a. **Cost** — n/a. **Confidence** — high.

### D5. `RunReport::of` and `RunReport::of_a_stored_cohort` repeat one constructor

- **Where** — `src/run/report.rs:120-170`.
- **What** — Two 20-line constructors that differ only in the `samples` variant and in whether
  `calling` is reached through `written`.
- **Proposal** — One private `fn new(calling, samples: WhatEachSampleDid, ...)` that both call.
- **Downside** — None worth naming.
- **Effect** — none. **Cost** — S. **Confidence** — high.

---

## Readability

### R1. Stale module and function documentation in the driver

Each of these describes code that has since changed. A reader who trusts them will reason wrongly
about which path runs.

| Where | What it says | What is true |
|---|---|---|
| `callers.rs:11-25` | the object returns every called locus at once; `pub` because no subcommand calls it yet | records are handed over one at a time; the subcommand landed; `merge_cohort`/`call_cohort` are test-only |
| `callers.rs:628-631` | every locus goes down the generic path; repeat-tract selection is unbuilt | STR and bundle loci are dispatched (`call_one_cohort_locus`, `:1418-1463`) |
| `callers.rs:531-534`, `:1997-2000` | both repeat-tract generator slots are refused as unbuilt | `generic_path_generators` fills all three (`walker.rs:1618-1619`) |
| `callers.rs:7424-7429` | the round driver is a prototype reached only through `NG_ROUND_CALL` | it is psp mode's production path; nothing reads the variable |
| `callers.rs:7431-7434` | tract counts leave on `Drop` because `map_init` gives no way back | the stream keeps scratch in a `Mutex<Vec<_>>` (`cohort_merge/stream.rs:206-270`); the `Drop` flush still works, but the reason given is out of date |
| `walker.rs:1622-1629` | bundles and tracts are set aside and counted by calling | both are called |
| `walker.rs:1643` | the function is named `generic_path_generators` | it builds all three paths' generators; rename it to `run_generators` |
| `mod.rs:9-28` | "Landed so far"; `psp_source`'s `PspObservationSource` is psp mode's source | `PspSummarySource` is the source; `PspObservationSource` is test-only (`psp_source.rs:621`) |
| `psp_source.rs:18-19` | "on the path a run takes (`PspObservationSource::over`)" | as above |
| `psp_caller.rs:811` | a test named `the_round_driver_writes_what_the_streaming_driver_writes` | no such test (B3) |

- **Proposal** — One documentation pass over these lines, plus the rename. Strip the history
  ("until 2026-09-02 …") where it no longer explains the present code.
- **Downside** — None; doc-only, plus a rename inside the crate.
- **Effect** — none. **Cost** — S. **Confidence** — high.

### R2. `CohortWalkTallies::of` pairs names with walkers by position and checks only the count

- **Where** — `src/run/callers.rs:2020-2059`.
- **What** — psp mode asserts that each name matches its source's own name
  (`psp_caller.rs:955-962`). Direct mode checks only that the two lists are the same length, though
  each walker carries its sample's name. A mis-pairing would put one sample's counts under another
  sample's name in the run report.
- **Proposal** — Make `AlignmentFilesWalker::sample_name` available outside tests and assert name
  equality, as psp mode does. Then `sample_names` need not be collected before `walkers()`
  consumes the run (`callers.rs:814`), because the names come back with the walkers.
- **Downside** — None.
- **Effect** — none. **Cost** — S. **Confidence** — high.

---

## Looked at, left alone

- **The merge keeps running after the output refuses a record** (`callers.rs:1086-1097`). It is
  already flagged in the code as a `cohort_merge` interface change (the sink cannot say *stop*).
  It wastes time but does not corrupt the output: the VCF writer goes through `<output>.tmp` and a
  rename (`vcf/writer.rs:163-177`).
- **The round driver's tract totals flushed through `Drop` into a mutex**
  (`callers.rs:7448-7459`). The scratches drop inside the stream's thread scope, before
  `tract_totals.into_inner()`, and addition does not depend on order, so the totals are the same
  at any thread count.
- **The round driver's `assert_eq!(displaced, 0)`** (`callers.rs:7779-7784`). A panic is right
  for a defect that would otherwise leave the VCF short with nothing saying so.
- **The psp source treats a record too large to hold as over the depth ceiling**
  (`psp_source.rs:546-553`), which direct mode never does. This is a deliberate memory guard, and
  it is counted and reported (`report.rs:659`).
- **`EdgeDeletionGathering`'s `HashMap`** (`edge_deletion_gathering.rs:449`) is used only for
  look-ups. The best pairing is chosen by a fold over a `Vec` in permutation order, so the result
  does not depend on hash order.
- **`EdgeDeletionGathering`'s order key has no contig** (`:346-351`). It is safe because any
  contig change flushes everything held first (`:265-278`).
- **The fragment keeper's thread split** (`fragment_keeper.rs:308-350`). Tasks are joined in
  sample order, so the result is the one-thread result.
- **`call_one_cohort_locus` and its three arms each take 11-13 arguments with
  `too_many_arguments` waivers** (`callers.rs:1391-1412` and below). D1's context struct would
  absorb most of them; on their own the waivers' reasons are fair.
- **psp mode's three per-sample lists, paired by position and asserted**
  (`psp_caller.rs:927-969`). They are verbose, but they guard a silent mis-pairing that would give
  a wrong copy number.
- **`calling_progress.rs` counting a failed locus in the streaming driver but not in the round
  driver.** This affects only a progress line on a run that is about to fail.
