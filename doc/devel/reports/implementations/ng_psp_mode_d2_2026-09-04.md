# ng psp mode — D2: one calling loop, and direct mode's VCF unchanged

**Date:** 2026-09-04
**Plan step:** [run_driver_psp_mode.md](../../ng/impl_plan/run_driver_psp_mode.md) Milestone D, step D2 (its own commit, deliberately not bundled)
**Spec:** [run_streaming.md](../../ng/spec/run_streaming.md) §3.1; arch [run_streaming.md](../../ng/arch/run_streaming.md) §3.1, §3.4
**Branch:** `ng-psp-mode`

## Plan

Spec §3.1 says both callers drive one merge and one calling loop. Until now that loop was a
method on `AlignedFilesVariantCaller`, so psp mode's caller would have had to copy it. D2 lifts
it out, generic over the source, and changes nothing else.

## Assumptions

- **The lift stops where the modes genuinely differ.** What the free function does not do is
  turn spent sources into per-sample tallies: a walker knows which regions it handled and what
  its generators counted, and a psp source knows none of that. So the loop hands the sources
  back in the run's sample order and each caller says what its own kind of source can.
- **The cache is taken by value**, so the sources can be handed back without the caller holding
  a borrow across the call.
- **The cohort's size is read off the cache** (`ObservationCache::sample_count`, new) rather
  than passed beside it. It was `self.samples.len()`; the two are the same number today, and
  the failure if they ever came apart is a genotyper told one cohort size while a different
  number of samples is drawn.

## Changes made

- **`call_cohort_from_sources_handing_each_record_over`** (`callers.rs`): the whole body from
  `parameters.view()` on, as a free function over `ObservationCache<Source>` +
  `CohortCallingInputs` + genotyper + sink.
- **`CohortCallingInputs<'a>`**: the six things a cohort is called with that are not its
  sources — the segmentation, the merge parameters, the model's numbers, the calling-loop and
  candidate-selection settings, and the reference accessor the padding base is read from.
  `AlignmentInputs` is the shape it follows.
- **`CohortCallingTallies`**: the five facts a finished run states that do not depend on where
  the observations came from. **`WrittenCohort` now holds one of these plus the walk tallies**,
  and so does the loop's own `CohortCallingOutcome<Source>`, so the five have one definition
  rather than two.
- `AlignedFilesVariantCaller::call_cohort_handing_each_record_over` is now: open the walkers,
  call the loop, turn the spent walkers into per-sample tallies. Its doc says what direct mode
  adds and points at the loop for the rest, rather than carrying a second copy of it.

## Tests added

One, and the lift is what made it writable: **a refused record is the run's answer even when
the source then fails.** Both failures are live at the same return — the sink refused a record,
and the source gave out afterwards — and the rule is that the refusal wins, because it happened
first. Staging that needs a source the run cannot open, so the test drives the lifted loop
directly over a `Vec` of the walker's own observations with a failure appended behind them. It
is also the crate's only instantiation of the loop over something that is not a walker, which
is the source-agnostic claim exercised rather than asserted.

`cargo test --lib 'ng::run'` goes from 474 to 475.

## Validation results

- **The oracle the plan asks for: direct mode's VCF is byte-identical across the lift.** Six
  tomato accessions (`SRR7279481`, `488`, `501`, `533`, `536`, `537`) over the first two 100 kb
  intervals of `benchmarks/tomato1/regions.bed`, with the shipped defaults at four threads:
  **598 records, sha256 `5f0903cf1319b0157dda8dfc0a81b8b70a22d452193bd423915ee018a2200bf8` on
  both sides.** The parameters file matches too (`21411ade…`), and so does the whole run report
  — which matters because the report is rendered from `WrittenCohort`, the type this step
  reshaped. `tmp/d2_oracle.sh before` / `… after`, run either side of the change.
- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` exit 0.
- `cargo test --lib` — **6,124 passed**, 0 failed, 14 ignored. `--lib 'ng::run'` 475 passed.
- `cargo test --tests` — all 21 integration binaries pass.
- **Mutation:** moving `merged?` ahead of the refusal check — the exact error-precedence swap
  that survived the review's own five mutations — now fails the new test and nothing else.

Standing red, untouched and predating the branch: the three locus-dump example tests, and
`cargo test --all-targets` aborting on the `ng_joint_fit_perf` bench, which needs
`--features bench-fixtures`.

## Tradeoffs and follow-ups

- **`WrittenCohort`'s five fields moved down one level** (`written.calling.records_written`),
  which touched `report.rs`, its tests, `callers.rs`'s own tests and one example. The
  alternative was two copies of five field definitions, and the review found two facts already
  recorded in only one of them on the day they were written.
- **The three new items are `pub(crate)`.** Milestone E's `PspVariantCaller` lands in this same
  file, so nothing needs them wider.
- **Nothing psp-shaped drives the loop yet.** That is E3.
