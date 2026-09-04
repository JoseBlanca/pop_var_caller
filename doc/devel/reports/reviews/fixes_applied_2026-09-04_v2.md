# Fixes applied — psp mode D2, the lifted calling loop

**Date:** 2026-09-04
**Review:** [ng_psp_mode_d2_2026-09-04.md](ng_psp_mode_d2_2026-09-04.md)
**Branch:** `ng-psp-mode`
**Outcome:** every finding applied; nothing deferred, nothing routed to the owner.

## The Major: the missing test, written where the lift made it possible

`a_refused_record_is_the_runs_answer_even_when_the_source_then_fails` stages both failures at
one return — the sink refuses the record at chr1:15, and the source gives out behind it — and
asserts the refusal is what comes back, because it happened first. It drives the lifted loop
directly over a `Vec` of the walker's own observations with a failure appended, since a walker
over real alignment files cannot be made to fail on cue. Building regions are narrowed to ten
bases so that the locus is built and refused in an earlier batch than the one that draws the
failure; without that, the cover returns before the sink is ever called and only one of the two
failures is live.

Re-running the review's surviving mutation — `merged?` moved ahead of the refusal check — now
fails that test and nothing else. It also closes the review's other Minor: this is the crate's
only instantiation of the loop over a source that is not a walker, so the source-agnostic claim
is exercised rather than asserted.

## The Minors

- **`sample_count` moved below `into_sources`.** It had landed between that method's doc
  comment and its signature, so rustdoc gave each the other's text.
- **The concurrency claim now names the test that exists.**
  `the_record_path_is_byte_identical_at_every_thread_count` compares this method's VCF bytes at
  pools of 1, 2, 4, 8 and 16 and has since 2026-09-01; the cover's fixpoint is pinned a layer
  down, at the merge, on fixtures that can hold a locus wider than one position. The doc had it
  the other way round and said the end-to-end fixture was unbuilt.
- **`CohortCallingOutcome::sources` says what is true of it**: the sources are spent, because
  the value is built only after both error returns.
- **The five shared tallies have one definition.** `CohortCallingTallies` holds them;
  `WrittenCohort` and `CohortCallingOutcome` each hold one. The alternative the review
  considered and rejected — `WrittenCohort` wrapping the generic outcome — would have made a
  run-report type generic over the observation source, which is the leak this step exists to
  prevent.
- **Renamed to `call_cohort_from_sources_handing_each_record_over`**, which carries the file's
  own distinction: `call_cohort` accumulates, this hands each record over and keeps none.

## The nits

`CohortCallingOutcome` derives `Debug`; the errors line names `TractGeneratorSettings` instead
of "its neighbours"; the three new items are `pub(crate)`.

## Validation

- **The byte-identity oracle re-run after every fix, including the `WrittenCohort` reshape**:
  direct mode's VCF over six tomato accessions and 200 kb is
  `5f0903cf1319b0157dda8dfc0a81b8b70a22d452193bd423915ee018a2200bf8` on both sides, 598
  records; the parameters file and the whole run report match too. The run report is the one
  that mattered here — it is rendered from `WrittenCohort`, the type the reshape touched.
- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` exit 0.
- `cargo test --lib` — 6,124 passed, 0 failed. `--lib 'ng::run'` 475 passed (474 before).
- `cargo test --tests` — 21 integration binaries, all ok.
