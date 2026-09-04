# ng psp mode — D1: a stored sample behind the merge's source interface

**Date:** 2026-09-04
**Plan step:** [run_driver_psp_mode.md](../../ng/impl_plan/run_driver_psp_mode.md) Milestone D, step D1
**Spec:** [run_streaming.md](../../ng/spec/run_streaming.md) §3.1, §3.4, §8; arch [run_streaming.md](../../ng/arch/run_streaming.md) §2, §5, §8
**Branch:** `ng-psp-mode`

## Plan

The merge asks a source one question — *what did this sample see next?* — and direct mode
answers it from alignment files. D1 is the other answer: one sample's observations decoded from
its psp, behind the same trait, so nothing above it can tell the two apart (spec §3.1).

## Assumptions

- **The adapter is `ObservationSource`, not an `Iterator`.** Both were possible and only one
  compiles: the blanket implementation at `observation_cache.rs:98` already makes every
  iterator of observations a source, so a type that was both would implement the trait twice.
  Taking the blanket implementation instead was the real alternative, and it cannot work —
  it passes the reader's error through untouched, where arch §5 requires an error that names
  the sample.
- **The sample's name comes from the file's own header** on the path a run takes
  (`PspObservationSource::over`), so the name a failure carries and the file it came from
  cannot come apart. The general constructor `new(sample, walk)` takes the name as an argument
  and is `pub(crate)` for exactly that reason.
- **It borrows the open reader rather than owning it.** `PspReader::records` hands back a walk
  that borrows the reader, so a source owning both would be self-referential. A calling run
  therefore holds its open psps and mints its sources at the call — which is the shape arch
  §3.4 already sketches for `PspVariantCaller` (`open psps, segmentation, params, pool`).
- **The spare record is dropped**, as the walker drops it and as the trait permits. A decoder
  is the reuse hook's best customer, but the store's walk hands records out rather than filling
  ones it is given, so taking the offer is a change to `PspReader` — deferred with the rest of
  psp-mode performance in the plan's Out.

## Changes made

- **`src/ng/run/psp_source.rs` (new).** `PspObservationSource<W>` — the sample's name, how far
  it has got, where the last record started, and the walk — implementing `ObservationSource`
  with `Error = RunError`. `PspObservationSource::reading(&mut PspReader)` is what a run
  builds; `over(sample, walk)` is the general constructor a fixture uses.
- **`PspSourceError`, three variants**, travelling as the cause under
  `RunError::SourceFailed`: observations out of coordinate order, an observation whose body was
  never built, and a draw made after this source already refused a record. **A refusal ends the
  source** — without that, the next draw hands back the record *after* the refused one and the
  stream goes on succeeding, one observation short with nothing to say so. Answering `None`
  instead would be the same silence in a different shape: the merge reads it as a sample that
  ran out.
- **`src/ng/run/mod.rs`:** the module, its two re-exports, and the module doc's landed list.

### One deviation from the plan's wording, absorbed and recorded

The plan asks for "an iterator adapter over `PspReader::records()`". What landed is generic
over the walk — any iterator of the store's `StreamedRecord` — for one reason that is a
correctness argument rather than a convenience: a psp can also be walked *selectively*
(`records_where`), which hands back records as a head with no body. **A merge fed such a walk
would build cohort loci with one sample's evidence silently missing** — wrong genotypes, no
error. The generic form is what makes that refusable, and `ObservationBodyNotBuilt` is the
refusal. It costs nothing at the call site: `reading` still fixes the walk to `RecordIter`.

## Tests added

Fifteen, all in the new module; `cargo test --lib 'ng::run'` goes from 459 to 474. (Eleven
landed first; the review added four — see the fix report.)

- What a psp holds comes back as the observations that were stored, field for field, over a
  file whose 40 records cross **8 blocks and 2 contigs** — a source that read one block and
  stopped fails here rather than at a cohort.
- A repeat tract arrives with its motif and both flanks (different lengths, so a source that
  swapped them would not pass).
- The sample a failure reports under is the one the header names.
- A psp holding no records is a source that is spent from the start — spec §8's
  analysed-and-empty, read as a sample that saw nothing rather than refused.
- A walk that skipped a body is refused rather than passed on without its evidence.
- Stored observations that go backwards are refused naming the sample, with `reached` the
  record *before* the refusal.
- Records starting on the same base are **not** backwards — the overlap a deletion produces,
  which a check against the previous record's end would refuse.
- A failing walk reports the sample and the last observation it handed over.
- A damaged block ends the source with the ground it had already covered (a real file, its
  third block's compressed frame flipped from its midpoint on).
- A walk that will not start fails before anything is reached.
- A source over a psp can cross a thread — the compile-time proof the merge's parallel cover
  needs, the psp half of `a_run_walker_can_cross_a_thread`.
- A source that refused a record refuses every later draw, rather than skipping it.
- `reached` is the last base the last observation covered, not the base it began on — the only
  fixture here with a multi-base record, and the only one that can tell the two apart.
- A spent source answers `None` for ever.
- A body declined part-way through reports the records already handed over.

### Mutation pass

Two passes. Before the review, five mutations and five killed (`tmp/d1_mutations/run.sh`).
After the review's fixes, eight — the first five plus the two the review found surviving and
one on the new latch — and **seven killed** (`tmp/d1_mutations/run2.sh`; the file's checksum was
compared before and after each pass to prove the restore landed):

| mutation | killed by |
|---|---|
| the order check removed | `stored_observations_that_go_backwards_are_refused_naming_the_sample`, and the latch test |
| `reached` takes the observation's first base, not its last | `reached_is_the_last_base_the_last_observation_covered` |
| the refusal latch removed | `a_source_that_refused_a_record_refuses_every_later_draw` |
| `reached` advanced before the order check | `stored_observations_that_go_backwards_…` (its `reached` assertion) |
| order compared against the previous record's reach, not its start | `records_that_start_on_the_same_base_are_not_backwards` |
| the out-of-order error names the previous record twice | `stored_observations_that_go_backwards_…` |
| the head-only refusal names a fabricated locus | 2 tests |
| **survived:** the constructor names a sample the header does not | nothing — no fixture reaches that arm through a `PspReader`, and the code says so where it sits |

The surviving one is `over`'s error arm: `PspReader::records` fails only on a seek inside
offsets `open` has already bounded, on a manifest `open` has already parsed, or on a buffer
ceiling refused where it is set. It is marked uncovered in the code, with the mutation named,
rather than left looking tested.

## Validation results

In the container, on this tree:

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` exit 0.
- `cargo test --lib 'ng::run'` — **474 passed**, 0 failed (459 at Checkpoint C).
- `cargo test --lib 'ng::run::psp_source'` — 15 passed.

Standing red, unchanged by this step and predating the branch: three locus-dump behaviour
tests (`ng_generic_loci_dump` ×2, `ng_ssr_loci_dump` ×1) and 11 unresolved intra-doc links
under `cargo doc --no-deps`.

## Tradeoffs and follow-ups

- **The out-of-order refusal is this source's, not the merge's.** The cache's release assertion
  (`observation_cache.rs`, `draw_next`) still stands and still aborts for any source that gets
  past it. Arch §8's item — *a source whose observations go backwards should return `RunError`*
  — is discharged for the psp path only; folding it into the cache is one of the four cleanups
  `cohort_merge.md` §8 gates on `RunError` landing in the merge, and none of those is this
  plan's.
- **Nothing constructs one yet.** E1 opens the cohort of readers and E3 drives them through the
  calling loop that D2 lifts out of `AlignedFilesVariantCaller`; until then the type's only
  callers are its own tests.
- **The source hands over every record its file holds**, and knows nothing about the run's
  analysed regions or the observation reach ceiling. Both belong to the calling stage — E1
  compares the analysed regions across the cohort, E4 reads the ceiling where the merge needs
  it — and the type says so.
