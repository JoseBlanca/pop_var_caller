# ng psp mode — E2: a stored record stops being one file's and becomes the run's

**Date:** 2026-09-04
**Plan step:** [run_driver_psp_mode.md](../../ng/impl_plan/run_driver_psp_mode.md) Milestone E, step E2
**Spec:** [run_streaming.md](../../ng/spec/run_streaming.md) §6.1, §6.2
**Branch:** `ng-psp-mode`

## Plan

Every psp is walked on its own, so every one numbers its read groups from zero and a cohort's
files collide by construction. E1 merged their tables into one run-wide numbering; E2 is where
each drawn record is renumbered into it.

## Assumptions

- **The map is a constructor argument, not an option.** `PspObservationSource::over` and `new`
  take the sample's map, so a source cannot be built without one. A source that sometimes
  renumbered would be a source somebody forgets to hand a map to, and the failure is silent: the
  numbers are all in range, and four libraries would be scored against two calibrations.
- **A single-sample run's map is the identity and is applied anyway.** The loop is wasted work
  there; making it conditional would buy a branch and cost the guarantee above. Named in the
  field's doc rather than optimised.
- **The renumbering happens at the draw**, which is the one place a record crosses from the
  file's terms into the run's. Doing it in the merge or at the call would mean every consumer
  knowing which sample a record came from in order to read its read group — the thing the
  run-wide numbering exists to remove.

## Changes made

- `PspObservationSource` gains the sample's map and applies it to every observation of every
  record it hands over.
- `PspSourceError::ReadGroupNotInThisFilesTable` — a record naming a group past the end of its
  own file's table is the file disagreeing with itself, so it is refused naming the record, the
  number and the table's size, rather than left to panic in the middle of a cohort where nothing
  would say which file to look at.

## Tests added

Two, and `cargo test --lib 'ng::run'` goes from 496 to 498.

- **Two samples whose read groups collide, renumbered apart** — the plan's own fixture. Both
  files call their groups 0 and 1; the first keeps its numbers and the second's start where the
  first's ended, and every stored group reaches the merge under the number this run gave it.
- **An observation naming a group the file does not declare is refused**, with the record, the
  number and the table's size in the error.

The module's other fifteen tests now pass the identity map, which is what a single-sample run
has, so each stays about what it was about before the numbering existed.

### Mutation pass

Three, all killed (`tmp/e1_mutations/e2.sh`):

| mutation | killed by |
|---|---|
| the renumbering computed and dropped | the collision fixture |
| an unknown group falls back to zero | the unknown-group refusal |
| every observation renumbered as if it were group zero | both |

## Validation results

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` exit 0.
- `cargo test --lib 'ng::run'` — **498 passed**, 0 failed.

## Tradeoffs and follow-ups

- **Nothing builds these sources from a cohort yet** — E3 does, and it is where
  `OpenPspCohort::read_group_remap` meets `PspObservationSource::over`. Until then the map is
  handed over by tests.
