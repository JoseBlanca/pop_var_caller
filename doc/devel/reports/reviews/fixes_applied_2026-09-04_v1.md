# Fixes applied — psp mode D1, the psp-backed observation source

**Date:** 2026-09-04
**Review:** [ng_psp_mode_d1_2026-09-04.md](ng_psp_mode_d1_2026-09-04.md)
**Branch:** `ng-psp-mode`
**Outcome:** every finding applied; nothing deferred, nothing routed to the owner.

## The Major, and what it cost to fix

**A refusal now ends the source.** `PspObservationSource` gains a `refused` flag, set by each
of the two refusals it mints itself, and every later draw returns the new
`PspSourceError::AlreadyRefused` — not `None`, and not the next record.

The measurement that forced it, from the review: after an out-of-order refusal the next draw
returned `Some(Ok(contig 0:201-201))`, so the refused record at position 1 was gone from the
stream and nothing said so. `None` would have been the same silence in a different shape — the
merge reads it as a sample that ran out.

The doc paragraph that was wrong now names all three failure paths separately: a failure inside
the walk fuses the walk (the walker's deviation, and the silent one), and the two refusals
latch. Pinned by `a_source_that_refused_a_record_refuses_every_later_draw`, which draws twice
past the refusal and asserts the cause both times.

## The Minors

- **`reached` pinned against a multi-base observation.** `reached_is_the_last_base_the_last_observation_covered`
  uses a six-base record: it starts at 11 and reaches 16. The mutation the review found
  surviving all eleven tests (`reach_position` → `start_position`) now fails it.
- **One place shapes a psp source's failure.** `source_failed(sample, reached, cause)` is called
  both by `refuse` and by the constructor, where a duplicated `RunError::SourceFailed` literal
  used to sit. The test named for that constructor never called it; it is renamed
  `a_source_that_fails_on_its_very_first_draw_reports_nothing_reached`, which is what it does
  test and a class the suite otherwise lacked. **The constructor's own error arm stays
  untested and now says so in the code**, with the reason (no fixture reaches it through a
  `PspReader`: the seek is inside offsets `open` already bounded, the manifest `open` has
  already parsed, and a bad buffer ceiling is refused where it is set) and the mutation that
  survives.
- **The enum's doc no longer contradicts its variants.** What the variants share is that the
  merge would otherwise meet them as an assertion or as silence — not whose mistake they are,
  which differs: an out-of-order file is damaged input, a head-only record is this crate's own
  wiring. The head-only message now names the walk rather than the file, so nobody is sent to
  rebuild a psp that is sound.
- **`StreamedRecord` is destructured with no `..`**, so a field added to it has to be considered
  in the function step E2's read-group remap lands in.
- **"Exhaustion is final" is pinned** by `a_spent_source_answers_none_for_ever`.
- **The out-of-order error carries two start positions**, rendered the same way, because start
  positions are what the check compares. `offered`'s doc and its value now agree.
- **A body declined part-way through** is a fourth new test: the existing head-only fixture
  refuses on the first record, where a `reached` that was never set and a `reached` that is
  right both read `NothingYet`.
- **`a_psp_of_under` inlined** into its one caller.

## The nits

`new`/`over` now name the same roles as the walker's `new`/`over`, and `new` is `pub(crate)` —
which is what makes the sample name trustworthy: it is whatever a caller passes, and three
documents asserted the name comes from the file. A hand-written `Debug` replaces the derive,
matching the sibling and dropping the `W: Debug` bound. The damage fixture guards its block
index. `run/mod.rs`'s doubled "and" is gone, as are "and this is that day", the "three things"
whose third is a discarded parameter, and the report's "the tail D2 lifts".

## Validation

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` exit 0.
- `cargo test --lib 'ng::run'` — **474 passed**, 0 failed (459 at Checkpoint C, 470 before these
  fixes).
- `cargo test --lib 'ng::run::psp_source'` — **15 passed** (11 before).
- Mutation pass re-run with the review's two survivors included: **8 mutations, 7 killed**
  (`tmp/d1_mutations/run2.sh`). The survivor is the constructor arm named above, and it is
  marked uncovered in the code rather than left to look tested.
