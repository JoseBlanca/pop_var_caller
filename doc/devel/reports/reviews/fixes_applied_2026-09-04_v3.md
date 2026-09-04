# Fixes applied — psp mode E1, opening a cohort of stored samples

**Date:** 2026-09-04
**Review:** [ng_psp_mode_e1_2026-09-04.md](ng_psp_mode_e1_2026-09-04.md)
**Branch:** `ng-psp-mode`
**Outcome:** every finding applied. Two are carried to Checkpoint E as owner questions rather than settled, and both are recorded where the code makes the choice.

## The three missing refusals

- **The descriptor budget, before the first `open(2)`.**
  `refuse_if_more_descriptors_are_needed_than_allowed_for_psps` and
  `RunError::NotEnoughFileDescriptorsForPsps`, whose message names the limit, the sample count
  and the arithmetic — one descriptor a psp, held for the whole run. `OpenPspCohort::open` now
  splits from `open_within_a_descriptor_limit`, so the refusal's *place* has a test and not only
  its message: the fixture's paths do not exist, so a run that opened first would name one of
  them.
- **The catalog against the run's reference**, which is direct mode's own refusal reused rather
  than re-implemented.
- **The psp's whole-assembly digest**, compared where both sides carry one — the field's one
  documented consumer, now consuming it, and the only comparison in this check that is exact.

## The tests that could not tell a loop from its first step

Every per-file refusal is now provoked on the **second** file of a two-file cohort, through a
shared `a_cohort_whose_second_file` helper that says why. The run-wide read-group order within a
file is pinned by asserting which `@RG` entry got which number and which library the remap's
second entry resolves to — not only how many groups there are.

## The rest

- **`of_merged_tables` checks the pairing**, not only the coverage, and has five tests of its own
  where it had none: the ordinary case and each of the four ways the two views can disagree.
- **The library is marked `Synthesized`** — the weaker of the two claims, and the one that cannot
  be false, since a psp records the library a walk resolved and not which of the two it was. The
  experiment is `Synthesized` too, which is direct mode's own rule.
- **The duplicated per-sample remap is gone**: the merged table already is the map, because a
  psp's header order is its walk-local numbering.
- `#[expect(dead_code, reason = …)]` replaces `#[allow]`, so E3 turns each line into a compile
  error rather than leaving a suppression behind.
- Six doc claims corrected: `sample_count`'s placement inside another method's doc comment; the
  segmentation check that *is* made; the header being the third thing a reader reaches; the
  layering constraint being about dependency direction rather than the compiler; §12.3 being the
  mode-equivalence oracle; and the memory figure for a stage that holds no cursor.

## Carried to Checkpoint E rather than settled

- **The duplicate-`@RG ID` refusal is not made**, and spec §6.2 names it. Recorded in full at the
  refusal, with a test — `a_sample_walked_across_two_files_sharing_a_read_group_id_calls` — that
  says why it must not be one.
- **The by-name parameters match is F1's**, not E1's, because `RunParameters` carries no names.
  Recorded at the call site.

## Validation

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` exit 0.
- `cargo test --lib` — **6,150 passed**, 0 failed. `--lib 'ng::run'` 496 passed (475 before E1).
- **18 mutations, 17 killed, 1 changed no reachable behaviour** — the last being when the
  read-group-table refusal is called, whose only remaining arm neither this store's writer nor
  its reader will produce.
