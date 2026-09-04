# ng psp mode — E1: a cohort of stored samples, opened and refused before a block is decoded

**Date:** 2026-09-04
**Plan step:** [run_driver_psp_mode.md](../../ng/impl_plan/run_driver_psp_mode.md) Milestone E, step E1
**Spec:** [run_streaming.md](../../ng/spec/run_streaming.md) §5.3, §6.1, §6.2, §7.1a; arch [run_streaming.md](../../ng/arch/run_streaming.md) §3.4, §5
**Branch:** `ng-psp-mode`

## Plan

Open every psp of a cohort, read every header, and run every refusal spec §6.2 asks for — all
of it before a single block is decoded, so a cohort of thousands is refused before any of it is
spent.

## Assumptions, and two departures from the architecture

- **Opening is two moves, not the one the architecture sketches.** `OpenPspCohort::open` reads
  every header and runs the checks that compare the files *with each other*; the caller then
  builds the run's segmentation over the ground those headers agree on; `PspVariantCaller::open`
  runs the checks that compare each file *against the run*. The sketch has one `open` that takes
  the catalog and rebuilds the segmentation itself, and it cannot: the one copy of the ground
  assembly — `run_ground::segments_over`, shared by `call-from-alignments` and `generate-psps` —
  lives in the command module, and a pipeline stage reaching down for the commands that drive it
  inverts the dependency direction. Keeping one copy of that assembly is what makes the
  file-against-run check mean anything. **Recorded in arch §3.4.**
- **`PspVariantCaller` went into `psp_caller.rs`, not `callers.rs`.** That file was 6,559 lines,
  and every other run object in the module has its own. **Recorded in arch §3's file tree.**
- **`IncompletePsp { sample }` became `PspNotRead { path, source }`.** A reader reaches a psp's
  header third — footer, then index, then header — so a file whose footer will not parse has no
  sample name to be reported under. The interrupted case is not lost: it arrives as
  `PspReadError::Incomplete` under the cause, rendering *"the writer did not finish"*.
  **Recorded in arch §5.**

## Changes made

- **`src/ng/run/psp_caller.rs` (new).** `OpenPspCohort` — the open readers, the ground they
  agree on, the run-wide read-group table, and the widest observation reach any of them
  declares. `PspVariantCaller` — that cohort checked against the run it is about to be called
  by, plus the parameters and settings the calling loop will be handed at E3.
- **Five `RunError` variants** beyond the three the architecture listed: `NoPsps`, `PspNotRead`,
  `PspAgainstAnotherReference`, `PspReadGroupsCannotBeMerged`, and
  `NotEnoughFileDescriptorsForPsps`.
- **`ReadGroups::of_merged_tables`** — a run-wide table assembled from tables built elsewhere,
  where the ordinary constructor reads alignment headers and mints identifiers itself. It
  asserts its two views agree, including that the sample claiming a read group is the sample
  that read group names.
- Three helpers in `callers.rs` widened to `pub(super)` so psp mode makes direct mode's checks
  rather than second copies of them.

## The refusals, against spec §6.2 clause by clause

| clause | where |
|---|---|
| across the cohort: analysed regions equal, naming both samples | `the_ground_every_file_agrees_on` |
| each file against the run: segmentation inputs match, naming the first field | the per-file loop, via `SegmentationInputs::first_difference` |
| across the cohort: two psps may not name one sample, naming both files | `refuse_a_sample_named_twice` |
| each file against the parameters | a **count**, and §6.2 asks for a match by **name** — see below |
| read-group tables merged, not compared | `merge_the_read_group_tables` |
| refused: a table that cannot be merged | the empty table only — see below |

Beyond §6.2, three checks direct mode makes and this now makes too: the descriptor budget
(§7.1a, before the first `open(2)`), the catalog against the run's reference, and the two views
of the reference being one genome. Plus one the psp header was designed for and nothing was
reading: the file's own whole-assembly digest.

### Two things for the owner at Checkpoint E

- **The duplicate-`@RG ID` refusal is deliberately not made**, and §6.2 names it. The reason the
  spec gives is that such a table "cannot be renumbered without guessing" — which is not this
  format's situation: identity is the walk-local *number*, which is the entry's own position,
  checked on both sides, and nothing in the merge reads the id. The psp format's own validator
  declares the case legal in as many words, because **a psp holds one sample, not one alignment
  file**, and a sample sequenced across lanes may carry two entries with one id and different
  libraries. Direct mode calls that cohort. Refusing it here would break goal 1 for every
  multi-lane sample whose lanes reuse an id. **Recommendation: §6.2's clause should say what it
  means, which is the empty table.**
- **The parameters match is a count, where §6.2 asks for a match by name.** It cannot be made
  here: `RunParameters` is assembled per sample by position and carries no names at all. The
  by-name match belongs where the parameters *file* is read against this cohort's sample list —
  `ParametersFile::to_run_parameters_for`, which the subcommand calls at F1. Direct mode has the
  same shape. **Recommendation: E1's plan entry should say the by-name half is F1's.**

## Tests added

Twenty-one in the new module, plus five for `of_merged_tables` in `read_groups.rs`;
`cargo test --lib 'ng::run'` goes from 475 to 496.

Every refusal has a test that provokes it and reads what it names. **The per-file refusals are
provoked on the *second* file of a two-file cohort**, which is what tells a loop over every psp
from one that checks the first and stops: measured, with them pinned on one-file cohorts a
mutant that checked only `psps[0]` passed all 490 tests. Also pinned: the run-wide read-group
order *within* a file (a merge that walked a table backwards would put every observation on the
wrong calibration, and still hand back the same remap); that a cohort of one sample opens and
calls; and that a sample walked across two files sharing an `@RG ID` calls.

### Mutation pass

**18 mutations, 17 killed, 1 changed no reachable behaviour.** Ten before the review, all
killed. Eight more from the review's findings: five killed, three survived — the descriptor
check's call site, the two-references check's call site, and guarding the read-group-table
refusal. The first two were real gaps and are now killed by
`the_descriptor_refusal_comes_before_any_psp_is_opened` (which uses paths that do not exist, so
a run that opened first would name one of them) and `a_run_checked_against_another_genome_is_refused`.
The third is the changed-no-behaviour one: with the duplicate-id arm gone, the only remaining
arm is the empty table, which neither this store's writer nor its reader will produce, so when
that refusal is called cannot be observed through any psp.

## Validation results

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` exit 0.
- `cargo test --lib` — **6,150 passed**, 0 failed, 14 ignored. `--lib 'ng::run'` 496 passed.

Standing red, untouched: the three locus-dump example tests, and `--all-targets` aborting on the
`ng_joint_fit_perf` bench, which needs `--features bench-fixtures`.

## Tradeoffs and follow-ups

- **A psp records the library a walk resolved, not whether the file declared it**, so the merged
  table marks every library `Synthesized` — the weaker of the two claims, because the field
  exists to tell "ours" from "the file's" and `run_report` states the rule: *a synthesized name
  reported as a declared one is a claim about the run that nobody made*. A header field would
  let it say the true one; nothing reads the origin today.
- **Four fields of `PspVariantCaller` are `#[expect(dead_code)]` until E3** spends them, so the
  first real reader turns each line into a compile error.
