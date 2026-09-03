# ng psp mode — C2+C3: what a walk says, and what it refuses to replace

**Date:** 2026-09-03
**Plan steps:** [run_driver_psp_mode.md](../../ng/impl_plan/run_driver_psp_mode.md) Milestone C, steps C2 and C3 — **bundled deliberately**: C2 is a report over the loop C1 already made sequential, and C3 is a flag and a refusal in that same loop, in one file. Two review fan-outs over one small coherent diff would have reviewed the same code twice.
**Spec:** [run_streaming.md](../../ng/spec/run_streaming.md) §5.2, §8; `psp_file_format.md` §10
**Branch:** `ng-psp-mode`

## Plan

C2: the run says what each sample stored. C3: a run never silently replaces a psp, and leaves
nothing that reads as whole when it stops.

## Assumptions

- **The report is a value, not a side effect.** `walk_every_sample` returns a `WalkReport`
  and `run_generate_psps` prints its `lines()`. This is the sibling command's own split, made
  for a reason its note records: the printing was the one part of that command a mutation
  could change with the whole suite still green.
- **Two moments, two audiences.** A progress line is printed as each sample finishes (a
  cohort of sixty is an hour; a command that says nothing until it is done cannot be told from
  a hung one), and the report proper at the end.
- **What the report carries** is what a psp cannot say for itself: how much of the analysed
  ground the walk actually spoke for. A segment with no generator is analysed ground this
  build cannot describe, and in a psp that is indistinguishable from ground where nothing was
  there — so the two kinds of uncovered ground are named separately (*not built yet*,
  temporary; *out of scope*, permanent).
- **The overwrite check runs before any sample is walked**, not per sample as it goes: a
  cohort whose second psp is already there must not spend the first sample's hours before
  saying so.
- **A stale `.partial` is this command's own scratch** and is overwritten freely; only the
  final `<sample>.psp` is protected.

## Changes made

- `generate_psps.rs`: `--force`; `GeneratePspsCliError::PspAlreadyThere`; the door check
  extended to every sample's psp; `walk_every_sample` split out of `run_generate_psps`;
  `SampleWalkOutcome` and `WalkReport` with `lines()`; the per-sample progress line;
  `write_psp`'s `(WriteStats, LocusCounts)` kept rather than discarded.

## Tests added

Seven, taking the command's own tests from 24 to 31: the report naming every sample with what
it stored and how much ground it covered (and distinguishing the sample with reads from the
one without); the per-sample line actually carrying its numbers rather than only the name; a
re-run refusing and leaving the psp untouched; the refusal arriving before the first sample is
walked (only the *second* sample's psp is in the way, and the first must not have been
written); `--force` replacing it; a truncated psp refused by a reader; and a stopped walk
leaving nothing at the sample's path or beside it.

## Validation results

In the container, on the committed tree:

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` exit 0.
- `cargo test --lib 'pop_var_caller_exp'` — **119 passed**, 0 failed, 1 ignored (was 112).
- **Run for real**, release build, one tomato accession over two 100 kb intervals:
  - first run: `SRS3394712: 193603 observations in 4 blocks, 914715 bytes -> …`, then
    `walked 1 sample(s) over 200000 analysed bases` /
    `SRS3394712: 193603 observations, 914715 bytes at …; spoke for 311 of 318 segments,
    199672 of 200000 analysed bases; not built yet 328 bases, out of scope 0 bases` /
    `1 psp(s), 914715 bytes in total`. Exit 0.
  - second run, no `--force`: `error: SRS3394712 already has a psp at …; pass --force to walk
    it again and replace it`. Exit 1, the psp untouched.
  - with `--force`: replaced, **914,723 bytes — 8 more than before**, which is the recorded
    command line being 8 characters longer (` --force`). A useful accident: it shows the
    header really does record the invocation.
- One defect found and fixed during the run rather than by a test: the per-sample report line
  was built from a format literal wrapped across source lines, which put eighteen spaces of
  its own indentation into the middle of the sentence. It is now assembled by pushing
  fragments, with a comment saying why.

## Tradeoffs and follow-ups

- **`LocusCounts` has eight fields and the report prints four.** The region counts and the two
  base counters are what answer *how much ground did this walk speak for*; `loci_emitted` is
  the record count under another name, and the two unhandled *region* counts are less useful
  than their base counterparts. Recorded rather than assumed — if a per-sample report needs
  the region counts too, they are one line away.
- The command still writes no census beside each psp (Milestone G), and the walk-stage
  carry-forwards from C1's review remain: the duplicated reference-open block, the
  `#[command(flatten)]` simplification for the routing flags, and the shared on-disk cohort
  fixture.
