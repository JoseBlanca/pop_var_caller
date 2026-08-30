# Code review — ng calling loop D1: the driver

**Scope:** the working-tree diff of step D1 of
[`calling_loop.md`](../../ng/impl_plan/calling_loop.md), on top of `5aa9c47d` — the
`SummariseConditionLoop` arm, the table build, the row map on `CallingScratch`, and the reworked
`run_frequency_loop` and `summarise_final_pass`.
**Date:** 2026-08-26. **Verdict: request changes** — **1 Blocker, 8 Majors, 12 Minors, 6 Nits**, and
**10 of 39 of the diff's own claims wrong**. All applied; see
[the fix report](fixes_applied_2026-08-26_d1.md).

**Three agents, each in its own worktree**, each detached at `5aa9c47d` with the diff applied as a
patch:

| agent | brief | outcome |
|---|---|---|
| reliability | tests and mutation testing | **21 mutations run, 6 survived, 2 changed no behaviour** — 1 Blocker, 5 Majors, 4 Minors, and eight test bodies with measured numbers |
| craft | naming, errors, idiomatic, smells, refactor-safety, module structure | 3 Majors, 8 Minors, 5 Nits; every proposed refactor compiled in its own worktree |
| numbers | re-derive every claim the diff makes about itself | **39 claims checked, 29 correct, 10 wrong** — every wrong one a mechanism |

---

## Blocker

**The join this step exists to get right was untested, and the mutation is silent.** The table
build reads `per_sample[scratch.run_sample_of_each_row()[row]]`; replacing that with
`per_sample[row]` — exactly the "row filled for sample *i*, read for sample *j*" defect the design
warns about — left the whole suite green. **Every set-aside fixture put the uncallable sample
last**, where a row's index and its sample's index are the same number and a table filled by row is
right by accident. On a fixture with the gap **first**, the clean code returns `1/1` with copies
`[0.0077, 1.9923]` and the mutant returns `0/0` with `[2.0000, 3.9e-6]` — a systematic permutation,
no panic. **Fixed** by that fixture, and by the final pass reading the same map the table build
reads rather than re-deriving it (see M2).

## Majors

**M1 — the record claimed a warrant it had not earned.** `call_locus` passed
`Provenance::FittedHere` unconditionally. The fixtures' calibration is
`ReadGroupCalibration::defaulted()`, whose own provenance is `Defaulted`, so the record was already
wrong on the code as it stood — the exact failure `LocusInference::weakest_provenance` exists to
prevent. **Fixed** by deriving it: the weakest warrant of the calibrations of the read groups whose
reads reached the locus, folded with `Provenance::weaker_of`. **And the field's doc was stale** —
it said `Provenance` defines no ordering and that the first step to compare two must settle it;
`parameter_estimation` states the ladder and implements it. Corrected.

**M2 — the row-to-sample map was derived two ways and only the counts compared**, which a
*permutation* satisfies. The final pass re-derived "which samples have rows" by inlining the
predicate and advancing a bare cursor. **Fixed**: it reads the stored map and asserts at every
sample that the map and the candidate step's ruling agree, in both directions. The count check
stays, for the one disagreement the per-sample one cannot see — rows claimed past the end of the
run — which now has its own test.

**M3 — three `sample_count()` accessors, and this step changed what one of them means.**
`LocusEvidence`'s and `FrozenParameters`' count the run's samples; `CallingScratch`'s now counts
only the callable ones, and two of the three appear forty lines apart in one function as
deliberately different numbers. **Fixed** by renaming the scratch's to `row_count()`, with the
field, the `prepare_for_locus` parameter and every message following.

**M4 — two doc comments stated the claim-then-prepare order that the code refuses.** A caller
following `claim_row_for`'s doc writes the bug; swapping the driver to that order fails two tests.
**Fixed.**

**M5 — the driver's `config` argument was unpinned.** Dropping it turns a capped locus
(`passes = 1, converged = false`) into a converged one. **Fixed** by a test at a cap of one.

**M6 — every driver fixture was outbred**, so the per-sample inbreeding lookup was invisible.
**Fixed** by a fixture at F = 0 beside one at F = 0.9: identical reads, qualities 32.32 and 20.75.

**M7 — the emission counter was asserted only on the shape its own doc names as the one that hides
the bug** — two samples of one observation each. **Fixed** by one observation beside two (6 against
a wrong 4), by a partial-carrying fixture, and by a four-pass locus that still builds once.

**M8 — the two refusals fired after the worker's scratch had been prepared and its rows claimed.**
The release profile aborts on a panic, so one such locus ends a whole cohort run — from three
frames down, having left a shared scratch prepared for a locus nobody scored. **Fixed** by moving
both to the seam's front door; the restatements inside the table build are now `debug_assert` and
`unreachable!`, because a release check no test can reach is one the suite cannot keep honest.

## Ten wrong claims, every one a mechanism

Every *counted* figure in the first draft was right — the test counts, the battery, the copies, the
`EmissionCost` triples, and the design-document quotes. The ten failures were all explanations:

- **Four doc comments told one stale story about a mis-sized row map**, and the two directions are
  each other's answer swapped: short panics with `index out of bounds` (loud, but naming a slice),
  long is the silent one. The prose was accurate about the *deleted* walk.
- **"A `RunnableCallingLoopConfig` cannot hold one and the bodies cannot be reached"** fails both
  halves: the private tuple field is visible to every descendant of `inference`, which is where the
  loops live, and the bodies run at every locus and hold all the work. What makes each run *once*
  is the unconditional `break`.
- **The `never_loop` reason's mutation claim**: `break` → `continue` is an `E0384`, and with the
  binding made `mut` it spins forever. Not "one build per round".
- **"Prepared first, claimed second … would hand this locus the last one's samples"**: the wrong
  order is caught immediately with zero claimed rows.
- **"Moved to where the two lengths are set"**: the check fires at the map's first read.
- Two counts: "four places" for the alternative shape is five, and "four fitted numbers" in an SSR
  scoring context is five fields of which two are fitted.

**One claim the reviewer singled out and confirmed:** the 0.015 of a copy in the driver fixture
really is the seed's pull toward the reference — a symmetric seed gives exactly 0.0, a mirrored one
−0.015067, and it grows monotonically to +0.443 at `α_ref = 100`.

## Minors and Nits, applied

`pooled_primary_alternative`'s successor already renamed in C3b; three `expect()` calls given the
repo's `// PANIC-FREE:` marker; `error_spreads` sized with `checked_mul` like its six siblings;
`#[expect]` rather than `#[allow]` on `never_loop`; `EmissionCost` made `pub(crate)`;
`called_inbreeding` → `inbreeding_coefficient_by_row` and `called_run_sample` →
`run_sample_of_each_row`; `GenericLocusSample::is_callable` as the one spelling of the predicate
four sites tested; the output vector renamed `calls` where `per_sample` named an input three lines
away; a 14-space run inside a panic message.

## Accepted without change

- **The long functions do not want splitting** — both are flat. The module question the craft agent
  did raise is for later: five helpers in the file are arm-independent, and a second arm of the
  seam would need all five verbatim.
- **`EmissionCost` counting from the shapes each row build was handed** rather than from inside the
  emission, which the SNP/indel path has no seam for. Its doc says which it is.

## Verification

After the fixes, in the container:

- `cargo fmt --all -- --check` — exit 0; `cargo clippy --all-targets --all-features -- -D warnings`
  — exit 0.
- `cargo test --lib` — `4691 passed; 0 failed; 14 ignored`.
- `cargo test --release --lib ng::calling --all-features` — `645 passed; 0 failed; 3 ignored`.
- The six release-held checks downgraded together: `638 passed; 7 failed`, every check reached.
