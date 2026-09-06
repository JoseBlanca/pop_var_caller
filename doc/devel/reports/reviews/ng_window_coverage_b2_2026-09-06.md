# Review — window coverage B2: the measurement behind the cheap rule

**Date:** 2026-09-06
**Reviewed:** commit `64c4b072` on `ng-window-coverage` (plan step B2 of
[window_coverage.md](../../ng/impl_plan/window_coverage.md))
**Impl report:** [ng_window_coverage_b2_2026-09-06.md](../implementations/ng_window_coverage_b2_2026-09-06.md)
**Per-category files (audit trail):** `tmp/review_2026-09-06_window-coverage-b2/`

## 1. Scope

One new example — a measuring program — and one doc comment.

- **In scope:** `examples/ng_window_coverage_probe.rs`, the doc-comment change in
  `src/ng/window_coverage/depth.rs`, and every number the step's prose claims.
- **Out of scope:** the rest of `src/`; production (`src/sample_summary/`, `src/paralog/`,
  `src/var_calling/`); the four checks red on `main` before this branch.
- **Categories dispatched:** two agents, each in its own worktree detached at `64c4b072` — one on
  `reliability` and `errors` plus the claims check, one on `naming`, `idiomatic`, `smells` and
  `refactor_safety`.

## 2. Verdict

**Approve with changes.** The probe measures the right thing — one reviewer verified the chain
rather than assuming it, from the rule's own branch through `LocusSummary::from(&streamed.head)`
on the psp path to `num_obs_along_locus` — and both reviewers reproduced the report's two tables
exactly against the real stores. What the review changed is the difference between a reading and a
check: as committed, the program could report the rule confirmed having checked nothing, and the
evidence offered that its comparison could fail did not exercise the mechanism the claim is about.

## 3. Top three

1. **"Zero disagreements" and "nothing was checked" were the same verdict.** The only assertion
   was that the disagreement count is zero, which holds vacuously over an empty set. A reviewer
   wrote a psp with no records and ran the shipped probe on it: every count printed `0` and the
   process **exited 0**. **Fixed:** it now asserts, per store, that a record covering one base was
   seen, and names the file when none was.
2. **The discrimination check the report offered did not exercise the shape the claim is about.**
   Pointing the comparison at the wrong head field fires on 197,051 records *none of which carry a
   partial witness*. Two mutations show the gap — disabling the partial-witness detector, and
   comparing the head against a recomputation of its own quantity (which the psp decoder already
   enforces, so it cannot fail on a readable file) — and **both leave every printed number
   identical on both tomato stores**. A reviewer built the missing control, round-tripped it
   through a psp, and confirmed the shipped probe catches it. **Fixed:** that fixture ships as a
   test, and it asserts both quantities off the record itself, so it pins that the shape is a
   genuine disagreement rather than one the probe merely labels.
3. **Two numbers in the report were wrong, both flattering the rule.** "About 1 record in 900" is
   1 in 871. "999 positions in every 1,000" was the *record* share relabelled: a wide generic
   record reports one position and a wide tract its whole span, so the positions figure is **992
   per 1,000 on the slice and 994 on the wider store**. **Fixed** in the prose, and the probe now
   prints the positions row itself so the claim comes from the output.

## 4. What's good

- **The chain was verified rather than assumed.** "One base" is tested equivalently on both sides
  (`len() == 1` ⟺ `start == end` ⟺ `min == max`); `records()` steps with `|_| true` so every body
  is built; every decode failure and every unreadable store aborts loudly — confirmed against a
  missing path, a non-psp file, and an older-layout ng store.
- **Both reviewers reproduced both tables exactly**, independently, on the real stores.
- **The mechanism the report explains is the one in the code**, checked by reading
  `num_obs_along_locus` against `non_reference_and_compared_reads` — and the same reading found
  that the report's two headline zeros are one fact rather than two, which the report now says.
- **The counter-example shape is reachable**, which is what stops the measurement being vacuous by
  construction: `ReadWitness::from_left`'s own doc says a saturating reach still answers `Partial`,
  and the psp codec preserves it — round-tripped and caught.

## 5. Findings, and what was done with each

### Major

- **M1 — the probe could pass having checked nothing.** Fixed with a per-store assertion naming
  the file. **Categories:** reliability, refactor_safety (convergent).
- **M2 — no positive control for the mechanism under test.** Fixed by shipping the fixture as a
  test; five tests now, where the file had none.
- **M3 — `add` and `print` enumerated nine fields by hand**, so a tenth would compile while going
  unsummed and unprinted — and plan step C3 adds fields here. Fixed by destructuring `Self` with
  no `..` in both, making a new field a compile error in each.

### Minor

- **Mi1 — the walk and its one measurement were welded together**, so C3's second measurement
  would have meant editing five sites, two of them silent if missed. Fixed with a
  `RecordMeasurement` trait: the walk decodes and shows each record to a slice of measurements,
  each owning its own arithmetic and its own output. C3 costs one type and one slice entry.
- **Mi2 — `first().copied().unwrap_or(0)` could only ever turn a missing answer into an
  agreement.** Now a slice pattern that panics with both numbers.
- **Mi3 — every disagreement was retained to print twenty.** Measured at 7,605,320 values, each
  with a heap string, on a run where the rule failed wholesale — allocated before the report the
  probe exists to print. Now capped at what is printed.
- **Mi4 — the wide records' bases were pooled across kinds**, so the positions claim could not be
  derived from the output. Split, and the positions row printed.
- **Mi5 — a one-base record's kind was never looked at**, so spec §3.1's "a single-base record is
  generic by construction" went unchecked. Now counted: zero in both stores.
- **Mi6 — the probe decodes bodies against a live set the run will not have.** Identical today
  because the encoder writes no chain ids; it stops being identical at the psp path's Milestone E.
  Recorded in the probe's module doc, so the day it changes the reader knows the measurement has
  to be retaken.
- **Mi7 — `!= ReadWitness::Complete` absorbed a future third witness kind.** Now an exhaustive
  match, which would be a compile error instead.
- **Mi8 — the output never said which file a row came from**, and two stores here share a
  basename. The path is printed once per store.
- **Mi9 — the totals line was printed only for two stores or more**, and a run that died part-way
  through several stores left a plausible shorter report. Printed unconditionally now, beside a
  `stores-walked` line.
- **Mi10 — the counter-example line re-read the head rather than carrying the value the
  comparison used**, which under a broken comparison printed two equal numbers beside the word
  `counter-example`. The compared value is passed in.
- **Mi11 — the measured share was quoted twice in two shapes, and named neither its stores nor
  their depth.** Quoted once now, with both stores and their mean head depth; the second site
  refers to it.

### Nits, applied

The output key now matches the field name it prints; the wildcard arm carries the
`// REVIEW ON UPGRADE:` marker; a duplicate disagreement counter that had to be kept equal to a
struct field by hand is gone; the module doc's "one block per store and one for the total" is now
true unconditionally rather than only for two stores or more.

### Not fixed, with reasons

- **The psp round-trip as a shipped test.** A reviewer verified it once — a one-base record with a
  full-span partial witness survives `PspWriter`/`PspReader` and the probe catches it — but
  building a header is about thirty lines of unrelated setup in an example, and the claim belongs
  to `src/ng/psp/`. Recorded in the report's follow-ups instead of pinned.
- **A borrowing form of `num_obs_along_locus`**, which allocates a vector per record where only
  the first entry is read at a one-base one. It is in `src/ng/locus_generation/` and out of this
  step's scope; plan step C3 needs every position's depth anyway and is where the question lands.
- **`samtools idxstats` on a bench CRAM**, to verify independently that the tomato CRAMs are cut
  to the benchmark's intervals rather than taking `bench.config.sh`'s word for it. The deviation
  the report records does not turn on it: what matters is that no readable whole-genome ng store
  exists here, which was checked directly.

## 6. The claims check

Every figure in the implementation report and in `PROJECT_STATUS.md` was re-derived, both tables
by re-running the commit's own probe against the real stores. **Three were wrong:**

1. "About 1 record in 900" — 1,164,208 / 1,336 is **1 in 871**.
2. "999 positions in every 1,000" — the record share relabelled; the positions figure is **992 per
   1,000 on the slice and 994 on the wider store**.
3. `PROJECT_STATUS.md`'s "every large `.psp` on this machine is production's format" — **wrong**:
   under `/Users/jose/devel` there are 330 ng-format `.psp` files above 1 MB. The fact the sentence
   was reaching for is that no ng store *this build can read* and none over a whole genome exists
   here; the one large ng store outside this branch's scratch is written to an older record layout
   that this build refuses. The report's scoped version was correct for
   `/Users/jose/devel/pop_var_caller` (1,568 of 1,568 production's) but inverted for the branch's
   own checkout, where all 13 are ng's.

Everything else checked out: both tables reproduced exactly, the totals, the wide-record splits,
the 40×/6.5× comparison, the BED's 80 intervals and 8 Mb over 12 contigs, the discrimination
count, and the mechanism explanation. One nuance the report now carries: the two headline zeros
are one fact, not two independent ones.

## 7. Validation after the fixes

In the container, on this worktree:

- `cargo test --release --example ng_window_coverage_probe` — **5 passed, 0 failed**.
- The probe re-run on both stores after every change — the report's tables, plus the positions and
  depth rows the review asked for.
- `cargo build --release --example ng_window_coverage_probe` clean; no clippy diagnostic for this
  example; `rustfmt --check --edition 2024` clean on both changed files.
- `cargo test --lib --all-features` — 6,343 passed, 0 failed, 15 ignored, unchanged by this step.

## 8. Deferred, with a home

- **The psp round-trip test** — `src/ng/psp/`, where the codec's preservation of a full-span
  partial witness is that module's claim to pin.
- **A borrowing `num_obs_along_locus`** — plan step C3, which needs every position's depth.
- **Retaking this measurement once chain ids are written** — the psp path's Milestone E; the
  probe's module doc names the dependency.
