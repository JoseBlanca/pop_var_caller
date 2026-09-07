# window coverage — C4 review: the pair on its way to the record, and two oracles that would have passed

**Date:** 2026-09-07
**Reviewed:** the working tree of plan step C4, at `bcc2ef2f wip: C4 for review`
**Implementation report:** [ng_window_coverage_c4_2026-09-07.md](../implementations/ng_window_coverage_c4_2026-09-07.md)
**Fixes applied in:** C4's own commit

Two agents on the step's diff, each detached at it in its own worktree: one on reliability, errors
and every number the step's prose claims, one on naming, idiom, module structure and refactor
safety. **1 Blocker, found by both. 4 Major. 12 Minor.** All applied.

## Verdict

**The measurement itself is right where it matters.** The sparse-to-dense translation maps the
merge's covering-sample index to the run's the correct way round; the pair is read at the locus's
first base, the same position the ownership rule uses and the same one the recorder keys its rows
by; and the scratch buffer that carries the pairs to the sink cannot leak a stale locus, because
the closure filling it is `FnOnce` and the only outcome reaching the sink is that closure's own
return value.

**What was wrong was the checking.** The step broke the build of two integration tests, and the
command this session was using to declare a step green — `cargo test --lib` — does not build them,
so twenty assertions had stopped running with nothing to say so. And **the extended
mode-equivalence oracle passed when neither mode measured anything**: two files of absent pairs
are byte-identical, so the comparison agreed loudest exactly where the measurement was off — which
is the failure the step's own design decision was taken to avoid, reappearing one level up.

## Findings

### 1. Blocker, CONFIRMED — two integration tests stopped compiling

`tests/ng_calling_loop_calls_genotypes.rs` and `tests/ng_candidate_selection_truth_recall.rs`.
`CohortObservation::over` gained a parameter and `CohortObservation` gained a field; neither call
site under `tests/` was updated, so `cargo check --tests` failed with `E0061` and `E0063`. **A
target that does not build reports no failing test**, so the twenty assertions in those two
binaries — including fifteen that pass — had silently stopped running.

**Both agents found it independently**, and it is the review's most useful result, because the
green criterion this branch has been using could not have caught it. Fixed, and `cargo check --lib
--tests` is now part of what "green" means here.

It also exposed the cost of the constructor's new parameter: `CohortObservation::over` is `pub`,
so an out-of-crate caller has to build a `WindowedCohort` — an observation-cache type with three
public fields and no constructor. **`WindowedCohort::with_nothing_measured()`** is now the one
documented definition of "a window that answers nothing", and the module's own fixture delegates
to it.

### 2. Major, CONFIRMED — the mode-equivalence oracle passes with the measurement off

`scripts/ng_mode_equivalence_oracle.sh`. It refused an *empty* file but never checked that a row
carried a number. **Measured**: with the cache's window read stubbed to `None`, all rows came back
absent and the script printed `identical, bit for bit` and exited 0.

**The fix is a second guard**, requiring at least one measured row — an absent pair is two `NaN`s,
whose `f32` bit patterns are all at or above `0x7FC00000`. The success line now reads
`13,866 rows, 13,589 of them a measurement`, which cross-checks against the whole-store
comparison's own count of loci that agree. Verified in both directions: with the measurement
stubbed out the guard trips and the script says why.

### 3. Major, CONFIRMED — "the locus's first base" had no test

Every fixture exercising `CohortObservation::over` with a real measurement uses **one-base loci**,
where the first base and the last are the same position. Changing `region.start` to `region.end`
**passed the whole library suite**. The real-data probe catches it; the mode-equivalence oracle
cannot, because both modes would read the same wrong position.

Fixed by `the_window_is_read_at_the_locus_first_base_and_not_its_last`: a five-base locus with a
different window at each end. Re-applying the mutation fails exactly that test.

### 4. Major, CONFIRMED — the two per-sample orders met at an unguarded index

`records.rs`. The merge's windows are one per *covering* sample and the record's evidence one per
*run* sample; this is the only place the two orders meet. Swapping them was caught — but by an
unlabelled `index out of bounds` inside two tests about read counts, which blames the reader for
what the builder did. There is now an `assert_eq!` naming both list lengths, a line in the
function's `# Panics`, and `a_covering_samples_window_lands_on_the_run_sample_it_names` — three run
samples with the merge covering the first and the **third**, so the two indices differ.

### 5. Major — three doc comments had stopped describing the code

- The recorder explained a `-` row **the writer can no longer produce** — measured, 0 of 13,866 —
  and still priced itself per *built locus* and justified its lock by "the parallel merge builds
  several regions at once", when it now fires per *written record* from the calling thread.
- The `Eq` removal was blamed on `NaN`, which is the one thing `WindowCoverage`'s hand-written
  comparison already handles correctly. As written it taught "floats mean no `Eq`", and a reader
  who believes that is one step from comparing the floats with `==` — the defect the bitwise
  comparison exists to prevent. The real reason is that a derive cannot reach past a field whose
  type withholds `Eq`, and `WindowCoverage` withholds it because bitwise equality separates `+0.0`
  from `-0.0`. The claim that "nothing used `Eq`" was also checkable and now is: nothing in the
  tree compares two of these, and the only type holding them never derived `Eq` either.
- `render`'s doc still told the reader that a field-by-field comparison is the thing to avoid,
  immediately above a comparison now written field by field — and its key sentence needed two
  readings. Both rewritten: the destructure is what lets one field be left out, and the tempting
  repair is named ("asserting only that both sides are absent") with what it would let through.

### 6. Major — the oracle script's usage message lost the usage

Its help text prints a fixed line range, and the paragraph this step added pushed the command line
past the end. Replaced with "every comment line from the shebang to the first line that is not
one", so the next paragraph cannot silently truncate it again. A missing-file check was added
beside it: under `set -eu` an older binary that ignores the variable died on `sort: cannot read`
rather than on the script's own message.

### 7. Minor, applied — the rest

The recorder's name took the wrong preposition (`record_the_windows_of` against a first argument
that is a position) and its parameter was a prepositional phrase; the new test shadowed its own
inputs and checked the length invariant on one side only; three fully-qualified paths sat in the
hottest loop in `callers.rs`; the probe's comment described a row shape a run no longer writes.

## The question the brief put first, answered

**Did collapsing "absent" and "missing" hide the failure C3 exists to prevent?** No, and it is
measured. The comparison's "the run had no window here" count went from 12 to 0, and the reason is
the *locus set*, not the collapse: C3's recorder fired at every locus the merge **built** and C4's
fires at every record it **writes** — 26,754 sample-loci against 13,866 here, and on the reviewer's
own smaller slice 15,632 against 4,862, a 3.2× fall that predicts about 10 against the 12 observed.
The two positions those twelve sat at establish no variant and become no record.

**Without the collapse a correct run reports 12 spurious disagreements** on the reviewer's slice
and 277 on the six-accession one, every one of them the run saying *absent* where the walk says
*none*. The collapse is necessary, not merely convenient.

**One thing it does hide, and it is not reachable today.** A window the *run* silenced against one
the walk measured would be counted as "the run had no window here" rather than as a disagreement —
the same bucket a truncated cover's centres land in. Nothing on this data silences a window
(`windows-under-the-floor` is 0), so the case does not arise; **plan step D1, which sets the floor
from a distribution, is where it could start to**, and the probe now says so where the collapse is
written.

## Mutations run

| mutation | outcome |
|---|---|
| the sparse and dense indices swapped | caught — 3 tests, one of them now naming the defect |
| the pair read at the locus's **last** base | **survived** the library suite; caught by the probe on real data, **not** by the oracle; now caught by a test |
| the scratch buffer not cleared | **survived** the library suite; the recorded file grew to 5.9 million rows and the probe panicked; the oracle would have passed |
| the recorder given the previous locus's slice | **survived** the library suite; the probe reported 4,824 of 4,862 disagreeing; **not** the oracle |
| `render` including the window coverage | caught — 9 driver-agreement tests |
| `a_usable_window` made the identity function | caught by its own test; on real data turns a correct run into 12 disagreements |
| the oracle's sort removed | **survived** — the two files are already identical unsorted, which is why the script's own claim about thread schedule was wrong |
| the cached driver stops measuring | caught by two library tests; **not** by the oracle before finding 2's guard |
| the measured-row guard removed, measurement off | caught by the guard added for it |

**Three of the four defects that survive the library suite are caught only by the real-data probe,
and none by the mode-equivalence oracle** — which compares two routes that share the defect. That
is the shape of this milestone's evidence and worth stating plainly: the oracle proves the two
modes agree, the probe proves what they agree on is right.

## Numbers checked

Every figure in the step's prose was run rather than recalled, and three were wrong — all three
about how many rows a run writes and where, all three corrected with the measurement in them. The
brief's own "12 → 0" was checked from both ends and is the locus set, above. Mode equivalence
independently reproduced on a different slice: 2,431 records, both VCFs identical, parameters
identical.

## Looked at and found sound

The translation and its absent-pair default; the three agreeing readings of "the locus's first
base", and that an inverted region gives absence rather than a panic; the `FnOnce` argument that
makes the scratch buffer safe on every outcome path; `render`'s exclusion balanced against the test
that states the asymmetry — breaking either side fails tests; the switch from `Debug`-on-the-outcome
to `render` in the builders' differential, which removes the last such comparison in the module;
withholding `Eq`; the sink's second parameter as the right shape given that the filter's other half
is per-run and leaves at C5; and the collapse itself, which is necessary rather than convenient.
