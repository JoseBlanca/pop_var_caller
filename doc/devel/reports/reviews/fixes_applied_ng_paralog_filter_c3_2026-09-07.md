# Fixes applied — ng hidden-duplication filter, step C3 (pass two)

**Date:** 2026-09-07
**Review:** [ng_paralog_filter_c3_2026-09-07.md](ng_paralog_filter_c3_2026-09-07.md) — 1 Blocker,
8 Major, 16 Minor, Request-changes; three sub-agents in isolated worktrees
**Reviewed commit:** `14c64563` · **Branch:** `ng-paralog-filter`

## The answer

**Applied 23 · Applied with adaptation 2 · Deferred 3 · Disputed 0 · Failed validation 0.**

The Blocker and all eight Majors are addressed. The step's tests go from 13 to 23; the module suite
from 121 to 130.

**And applying the fixes found a tenth defect, in one of the fixes.** The test written for M5 — that
the range the run reports is the range its ratios were folded into — compared the recorded shape
against the constants and never against the histogram. A mutation that built the histogram over a
different range survived it. That is now closed by construction rather than by a test: the shape
builds its own histogram, so there is no second place where a range is written down, and the test
asserts the shape's own histogram is the shipped one.

## What was applied

| # | severity | finding | outcome |
|---|---|---|---|
| B1 | Blocker | nothing asserts an unscored record is never flagged | **Applied** |
| M1 | Major | `target_fdr` unvalidated; both wrong ends silent | **Applied** — `TargetFdr` newtype |
| M2 | Major | the fallback rate and the EM's seed are both `0.03` | **Applied** |
| M3 | Major | the warning's count checked only where both counts agree | **Applied** |
| M4 | Major | every fixture on contig 0 | **Applied** |
| M5 | Major | the verdicts record none of the constants that produced them | **Applied** |
| M6 | Major | "in spill order" is prose | **Adapted** — assertion + doc now; the cursor is C4's |
| M7 | Major | the calibration and its differential are in the run stage | **Applied** — both moved |
| M8 | Major | one held-constant dimension named where there are four | **Applied** |
| Mi1 | Minor | the `spill.read()` failure path unreachable from any test | **Applied** |
| Mi2 | Minor | the spill only tested short, never long | **Applied** |
| Mi3 | Minor | the single-sample test asserts only finiteness | **Applied** |
| Mi4 | Minor | a comment promising a count no caller keeps | **Applied** |
| Mi5 | Minor | spill order and sorted order coincide in the order test | **Applied** |
| Mi6 | Minor | `records_scored` names the wrong set | **Applied** — `records_in_the_fit` |
| Mi7 | Minor | two private one-line getters after the narrowing | **Applied** — inlined |
| Mi8 | Minor | the histogram's fixed ±100 against a cohort-growing ratio | **Adapted** — counted, not changed |
| Mi9 | Minor | a test passing on an exact floating-point zero | **Applied** — commented |
| Mi10 | Minor | a wildcard arm over a `#[non_exhaustive]` enum in-crate | **Applied** |
| Mi11 | Minor | the capacity reservation discards its conversion failure | **Applied** — commented |
| Mi12 | Minor | `mod.rs`'s six-clause sentence | **Applied** — a list |
| Mi13 | Minor | `Clone` on a type whose main field is one `f64` a record | **Applied** — dropped |
| Mi14 | Minor | `pass_one.rs`'s stale claim about C3's check | **Applied** |
| Mi15 | Minor | "three lines" against "four lines" | **Applied** |
| Mi16 | Minor | `ParalogVerdicts` holds no verdict | **Deferred** — spec §3.7 names it |
| 6a | numbers | "eleven behavioural tests" | **Applied** — twelve |
| 6a | numbers | the twenty-record fixture's mechanism | **Applied** — rewritten |
| 6a | numbers | "12 files fmt-dirty before, 9 after" | **Applied** — reworded |
| nits | — | test import order, doubled `pub mod`/`pub use` | **Deferred** |

## The three that carried the most

**B1 — the third of the contract nobody tested.** `an_unscored_record_is_never_flagged_however_loose_the_target`
runs the whole pass at four targets up to the loosest the type admits, asserts the unscored record
is never removed, and — so the test cannot pass by removing nothing — asserts that both *scored*
records are removed at that loosest target. The measured number that makes it matter is in the
test's own doc: the curve's answer for a value that is not a number is `0.4999999999999852`, inside
a target of one in two by about 1.5 × 10⁻¹⁴, so the finiteness screen in the frozen verdict is the
only thing holding the record back.

**M2 — two constants that are the same number.** `DEFAULT_EM_START` and
`DEFAULT_FALLBACK_PARALOG_PRIOR` are both `0.03`, and an empty histogram's estimate *is* the seed —
so no test could tell "the fallback was substituted" from "nothing happened", nor "the fallback was
read" from "the seed was read". Every test that asserts the fallback now runs on a configuration
where the two are `0.41` and `0.17`, and a new `the_shipped_fallback_rate_is_the_documented_one`
pins the number an ordinary run uses. The differential against production was widened from two
configurations to three for the same reason, and **production's configuration is now built from
ng's field by field** — the two sides drifted into comparing different questions the moment ng's
fallback was separated, which the differential caught on its first run.

**M7 — a placement whose stated reason was wrong.** `calibrate_from_the_ratio_histogram` now lives
in [src/ng/paralog/calibrate.rs](../../../../src/ng/paralog/calibrate.rs) and its differential in
[production_parity.rs](../../../../src/ng/paralog/production_parity.rs) beside its nine siblings.
The original reason — that `calibration.rs` may not gain a line — is true of that file and not of
the module, which already held three files of ng's own. The copy guard's inventory now names four,
with the reason `calibrate.rs` cannot be a guarded copy: production keeps its counterpart below
four items ng deliberately does not port, so no span can express it, and a differential checks it
instead. **The guard caught the new file itself** — its directory-listing test failed until the
file was declared, which is what that test exists for.

## The mutation ledger, before and after

**First sweep, on the committed step: 8 run, 8 killed.**

**The review's own sweep: 6 run, 4 survived.** Reading the seed instead of the fallback; skipping
the substitution on an empty histogram; reporting contig 0; and counting parked records instead of
records in the fit.

**Second sweep, after the fixes: 8 run, 7 killed, 1 survived.**

| # | the defect | result |
|---|---|---|
| R-A | the fallback reads the iteration's seed | killed by 4 |
| R-B | no substitution when the histogram is empty | killed by 2 |
| R-C | the mismatch reports contig 0 | killed by 1 |
| R-L | the warning counts what was parked | killed by 1 |
| N-1 | an unscored ratio becomes `0.0` before the verdict | killed by 4 |
| N-2 | a target above one is admitted | killed by 1 |
| N-3 | the outside-the-range count is always zero | killed by 1 |
| N-4 | the histogram is built over a different range than the one recorded | **SURVIVED** |

**N-4 is the tenth defect, and it was in a fix.** The shape and the histogram were two expressions a
few lines apart, and the test compared the shape against the constants rather than against the
histogram. Closed by making `LrHistogramShape::histogram` the only way pass two builds one, and by
having the test call that same constructor. Re-run afterwards:

| # | the defect | result |
|---|---|---|
| N-4b | the shape's own histogram uses a different range | killed by `the_recorded_shape_is_the_histogram_the_ratios_were_folded_into` |

**Nine mutations, nine killed.** Every file was restored afterwards and its sha256 compared with
the original, both sweeps.

## Deferred, with reasons

- **M6 — the consuming `RatiosInSpillOrder` cursor.** The reviewer is right that "in spill order"
  is the only thing tying a ratio to its record and that prose does not enforce it, and right that a
  cursor costs nothing: the review's open question 1 is answered — spec §3.5 walks the spill "per
  entry, in spill order, with its ratio from pass two", so pass three streams and never needs an
  index. **Deferred because spec §3.7's type block declares `ratios: Vec<f64>`**, and replacing it
  is a design change rather than a fix. Applied instead: a `debug_assert_eq!` against
  `entries_written()` and a field doc that says plainly what a later change would break. **C4 is
  where the cursor lands**, together with the consumer that makes it checkable.
- **Mi8 — the histogram's fixed ±100 range.** Measured, on one duplication-shaped record: 24.2 at
  one sample, 156.3 at six, **1,663 at 63**. The score is a sum over the samples weighed, so it
  grows with the cohort while the range does not. Harmless while only one class saturates — its
  posterior is saturated anyway — and not harmless if a cohort is large enough to push both classes
  past the same edge, where they share one bin and the target has nothing to move. **Not fixed,
  made countable**: `ratios_outside_the_histogram` on the verdicts. Whether the range should scale
  with the cohort is a spec question and is raised at Checkpoint C; D2 and D3 are the runs that
  would answer it.
- **Mi16 — renaming `ParalogVerdicts`.** It holds a calibration, a vector of ratios and three
  counts, and no verdict. Spec §3.7 gives it that name, so renaming is a spec change.

## What the numbers said, re-measured

The review's §6a found one wrong number and three claims that needed rewording. All corrected in
the implementation report, the `PROJECT_STATUS.md` entry and this branch's prose:

| claim | was | is |
|---|---|---|
| how many behavioural tests are blind to mutation 7 | eleven | **twelve** — 13 tests, one of which bypassed the pass |
| why the twenty-record fixture gave one cut at every target | "the classes separate completely" | the duplicated class's tail value **underflows to exactly zero**, and the other class sits at **0.75 — above every target tried**. Separation gives two values; it does not say where the second falls |
| the fallback's size | "three lines" in the code, "four" in the report and commit | **four**, everywhere |
| the fmt file counts before and after | "12 files before, 9 after" | reworded — a pre-commit working-tree state git does not hold |
