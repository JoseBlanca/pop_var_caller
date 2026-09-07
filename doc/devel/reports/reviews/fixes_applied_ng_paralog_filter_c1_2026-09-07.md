# Fix Application Report: ng_paralog_filter_c1_2026-09-07.md

**Date:** 2026-09-07
**Source review:** [ng_paralog_filter_c1_2026-09-07.md](ng_paralog_filter_c1_2026-09-07.md)
**Source state reviewed against:** `0ceb463b` on `ng-paralog-filter`; **fixed forward in its own
commit**, not folded into that one — unlike B1, B2 and B3, whose reviews ran before anything else
landed on top. C1 had already been reported as committed and a `PROJECT_STATUS` commit sits above
it covering more than this step, so amending would have misattributed that and tidied a record the
project's rule says to leave alone.
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 4
- Minors: 13
- Nits: 2

### Outcome totals
- Applied: 16
- Applied with adaptation: 2
- Deferred: 2
- Disputed: 0
- Already fixed: 0
- Failed validation: 0

### Validation summary

Run in the dev container from this worktree, on the tree being committed.

- `cargo test --all-features --lib "ng::run::paralog_filter::scoring_context"` → `ok. 19 passed; 0 failed`
- `cargo test --all-features --lib "ng::run::paralog_filter"` → `ok. 94 passed; 0 failed`
- `cargo test --all-features --lib --bins --tests` → **6,563 lib tests passed, 0 failed, 15
  ignored**; one integration test failed
  (`a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele`), which is `main`'s
- `cargo clippy --lib --bins --tests --all-features -- -D warnings` → 9 errors, **none in
  `src/ng/run/paralog_filter/`** — unchanged from the merge base
- `cargo fmt --check` → dirty on the same 9 files as the merge base, none of them this step's
- Performance check → not applicable: no `Apply` touched code reachable from a bench harness, and
  the review measured the one hot-path question directly (below)

**The lib count moves 6,554 → 6,563**: nine tests added to `scoring_context`.

### Unresolved high-priority findings

None. The Blocker and all four Majors are Applied or Applied with adaptation.

## 2. Findings table

| ID | Severity | Title | Final status | Files changed |
|---|---|---|---|---|
| B1 | Blocker | the fixture has one GC bin, so a wrong GC passes all ten tests | **Applied** | `scoring_context/tests.rs` |
| M1 | Major | neither half of the finiteness guard has a test that fails when deleted | **Applied** | `scoring_context.rs`, its tests |
| M2 | Major | a record narrower than the cohort is scored on a prefix | **Applied** | `scoring_context.rs`, `mod.rs`, its tests |
| M3 | Major | a bad fit configuration is reported once per sample as a data verdict | **Applied with adaptation** | `scoring_context.rs`, `mod.rs`, its tests |
| M4 | Major | the sample's row is read twice and by field, so a new field is dropped silently | **Applied** | `scoring_context.rs` |
| Mi1 | Minor | four reasons for `None`, one of them corruption | **Applied with adaptation** | `scoring_context.rs` |
| Mi2 | Minor | the guard diverges from `is_absent` and nothing says so | **Applied** | `scoring_context.rs`, its tests |
| Mi3 | Minor | the tolerance is 0.1 where the fixture is exact | **Applied** | `scoring_context/tests.rs` |
| Mi4 | Minor | two parallel `Option` vectors spell one sum type | **Deferred** | — |
| Mi5 | Minor | `new` takes `&[f64]` where the parameters file exposes `&[InbreedingF]` | **Applied** | `scoring_context.rs`, its tests |
| Mi6 | Minor | `WhyNoCoverageModel` has no `Display` | **Applied** | `scoring_context.rs` |
| Mi7 | Minor | `samples_with_a_coverage_model` reads as a collection | **Applied** | `scoring_context.rs`, its tests |
| Mi8 | Minor | `precompute()` is a verb naming an accessor | **Applied** | `scoring_context.rs` |
| Mi9 | Minor | the by-value `Vec` is correct and unexplained | **Applied** | `scoring_context.rs` |
| Mi10 | Minor | the fit configuration is dropped, so it cannot be reported | **Deferred** | — |
| Mi11 | Minor | five `pub` accessors; production keeps three private behind one `score` | **Deferred** | — |
| Mi12 | Minor | `mod.rs` still lists the scoring context among "later steps" | **Applied** | `mod.rs` |
| Mi13 | Minor | the report defers a cost that is now measured | **Applied** | the impl report |
| §6a | — | three wrong counts and one wrong list | **Applied** | four files |
| §7 | — | a third `SpilledSamples` variant would pass as "not a tract" in `spill.rs` | **Applied with adaptation** | `spill.rs` |

## 3. Questions asked and answers

None asked. The review's two open questions were both answerable from the code and are resolved in
§4 rather than sent to the owner: a half-absent window pair *is* representable (the owning type's
own doc says so), and the record-level entry point *does* get its own error type.

## 4. Adaptations — where the fix differed from the review's suggestion

- **M3** — the review proposed validating the configuration once before the loop. `validate_config`
  is private to `coverage_model.rs`, which is a byte-for-byte copy of production behind a textual
  guard, so it cannot be called from outside and must not be edited to make it public. **Detected
  instead**: the first sample whose fit fails with `InvalidConfig` short-circuits and `new` returns
  `Err(CoverageFitConfigRefused)`. Equivalent outcome, and it touches nothing frozen. The
  consequence is that `new` now returns a `Result` — which is the right shape anyway, and free,
  since C2 is unwritten.
- **Mi1** — the review asked for either a counted total of unsummable read pairs or a doc sentence.
  **Took the sentence**, at the point the `checked_add` happens, and the record-level entry point
  removes the index-past-the-cohort case entirely by turning it into an error once per record. A
  counter with no reader would be surface for a step that does not exist; C3 can add one when it
  has somewhere to report it.
- **§7** — the writer's `matches!` became a `match`, so a third row shape cannot pass as "not a
  tract". **The decoder cannot be fixed the same way and is not**: the file stores one bit, which
  distinguishes two shapes and no more, so a third variant needs a wider discriminant in the format
  itself. That is now a marked comment at the branch rather than a silent constraint.

## 5. Deferred, with reasons

- **Mi4 — collapsing `coverage_models` and `why_no_model` into one `Vec<Result<…>>`.** Real, and
  the loop already builds exactly that `Result` before taking it apart. Deferred because the review
  itself found the weaker half harmless — there is no `pub` field and no `&mut self` method, so all
  four vectors are filled in one function and cannot drift afterwards — and because `precompute`'s
  argument and `single_copy_depth_sd` both want contiguous slices, so the change is less local than
  it looks. Worth doing when C3 shows which shape it actually reads.
- **Mi10 — keeping the fit configuration so the run report can state it.** Agreed in principle:
  spec §3.1 flags those constants as prototype-tuned and never re-measured, which is exactly what a
  run report should say. Deferred to C4, which is the step that writes the report and will know
  what it wants to print.
- **Mi11 — narrowing the public accessors behind one `score` method, as production does.** This is
  a recommendation about C3's shape, filed at Medium confidence by the reviewer, and C3 is where it
  can be judged against real code. Recorded so it is not lost.

## 6. What the numbers said, re-measured

The review's §6a found three wrong counts and one wrong list, all mine, all corrected:

| claim | was | is |
|---|---|---|
| the guarded span of `coverage_model.rs` | 694 lines, in two places | **644** — the header is 54 `//!` lines and `#[cfg(test)]` is at line 700, less the trailing blank |
| the same file's length past its header, in two module docs | 1,158 | **1,157** |
| the module header's list of the fit's refusal reasons | three, and omitting the one the file's own test asserts on | **all six**, named |
| `ONE_COPY_DEPTH` is one copy "within 0.1" | true but weak — it is exactly 1.0 | tolerance tightened to `1e-9`, with the reason the fixture is exact |

**And one measurement replaced a deferral.** The implementation report deferred the duplicate row
lookup's cost to step D2. The review measured it: `observation_of` costs **1.82 ns a sample**
against **2,828 ns** to score the same sample, about **1,550×**, so it could never have appeared in
D2's pass shares. The lookup is gone anyway — but because of M2's missing length check, not the
cost, and the report now says so.

## 7. The mutation ledger, before and after

The review ran five mutations against the ten committed tests: **three survived**. All three are
now killed, verified on the tree being committed, each restored from a byte-compared backup:

| mutation | before | after |
|---|---|---|
| the GC argument replaced by the literal `0.9` | survived all 10 | **1 of 19 fails** |
| the GC half of the finiteness guard deleted | survived all 10 | **2 of 19 fail** |
| the depth half deleted | survived all 10 | **2 of 19 fail** |
| the tract arm returns `None` (trap 1) | killed | killed, 1 of 19 |
| σ₀ pushed as `0.0` rather than `NaN` | killed | killed |

## 8. Commands to re-verify

    ./scripts/dev.sh cargo test --all-features --lib "ng::run::paralog_filter"
    ./scripts/dev.sh cargo test --all-features --lib --bins --tests
    ./scripts/dev.sh cargo clippy --lib --bins --tests --all-features -- -D warnings
    ./scripts/dev.sh cargo fmt --check
