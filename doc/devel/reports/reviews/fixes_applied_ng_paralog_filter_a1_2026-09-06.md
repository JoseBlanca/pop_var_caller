# Fixes applied: ng_paralog_filter_a1

**Date:** 2026-09-06
**Review:** [ng_paralog_filter_a1_2026-09-06.md](ng_paralog_filter_a1_2026-09-06.md)
**Branch:** `ng-paralog-filter` · reviewed at `1823bba1`, fixes folded into that commit
**Outcome:** 4 Major and 8 Minor findings; **10 applied, 2 deferred**, all nits applied.

## The answer

**The guard now catches every divergence the review proved it missed, and each is proved
caught by re-running the mutation that used to pass.** Four mutations were run in the working
tree, each with the file restored from a byte-compared backup afterwards:

| what was changed | before the fixes | after |
|---|---|---|
| a doc line in the copied constants (`/// \`EPS\`.` → `/// \`EPS2\`.`) | all 6,303 lib tests passed | fails, naming `model_params.rs, line 5 of the copied content` |
| a stray `src/ng/paralog/likelihood.rs` added | nothing noticed | fails: `likelihood.rs … is neither guarded, ng's own, nor released` |
| the guarded file list emptied | the guard printed `ok` having compared nothing | fails, naming `coverage_model.rs` as unclassified |
| the copy's final newline stripped | passed | fails: `a byte-level difference the line comparison cannot see … 10392 bytes to production's 10393` |

The suite goes from 28 tests to **31**, and the whole scoped suite from 6,303 passing to
**6,306** — the three additions being ng's own finiteness test, the comparison's own fixtures,
and the directory-coverage check.

## Findings table

| ID | Severity | Decision | Status |
|---|---|---|---|
| M1 — "the swap is one `use` line" is wrong | Major | Apply (prose) + Defer (the design half) | **Applied with adaptation** |
| M2 — the guard's file list could go empty unnoticed | Major | Apply | **Applied** |
| M3 — 236 copied lines in the one file the guard exempted | Major | Apply | **Applied with adaptation** |
| M4 — the guard's comparison had no test | Major | Apply | **Applied** |
| Mi1 — `crate::ng::run::paralog_filter` in the present tense | Minor | Apply | **Applied** |
| Mi2 — "every file beside this one is production's" is false | Minor | Apply | **Applied** |
| Mi3 — nothing makes the temporary import temporary | Minor | Defer | **Deferred** |
| Mi4 — the note cites a spec file not in the tree | Minor | Apply | **Applied** |
| Mi5 — a 3-tuple named `pairs` | Minor | Apply | **Applied** |
| Mi6 — failure messages name no location, and one misdirects | Minor | Apply | **Applied** |
| Mi7 — "byte for byte" promised, line-for-line checked | Minor | Apply | **Applied with adaptation** |
| Mi8 — the `// non-finite` case does not test finiteness | Minor | Apply | **Applied with adaptation** |
| Nits (6 groups) | Nit | Apply | **Applied** |

## What was done, finding by finding

### M1 — the histogram swap costs four lines, not one · Applied with adaptation

The review recommended settling the field set in prose *and* editing
`doc/devel/ng/spec/window_coverage.md` §7. **The prose half was applied; the spec edit was
not, and deliberately.** That document belongs to the window-coverage plan, being built in a
parallel session on its own branch, and the plan-driven skill does not let this loop edit
another plan's design.

[coverage_model.rs](../../../../src/ng/paralog/coverage_model.rs)'s appended note now says
what is true: the fit reads five of the histogram's fields and ng's copy carries all five,
**but the transcribed test fixture builds it with all eight**, and ng's copy drops
`n_skipped_tiles` and `callable_positions` and renames `n_positions` to `windows_folded`. So
the swap is the `use` line plus three lines of that fixture — inside a file the guard forbids
editing, which is why it must land as its own commit that releases the file.

**Verified against the other branch's code, not its spec.** `ng-window-coverage` committed
`09bad846` while this step was under review; its `CoverageByGcHistogram` carries `window_bp`,
`gc_bins`, `depth_bin_width`, `depth_bins`, `windows_folded`, `counts`.

**The design half is deferred to Checkpoint A**: `window_coverage.md` §7's instruction to drop
the fields the fit does not read cannot stand alongside a verbatim port of production's tests,
and which way that resolves is the other spec's call.

### M2 — the guard can no longer go quietly empty · Applied

`every_file_in_the_module_is_guarded_ng_s_own_or_released` enumerates `src/ng/paralog/` at test
time and requires every `.rs` file to appear in one of three lists — guarded, ng's own, or
explicitly released — and requires the guarded list and the file list to name the same set.
The release protocol *deletes* an entry, so a dropped entry and a deliberate release look
identical; the directory listing is now the authority instead.

### M3 — the copied constants moved into a file the guard can reach · Applied with adaptation

Two fixes were proposed: extract the block to its own file (module_structure) or guard the
region by substring inside `mod.rs` (reliability). **The extraction was chosen.** It makes the
header's claim — "this file is ng's own" — true as written rather than true of part of the
file; it uses the guard machinery already there instead of a second, weaker one; and it makes
the directory uniform for A2 and A3, which each bring one more copied file.

[model_params.rs](../../../../src/ng/paralog/model_params.rs) is production's `src/paralog/mod.rs`
from its first item to its end, 236 lines, under an ng header.
[mod.rs](../../../../src/ng/paralog/mod.rs) keeps only what ng wrote.

Production keeps those items in a `mod.rs` whose declarations ng cannot share, so this copy
cannot be compared from the end of a shared module header. `CopyBegins` names the two shapes:
`AfterTheModuleHeader` for a whole-file copy, `AtTheFirstItemDoc` for this one.

### M4 — the comparison is now driven by fixtures · Applied

`how_the_copy_differs` is extracted from the `#[test]` and returns `Option<String>` — the
message a failing guard prints — so
`the_comparison_accepts_an_appended_note_and_rejects_every_other_edit` can drive every
rejecting branch with synthetic file pairs: a rewritten production header line, a note
prepended rather than appended, a dropped header line, a drifted body line, a deleted one, an
added one, a stripped final newline, and a CRLF ending on a body line and on a header line.

### Mi6 — the messages · Applied

The per-line loop now runs **before** the length check, so a deletion is reported with the line
where the two sides part instead of a bare count. Both sides are labelled `ng's copy:` and
`src/paralog/:`. The remedy is phrased symmetrically — *"whichever side moved, revert it"* —
because the guard fires just as correctly when production is what was edited, and the old
message told the reader to release ng's copy in that case. The header check names the
module-header line it failed on.

### Mi7 — byte for byte, now actually · Applied with adaptation

The review offered the cheap fix (weaken the prose to "line for line") or the strong one (add
a byte comparison). **The strong one was taken**, because fidelity is what the commit is for.
The copied content is now compared line for line *and* as bytes, and the module-header prefix
is compared with `split_inclusive('\n')` rather than `lines()`.

**The review's fixture for the CRLF case put the `\r` on a header line, and that case failed
when first run** — the header was still being compared with `lines()`, so nothing saw it, and
the body comparison does not reach the header. Fixed by comparing the header prefix by bytes
too; both the header-line and body-line CRLF cases are now pinned.

### Mi8 — ng's own test, beside production's · Applied with adaptation

The finding is real — `NAN < 0.6` is `false`, so the transcribed `// non-finite` case is
refused by the `lo < hi` clause with or without `is_finite`. But the test is production's and
lives in a file that may not gain a line. `grid_spec_new_refuses_an_infinite_endpoint` is
therefore ng's own, in `mod.rs`, asserting the three inputs that separate the two: an infinite
`hi`, a negative-infinite `lo`, and a `NaN` in `hi` rather than `lo`.

### Mi1, Mi2, Mi4, Mi5 and the nits · Applied

The header now says the wiring *will* live in `crate::ng::run::paralog_filter` at plan step B1
and has not landed; names which files are copies and which is ng's own, with their sizes;
names the `ng-window-coverage` branch wherever it cites something only that branch has; and
opens by saying what the four quantities *are* in plain terms rather than stacking their names.
`pairs` is now `guarded_copies`, a `Vec<GuardedCopy>` with named fields — the transposition the
tuple allowed cannot be written. The five identifier slips are renamed, the `#[cfg(test)] mod
copy_fidelity;` declaration carries the "ng's, not a copy" marker its siblings carry, the
`pub use` block says why it reproduces production's surface name for name, the explicit array
type is gone, and the header-length closure is a function.

## What was deferred, and why

**Mi3 — nothing makes the temporary import temporary.** The review's fix is to add the swap as
a step to `doc/devel/ng/impl_plan/hidden_paralog_filter.md` with the window-coverage plan's
first step as its precondition. The plan-driven-implementation skill explicitly forbids this
loop from editing the plan's design — only its step checkboxes — so adding a step is a
stop-and-ask, not a fix. It is recorded in `PROJECT_STATUS.md`'s open items for this feature
and raised at Checkpoint A. The alternative the review offers, copying the eight-field struct
into `src/ng/window_coverage/` now, would write into a directory another branch owns.

**The nit on `grid_specs_have_expected_defaults`** — a sentence saying which case it adds over
its neighbour, so a later tidy does not delete it as duplication — could not be applied: the
test is production's, inside `model_params.rs`, which is guarded byte for byte. It is recorded
here instead.

## Validation

Run in the dev container, on the fixed tree.

    cargo test --all-features --lib "ng::paralog::"
    → test result: ok. 31 passed; 0 failed; 0 ignored; 6290 filtered out

    cargo test --all-features --lib --bins --tests
    → 6306 lib tests passed, 0 failed, 15 ignored
    → one integration test failed: a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele

    cargo clippy --lib --bins --tests --all-features -- -D warnings
    → 3 errors, all needless_lifetimes in cohort_merge/build.rs and serial.rs

    cargo fmt --check
    → dirty on 9 files, none of them src/ng/paralog/

The integration failure, the three lints and the nine files are `main`'s at `a33ada0f` and were
proved so before this step began: with `pub mod paralog;` removed, the same suite gives 6,275
passing lib tests and the same single failure. The gate — no failure `main` does not already
have — holds, and the lib count moves 6,275 → 6,306.

Every mutation run above was reverted from a backup and the restoration confirmed by `diff`
before the next command; `src/ng/paralog/` ends with exactly four files.
