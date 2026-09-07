# Fixes applied: ng_paralog_filter_a3

**Date:** 2026-09-06
**Review:** [ng_paralog_filter_a3_2026-09-06.md](ng_paralog_filter_a3_2026-09-06.md)
**Branch:** `ng-paralog-filter` · reviewed at `2f6a6c32`, fixes folded into that commit
**Outcome:** 1 Blocker, 6 Major, 7 Minor and 2 nits; **all applied**, none deferred.

## The answer

**The guard could be told to sanction an arbitrary edit, and now cannot.** Two of its three
checks tested only that the two declared lines *contained* certain substrings; everything else on
the line was free. A declaration reading

    ng_line: "use crate::ng::paralog::{GridSpec, SfsPriorSpec}; const SMUGGLED: f64 = 0.05;"

together with that same line in ng's copy left the whole guard green. **Every kind is now exact**
— ng's line must be production's with one named substring replaced and nothing else changed — and
each of the three has a test in its *rejecting* direction, which none had before.

**The differential grew four blind spots' worth of coverage**, and one of them turned out to be
unreachable by the route the review proposed:

| what was uncompared | how it is compared now |
|---|---|
| the verdict itself — `flags` and `posterior` | against production's, over four stream shapes × five targets × sixteen ratios |
| an empty histogram, and an EM that cannot converge | its own test; neither is reachable from a 5,000-ratio stream that settles in seven steps |
| both of `bin_index`'s clamps | **a narrow histogram** — see below |
| what the ratio streams contain | a fifth shape that reaches past the histogram's range |

**The clamps could not be caught the way the review suggested, and finding out why is the useful
part.** The review's diagnosis was that no drawn ratio left the histogram's `[-100, 100]`. Adding
a shape that does — measured, 1,251 of 5,000 below `-100` and 1,250 above `+100` — left both
mutations alive. The reason is that the shipped histogram's end bins sit at `±99.95`, where the
logistic that turns a ratio into a probability has saturated: `σ(-99.95 + logit π)` is exactly `0`
in `f64` and `σ(99.95 + logit π)` exactly `1`, so moving that mass one bin changes nothing π or
the curve can see. **A bin-index error at either end is invisible through the shipped range,
however far outside it the ratios reach.** A histogram over `[-6, 6]` puts its end bins where the
logistic is still moving; both mutations then die. That case is also the only place the
differential runs at a model setting other than the shipped one.

## Findings table

| ID | Severity | Decision | Status |
|---|---|---|---|
| B1 — two sanction checks were containment tests, with no rejecting-direction test | Blocker | Apply | **Applied** |
| M1 — the ratio streams never left the range, ran dry, or failed to converge | Major | Apply | **Applied with adaptation** |
| M2 — the differential stopped before the verdict | Major | Apply | **Applied** |
| M3 — `flags`'s oracle was its own body, at one target | Major | Apply | **Applied** |
| M4 — `lr_threshold` unread, and its documented equivalence false | Major | Apply | **Applied with adaptation** |
| M5 — the published recipe cannot see an `include_str!` | Major | Apply | **Applied with adaptation** |
| M6 — the guard's header stated the pre-A3 rules | Major | Apply | **Applied** |
| Mi1 — the end normalisation swallowed a CRLF on the last content line | Minor | Apply | **Applied** |
| Mi2 — `ends_before` had no "occurs exactly once" check | Minor | Apply | **Applied** |
| Mi3 — three stale counts in doc comments | Minor | Apply | **Applied** |
| Mi4 — `src/ng/mod.rs`'s inventory under-reported the module | Minor | Apply | **Applied** |
| Mi5 — the mirror visibility rule lived only in a test module | Minor | Apply | **Applied** |
| Mi6 — what is deliberately not copied was not in the guard's own tables | Minor | Apply | **Applied** |
| Mi7 — `CalibrationConfig` and the fallback constant untested in ng | Minor | Apply | **Applied** |
| Nits (2) | Nit | 1 Apply, 1 informational | **Applied** |

## What was done

### B1 — every sanction is exact, and every one is tested in its rejecting direction · Applied

`PathIntoProduction` now requires `production_line.replace("crate::paralog", "crate::ng::paralog")
== ng_line`. `DocLinkToAnItemNgDoesNotCopy` carries the link it removes and the plain words that
replace it, and requires the same one-substring equality. The visibility check was already exact.
A failed sanction returns through the enclosing `Err(String)` rather than an `assert!`, so it
carries the `THE GUARD COULD NOT RUN` banner and says the declaration is at fault, not ng's copy.

`a_declaration_that_smuggles_an_edit_is_refused` drives all three kinds both ways, building ng's
side *from the declaration itself* so that only the declaration can be what refuses it.

**Measured, after:** the review's own mutation — the declaration's `ng_line` changing a constant
while production's stays right — is refused with the reason named.

### M1, and why the clamps needed something else · Applied with adaptation

A fifth stream shape, `BeyondTheHistogramsRange`, and a narrow-histogram case. The review proposed
pinning the ratio streams' contents the way the scorer's stream is pinned; that was not applied,
because the narrow-histogram case makes the property it would protect observable directly — every
bin index is read at least once there, past both ends — and a second table of pinned counts is a
second thing to maintain for the same guarantee.

### M2 — the verdict compared with production's · Applied

`the_copied_verdict_agrees_with_productions_bit_for_bit` names
`crate::var_calling::paralog_filter::calibrate::ParalogCalibration` and compares `flags` and
`posterior` at sixteen ratios — the saturating ends and the three unscorable values among them —
for five targets across four stream shapes.

### M3, M4 — `flags` and the recorded cut · Applied with adaptation

The target is now swept over five values rather than fixed, which kills an implementation that
reads a literal `0.01` instead of the operator's knob; and a target set exactly to a q-value the
curve produces separates `<=` from `<`. **Both mutations were run and both now fail.**

For M4 the review proposed replacing the false equivalence with a "within one bin" relation, and
that is what was done — but stated as what it costs rather than as a tolerance:
`the_recorded_cut_and_the_flag_agree_to_within_one_bin` asserts that the two disagree **only**
inside the crossing bin, that they *do* disagree there (so the fixture cannot silently stop
exercising it), and that the recorded cut is itself flagged. The doc comment says why it matters:
`lr_threshold` is what goes in the VCF header, and an operator reading *"records at or above this
ratio were dropped"* is wrong for every record inside a half-bin — `0.05` on the likelihood-ratio
axis at the shipped resolution.

### M5 — the recipe, corrected but not as proposed · Applied with adaptation

The recipe in `src/ng/mod.rs` now matches compile-time file inclusions as well as `use` lines. The
review's wording said the new `include_str!` "is the only place `src/ng/` names `src/var_calling/`
at all". **That is false and was not written down**: `calling/genotype_table_parity.rs`,
`loop_parity.rs`, `quality_parity.rs` and `calling/likelihood/generic.rs` all name it through
ordinary `use` lines. What is true, and what the correction says, is that the two copy guards read
production's source as *text* and a `use`-only sweep sees neither.

### M6, Mi1–Mi7 · Applied

The guard's header gains a section on where a span copy ends and what the one normalisation is,
and its contract section now names all three sanctioned kinds and says each is exact. The
normalisation drops **whole blank lines only**, so a trailing space or a `\r` on the last content
line still fails — extracted into `without_trailing_blank_lines` with the rule stated once, in
place of the inline trim that had grown unreadable. `ends_before` must occur exactly once past the
copy's start. Three stale counts corrected. `src/ng/mod.rs` names all four statistics and states
the mirror visibility rule beside the existing one, with both subjects named. `copy_fidelity.rs`
gains a *"What is deliberately not copied"* section and
`the_item_the_span_stops_before_is_still_productions`, which asserts production still declares
`cohort_inbreeding`.

**One naming error of the author's was found in the process.** The commit message, the report and
three doc comments called it `CohortInbreeding`, a type. It is `cohort_inbreeding`, a function
that fills a slice with one coefficient repeated per sample — and the first version of the new
test asserted `struct CohortInbreeding` and failed, which is how it surfaced. Corrected
everywhere.

## Validation

Run in the dev container, on the fixed tree.

    cargo test --all-features --lib "ng::paralog::"
    → test result: ok. 83 passed; 0 failed; 0 ignored; 6290 filtered out

    cargo test --all-features --lib --bins --tests
    → 6358 lib tests passed, 0 failed, 15 ignored
    → one integration test failed: a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele

    cargo clippy --lib --bins --tests --all-features -- -D warnings
    → 3 errors, all needless_lifetimes in cohort_merge/build.rs and serial.rs

    cargo fmt --check
    → dirty on 9 files, none of them src/ng/paralog/

The integration failure, the three lints and the nine files are `main`'s at `a33ada0f`. The gate —
no failure `main` does not already have — holds, and the lib count moves 6,352 → 6,358.

**Five clippy errors of ng's own appeared during this work and were fixed rather than silenced**:
`enum_variant_names` on `WhyRepointed`, and four `clone_on_copy` on `ParalogPrior`.

Every mutation quoted above was applied to a file backed up first, and each restoration was
confirmed by `diff` before the next command ran.
