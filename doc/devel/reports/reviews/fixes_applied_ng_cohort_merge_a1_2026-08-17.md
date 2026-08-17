# Fixes applied — ng cohort merge, step A1

*2026-08-17, branch `ng-cohort-merge`. Input:
[the A1 review](ng_cohort_merge_a1_2026-08-17.md) — 2 Major, 9 Minor, 7 Nits over six
category checklists. Every finding is accounted for below.*

## Findings table

| ID | Title | Severity | Decision | Status |
|---|---|---|---|---|
| M1 | `get()` exercised only at its own default | Major | Apply | **Applied** |
| M2 | production's number over a rule that is not production's | Major | Apply (doc) | **Applied** |
| Mi1 | non-zero guard detached from the calls it protects | Minor | Apply | **Applied** |
| Mi2 | `the_region_width_is_not_the_locus_bound` cannot fail | Minor | Apply | **Applied** (deleted) |
| Mi3 | the `pub` constants were never pinned, only `Default` | Minor | Apply | **Applied** |
| Mi4 | bare `u32` defaults interchangeable below the newtypes | Minor | Apply | **Applied** |
| Mi5 | arch §1's "recorded in the run's output" clause dropped | Minor | Apply | **Applied** |
| Mi6 | `MinAltObs`'s doc did not name its default constant | Minor | Apply | **Applied** |
| Mi7 | `MinAltObs` abbreviates the term this module redefines | Minor | Apply (doc, not rename) | **Applied with adaptation** |
| Mi8 | two parent docs announced a stage not yet written | Minor | Apply | **Applied** |
| Mi9 | the `pub`/crate-private mismatch recorded nowhere | Minor | Apply (code half) | **Applied with adaptation** |
| Nit | `nonzero` named as an adjective | Nit | Apply | **Applied** |
| Nit | `get()` not `const fn` | Nit | Apply | **Applied** |
| Nit | line-pinned citations into production rot silently | Nit | Apply | **Applied** |
| Nit | threefold newtype boilerplate | Nit | Dispute | **Disputed** |
| Nit | derived `Ord`/`Hash` untested | Nit | Defer | **Deferred** to A4 |
| Nit | module name `run` needs an unmet two-word phrase | Nit | Dispute | **Disputed** |
| Nit | a `parameters.rs` would keep `mod.rs` a front door | Nit | Defer | **Deferred** |
| — | rename `CohortLocusBuilderRegionsLen` (plural) | Minor | Defer | **Deferred** to the owner |
| — | five design-document defects | — | Ask | **Raised at Checkpoint A** |

## What changed

All in [`src/ng/run/cohort_merge/mod.rs`](../../../../src/ng/run/cohort_merge/mod.rs)
unless named otherwise.

**M1 — the new test `get_returns_the_wrapped_value_not_the_default`.** It builds each
newtype from a value that is *not* its default (200, 7, 100) and reads it back, plus
`NonZeroU32::MIN` and `NonZeroU32::MAX`. This is the value class every parameter exists
for — all three are set from a command line — and it was the one class no test built.

**M2 and Mi5, Mi6, Mi7 — four doc corrections, no behaviour.** `DEFAULT_MIN_ALT_OBS`
now says that ng's keep rule sums every sample's non-reference reads where production
sums the per-position maximum across samples, that a maximum is never larger than a sum,
and that the performance claim it inherits was measured under the other rule.
`MaxCohortLocusSpan` regains the arch clause requiring its effective value to be recorded
in the run's output, marked **Owed** because the summary's surface is the emission step's.
`MinAltObs` now names its default constant, as the other two already did, and says
plainly that its `Obs` counts reads rather than this module's *observations*.

**Mi1 and Mi4 — one structural change that answers both.** Each newtype gained
`pub const DEFAULT: Self`, and `Default` returns it. A `const` item is evaluated when the
crate is compiled, so a zeroed default is now a build error *because the default is
used*, not because someone remembered to write a matching assertion. The three detached
`const _: () = assert!(… > 0)` lines are gone. The same associated constant is the typed
spelling a call site can name, so the three `u32` constants are no longer the only route
to a default.

**Mi2 — `the_region_width_is_not_the_locus_bound` deleted.** Its two operands were
constants the preceding test already pinned, so it could not fail on its own; and
`assert_ne!` demanded the two defaults differ, which no document requires — the region
width's own doc says a sweep will settle it, and a sweep landing on 50 would have turned
the test red with nothing wrong.

**Mi3 — `the_defaults_are_the_documented_values` now asserts on both spellings**: the
three `pub` constants a command line will advertise, and the three typed defaults a run
will use.

**Mi8 — [`src/ng/run/mod.rs`](../../../../src/ng/run/mod.rs) and
[`src/ng/mod.rs`](../../../../src/ng/mod.rs)** now say that what landed is the cohort
merge's parameters, not the stage.

**Mi9 — the module doc records why it is `pub`** and that the intent is the
architecture's crate-private one, to be narrowed when the caller objects land.

**Nits.** `nonzero` → `non_zero_default(default_value: u32)`; all three `get()` are
`const fn`; the production citations keep their symbol names and drop their line numbers.

## Applied with adaptation

- **Mi7** — the reviewer offered a rename (`MinAltReads`) or a doc clause. Took the
  clause: the name is fixed verbatim by the plan, the arch and the spec, and this loop
  does not edit design documents.
- **Mi9** — the reviewer's fix had two halves, one in the module doc and one in
  `doc/devel/ng/arch/cohort_merge.md`. Only the code half is applied; the arch half is
  the owner's, and is question 4 of the review's open questions.

## Disputed

- **Threefold newtype boilerplate.** `src/ng/types.rs` repeats the identical shape for a
  dozen newtypes, each carrying its own doc comment. A `macro_rules!` would hide exactly
  the prose that makes these types worth having.
- **The module name `run`.** `read` sets the precedent for a one-word domain noun in
  `src/ng/`, and the arch fixes `src/ng/run/` as the path.

## Deferred

- **Renaming `CohortLocusBuilderRegionsLen`** to the singular, and/or to spec §1.3's
  *building region*. Correct as criticism — the type is the only place that says
  *regions*, plural, for one region's width — but the name is written verbatim in the
  spec, the arch, the plan and a future flag, so the rename is the owner's call.
- **Testing the derived `Ord`/`Hash`** — to A4, where a locus span is first compared
  against `MaxCohortLocusSpan` and the comparison can carry its own test.
- **A `parameters.rs` split** — the arch's *Module home* puts the constants in `mod.rs`.

## Validation

Re-run in the container after every fix above; see the step's implementation report for
the verbatim output.
