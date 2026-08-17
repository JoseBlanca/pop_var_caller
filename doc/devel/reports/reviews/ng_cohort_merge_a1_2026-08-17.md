# Review — ng cohort merge, step A1 (the three parameters)

*2026-08-17, branch `ng-cohort-merge`, working-tree diff (uncommitted at review time).
Six category checklists, five sub-agents. Per-category audit trail:
`tmp/review_2026-08-17_ng-cohort-merge-a1/`.*

## 1. Scope

- **Reviewed:** the working-tree diff of plan step A1 — the module home
  `src/ng/run/cohort_merge/` and the three parameters a calling run sets.
- **In-scope files:** [`src/ng/run/cohort_merge/mod.rs`](../../../../src/ng/run/cohort_merge/mod.rs)
  (new), [`src/ng/run/mod.rs`](../../../../src/ng/run/mod.rs) (new),
  [`src/ng/mod.rs`](../../../../src/ng/mod.rs) (one `pub mod` line and one doc clause).
- **Out of scope:** the rest of the crate; the pre-existing lint state of `examples/`,
  `benches/` and other modules' test code.
- **Categories dispatched:** `naming`, `defaults` (this step *is* configuration),
  `module_structure` (a new module tree), `reliability` (the tests), `idiomatic`,
  `smells`. **Not dispatched:** `errors` (the step adds no error path), `tooling` (no
  `Cargo.toml` change), `unsafe_concurrency` (no `unsafe`, no threads, no shared state),
  `extras` (no parser, no untrusted input, no hot path), `refactor_safety` (folded into
  `smells`, which covered the line-pinned citations).

**One deviation from the review skill, recorded:** the sub-agents reviewed by reading in
the author's checkout rather than each taking a git worktree. A1 adds no executable
behaviour beyond `Default` and `get()`, so mutation testing had almost nothing to
mutate, and a worktree apiece would have cost five full container builds to prove it.
The agents were told to reason about mutants instead and to touch nothing; the
`reliability` agent reports **0 mutations applied, 0 run, 0 survivors observed**, and
labels its mutant claims as reasoned rather than measured. A3 and A4 have real
behaviour and will get isolated worktrees.

## 2. Verdict

**Approve-with-changes.** No Blocker. Two Major, both real, both cheap: one gap in the
tests, one claim in a doc comment that does not survive checking.

## 3. Execution status

Run in the container (`./scripts/dev.sh`), verbatim results:

| command | result |
|---|---|
| `cargo fmt --check` | clean, exit 0 (after a `cargo fmt` pass over three files the gate already flagged on this branch's base) |
| `cargo clippy --lib --all-features -- -D warnings` | clean — `Finished \`dev\` profile` |
| `cargo test --lib` | `test result: ok. 3613 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out; finished in 577.64s` |
| `cargo clippy --all-targets --all-features -- -D warnings` | **fails, 49 errors** — every one pre-existing on this branch's base, in `examples/`, `benches/`, and the test code of `census_file.rs`, `ssr_fit.rs`, `open_bam.rs` and `src/ssr/cohort/sim.rs` (which ng may not edit). None in the in-scope files. Confirmed pre-existing by stashing the step's changes and re-running. |

Findings labelled "needs verification": 0.

## 4. Open questions — for the owner, not for this step

Five of these are defects in the **design documents**, which this loop may not edit.
Two of them reach A4, the next step but one.

1. **Arch §6's keep-rule bullet reads two ways, and one of them makes `MinAltObs`
   inert.** It says the keep rule is *exact "any non-reference observation in any
   sample"*, which as a statement about the threshold means "keep on one read" — i.e.
   `min_alt_obs = 1`, and the default of 2 is then wrong. Spec §4.3 forecloses that
   reading explicitly ("Not 'any non-reference read at all'"), and the code follows the
   spec. **A4 is coded from that bullet**, so it wants one line of correction.
2. **Spec §6.1 says a region is "about a hundred" bases while §6.4 and the arch say
   20** — residue of the retired "twice `max_cohort_locus_span`" derivation. Which is
   the intended default?
3. **Arch §1 declares `DEFAULT_MAX_COHORT_LOCUS_SPAN` twice**, in two adjacent doc
   comments. Harmless, but it is the kind of thing a reader trips on.
4. **The arch says this module is crate-private; the code is `pub`.** The code now
   records why (no caller objects yet, so `pub(crate)` is dead code under `-D warnings`;
   ng's probes live in `examples/`, a separate crate target). The arch sentence still
   says the opposite.
5. **Two names would read better but are fixed verbatim by the plan and arch**:
   `CohortLocusBuilderRegionsLen` is plural for one region's width, and `MinAltObs`
   abbreviates *observation* for a quantity counted in *reads* — the word this module
   redefines. Renaming either touches spec, arch, plan and a future flag name, so it is
   the owner's call; the code carries a disambiguating doc clause instead.

## 5. Top 3 priorities

1. **M1** — nothing tests `get()` on a value that is not the default, so an accessor
   ignoring its argument passes the suite while silently running every cohort at
   50/2/20 whatever the operator set.
2. **M2** — the `DEFAULT_MIN_ALT_OBS` doc borrows production's measured justification,
   but ng's keep rule is not production's, and the two diverge more as the cohort grows.
3. **Mi1** — the "a zero default is a build error" guarantee was maintained by hand, in
   three assertions detached from the three calls they protect.

## 6. Findings

### Major

**M1: `src/ng/run/cohort_merge/mod.rs` — `get()` is exercised at exactly one input per
type, its own default.**
**Categories:** reliability. **Confidence:** High.
Every `get()` call in the crate went through `Default`, so replacing
`MaxCohortLocusSpan::get`'s body with the literal `DEFAULT_MAX_COHORT_LOCUS_SPAN` passed
the suite unchanged. That is not a no-op mutant: the tuple field is `pub`, so
`MaxCohortLocusSpan(NonZeroU32::new(200).unwrap()).get()` is constructible today and
returns 200 under the real code, 50 under the mutant. The untested value class is the
operator-set one — the case all three doc comments call "a command-line parameter", and
the one spec §3.1 says a long-read run will use.
**Fix:** the new test `get_returns_the_wrapped_value_not_the_default`, which also covers
both ends of `NonZeroU32`.

**M2: `src/ng/run/cohort_merge/mod.rs` — production's number, carried over a rule that
is not production's.**
**Categories:** defaults. **Confidence:** High on the mechanism.
The constant's doc said "production's value, carried over", and `MinAltObs` imported
production's measured justification with it. Production's `derive_is_kept` sums, across
a group's positions, the *maximum over samples* of that sample's non-reference
observations; ng sums every sample's non-reference reads (spec §4.3, and spec §15 pins
the difference with a test that inverts production's). A maximum is never larger than a
sum, so at the same threshold ng keeps everything production keeps and more —
identically at one sample, by a widening margin as the cohort grows. At 63 samples, one
non-reference read in each of two samples at one position reaches ng's 2 and never
reaches production's.
**Not a design defect** — spec §4.3 chose this deliberately. It is the doc comment that
implied an inherited measurement transfers at full strength.
**Fix:** the doc now states the rule difference and marks the performance claim as
inherited evidence rather than evidence about ng.

### Minor

**Mi1: the compile-time non-zero guard was detached from the calls it protects.**
**Categories:** idiomatic, reliability, defaults — convergent, all three found it
independently. **Confidence:** High.
`nonzero` was `const fn`, but its three callers were `fn default()` bodies — ordinary
functions — so the call was a run-time call and its `panic!` arm a run-time path. What
actually made a zero a build error was three separate `const _: () = assert!(…)` lines,
each paired by hand with one constant. Arch §1 has a fourth parameter of the same shape
still to land (`ObservationReachCeiling`); one added without its assertion would compile
clean and panic on first use, with the doc still claiming otherwise.
**Fix:** each newtype gained `pub const DEFAULT: Self`, which is const-evaluated by
definition, and `Default` returns it. The guarantee is now made by the same line that
needs it. The three detached assertions are gone.

**Mi2: `the_region_width_is_not_the_locus_bound` could not fail, and would have failed
on a legitimate retune.**
**Categories:** reliability, smells, defaults — convergent. **Confidence:** High.
Both operands were `const u32` literals already pinned by the preceding test, so nothing
that left that test green could turn this one red. Worse, `assert_ne!` requires the two
defaults to *differ*, which no document demands: both are documented in the same file as
unmeasured, and a sweep landing the region width on 50 would have turned the test red
with nothing wrong. The property it claimed to pin — that the width is not *derived*
from the bound — is structural, and the swap it worried about is already prevented by
the two being distinct types.
**Fix:** deleted. The argument lives on `CohortLocusBuilderRegionsLen`'s doc, where it
already was.

**Mi3: the test pinned the `Default` impls but never the `pub` constants.**
**Categories:** reliability. **Confidence:** High.
A `clap` `default_value_t` reads the constant; a run without the flag reads `Default`.
Two advertised defaults per parameter, only one of them checked — so an impl
open-coding a literal while its constant was retuned would print one number and run
another, with nothing failing.
**Fix:** three constant assertions added to `the_defaults_are_the_documented_values`.

**Mi4: the three bare `u32` defaults are interchangeable one layer below the newtypes.**
**Categories:** naming, defaults — convergent. **Confidence:** Medium.
The newtypes stop the two spans being swapped, but their `u32` defaults are freely
transposable: `MaxCohortLocusSpan(nonzero(DEFAULT_MIN_ALT_OBS))` compiled. The
command-line layer is where all three will first be handled together.
**Fix:** the typed `pub const DEFAULT: Self` of Mi1 gives each default a spelling that
is both explicit and type-checked. The readable `u32` constants stay, for help text.

**Mi5: arch §1's clause requiring the effective bound to be recorded in the run's output
was dropped from an otherwise faithful copy.**
**Categories:** defaults. **Confidence:** High.
It is the one place the design says a resolved default must be observable from a
finished run, and it was written down nowhere in the code.
**Fix:** restored on `MaxCohortLocusSpan`, marked **Owed** — the summary's surface is
the emission step's (spec §13), so nothing here writes it yet.

**Mi6: `MinAltObs` was the one parameter whose doc did not name its default constant**,
where the other two did. **Categories:** defaults. **Fix:** applied.

**Mi7: `MinAltObs` abbreviates the one term this module redefines.** The module's first
sentence is about *observations*, which spec §1.3 fixes as one sample's whole record
over a stretch of genome; `MinAltObs` counts *reads*, and one observation carries many.
Misreading it gives a different number, not a vaguer one.
**Categories:** naming. **Fix:** the reviewer's cheaper option — a disambiguating clause
on the type. The rename is open question 5 above.

**Mi8: two parent module docs announced a stage that has not been written.**
`src/ng/run/mod.rs` and `src/ng/mod.rs`'s "Landed so far" index both described the merge
stage as landed when three parameter newtypes had; the module's own doc was honest, so
the parents contradicted the child.
**Categories:** smells, reliability. **Fix:** both now say what landed.

**Mi9: the `pub`/crate-private mismatch was recorded nowhere.**
**Categories:** module_structure, idiomatic. **Fix:** the module doc now states why it
is `pub` and that the intent is the architecture's. The arch half is open question 4.

### Nits

- `nonzero(default: u32)` was a converter named as an adjective with a bare-adjective
  parameter — now `non_zero_default(default_value: u32)`, which is also how std spells
  it. *(naming, smells)*
- `get()` was not `const fn`, so a default could not be used in a const context. Now it
  is. *(defaults)*
- Line-pinned citations into production (`cohort_integration.rs:166`,
  `variant_caller.rs:380`) rot silently. All three were correct at review time — verified
  by two agents independently. The line numbers are dropped; the symbol names, which
  survive an edit, stay. *(smells, module_structure)*
- The threefold newtype boilerplate trips the duplication heuristic. **Kept**:
  `src/ng/types.rs` repeats the same shape a dozen times, each with its own doc comment,
  and a macro would hide those. *(smells)*
- The derived `PartialOrd`/`Ord`/`Hash` are public behaviour with no test. **Deferred to
  A4**, which is where a locus span first gets compared against `MaxCohortLocusSpan` —
  the comparison should carry its own test rather than riding in on `derive`.
  *(reliability)*
- The module name `run` needs a two-word phrase (*calling run*) the reader has not met.
  **Won't fix**: `read` sets the precedent and the arch fixes the path. *(naming)*
- A `parameters.rs` beside `close.rs`/`build.rs`/`organise.rs` would keep `mod.rs` a
  front door. **Deferred**: arch's *Module home* puts the constants in `mod.rs`.
  *(module_structure)*

## 7. Out of scope observations

- **49 pre-existing clippy errors** under `--all-targets --all-features`, across 20
  files in `examples/`, `benches/`, and the test code of `census_file.rs`, `ssr_fit.rs`,
  `open_bam.rs` and `src/ssr/cohort/sim.rs`. All are mechanical lints
  (`needless_range_loop`, byte-string literals, `useless_vec`, unused `Result`). Two of
  them sit in `src/ssr/`, which ng is forbidden to edit, so greening the gate fully is
  not this milestone's to do. Suggested follow-up: its own commit, as `ce3f0b4` did for
  the same class of red.

## 8. Missing tests added now

- `get_returns_the_wrapped_value_not_the_default` — the operator-set value class, plus
  `NonZeroU32::MIN` and `NonZeroU32::MAX`. Catches an accessor that ignores its argument
  (M1).
- Three constant assertions inside `the_defaults_are_the_documented_values` — the
  `pub` constants a command line will advertise, read directly (Mi3).

## 9. What's good

- Copying production's *value* while refusing its *name* is the right call, and the
  reasoning holds against both production call sites: one constant there feeds two
  different tests, and binding to it would let a retune of the per-sample filter move
  ng's cohort keep silently. *(defaults agent, checked at both sites)*
- The module has no code dependency on frozen production at all — its only `use` is
  `std::num::NonZeroU32`, and the single `crate::` mention is a rustdoc link. That is
  exactly the copy-don't-depend shape the freeze rule asks for. *(module_structure)*
- Keeping the defaults as `u32` at their declaration, where an operator looks for them,
  while the typed route is the associated constant.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --lib --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib`
- `./scripts/dev.sh cargo test --lib ng::run::cohort_merge` — the two tests of this step.
