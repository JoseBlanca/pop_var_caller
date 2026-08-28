# Code Review: ng_parameters_file_b1
**Date:** 2026-08-28
**Reviewer:** rust-code-review skill (orchestrator), five category agents in isolated worktrees
**Scope:** step B1 of the parameters-file plan — `RunParameters` → the file shape
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** the uncommitted working-tree diff of step B1, handed to each agent as a
  patch applied over `90885e48` in its own worktree.
- **In-scope files:**
  - [from_run_parameters.rs](../../../../src/ng/calling/parameters_file/from_run_parameters.rs) — new
  - [mod.rs](../../../../src/ng/calling/parameters_file/mod.rs) — `SeedRung` gained a fourth variant
  - [run_parameters.rs](../../../../src/ng/calling/run_parameters.rs) — seven read-only accessors
  - [stratum_fits.rs](../../../../src/ng/parameter_estimation/joint/stratum_fits.rs) — three iterators
- **Out of scope:** the TOML writer (B2), the reader and its `validate` (C1, C2), the bindings'
  refusals (D2, D3), and the pre-existing `benches/psp_writer_perf.rs` panic.
- **Categories dispatched:** reliability (always; owns the mutation pass), errors (always; this
  step chose release panics over a typed error), naming (always; and the numbers in changed prose
  are its to re-derive), idiomatic + smells (one agent, both checklists), and a design-fidelity
  pass reading the code against spec §3 and §5 section by section.

## 2. Verdict

**Request-changes.** Two Blockers, both measured rather than argued, and both about the tests
rather than the projection: one of the step's own `#[should_panic]` tests could not fail, and the
one join the code documents as checked was not.

## 3. Execution status

Run by the orchestrator in the container, on the reviewed tree:

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | clean |
| `cargo test --lib ng::calling::parameters_file` | 0 | 30 passed, 0 failed, 1 ignored |
| `cargo test --lib` | 0 | 4,950 passed, 0 failed, 12 ignored |

`cargo test --all-targets --all-features` was **not run as a gate**: it exits 101 on a pre-existing
panic in `benches/psp_writer_perf.rs:386`, verified on clean `main` and standing in
`PROJECT_STATUS.md`. `cargo doc --no-deps` likewise fails on 25 pre-existing unresolved intra-doc
links elsewhere in the tree.

Findings labelled "Needs verification": **none**. Every Blocker and Major below was proved by
running a mutation or a probe in the reporting agent's own worktree, with the output quoted in the
per-category file.

## 4. Open questions and assumptions

1. **Should the projection return a `Result` rather than panic?** (affects M1.) It runs after the
   last locus, so a panic discards a cohort's calling work where `RunParameters::assemble`'s
   equivalent checks discard a startup. No call site exists yet.
2. **Is a slippage number's *origin* a substitute for a *warrant*?** (affects M9.) Spec §3.7 asks
   for "the slippage numbers, with the warrant on each"; the shape carries an origin instead, and
   §2.1's wholesale demotion to `Supplied` then has nowhere to write itself for those numbers.
   This is the same open item Checkpoint A raised and the owner has not yet ruled on.
3. **Does a median over a run's other strata read as `fitted_here` or as `borrowed`?** (affects
   M10.) The spec settles neither.

## 5. Top 3 priorities

1. **B1** — `a_read_group_table_that_is_not_the_runs_is_refused` cannot fail: deleting the guard it
   names leaves all 30 tests green.
2. **B2** — the calibration row's evidence count is joined to an unchecked argument, and no part
   of that join is tested; a rate map from another fit is accepted and writes its counts.
3. **M2** — swapping `intercept` and `slope` in both curve conversions survives the whole suite.

## 6. Findings

### Blocker

**B1: from_run_parameters.rs — `a_read_group_table_that_is_not_the_runs_is_refused` cannot fail.**
*Categories: reliability.* Its `#[should_panic(expected = "the read-group table names")]` matches
both of `of_run`'s guards, and its one-lane fixture trips both. With the read-group guard removed,
the suite stays green and `of_run` accepts a four-lane table against a three-read-group run,
returning a file with four identity rows, three calibration rows and three batching rows — `lib6`
named and never calibrated. **Fix:** give the two guards messages that share no opening, and split
the test into two fixtures, each tripping only its own guard.

**B2: from_run_parameters.rs — the calibration row's count comes from an unchecked argument.**
*Categories: reliability, design-fidelity (convergent).* The value and warrant come from the run's
own `ReadGroupCalibration`; the count comes from the passed-in rate map, checked only for the id's
presence. A rate map with the run's three ids and another fit's numbers was accepted, writing
`Reads(7)` beside a multiplier built from 812,344 reads. The `# Panics` doc claims more than the
code does. The inbreeding path, two functions down, does check its equivalent. **Fix:** compare the
calibration's warrant against the rate's, allowing the one legitimate disagreement (assembly
substitutes `Defaulted` when it refuses a rate), and test both halves.

### Major

**M1: from_run_parameters.rs — seven panics abort a finished run to protect a record, not a
result.** *errors.* The doc justifies them by pointing at `RunParameters::assemble`, whose checks
prevent *wrong genotypes* and fire before the first locus. Here every genotype is already computed
and a mismatch corrupts a provenance record beside the VCF. **Fix (or the deferral's terms):**
return a `Result` and let the driver log and keep its VCF; or keep the panics and argue *this*
function's case in the doc.

**M2: from_run_parameters.rs — the two curve conversions carry sixteen fields and the tests assert
eight.** *reliability.* Swapping `intercept` and `slope` in both survives; so do zeroing
`held_out_error` on both, `bend` on the share curve, and the share curve's fitted range. These are
field-by-field renames between structs whose names differ, which is where a transposition happens.
**Fix:** assert every field of both curves once, on a fixture where the two curves' fitted ranges
differ.

**M3: from_run_parameters.rs — `was_declared_by_the_run` is only ever tested `true`.** *reliability.*
Hard-coding it `true` survives. It is the only thing separating a run that declared nothing from
one that declared a single batch holding every library. **Fix:** a test on a defaulted batching.

**M4: from_run_parameters.rs, stratum_fits.rs — eight enum-mirror match arms are entered by no
test.** *reliability.* `LevelSource::Cell`, `ShareSource::Blend`, `Provenance::Supplied`, three
`ShareCurveSource` fallbacks and two `ShareShape` variants — among them `ShareShape::Flat`, which
the project's own measurement makes the commonest real answer. Exhaustiveness catches a variant
*added* upstream, not one *crossed* here: `Flat => Sloping` compiles and spells correctly. **Fix:**
a table naming every pair, plus fixtures reaching the two arms that are not `From` impls.

**M5: from_run_parameters.rs — the documented "exactly the run's read groups" check is
half-enforced.** *errors.* A rate set covering *more* read groups is accepted in silence (probed: a
five-read-group rate set accepted for a three-read-group run). **Fix:** one size assertion.

**M6: stratum_fits.rs — three `expect()`s with no `// PANIC-FREE:` comment, one naming neither the
stratum nor the group.** *errors.* The level-provenance `expect` fires inside a walk over every
stratum of a run. **Fix:** name the stratum and the group; remove the other two by having the
length-spectrum iterators hand out the borrowed shape rather than a `LengthSpectrum` whose
`fitted_weights` is an `Option` the producer can never leave empty.

**M7: from_run_parameters.rs, stratum_fits.rs — the three new `pub` iterators have no tests of
their own, and one stamps a rung nobody observes.** *reliability.* Changing
`FittedFrom::ThisStratum` to `ItsPeriodsPooledTracts` survives the whole 4,950-test library suite.
**Fix:** test them where they live, asserting the rung.

**M8: from_run_parameters.rs — five documented guards have no test.** *reliability.* The
sample-count assert, the missing-rate panic, the inbreeding length assert, `reach_of`, and the
share path's `curve_of`.

**M9: mod.rs — a slippage number carries an origin and no warrant.** *design-fidelity.* Spec §3.7
asks for a warrant on each; `SlippageRow`, `OrdinarySitePrior` and both length-spectrum rows have
none, so §2.1's wholesale demotion has nowhere to write itself — which is the argument the code
itself makes for why the outlier weight *does* carry one. **Convergent with the standing open item
from Checkpoint A**, and the owner's to rule on.

### Minor

- **Mi1** *design-fidelity, errors (convergent).* A substitution rate keyed past the run's
  read-group axis is written, naming a read group the identity block does not list (probed:
  identity names 0, 1, 2 and the substitution table names 0, 2, 7).
- **Mi2** *design-fidelity.* The tract ladder's stated concentration is stamped `FittedHere` when
  it is a median over *other* strata — §2's definition of `Borrowed`, generalised from read groups
  to strata. The spec settles neither.
- **Mi3** *design-fidelity.* The projection can write a contamination section in which no row has
  a measurement — the uncontaminated run longhand, which C2 is told to refuse. Not reachable from
  the estimator (it refuses below 100 markers), only by hand.
- **Mi4** *design-fidelity.* An unmeasured contamination row drops a fraction the run may have
  scored with. Safe only because `fit_alpha` returns exactly zero where no marker carries a read —
  a coupling stated nowhere.
- **Mi5** *reliability.* `ploidy` is written from a fixture where four candidate sources all give 2.
- **Mi6** *reliability.* The shares are never looked up at a slippage group other than 0, so
  `row.shares[0]` is byte-identical on the fixture.
- **Mi7** *idiomatic.* Two helpers rebuild a `ReadGroupId` from an index where `ReadGroups::iter()`
  hands one over with the crate's own widening idiom.
- **Mi8** *idiomatic.* `warranted` takes the count and its unit as two parameters, one a function
  pointer, where one `EvidenceCount` would do.
- **Mi9** *idiomatic.* `ssr_substitution_rate` hands out a `&BTreeMap` where the same change added
  iterators everywhere else.
- **Mi10** *idiomatic.* Three helpers take the whole `&RunParameters` to read one field.
- **Mi11** *smells.* The `SeedRung` doc opens with the history of a bug, told in three places.
- **Mi12** *smells.* The module's design rationale lives in a private module's `//!`, which
  `cargo doc` never renders, and `of_run` points at it.
- **Mi13** *smells.* `each_stratum_and_group_with_numbers` re-implemented the tail of `at`, down to
  a duplicated `expect` string.
- **Mi14** *smells.* The six-argument `of_run` call is spelled out five times in the test module.
- **Mi15–Mi24** *naming.* `of_each_…` promises a dense walk for two iterators that filter;
  `parameters` is the one word this module cannot use unqualified; `|of|`, `held`, `views`,
  `batching`, `what` are half-names; `warranted` is a bare participle; and **four quantitative
  claims in changed prose are wrong** (below).

### Wrong numbers in the diff's own prose (review step 8a)

Ten claims were re-derived and correct — the "fourteen of twenty-two" spelling count, the
`0.4999999998066437` multiplier and its 3.9 in 10¹⁰ shortfall, the 4.8 in 10⁷ bound, the 4,242
reads, the 0.25 scale, the six-cells-two-filled count, and every spec citation. **Four were
wrong**, all of them the author's own claims about the author's own fixture:

| claim | truth |
|---|---|
| the two inbreeding counts are "five orders of magnitude apart" | 180,600,412 against 9,411,027 is a factor of **19** |
| the module doc calls the stated concentration "defaulted by construction" | the code marks it `FittedHere` on any run that fitted a stratum, and a test asserts it |
| `ZeroDiversity` "sits between two long comment blocks" | there is **one**, above it |
| "eleven of twenty-one was wrong on both halves" | the **denominator was right**; twenty-one became twenty-two because B1 added a variant |

### Nits

Redundant `+ '_` on RPITs; `SsrPeriod::try_new(period as usize)` where `usize::from` is
infallible; a repeated digest literal; an assertion message whose line continuation renders with
35 embedded spaces; `_of` doing two unrelated jobs in one file; a `?` inside a struct literal
inside a `Some`; and a fixture doc that claims more than it should about which batchings
`declared` accepts.

## 7. Out of scope observations

- **Spec §3.6 is inaccurate about a value this code reads.** "An alternative total of exactly zero
  … is not floored on the way in or out" — the seed builder floors it at `MIN_ALT_CONCENTRATION`
  before `RunParameters` holds it, and the new `SeedRung::ZeroDiversity` is what carries that
  state. Offered for the spec, not the code.
- **`StratumFits::over` checks that its three per-group vectors are the same length and not that
  they agree cell by cell**, which is what makes the level-provenance `expect` reachable from a
  hand-built outcome.
- **A `SlippageGroup(u32)` newtype** would stop a bare `u32` crossing five new signatures
  positionally; `slippage_group_of` set that precedent before this step.

## 8. Missing tests to add now

The eight the reliability agent specified, plus the two split refusal tests from B1 and the
provenance-agreement test from B2: `a_read_group_table_with_another_runs_library_count_is_refused`,
`a_read_group_table_with_another_runs_sample_count_is_refused`,
`a_rate_set_missing_one_of_the_runs_read_groups_is_refused`,
`a_rate_set_over_more_read_groups_than_the_run_has_is_refused`,
`a_rate_set_from_another_fit_is_refused`, `an_inbreeding_estimate_list_of_another_length_is_refused`,
`a_level_off_a_curve_that_does_not_say_how_far_it_reached_is_refused`,
`a_share_that_came_off_a_curve_with_no_curve_recorded_is_refused`,
`a_level_that_is_the_stratums_own_carries_no_curve`,
`every_field_of_both_curves_reaches_the_file_under_its_own_name`,
`a_run_that_declared_no_batching_writes_the_flag_false`,
`every_pre_pass_word_maps_to_its_own_word_in_the_file`,
`of_run_writes_a_single_sample_single_read_group_run`,
`an_evidence_count_of_zero_is_written_and_is_not_absence`,
`a_substitution_rate_keyed_past_the_runs_read_groups_is_refused`,
`a_read_group_measured_at_nothing_may_not_carry_a_fraction`, and three in `stratum_fits.rs` for the
new iterators.

## 9. What's good

- **`a_defaulted_number_carries_no_evidence_count` is the strongest test in the diff** — read group
  1's rate is exactly zero, which assembly refuses, and it still carries 4,242 reads, so the
  absence cannot come from a zero count.
- **The contamination projection gates on `was_measured()` and never on the fraction**, which is
  precisely the distinction spec §5's second row exists to protect.
- **The design-fidelity pass found §3 covered in full**, at the right grain and unit in every
  section, including §9's correction of the substitution rate's axis.
- **`every_seed_regime_has_a_rung_in_the_file`** is the shape the other three enum mirrors should
  have copied, and the reliability agent said so.
- **The `# Panics` sections and the assertion messages name both differing values**, in the shape
  the census's own refusals use.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib ng::calling::parameters_file`
- `./scripts/dev.sh cargo test --lib ng::parameter_estimation::joint::stratum_fits`
- `./scripts/dev.sh cargo test --lib`

Audit trail: the five per-category files in the gitignored
`tmp/review_2026-08-28_parameters-file-b1/`.
