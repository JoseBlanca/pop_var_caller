# ng parameters file — A2: a value and its warrant, as one serialised shape

**Date:** 2026-08-28
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone A, step A2
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §2
**Code:** [src/ng/calling/parameters_file/mod.rs](../../../../src/ng/calling/parameters_file/mod.rs)

---

## 1. Plan

Spec §2 says every number in the file is a value plus a warrant plus a count of what was behind
it. A1 shipped three spellings of that idea and one number with no warrant at all. A2 gives it one
shape, `WarrantedValue`, so that a reader who has understood one warranted number has understood
all of them.

## 2. What it is, and the five numbers that use it

```rust
pub struct WarrantedValue {
    pub value: f64,
    pub warrant: Warrant,
    pub observations: Option<EvidenceCount>,
}
```

Five numbers carry it: the base-quality calibration multiplier, the per-sample inbreeding
coefficient, the tract ladder's stated concentration, the repeat-tract substitution rate, and the
repeat-tract outlier weight. **The last of those was not in the first draft** and is the review's
find: spec §8 names three parameters with an honest default that must all be marked `Defaulted`,
and A2 had moved two of the three.

**The outlier weight has two reachable states, which is why it needs the shape at all.** The run's
own value is `defaulted` — 0.01 from `likelihood/ssr.rs`. A value a person typed into the file is
`supplied`, and spec §3.8 says that a person editing it is the whole point of writing it down.
Without a warrant, a run reports an edited guess as the project's own constant, and spec §2.1's
wholesale demotion of a mismatched file has nowhere to write itself for this one number.

## 3. The evidence count names its unit, in the file

`EvidenceCount` is `reads` / `covered_positions` / `bases_compared` rather than a bare integer, so
the file says what it counted:

```toml
observations = { reads = 812344 }
observations = { covered_positions = 180600412 }
observations = { bases_compared = 40122 }
```

**This is not a hypothetical improvement.** The first draft left the unit out of the key on the
grounds that it follows the quantity and belongs in the doc comment — and then the doc comment
gave the wrong unit for two of the three. An inbreeding coefficient is fitted over **covered
reference positions** (`generic/runs.rs:698`, summed across windows of 100,000 bases), not over
windows; a repeat-tract substitution rate is fitted over **bases compared**
(`ssr/mod.rs:1093`), not over reads. The three units differ by orders of magnitude on one cohort,
so a reader comparing two numbers' evidence without knowing which unit is which is not comparing
anything.

## 4. Assumptions, and what B1 inherits

- **Two of the five numbers need data `RunParameters` does not hold, and the projection has to be
  handed it.** The calibration's read count is on the `Estimate<ErrorRate>` that
  `RunParameters::assemble` reads and does not store; the inbreeding coefficient's warrant *and*
  count are on the `Estimate<InbreedingF>` the same seam reduces to a bare `Vec<InbreedingF>`.
  Spec §3.3 asks for the first by name. **So step B1 projects from the pre-pass's estimates, not
  from `RunParameters` alone**, or the file says *supplied* where the run fitted and drops two
  counts it has.
- **A deviation from the plan's letter, recorded.** The plan says the shared shape is "the value,
  the warrant, and the observation count behind it". The count is `Option`, because two of the
  five have none: the stated concentration is a median over strata rather than an estimate with a
  sample size, and the outlier weight is a stated constant. Absent, never zero.
- **`Warrant` is not a synonym for "whatever stands behind a number".** Five other quantities in
  the file carry a different kind of warrant and each keeps its own word — a contamination
  fraction's two evidence counts, the ordinary-site prior's `rung`, a slippage number's *origin*,
  a length spectrum's *placement*, and a read group's slippage group, which is declared rather
  than estimated. `WarrantedValue`'s doc lists all five, so that a reader looking for the sixth
  finds out why there isn't one.

## 5. Changes made

One file, [src/ng/calling/parameters_file/mod.rs](../../../../src/ng/calling/parameters_file/mod.rs),
and its golden copy. `WarrantedValue` and `EvidenceCount` added; five fields moved into the first;
the module's two `warrant`-adjacent doc sections rewritten. No other module is touched, and nothing
outside `src/ng/calling/mod.rs`'s one `pub mod` line references any of it.

## 6. Tests

Eleven became thirteen, one of them ignored:

| test | what it holds that nothing else does |
|---|---|
| `every_warranted_number_is_written_the_same_way` | **which five fields are warranted, by path** — a sixth, or a fifth that lost its warrant, fails here; so does a number written flat |
| `an_evidence_count_names_its_unit_and_is_absent_where_there_is_none` | the three units, by name, and the two absences |
| `a_stated_concentration_of_one_reads_back_saying_whether_it_was_fitted` | re-pointed at what a *reader* recovers; asserting only that two documents differ passes for any shape |

The first asserts a list of field paths rather than a count of rows, so adding a read group to the
fixture for an unrelated reason no longer breaks it with a message about warranted numbers.

## 7. Validation

All in the container.

| command | result |
|---|---|
| `cargo fmt -- --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib ng::calling::parameters_file` | **12 passed, 0 failed, 1 ignored** |
| `cargo test --lib` | **4,932 passed, 0 failed, 12 ignored** |
| `cargo doc --no-deps` | zero unresolved links in this module |

## 8. Tradeoffs and follow-ups

- **`serde` gives a section no header line when every one of its fields is a table.** Both
  `[repeat_tracts]` and `[stated_constants]` have vanished from the golden file, so the largest
  section of it opens unnamed. That is `serde`'s rendering, not the artefact's design — step B2's
  writer emits its own headers — and the module doc now says to read the golden file as a record
  of key names rather than as the file a run will produce.
- **The nesting costs a level.** `error_probability_multiplier = { value = …, warrant = …,
  observations = { reads = … } }` is three levels where the first draft had one. In B2's inline
  form it is still one line a row; under `serde`'s array-of-tables rendering it is four header
  lines. The alternative — the count on the row under a unit-naming key like
  `reads_behind_the_rate` — keeps the shape to two keys and was argued for in review; it was not
  taken because three of the five numbers do carry a count and splitting a value from its count is
  the thing spec §2's spine exists to prevent.
