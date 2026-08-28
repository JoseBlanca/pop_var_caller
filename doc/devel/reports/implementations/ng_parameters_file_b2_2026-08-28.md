# ng parameters file — B2: the file shape → TOML text

**Date:** 2026-08-28
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone B, step B2
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §4
**Code:** [src/ng/calling/parameters_file/to_toml.rs](../../../../src/ng/calling/parameters_file/to_toml.rs), with the golden file it produces at [testdata/every_shape_as_written.toml](../../../../src/ng/calling/parameters_file/testdata/every_shape_as_written.toml)

---

## 1. Plan

`ParametersFile::to_toml` — the text a run writes, emitted by hand rather than by `serde`,
because spec §4 chose TOML for the comments and no serde serializer can emit one. Step B3 attaches
the comments; this step is the text they attach to.

## 2. The layout, and the one place it departs from §4

**One row a line, as an inline table**, for every table keyed by an axis: the read groups, the
census terms, the calibration, the contamination, both batching axes, the samples, the
`(stratum × slippage group)` pairs, both length-spectrum rungs, and the substitution rates. A row
names its own axis value, so a line reads on its own and a person editing one sample's coefficient
edits one line. An empty table is `key = []`.

**§4 suggests the numeric rows of §3.7 be "arrays of arrays rather than arrays of tables", and
this writer emits inline tables there too.** The plan's preconditions put the choice here — *"the
spelling is the coder's: key names, nesting, where a table starts, **whether a row is inline**"* —
and two things decide it:

- **a slippage row is not flat.** It carries where its level and its two shares came from, and a
  number off a curve carries the whole curve — seven fields for a slippage curve, ten for a share
  curve, four braces inside the row's own fields. Flattened to a bare array the row becomes a sequence of
  numbers whose meaning is a column legend the reader has to hold, and the nesting has nowhere to
  go.
- **the reason §4 gave for the suggestion no longer applies.** It says arrays of arrays "neither
  needs a custom encoder" — and this *is* the custom encoder, written because the comments
  demanded one.

**Arrays are one element a line where they are a table, and inline where they are a number.** A
length spectrum's `shares_by_repeat_offset` is a handful of numbers describing one thing and stays
on its row.

## 3. Assumptions and choices

- **Every section opens with its own `[header]`, including the two whose every field is a table.**
  `serde` gives those none, so the largest section of the file — `[repeat_tracts]` — opens unnamed
  under the derived writer. That is one crate's rendering, not the file's design.
- **A whole number is written with its decimal point.** In TOML `1` is an integer and would not
  deserialise into an `f64`; `1.0` does. Rust's `Debug` formatting for a float always carries a
  point or an exponent and emits the shortest decimal that round-trips, which is exactly what is
  wanted here. **Whether it round-trips every value is step C3's**, on adversarial ones; this step
  owes only that what it writes is a TOML float.
- **`nan`, `inf` and `-inf` are written in TOML's spelling**, where Rust's own formatting gives
  `NaN` and `inf`. A run should produce none of them; a file that carries one should say so in a
  form a reader can parse rather than one it fails at.
- **A string is escaped only where TOML requires it.** Backslash, quote and the control
  characters; everything at or above a space passes through as itself, including non-ASCII —
  escaping a sample name would make the one field a person is most likely to search for
  unfindable.
- **The writer spells the enums itself** rather than deriving them from `Serialize`. Deriving
  would make the golden file a tautology: a renamed variant would move the writer, the golden file
  and the reader together. `the_hand_written_words_are_serdes_words` compares the two lists, so
  they agree only if both are right.
- **A second golden file rather than a shared one.** `testdata/every_shape.toml` is serde's output
  and stays as the key surface under that writer; `testdata/every_shape_as_written.toml` is this
  one's, from the same fixture. **269 lines against 81** for the same content, which is the size
  of the layout difference.

## 4. Changes made

| file | change |
|---|---|
| `parameters_file/to_toml.rs` | new, 982 lines — `ParametersFile::to_toml`, three layout primitives, two scalar formatters, one function a row, and the file's word for each unit-variant enum |
| `parameters_file/mod.rs` | `mod to_toml;`, and `a_file_using_every_shape` widened to `pub(super)` so both writers' golden files come from one fixture |
| `parameters_file/testdata/every_shape_as_written.toml` | new — the golden file this writer produces |
| `parameters_file/from_run_parameters.rs` | one test: a fitted run's parameters written and read back |

## 5. Tests

**Fourteen new, one of them ignored** — thirteen in `to_toml.rs` and one in `from_run_parameters.rs`. The module's suite went 45 → 58 passing; ten landed with the step and **four came out of the review**, which found that two of spec §5's absences and both empty-array paths were written by code no test exercised.

- `the_written_text_reads_back_as_the_same_file` — the half of goal 1 step B can hold: what the
  writer puts on disk is the value it was handed. C4 owns the other half.
- `the_whole_shape_writes_the_documented_toml` and `regenerate_the_written_golden_file` — the key
  surface and layout, pinned against a checked-in copy, with the regeneration ignored so it never
  makes its own test pass.
- `every_row_of_every_table_is_one_line` — the layout claim, checked structurally: between a
  `key = [` and its bracket there are exactly as many lines as the table has rows, each indented
  and each ending in a comma.
- `a_whole_number_is_written_as_a_float` — including a parse-back proving `1.0` comes out a float
  and not an integer.
- `the_three_values_that_are_not_numbers_are_written_in_tomls_words`.
- `a_name_that_needs_escaping_is_escaped_and_reads_back` — a quote, a backslash, a tab, a newline,
  a control character and a non-ASCII name, each round-tripped through the `toml` parser.
- `the_hand_written_words_are_serdes_words` — all eighteen unit variants (four warrants, four
  seed rungs, three reaches, four share-curve rungs, three share shapes) plus the two
  contamination sources.
- `every_absence_is_a_missing_section_or_a_missing_key` — no `[contamination]` for an
  uncontaminated run, no `measurement` key for an unmeasured read group, no `observations` beside
  a defaulted value.
- `a_fitted_run_writes_a_file_that_reads_back` (in `from_run_parameters.rs`) — Checkpoint B's own
  claim, end to end from a `RunParameters`.

## 6. Validation

Run in the dev container, `./scripts/dev.sh`:

| command | result |
|---|---|
| `cargo fmt --check` | clean, exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean, exit 0 |
| `cargo test --lib ng::calling::parameters_file` | **58 passed, 0 failed, 2 ignored** (54 before the review's fixes) |
| `cargo test --lib` | **4,981 passed, 0 failed, 13 ignored** |

`cargo test --all-targets --all-features` is not the gate: pre-existing panic in
`benches/psp_writer_perf.rs:386`, verified on clean `main`.

## 7. What the first produced file shows, and it is the trigger the plan named

The plan says the key names and the layout are revised "the first time a person reads a file this
writer produced". There is now one, and one thing in it is worth a decision.

**A slippage row is 758 characters on the fixture.** Measured on
`every_shape_as_written.toml`: of 81 lines, 11 are over 120 characters, and the three longest are
758, 561 and 365 — all of them `slippage_by_stratum_and_group` rows. The cause is A1's shape
decision that a number off a curve carries the whole curve, so that a reach cannot be written
without the curve it is a reach into: a row that took a level curve *and* a share curve inlines
seventeen curve fields, and **the same curve is repeated on every stratum of its period**.

Everything else is short. The longest line that is not a slippage row is **183 characters** — a
contamination row carrying a measurement — and the two per-sample inbreeding rows are **142 and 158
bytes** — 142 and 153 characters, the second carrying the fixture's non-ASCII name — which
brackets the 146 bytes spec §9 prices the per-sample axis at.

This is not B2's to change — the curve rides on the smoothing variant by A1's design, and moving
it to a table of its own keyed by (period, slippage group) is a change to the shape. **Raised at
Checkpoint B** with the numbers.
