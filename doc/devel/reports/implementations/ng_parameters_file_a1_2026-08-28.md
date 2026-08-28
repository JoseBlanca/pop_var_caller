# ng parameters file — A1: the file's Rust shape, and the provisional TOML tree

**Date:** 2026-08-28
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone A, step A1
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §3
**Code:** [src/ng/calling/parameters_file/mod.rs](../../../../src/ng/calling/parameters_file/mod.rs)

---

## 1. Plan

One type per section of spec §3 — identity and binding, ploidy, the per-read-group calibration
rows, the contamination table and its batching, the per-sample inbreeding rows, the prior seed,
the repeat-tract tables, and the compiled-in constants — with `serde` derives and no reading or
writing. **The key names and the shape of the tree are the coder's proposal** by the owner's
decision of 2026-08-28, revised later against a file a fitted run actually produced; §2 below is
what that revision has to argue with.

## 2. The TOML tree that was chosen, and what was rejected

### 2.1 The tree

**As it shipped, after the review.** The first draft's tree and the changes the review made to it
are in [the fixes report](../reviews/fixes_applied_ng_parameters_file_a1_2026-08-28.md) §4.

```toml
format_version = 1
ploidy = 2

[fitted_from]                    # §3.1
reference_digest = "…"
samples = ["TS-1", "TS-2"]
read_groups = [ { read_group = 0, declared_id = "HWI.3", library = "lib3", sample = "TS-1" } ]
[fitted_from.census]
terms = [ { term = "the loci actually kept", digest = "…" } ]

[base_quality_calibration]       # §3.3
by_read_group = [ { read_group = 0, error_probability_multiplier = 1.0324,
                    warrant = "fitted_here" } ]

[contamination]                  # §3.4
by_read_group = [ { read_group = 0, library = "lib3", fraction = 0.031,
                    markers_with_reads = 4211, reads_on_markers = 90233,
                    fitted_from_reads_of = "this_read_groups_own_reads" } ]

[sequencing_batches]             # §3.4
was_declared_by_the_run = false
by_read_group = [ { read_group = 0, batch = 0 } ]
by_sample = [ { sample = "TS-1", batch = 0 } ]

[inbreeding]                     # §3.5
by_sample = [ { sample = "TS-1", inbreeding_coefficient = 0.42 } ]

[ordinary_site_prior]            # §3.6
reference_concentration = 1.0
alternative_concentration_total = 0.0006
rung = "fitted_curve"

[repeat_tracts]                  # §3.7
stated_length_spectrum_concentration = 1.0
stated_length_spectrum_warrant = "defaulted"
slippage_group_by_read_group = [ { read_group = 0, slippage_group = 0 } ]
slippage_by_stratum_and_group = [ … one row a (stratum × slippage group) … ]
length_spectrum_by_stratum = [ … ]
length_spectrum_by_period = [ … ]
substitution_rate_by_stratum = [ … ]

[stated_constants]               # §3.8
repeat_tract_outlier_weight = 0.01
```

`src/ng/calling/parameters_file/testdata/every_shape.toml` is the file the fixture actually
produces, 237 lines, and it is checked in as the golden copy every key is pinned against.

### 2.2 The four choices that shape everything else

**`format_version` and `ploidy` are bare keys at the top, not a `[run]` table.** They are the two
questions a person opens the file to answer first, and TOML puts bare keys before the first table
header anyway. *Rejected:* a `[run]` table holding one integer, which buys a header line and
costs a level of nesting on the one value a user is most likely to check.

**Every row that names a read group uses the key `read_group`, and it is always the run's dense
index `0..n`.** That index is what
[`ReadGroupParameters::calibration_of`](../../../../src/ng/calling/likelihood/mod.rs) indexes by,
so it is the join five sections are written against. The `@RG ID` string appears exactly once, as
`name` in the read-group table of `[fitted_from]`. *Rejected:* `id` for the index, which collides
with the `@RG ID` a reader is thinking of. **The string is `declared_id`** — declared by the
alignment file, as against the run's own index. The first draft called it `name`, which says less
than either, in a row that already carries two other names.

**`reference_repeats` is the repeat count a stratum is the bin for, everywhere.** The two
`Stratum` types in the tree spell that field differently —
[`census::Stratum`](../../../../src/ng/parameter_estimation/joint/census.rs) calls it
`reference_repeats` and [`ssr::Stratum`](../../../../src/ng/parameter_estimation/ssr/mod.rs)
calls it `repeats` — and they are the same quantity. *Rejected:* following each type's own
spelling, which would put two words for one number in one file.

**Sections are named for the quantity, not for the code that produces them.** `[fitted_from]`
rather than `[identity]` or `[binding]`: what a person wants to know is what these numbers were
fitted from, and both alternatives name the mechanism instead. `[base_quality_calibration]`
rather than `[calibration]`: there are three different things in this file a run calibrates
against, and the long name says which — **and the row type is `BaseQualityCalibrationRow` for the
same reason**, where the first draft called it `CalibrationRow`, the bare word this paragraph had
just rejected. `[stated_constants]` rather than `[defaults]`: they are
constants the project has stated, and `[defaults]` would read as *values a reader may omit*,
which they are not.

### 2.3 One departure from spec §4's letter, keeping its intent

Spec §4 asks for "each sample's row as a single inline table on one line, and the numeric rows of
§3.7 as arrays of arrays rather than arrays of tables". **Every row in the tree is an inline
table, including §3.7's** — one row a line, no `[[array of tables]]` headers, which is the
readability §4 is buying. An array of arrays is not workable for the slippage rows: each carries
two warrants, and a warrant is a nested structure (a source, a curve, a reach, a read count), so
a bare array would need a column legend beside it to be read at all. The length-spectrum
`weights` are arrays of floats, which is what an array of arrays would have given them anyway.

### 2.4 Two shapes chosen to make an illegal state unwritable

**The blend's weight rides on the variant, not beside it.** `BlendedSource` is
`this_stratum` / `its_periods_curve` / `blend { curve_weight }`, so *a blend with no weight* and
*a weight with no blend* are both unwritable. The first draft had `source` and a separate
`curve_weight: Option<f64>` and could express both. In the file that is serde's own enum
spelling: `source = "this_stratum"`, or `source = { blend = { curve_weight = 0.31 } }`.

**One `BlendedSource` for the level and for the two shares.** In memory these are two enums —
[`LevelSource`](../../../../src/ng/parameter_estimation/joint/slippage_curve.rs) and
[`ShareSource`](../../../../src/ng/parameter_estimation/joint/share_curve.rs) — whose variants
differ only in that one calls the stratum's own answer `Cell` and the other calls it `Stratum`.
They are the same three states, and a reader of the file should not have to learn two spellings
of them.

### 2.5 The census binding is twelve named digests, not the terms themselves

[`RecordingTerms`](../../../../src/ng/parameter_estimation/joint/census.rs) holds a per-stratum
locus-count table among its twelve values; written out, it would be the largest thing in the file
for a binding whose only use is an equality. Kept term by term rather than as one digest over all
of them, so that whatever reports the demotion (step D3) can say *which* term differed, in the
words `RecordingTerms::first_disagreement` already uses. *Rejected:* holding the values, and one
digest over all twelve.

## 3. Assumptions, and what the next steps inherit

- **The calibration row carries no observation count.**
  [`ReadGroupCalibration`](../../../../src/ng/calling/likelihood/mod.rs) is a scale and a
  provenance and keeps no count, so there is none to write. Spec §2 asks for "a value plus a
  warrant plus a count of what was behind it"; step A2's shared shape therefore has to make the
  count optional rather than required, and this row is the one that needs it absent.
- **The `Option` fields already here mirror `Option`s that exist in memory** — an absent curve,
  an absent reach, an absent slipped-read count, an absent shares warrant. **They are not step
  A3's work.** A3's five states are the ones spec §5 lists, and none of them is expressed yet:
  `contamination` is still a required section, not an `Option`.
- **The writer cannot be serde's, and this is now visible rather than suspected.** Spec §4 makes
  the format TOML *for the comments*, and no serde serializer emits a comment; the `toml` crate's
  own serializer offers one style knob (`multiline_array`) and no control over inline tables. So
  step B2 writes the emitter by hand. The `Serialize` derives stay, because they are what lets a
  test cross-check the hand-written writer against an independent one.
- **The reader can be serde's.** `toml::from_str` reads an inline table into a struct whatever
  wrote it, so the hand-written writer and a derived reader meet in the middle at step C1.
- **Unknown keys are refused, and that could not have waited for step C1.**
  `#[serde(deny_unknown_fields)]` is a *type* attribute with no knob on `toml::from_str`, on
  `toml::Deserializer` or on `serde::Deserialize`, so "decide it at the reader" was not an
  available option. The review measured what its absence costs: misspelling the `curve` table by
  one letter parses and yields `curve: None`, which this module's own documentation defines as
  *this stratum's period had no curve* — **a typo read back as a fitted fact**. The trade it buys
  is that a file from a later build carrying a key this one does not know is refused rather than
  partly read; for a file where a dropped key changes what a number means, that is the safer half,
  and it is written down beside `FORMAT_VERSION` rather than left to be inferred.
- **No step of this plan owns range checking**, and after step C ships, a file with an inbreeding
  coefficient of 1.7 or a length spectrum summing to 1.4 will be accepted. A3 owns `Option`, C
  owns reading; neither owns refusing a value outside its stated range, and spec §9 promises "a
  malformed file fails at read with a line number". **This is Checkpoint A's question, not a
  coding choice.**

## 4. Changes made

- **New** [src/ng/calling/parameters_file/mod.rs](../../../../src/ng/calling/parameters_file/mod.rs),
  960 lines: 25 structs and 7 enums, every one deriving `Debug`, `Clone`, `PartialEq`,
  `Serialize` and `Deserialize`. Directory form rather than a flat `.rs` because the writer and
  the reader are siblings in steps B and C.
- **Modified** [src/ng/calling/mod.rs](../../../../src/ng/calling/mod.rs): one line, `pub mod
  parameters_file;`.

Nothing else in the tree is touched, and nothing reads or writes the new types.

## 5. Tests added

Four, all in the module's own `#[cfg(test)]` block, over one fixture built so that **every
section is non-empty and every enum variant is used at least once**, and every numeric value is
distinct so a field-for-field comparison cannot pass on two fields that were swapped.

| test | what it would catch |
|---|---|
| `every_section_of_the_shape_survives_a_serde_round_trip` | a field ordering TOML refuses (a scalar after a table), an enum shape it cannot spell, a field that reads back as a different one |
| `each_warrant_spells_as_snake_case` | a renamed Rust variant silently re-interpreting every file on disk |
| `a_blended_source_carries_its_weight_and_the_others_carry_none` | the two illegal states §2.4 removed coming back |
| `an_absent_curve_writes_no_key` | an absence written as something a reader could mistake for a measurement |

**They use serde's own serializer, not the file's writer, which does not exist.** What they pin
is the shape.

## 6. Validation

All in the container, by absolute path (`CLAUDE.md`).

| command | result |
|---|---|
| `cargo fmt -- --check` | clean (after one reorder of the `mod` line) |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib ng::calling::parameters_file` | **4 passed, 0 failed** |
| `cargo doc --no-deps` | fails on the tree, on 25 pre-existing unresolved intra-doc links in other modules; **zero of them are in this module** after one fix (`ShareProvenance` lives in `ssr_fit`, not `share_curve`) |

## 7. Tradeoffs and follow-ups

- **The key names are provisional by decision.** The revision has a trigger rather than a date:
  the first time a person reads a file this writer produced and has to ask what a key means.
- **A curve is written on every row that used it**, so a period's curve is repeated across that
  period's strata. On a tomato run that is about 141 rows against one curve a period. Hoisting
  the curves into their own table keyed by period would remove the duplication and add a join;
  it is worth doing only if a real file looks bad, which is exactly what the naming revision is
  for.
- **`cargo doc` does not pass on this tree** and did not before this step. Out of scope here;
  worth a sweep of its own.
