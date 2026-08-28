# Fixes applied — ng parameters file, A1

**Date:** 2026-08-28
**Review:** [ng_parameters_file_a1_2026-08-28.md](ng_parameters_file_a1_2026-08-28.md)
**Code:** [src/ng/calling/parameters_file/](../../../../src/ng/calling/parameters_file/)

---

## 1. What changed, in one paragraph

Every Blocker and Major was applied, and the shape moved enough that the module was rewritten
rather than patched: it is **1,508 lines against 960**, and **34 types against 32** — 26 structs and 8
enums, where `ShareWarrant` and `BlendedSource` went and `LevelSmoothing`, `ShareSmoothing`,
`ReadGroupBatchRow` and `SampleBatchRow` arrived — with **eleven tests against four**, one of them
a golden file and one an ignored regenerator for it. The two changes with the largest reach are that the curve and
the reach now ride on the smoothing variant, so three meaningless states stopped being writable;
and that every struct refuses an unknown key, so a typo in a hand-edited file is an error rather
than an absence.

## 2. The two mutations that prove the tests changed

Both were survivors in the review — every test passed while the file on disk changed. Run in the
container against the fixed tree, from a checksummed copy of the file and restored to the same
checksum afterwards (`8ed0427e…`, verified before and after):

| mutation | before the fixes | after |
|---|---|---|
| drop `#[serde(rename_all)]` from `ShareShape` | all 4 tests passed | **2 fail** — `every_enum_variant_spells_as_the_file_says` and `the_whole_shape_emits_the_documented_toml` |
| rename the field `markers_with_reads` to `marker_positions` | all 4 tests passed | **1 fails** — `the_whole_shape_emits_the_documented_toml` |

The second is the one worth reading twice: **a `serde` round trip can never pin a key name**,
because it moves both sides of a rename together. Only a golden file can.

## 3. Findings table

| # | severity | status | note |
|---|---|---|---|
| B1 | Blocker | **Applied** | the round-trip test's doc corrected to what it holds, and two tests added that hold the rest |
| M1 | Major | **Applied** | batching axes now three read groups over two samples, and named rows rather than positional arrays |
| M2 | Major | **Applied** | `every_enum_variant_spells_as_the_file_says` — all 21 unit variants |
| M3 | Major | **Applied** | `stated_length_spectrum_warrant` beside the concentration, with its own test |
| M4 | Major | **Deferred → Checkpoint A** | no step of the plan owns range checking; open question 1 |
| M5 | Major | **Deferred → step C** | the property-based round trip belongs with the reader |
| M6 | Major | **Applied** | `the_documented_inline_form_parses` |
| M7 | Major | **Applied with adaptation** | rewritten against the parsed value's keys, both directions |
| M8 | Major | **Deferred → A2** | see §5 |
| M9 | Major | **Applied** | curve and reach ride on the smoothing variant |
| M10 | Major | **Deferred → C4, with the check named in the field's own doc** | see §5 |
| M11 | Major | **Applied** | `deny_unknown_fields` on every struct, with the trade recorded |
| M12 | Major | **Applied** | `fitted_from_repeats` / `fitted_to_repeats` / `fitted_from_reads_of` |
| M13 | Major | **Applied** | `LevelOrigin` / `SharesOrigin`, and *smoothing* for what they hold |
| Minors — widths | Minor | **Applied** | `reference_repeats` `u64` everywhere; `cells` and `strata` `u64` |
| Minors — newtypes | Minor | **Deferred → A2** | see §5 |
| Minors — `Copy` | Minor | **Applied** | dropped from the five composite types, which now nest a curve anyway |
| Minors — names | Minor | **Applied**, 14 of them | see §4 |
| Minors — doc accuracy | Minor | **Applied**, all 3 | see §4 |
| Minors — the 32 `pub` types | Minor | **Applied** as a recorded decision | narrowing is not available at A1; the module doc now says so with its revisit point |
| Nits | Nit | **Applied** where a name changed anyway; the rest folded into §4 | |

## 4. What was applied

**The shape.**

- **The curve and the reach ride on the smoothing variant.** `LevelSmoothing` and `ShareSmoothing`
  are `this_stratum` / `this_periods_curve { curve, reach }` / `blend { curve_weight, curve,
  reach }`. Three states stopped being writable: a reach with no curve, a curve with no reach, and
  "its period's curve, whole" with no curve recorded. **The module's own fixture had written the
  first of the three.** The cost is that the one shared `BlendedSource` becomes two enums, because
  a level's curve and a share's curve are different types; that is recorded in the type's doc.
- **`deny_unknown_fields` on all 26 structs**, with the trade written beside `ParametersFile`: a
  file from a later build carrying an unknown key is now refused rather than partly read, and for
  a file where a dropped key changes what a number means, refusing loudly is the safer half.
- **`stated_length_spectrum_warrant`**, so the run's own fitted median and the stated constant —
  both able to be exactly 1.0 — are not the same line.
- **`reference_repeats` is `u64` in all three rows**; `cells` and `strata` are `u64` rather than
  `usize`, which is platform-width and has no business in a file format.
- **The sequencing batching is named rows**, `{ read_group, batch }` and `{ sample, batch }`,
  rather than two positional integer arrays — which is the module doc's own stated rule, and the
  one place the first draft broke it.
- **`ContaminationRow` names its library**, which spec §3.4 asks for and the first draft omitted.
- **`Copy` dropped** from `LevelOrigin`, `SharesOrigin`, `SlippageCurve` and `ShareCurve`.

**The names.** `FittedFromInputs` → `InputsFittedFrom`; `CensusBinding` → `CensusIdentity`;
`ReadGroupRow::name` → `declared_id`; `CalibrationRow` → `BaseQualityCalibrationRow` and its
`scale` → `error_probability_multiplier`, whose doc now says which way it moves;
`ContaminationRow::fitted_from` → `fitted_from_reads_of` and its two variants →
`this_read_groups_own_reads` / `every_read_of_this_sample`; `SequencingBatching` →
`SequencingBatches` and `declared` → `was_declared_by_the_run`; `InbreedingRow::coefficient` →
`inbreeding_coefficient`; `SeedRegimeInFile` → `SeedRung`, `regime` → `rung`, and
`FallbackDiversity` → `StatedHeterozygosity`; `CurveReachInFile` → `CurveReach`;
`ShareShapeInFile` → `ShareShape`; `SlippageCurveRow`/`ShareCurveRow` → `SlippageCurve`/
`ShareCurve` (neither was a row); `ShareCurve::source` → `rung`, and `centre` →
`centre_repeats`; `weights` → `shares_by_repeat_offset`; the three `by_`-shaped tables of
`[repeat_tracts]` → `slippage_group_by_read_group`, `slippage_by_stratum_and_group`,
`substitution_rate_by_stratum`; `ItsPeriodsCurve` → `ThisPeriodsCurve`, symmetric with
`ThisStratum`.

**The documentation.**

- The round-trip test no longer claims to pin a field ordering `toml` refuses. **It does not
  refuse one:** the crate emits a struct's scalar fields before its table-valued ones whatever the
  declared order, so moving `ploidy` to the last field of `ParametersFile` leaves the emitted
  bytes identical. The doc now says what the test holds, what it cannot, and which test holds the
  rest.
- The module header records **why the module lives under `calling/`** — `calling` already imports
  from `parameter_estimation` in thirteen files and the reverse is zero, so this placement adds no
  edge where a top-level peer would add one into each.
- It records that **`serde` writes a struct's tables last**, so a `Blend` row's `smoothing` is
  emitted after its siblings — which a hand-written writer meaning to match serde's bytes has to
  reproduce, and which the module's inline sketch does not show because inline tables do not have
  the problem.
- It records that the 32 `pub` types are `pub` **because narrowing is not available at A1** —
  per-type `pub(crate)` emits `private_interfaces` warnings that `-D warnings` makes fatal, and
  `pub(crate) mod` fails clippy with 32 dead-code errors — with B1/C2 as the revisit point.

**The tests.** Four became eleven, one of them ignored:

| test | what it holds that nothing else does |
|---|---|
| `every_section_of_the_shape_survives_a_serde_round_trip` | the tree holds no enum shape TOML cannot spell and no field that reads back as another |
| `the_whole_shape_emits_the_documented_toml` | **every key and every spelling**, against a checked-in golden file |
| `every_enum_variant_spells_as_the_file_says` | all 21 unit variants, of which the golden file sees 11 |
| `a_smoothing_that_used_a_curve_carries_it_and_a_plain_one_carries_neither` | the three states the fix removed stay removed |
| `an_absent_shares_origin_writes_no_key_and_a_present_one_does` | absence *and* presence, on the parsed value's keys rather than a substring |
| `a_file_with_every_table_empty_round_trips` | the empty boundary, which is inside the committed input range |
| `a_stated_concentration_of_one_says_whether_it_was_fitted` | a fitted 1.0 and a stated 1.0 are not the same file |
| `a_mistyped_key_is_refused_rather_than_absorbed` | both an extra key and a misspelled optional one |
| `the_documented_inline_form_parses` | the shape the hand-written writer will emit and serde never does, four levels deep |
| `the_format_version_this_build_writes_is_one` | bumping the version is a deliberate act with a test to change |
| `regenerate_the_golden_file` *(ignored)* | the only thing that may rewrite the golden file, run deliberately |

**The fixture.** Three read groups over two samples, so the two batching axes have different
lengths and exchanging them cannot produce a file that parses; two libraries of one plant, which
is the grain the contamination fraction exists at; a slippage-group map that is not the identity;
a sample name that needs escaping and is not ASCII (`Ailsa ‘Craig’ "×2"`); and one float that
needs all its digits (`1/3`), because with every float a short decimal, narrowing one from `f64`
to `f32` had left the emitted file byte-identical. **The doc comment claiming "every numeric value
is distinct" was false and is replaced by the three properties the fixture actually has**, with
its one deliberate exception named.

## 5. What was deferred, and to whom

- **M4, the seven unenforced invariants → Checkpoint A**, as open question 1. No step of the plan
  owns range checking: A3 owns `Option`, C owns reading, and neither owns refusing a value outside
  its stated range. This is a gap in the plan rather than a coding choice, so it goes to the owner
  rather than being invented here.
- **M5, the property-based round trip → step C**, whose heart is the round trip.
- **M8, drift between the twelve mirrored types and their upstream originals → A2.** The fix is a
  `#[cfg(test)]` exhaustive-match witness per type — about 90 lines that turn a future upstream
  addition into a build failure here. It belongs beside A2's shared value+warrant shape, which is
  where the last of the mirroring lands.
- **M10, the slipped-read count written twice → C4, with the check named in the field's own doc.**
  Nothing available here can settle whether the two upstream fields hold one number: their docs
  are near-identical and their absence conditions are worded differently. `SharesOrigin::
  slipped_reads` now says so and names C4's real fit as where to compare them across every
  stratum.
- **File-local scalar newtypes for the read-group index, the repeat count and the motif period →
  A2.** Two categories asked for them and the argument is good — a newtype defined in this module
  changes none of the file's bytes and would have made the width split a compile error. The
  reason to wait is that A2 introduces the file's own value shapes and doing both at once would
  spell the same quantity twice; the observed hazard (a review swapped a row's two `u32` field
  names and every test passed) is closed in the meantime by the golden file.
- **Full-precision floats, non-finite values and counts past `i64::MAX` → C3**, which is the step
  the spec assigns the float question to. The fixture now carries one full-precision value so that
  narrowing a float is no longer invisible.
- **`SubstitutionRateRow::observations` keeps its generic name → A2**, which gives every
  value+warrant+count row one spelling; renaming it here and again there is churn.

## 6. Validation

All in the container, by absolute path.

| command | result |
|---|---|
| `cargo fmt -- --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --lib ng::calling::parameters_file` | **10 passed, 0 failed, 1 ignored** |
| `cargo test --lib` | **4,930 passed, 0 failed, 12 ignored** |
| `cargo test --all-targets --all-features` | every test target passes — **4,930 in the library and 92 across the fifteen other test targets, 0 failed**. The run exits 101 on a **benchmark**, not a test: `benches/psp_writer_perf.rs:386` panics with `index out of bounds: the len is 3300000 but the index is 3300000` in its own priming loop. **Verified pre-existing**, by running that bench in a worktree at `main` with none of this change: same line, same numbers, same exit code. Nothing outside `src/ng/calling/mod.rs`'s one `pub mod` line references this feature's module |
| `cargo doc --no-deps` | zero unresolved links in this module (25 remain elsewhere in the tree, pre-existing) |

**One environment note, because it cost two aborted runs.** `cargo test --all-targets` was killed
twice with exit 137 while four containers from a concurrent session held 56 GB of the host's 64.
Neither was a failure; both re-ran clean once the contention passed.
