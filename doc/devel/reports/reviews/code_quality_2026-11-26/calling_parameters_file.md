# Code-quality map: `src/calling/parameters_file/` (2026-11-26)

Read-only review. Nothing was built or run; every finding is reasoned from the code at
`f84129cf` on branch `code-quality`.

## Summary

1. **Bug, likely: a direct-mode run that numbers its read groups differently from the run that
   wrote the file gives four tables to the wrong read group.** Each read group's repeat-tract
   substitution rates, stutter multipliers, slippage group and misread-after-runs table are joined
   on the file's own read-group number, while the run numbers its read groups in the order its BAM
   headers list them. `estimate-parameters` numbers them by sample name and `@RG ID`, so the two
   orders differ whenever the BAMs are not given in that order. A single sample with two lanes is
   enough. The result is wrong repeat-tract (STR) likelihoods, with no message.
2. **Bug, likely: a parameters file whose `fitted_from.read_groups` rows are out of order
   swaps calibrations and contamination between lanes** in a calling run. Only a hand-written or
   hand-edited file can have this, but the spec invites both, and the test that claims such a
   file "is the same file" never projects it.
3. **On the size question: the schema is not spelled out five times.** About 10,700 of the 24,500
   lines are production code, and 4,000 of those are comments. Reading already uses serde. The
   writer is hand-written, which spec §4 requires, because a serde serializer cannot emit the
   comments that are the reason the format is TOML. The one mechanical second copy of the
   schema is the writer's ~500 lines of one-function-per-row formatters. A custom serde serializer
   could replace them for a net saving of about 300 lines. I recommend leaving them.
   Validation and the two projections between the file and `RunParameters` encode rules, not
   field lists, so no derive can replace them.
4. **The `crate::cli` import is test-only.** Three tests in `from_run_parameters.rs` call
   `estimate-parameters`' wrapper around `validate()`. Calling `validate()` directly removes the
   edge in about 15 minutes.
5. **The documentation has gone stale in places that matter to a fixer**: false statements about
   what calls `validate`, the default multiplier, how many numbers carry a warrant, and how
   read-group rows are indexed (the comment at the heart of bug 2). One function also carries
   another function's doc comment.
6. **There are small duplications worth merging**: the same "format-N table in an older
   file" refusal written seven times, and row key paths spelled two ways in error messages.

---

## Bugs

### 1. Read-group-keyed repeat-tract tables are joined on the file's numbering, not the run's

- **Where** — `src/calling/parameters_file/to_run_parameters.rs:376-380` (slippage group by read
  group), `:514-529` (stutter multipliers), `:543-544` (substitution-rate key),
  `:643-644` (misread-after-runs table). Compare `:233`, `:295`, `:315`, which correctly map
  through `order.read_group[...]`.
- **What** — `project_onto` receives `order: &WhereTheFilesRowsBelong`, which maps each of the
  file's read-group numbers to the run's position. This follows the owner's ruling of 2026-09-04
  that a file is matched to a run by sample and `@RG ID`, never by number. The calibration,
  contamination and batch tables use `order`. But `slippage_fits()`, `substitution_rates()` and
  `misread_after_runs()` take only `&self` and key their maps with `ReadGroupId(row.read_group)`,
  the file's own number. Calling then looks these up with the run's number.

  **The failure.** `estimate-parameters` builds its read-group table with samples in name order
  and each sample's read groups in `@RG ID` order (`src/cli/estimate_parameters.rs:714-716`).
  `call-from-alignments` numbers read groups in the order of its BAM paths and each header's
  `@RG` lines (`src/read/input/read_groups.rs:288-301`). Take one sample whose BAM header declares
  lane `L2` before `L1`. The file numbers `L1 = 0, L2 = 1`, and the run numbers `L2 = 0, L1 = 1`.
  `refuse_if_not_this_runs_inputs` accepts the file, as the ruling says it must. `L1`'s stutter
  multipliers, substitution rates and misread table are then applied to `L2`'s reads, and the
  other way round. With one slippage group and similar lanes the effect is small. With a
  PCR library beside a PCR-free one, stutter multipliers differ several-fold, so STR genotype
  likelihoods change. Every multi-lane or multi-sample direct-mode run whose input order is not
  name order is exposed. A run with one read group is not.
- **Why the tests miss it** — `the_same_samples_in_another_order_are_the_same_cohort`
  (`bindings.rs:1416`) says "every value lands on the plant it was fitted for", but it
  asserts only the inbreeding coefficients and the calibration counts.
- **Proposal** — Pass `order` into the three functions and key with
  `ReadGroupId(order.read_group[row.read_group as usize] as u32)`. `misread_after_runs` checks
  `row.read_group < read_group_count` before indexing, so the lookup is safe. Extend the reorder
  test to assert a substitution rate, a stutter multiplier and a misread entry for a moved lane.
- **Downside** — None in the code. It is a bug fix, so it needs approval.
- **Effect** — **Output changes**, only for direct-mode runs whose read-group order differs from
  the file's. psp mode and identically ordered runs are unchanged, because there
  `of_the_files_own_order` and `of_this_run` are both the identity. No speed or memory change.
- **Cost** — S. **Confidence** — high. **Status: likely** (traced from `calling_run.rs:1433`
  through `to_run_parameters_for` to `project_onto`).

### 2. `WhereTheFilesRowsBelong::of_this_run` indexes by row position, and its callers index by read-group number

- **Where** — `src/calling/parameters_file/to_run_parameters.rs:165-174`. The comment at `:176-178`
  states the opposite of what the code does. Callers are at `:233`, `:295`, `:315`. The related
  claim is at `bindings.rs:355-358`, and `validate.rs:588-607` sorts the ids before checking them.
- **What** — `read_group` is built by `file.fitted_from.read_groups.iter().map(...).collect()`, so
  entry *i* belongs to the *i*-th **row**. Every caller indexes it with `row.read_group`, the
  read-group **number**. `validate` accepts rows in any order, because it sorts the ids before its
  density check. So a file listing `fitted_from.read_groups` as `[rg 1, rg 0, rg 2]` validates
  and passes the binding check, then puts read group 0's calibration multiplier, contamination
  measurement and batch onto the run's lane for read group 1, and the other way round.
  `of_the_files_own_order` is unaffected, because it uses `0..n`.
- **Why the tests miss it** — `a_file_whose_rows_are_written_in_another_order_is_the_same_file`
  (`bindings.rs:1648-1658`) swaps two rows and asserts only that `validate` and the binding check
  pass. It never projects, so "is the same file" is untested.
- **Proposal** — Build the vector by number: `let mut v = vec![0; n]; for row in rows {
  v[row.read_group as usize] = run_position; }`. Extend that test to call
  `to_run_parameters_for` and compare against the unswapped file. Fix the comment at `:176-178`.
- **Downside** — None.
- **Effect** — Output changes only for a file with out-of-order rows. The writer never produces
  one, so only hand-written or hand-edited files are affected.
- **Cost** — S. **Confidence** — high. **Status: likely.**

### 3. The misread-after-runs table skips `validate`, so `estimate-parameters`' self-check does not cover it

- **Where** — `validate.rs:120` and `:128-145` (only a format-version gate);
  `to_run_parameters.rs:580-679` (the real checks); `src/cli/estimate_parameters.rs:582` and
  `:693-701` (`a_file_the_next_run_can_read`, which runs `validate()` only).
- **What** — The checks on the misread table live in the projection, not in `validate`: a row
  naming no table entry, a duplicate entry, an excess outside [0, 1), and a non-zero `as_the_run`
  further on. Two consequences follow:
  1. `estimate-parameters` checks the file it writes with `validate()` alone, so it can write a
     misread table the next run's reader refuses. That is the failure the check exists to stop.
     It is reachable only if the fit itself emits a bad row, which I could not show.
  2. The table's two warranted numbers never go through `validate.rs`'s `a_warranted_value`.
     A NaN or negative `standard_error`, or `observations` in the wrong unit, is accepted. The
     error is used only for reporting, so no call changes. A `fitted_here` row with no
     `standard_error` reads back as 0.0 (`to_run_parameters.rs:655-656`) and is written out as
     `standard_error = 0.0`, which the file's own header says never means zero.
- **Proposal** — Move the row checks into a `validate` function, as every other table has, and
  run `a_warranted_value` on `as_the_run` and `other_wrong`. Keep the projection's
  `RunEntry::of` lookup, which `validate` can call too.
- **Downside** — `validate` would then need `run_context::RunEntry`, one more import. The refusal
  messages and paths stay the same.
- **Effect** — No output change on any file a fit writes. Hand-edited files with a bad standard
  error are refused where they were accepted.
- **Cost** — S–M. **Confidence** — medium. **Status: suspected** for part 1, likely for part 2.

---

## Architecture

### 4. The `crate::cli` dependency is test-only, and it can go

- **Where** — `from_run_parameters.rs:3281-3282`, `:3598`, `:3651` (all inside `mod tests`, which
  starts at `:1266`). The target is `src/cli/estimate_parameters.rs:693-701`.
- **What** — Production code in this module imports nothing from `crate::cli` or `crate::run`.
  Three tests call `a_file_the_next_run_can_read`, which is `file.validate()` with its error
  wrapped in `EstimateParametersCliError::RefusedByItsOwnReader`. One test
  (`a_fit_with_an_error_of_zero_beside_a_penalised_number_stops_before_writing`) also matches on
  that variant. The test's own doc says the check "is `estimate-parameters`'s own; it is tested
  here because the fixtures … are here."
- **Proposal** — In all three tests, call `projected(&run, &read_groups).validate()` and assert on
  the `ParametersFileError` directly. The property under test, that the writer's output passes
  the reader's check, is unchanged. If the wrapping into `RefusedByItsOwnReader` needs pinning,
  add a three-line test in `src/cli/estimate_parameters/tests.rs` on any refused file.
- **Downside** — The test no longer exercises the exact function `estimate-parameters` calls. That
  function is one line, so little is lost.
- **Effect** — none. **Cost** — S. **Confidence** — high.

### 5. The module's shape is justified, and I would not restructure it

- **What was checked** — I read spec `parameters_file.md` §4 and §8, the implementation plan
  `doc/devel/ng/impl_plan/parameters_file.md` (milestones A–C), and the module headers that
  record the reasons. Four decisions explain the size, and each was taken on purpose:
  - **The file's types are separate from `RunParameters`** (`mod.rs:11-31`). The file names its
    rows by sample and `@RG ID` and expresses absence as a missing key. `RunParameters` is dense
    and cannot say "absent". So two projections, `from_run_parameters` (818 code lines) and
    `to_run_parameters` (721), are unavoidable. They encode rules such as "a defaulted value
    carries no count" and "an absent contamination section is the plain formula". They do not
    copy fields.
  - **Mirror enums with `From` impls in both directions** (`from_run_parameters.rs:1130-1225`,
    `to_run_parameters.rs:910-970`; about 120 lines). The reason is at `mod.rs:443-448`: the
    file's spelling is a compatibility surface, and renaming an in-memory variant must not
    silently re-read every file on disk.
  - **The writer is hand-written**, because no serde serializer emits comments (spec §4;
    `to_toml.rs:3-21`). See finding 7.
  - **`validate` is separate from parsing** (`validate.rs:11-16`), so a caller can read and
    inspect a file without meeting the run's constraints.
- **Proposal** — Leave the module structure as it is. The cheap structural cleanup is finding 4.

---

## Duplication

### 6. The answer to "is the schema spelled out four or five times?"

The numbers, from counting non-blank lines outside `#[cfg(test)]` blocks:

| file | code | comments | what the code is |
|---|---|---|---|
| `validate.rs` | 1,564 | 462 | meaning checks: ranges, density, warrant rules, format gates |
| `to_toml.rs` | 1,473 | 401 | section layout, ~450 lines of prose written into the file, ~500 lines of row formatters |
| `from_run_parameters.rs` | 818 | 406 | projection from `RunParameters` to the file |
| `to_run_parameters.rs` | 721 | 212 | projection from the file to `RunParameters` |
| `mod.rs` | 571 | 1,222 | the file's types, with serde derives |
| `bindings.rs` | 520 | 429 | refusing a file fitted on other inputs, and demotion |
| `from_toml.rs` | 339 | 244 | parse, migrate format 1–2, reword two parser messages |
| other three | 361 | 617 | defaults, what-was-fitted, between-locus share |
| **total** | **~6,370** | **~4,000** | |

Tests hold the other ~13,800 lines, in 304 `#[test]` functions.

Where a key name or variant word appears in production code:

1. **The types in `mod.rs`, with serde attributes.** This is the single declarative description,
   and the reader already uses it (`from_toml.rs:125`). Only `EvidenceCount` is read by hand,
   for one owner-ruled message (`from_toml.rs:538-607`).
2. **The row formatters in `to_toml.rs:1457-1876`**, one function per row type (about 25),
   restating each row's keys and order.
3. **The variant words**: `the_word_for_*` in `to_toml.rs:1888-1964`, plus `the_word_for` and
   `the_unit_of` in `validate.rs:1408-1415` and `:1555-1562`. This copy is deliberate
   (`to_toml.rs:1878-1886`), so that a rename cannot move the writer and the reader together,
   and a test compares it with serde's words.
4. **Key paths inside error messages** in `validate.rs` and `to_run_parameters.rs`. Any
   path-naming error message needs these. See finding 9.

So there are two real spellings of the schema (1 and 2), one deliberate guard copy (3), and
error-message paths (4). The projections and `validate` are not schema copies.

### 7. Row formatters could become a custom serde serializer. Recommendation: leave them

- **Where** — `to_toml.rs:1329-1343` (`an_inline_table`) and `:1457-1876` (row functions);
  `:1888-1964` (variant words); the test `the_hand_written_words_are_serdes_words`
  (`to_toml.rs:2293`).
- **What** — The rows could be rendered by a small `serde::Serializer` that writes any
  `Serialize` value as a one-line TOML inline table. It would use the module's own
  `a_toml_float`, `a_toml_integer` (which saturates) and `a_toml_string`, skip `None` fields, and
  spell unit variants as strings and data-carrying variants as `{ variant = {...} }`. Sections and
  their notes would stay hand-written. About 500 lines of row code and 80 lines of variant words
  would become a serializer of roughly 200–250 lines, a **net saving of about 300 lines**.
  Byte identity is checkable: `testdata/every_shape_as_written.toml` pins the writer's output.
  One struct would need its fields re-declared in writer order: `WarrantedValue` (`mod.rs:508`)
  is declared value, warrant, observations, standard_error, and written value, standard_error,
  warrant, observations (`to_toml.rs:1850-1864`). That reorder changes only the test golden
  file `every_shape.toml`.
- **Downside** — Each row function today is ten obvious lines. A `Serializer` implementation is
  ~30 trait methods that a geneticist reader will not follow, so readability goes down. Omissions
  are already caught: the round-trip test fails if the writer drops a field. The saving is about
  3% of production lines. The deliberate separation of the writer's words from serde's
  (point 3 above) would also go, though the as-written golden file still catches a rename.
- **Effect** — Output must be byte-identical. That is checkable, but the risk is not zero.
  **Cost** — M. **Confidence** — medium. **Recommendation: leave it**, unless the format keeps
  adding a table a month, in which case it removes one of the four places each new table
  touches.

### 8. Seven copies of the "format-N table in a file stamped older" refusal

- **Where** — `validate.rs:130-145`, `:157-167`, `:245-255`, `:260-271`, `:307-320`,
  `:404-423`, `:1165-1176`; seven copies in all. Each has the same message text: "is a format-{N} table in a file whose
  `format_version` says {}, so it is not clear which format the file was written in; set
  `format_version` to …".
- **Proposal** — `fn not_older_than_its_format(at: &str, added_in: u32, found: u32, kind: &str)`,
  with `kind` being "table" or "key". This saves about 50 lines and gives the message one source.
- **Downside** — None of substance.
- **Effect** — none (identical messages). **Cost** — S. **Confidence** — high.

### 9. Row paths in error messages are built in two places, and spelled two ways for one key

- **Where** — `validate.rs:1097-1100` writes `inbreeding.by_sample[sample = "TS-1"]…`.
  `to_run_parameters.rs:346-349` writes `inbreeding.by_sample["TS-1"]…`, and so does the doc at
  `mod.rs:360` and `validate.rs:26`. The slippage-row path is at `validate.rs:1197` and
  `to_run_parameters.rs:408`. The substitution-rate path is at `validate.rs:1316` and
  `to_run_parameters.rs:540`.
- **What** — The same row's path is formatted independently in both files. For inbreeding the two
  spellings already differ, so a reader meets two forms for one key. The `to_run_parameters`
  copy is unreachable today, because `validate` refuses first.
- **Proposal** — Give each keyed row type a `fn path(&self) -> String` in `mod.rs`, using the
  `[key = value]` form `validate` uses everywhere else, and call it from both files.
- **Downside** — It adds a method to each row type, which spreads the message vocabulary into
  `mod.rs`.
- **Effect** — Only the text of an unreachable message changes. **Cost** — S.
  **Confidence** — high.

### 10. Small one-liners

- `to_run_parameters.rs:720-735` (`an_evidence_count`) is `count.map_or(0, EvidenceCount::count)`,
  using `mod.rs:567-573`. Effect none, cost S.
- `validate.rs:1408-1415` (`the_word_for`) duplicates `to_toml.rs:1923-1930`
  (`the_word_for_warrant`). Making the latter `pub(super)` and using it in both places is cost S.
  It is a different copy from the deliberate guard copy in finding 6, point 3, because both
  of these are hand-written.

---

## Readability

### 11. Stale or false statements a fixer will trip over

Each item is S and has no effect on output. **Confidence** — high (each was checked against the
code).

- `validate.rs:17-19`: "**Nothing in a run calls this yet.**" This is false. `bindings.rs:583` and
  `to_run_parameters.rs:206`, `:221` call it.
- `mod.rs:5-9`: "Today those numbers exist only in memory … Nothing writes them down and nothing
  can read them in." The module does both.
- `mod.rs:25`: the table says an unfittable error rate is "a multiplier of 1.0 whose `warrant` is
  `defaulted`". Spec §5 and the file's own note (`to_toml.rs:196`) say it is 2.0 or a stated rate
  divided by the library's mean.
- `mod.rs:53-56` and `:162-187`: these describe the hand-written writer as future ("step B2 is
  meant to emit … today"). The example tree at `:59` says `format_version = 3`, and the current
  version is 10.
- `bindings.rs:679`: "Five numbers in this file carry a `Warrant`". The demotion below it moves
  nine keys: the multiplier, inbreeding, two fallbacks, the substitution rate, the stutter
  multiplier, two misread excesses, and the outlier weight.
- `bindings.rs:667`: "`defaulted`-at-the-project's-own-0.01". The built-in is 0.07
  (`calling/likelihood/ssr.rs:181`).
- `bindings.rs:704-706` promises that rows are destructured, so a new warranted key stops
  compilation. The misread rows at `:779-783` are not destructured, so a third warranted field
  there would be silently left undemoted. Destructure them like the others.
- `bindings.rs:355-358`: "the projection reads the table only for its length". `of_this_run`
  reads it positionally (bug 2).
- `to_run_parameters.rs:176-178`: "The vector is indexed by the file's own number." This is false,
  and it is bug 2.

### 12. A doc comment sits on the wrong function

- **Where** — `validate.rs:746-765`.
- **What** — The doc for `no_repeat_tract_table_names_one_thing_twice` ("No repeat-tract table
  names one thing twice…") is joined onto the doc of `every_group_number_is_below_the_read_group_count`.
  The function it describes, at `:824`, has none.
- **Proposal** — Move `:746-757` to above `:824`.
- **Effect** — none. **Cost** — S. **Confidence** — high.

### 13. The production writer's documentation lives on a test-only function

- **Where** — `from_run_parameters.rs:203-355` (`of_run`, `#[cfg(test)]`, about 130 lines of
  doc). The production function is `:357-373` (`of_run_under`), whose doc says only
  "As `of_run`" for both its behaviour and its panics.
- **What** — Rendered docs and anyone reading the production path find the explanation on a
  function that does not exist outside tests.
- **Proposal** — Move the doc onto `of_run_under`. Possibly also rename it `of_run`, and have the
  test helper take a plain `&StrRepeatCriteria`.
- **Effect** — none. **Cost** — S. **Confidence** — high.

### 14. The comment volume, and internal step labels

- **What** — Comments are about 4,000 of 10,700 production lines. In `mod.rs` they outnumber code
  two to one (1,222 against 571). Much of it is design rationale the owner has asked to keep.
  A noticeable share is history: "until 2026-08-31", "Measured before this type existed",
  "the first draft of this step's own test", and plan step codes such as B2, C2, D3, E3 and F1,
  which mean nothing to a reader without the plan (writing skill, Rule 5).
- **Proposal** — Fix the false statements in finding 11 first. Then, separately, replace step codes
  with what they name ("the reader", "the bindings check"), and move change history that no
  longer explains the current code into the commit log. Do not cut rationale to reach a size.
- **Downside** — Judging which history still explains the code is a call the owner may make
  differently. Prose edits on this scale are slow to review.
- **Effect** — none. **Cost** — M. **Confidence** — medium.

### Tests: duplicated setup and claims not tested

- The three-lane `ReadGroups` fixture (HWI.3, HWI.4, HWI.5) is built separately in
  `defaults.rs:530`, `from_run_parameters.rs:1299` and `bindings.rs:1317` (`the_run_the_fixture_was_fitted_by`),
  and they differ only in the second sample's name. One shared builder next to
  `a_file_using_every_shape` (`mod.rs:1979`) would do. This is low value, and each copy's comment
  explains its shape.
- Two tests claim more than they check, and both hide a bug above:
  `the_same_samples_in_another_order_are_the_same_cohort` (bug 1) and
  `a_file_whose_rows_are_written_in_another_order_is_the_same_file` (bug 2).

---

## Looked at, left alone

- **`validate` runs twice when a run reads a file** (`bindings.rs:583` and again inside
  `project_onto`). It runs once per run, on a file of at most tens of MB, and the reason is
  documented at `:598-603`.
- **Demotion clones the whole file** (`bindings.rs:696`). This happens once per run and is
  dropped after projection. Demoting the file rather than the parameters is what makes demotion
  per-file, as spec §2.1 requires.
- **`ReadsBehindEachCalibration::nothing_was_fitted` and `no_count_in_reads` have identical
  bodies** (`from_run_parameters.rs:88-99`). They name two different sources, so keep both.
- **The mirror enums and their `From` impls in both directions** keep the file's spelling
  independent of in-memory variant names (`mod.rs:443-448`).
- **`testdata/every_shape.toml` (serde's own layout) beside `every_shape_as_written.toml`.** The
  first also proves that the reader accepts a layout the writer never produces
  (`from_toml.rs:774`), which is spec goal 3.
- **`EvidenceCount`'s hand-written `Deserialize`** (`from_toml.rs:551-607`) exists by owner ruling
  for one message, and it refuses exactly what the derive refuses.
- **A derive-based validator crate (such as `garde` or `validator`) instead of `validate.rs`.**
  It would cover only the per-field range checks, about a third of the file. It would name rows
  by index (`by_sample[3]`) instead of the file's `[sample = "TS-1"]`, and the cross-field rules
  (density, warrant-by-key rules, format gates, sample order) would stay hand-written anyway.
- **Library dependencies**: `calling::parameters_file` imports from `parameter_estimation`,
  `read::input`, `repeat_catalog`, `alignment` and `region_typing`. All of these are in the
  allowed direction (`mod.rs:33-42`). Nothing in production imports `cli` or `run`.
- **Integer formatting**: some `u64` keys are written with `to_string()`
  (`to_toml.rs:820-836`, `:539-547`) and others with the saturating `a_toml_integer`. They differ
  only above 2^63, which no count or length reaches.
