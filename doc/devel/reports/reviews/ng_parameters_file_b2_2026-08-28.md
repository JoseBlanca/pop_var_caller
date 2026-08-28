# Code Review: ng_parameters_file_b2
**Date:** 2026-08-28
**Reviewer:** rust-code-review skill (orchestrator), three category agents in isolated worktrees
**Scope:** step B2 of the parameters-file plan — the file shape → TOML text
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** the uncommitted diff of step B2, handed to each agent as a patch applied
  over `021f7d8a` in its own worktree.
- **In-scope files:**
  - [to_toml.rs](../../../../src/ng/calling/parameters_file/to_toml.rs) — new, the writer
  - [testdata/every_shape_as_written.toml](../../../../src/ng/calling/parameters_file/testdata/every_shape_as_written.toml) — new, its golden file
  - [mod.rs](../../../../src/ng/calling/parameters_file/mod.rs) — two lines
  - [from_run_parameters.rs](../../../../src/ng/calling/parameters_file/from_run_parameters.rs) — one test
- **Out of scope:** the projection (B1, committed); the reader (C1); float round-trip fidelity as a
  proven property (C3) — though a float this writer emits that is not a legal TOML float at all was
  in scope.
- **Categories dispatched, and one of them is not a category:** reliability (owns the mutation
  pass); a **TOML-correctness** pass, whose only question was whether every document this writer
  can produce is valid TOML that reads back as the value it was given, worked against the v1.0.0
  specification fetched during the review; and naming, extended to **read the produced file as the
  geneticist would** — the plan names the first person to read a produced file as the trigger for
  revising the key names, and there had never been such a file before.

## 2. Verdict

**Request-changes.** Two Blockers, both from the mutation pass, both the same shape: spec §5's
absences are written by code no test exercised, so a writer that turned an absence into a sentinel
left all 54 tests green.

## 3. Execution status

| command | exit | result |
|---|---|---|
| `cargo fmt --check` | 0 | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | clean |
| `cargo test --lib ng::calling::parameters_file` | 0 | 54 passed, 0 failed, 2 ignored |
| `cargo test --lib` | 0 | 4,977 passed, 0 failed, 13 ignored |

**Mutation numbers:** 11 run, **4 survived the whole suite**, 1 survived only the test whose claim
it violates, 1 changed no behaviour, 5 killed. The five killed are listed in the reliability file,
so the report is not one-sided.

Findings labelled "Needs verification": **none**. Every finding below was run.

## 4. Open questions and assumptions

1. **Should the slippage table alone be emitted as `[[array-of-tables]]` headers?** It would trade
   three 758-character lines for roughly forty short ones with no shape change. It is an
   alternative to fixing the shape, not a complement — raised at Checkpoint B.
2. **Which key names does the file's first reader want changed?** The readability pass is that
   reader and named twenty-two keys it had to ask about. Six of them it guessed *wrongly*. The
   plan makes this the owner's revision.

## 5. Top 3 priorities

1. **B1** — two of spec §5's absences are written by untested code; turning either into a sentinel
   is invisible.
2. **B2** — both empty-array paths are untested, and breaking either writes a file that cannot be
   read back at all.
3. **M1** — a count above `i64::MAX` writes digits TOML does not define, and three readers give
   three answers.

## 6. Findings

### Blocker

**B1: to_toml.rs — two of spec §5's absences are written by code no test exercises.**
*reliability.* Hard-coding `was_declared_by_the_run` to `true`, and writing `SharesOrigin::slipped_reads`
as `unwrap_or(0.0)` rather than omitting it, each leave all 54 tests green. Both mutants emit valid
TOML that **reads back as a different value** — an undeclared batching recorded as declared, and a
zero sentinel where the rule is "absence, never a sentinel". Neither state exists in either
fixture. The derived writer has a dedicated test for the second; B2 replaced that writer without
carrying the witness across.

**B2: to_toml.rs — both empty-array paths are untested, and breaking either writes an unreadable
file.** *reliability.* Deleting `one_a_line`'s `key = []`, or `a_float_array`'s `[]`, survives
everything and produces a document that fails to parse — none of the shape's `Vec` fields carries a
serde default, so a missing key is a hard error. The empty table is the single-sample, no-tract end
of the committed range.

### Major

**M1: to_toml.rs — a `u64` above `i64::MAX` writes a document TOML does not define.**
*TOML-correctness.* Eight sites write a `u64` with `to_string()`. The specification caps integers
at 2⁶³−1 and says one that cannot be represented losslessly must throw. On `u64::MAX` in an
evidence count, three readers gave three answers: the shape's derived reader accepted it and the
whole file compared equal; `toml::Value` refused with *u64 value was too large*; Python's `tomllib`
accepted it because Python's integers are unbounded. **This is the "parses as an equal value only
by luck" case** — the round-trip test passes while the artefact is unreadable by anything modelling
integers as `i64`. Out of reach for a real count; an unsigned underflow upstream lands on it.

**M2: to_toml.rs — `every_row_of_every_table_is_one_line` inspects the wrong `by_sample`.**
*reliability.* Two sections hold a `by_sample` and the search runs from the top of the file, so it
finds the batching table; both happen to hold two rows. Emitting the inbreeding rows doubled left
the test green.

**M3: to_toml.rs — the stated reason for writing floats with a decimal point is false.**
*reliability.* The doc says `1` "would not deserialise into an `f64`". Measured: the `toml` crate's
derived deserialiser accepts it, so formatting with `Display` instead of `Debug` leaves the round
trip green. The rule stays right; the reason does not, and **C2's reader will be designed against
it**.

**M4: every_shape_as_written.toml — the one array a user can edit wrongly in silence does not say
where it starts.** *readability.* `shares_by_repeat_offset = [0.1, 0.8, 0.1]` runs `-span ..= +span`
with the middle entry the reference length, and that convention exists only in the Rust. A user who
reads it as `[0, +1, +2]` — the natural reading of an array called "by offset" — produces a file
that parses, deserialises, and shifts every length prior in that stratum by one repeat unit. **The
only edit in the file that is wrong without being invalid.**

**M5: the artefact — spec §1.2 goal 3's own worked example does not work.** *readability.* Raising
one library's error rate: the user knows the library as `lib3`, which appears in the read-group
table and in the contamination table but **not** in the calibration table, which is keyed by
`read_group` alone. So it is a two-line read for a one-line edit. And after the edit,
`warrant = "fitted_here", observations = { reads = 812344 }` still stands, so the file asserts that
a hand-typed number was fitted from 812,344 reads. **Goal 2 is what a goal-3 edit silently
breaks.**

### Minor

- **Mi1** *readability.* Twenty-two keys the file's first reader had to ask about, six of them
  guessed **wrongly**: `stated_length_spectrum_concentration` (read as *user-stated*, beside
  `warrant = "defaulted"`), `fraction` (of DNA or of reads?), `was_declared_by_the_run` (declared by
  the *sequencing* run?), `slipped_reads = 8000.5` (a fractional read count reads as a typo),
  `rung` (used for two different ladders), and `level` — "the softest number in the file, and it has
  the emptiest name in the file".
- **Mi2** *readability.* `reach = "inside"` — inside what? The answer is three keys back on the same
  line. `inside_the_fitted_range` costs 19 characters on a line already 758 long.
- **Mi3** *readability.* "Stratum" names three tables and appears in no row; the rows spell it
  `period = 2, reference_repeats = 6`.
- **Mi4** *naming.* The article convention is applied consistently — all 30 `String`-returning
  helpers carry an article, and the only three that do not are exactly the three that append to
  `out` — but the rule is never stated, and the sibling projection spends its articles on a
  different distinction.
- **Mi5** *naming.* The `*_word` family are dropped possessives that do not survive being read
  aloud: `a_reachs_word`, `a_share_shapes_word`.
- **Mi6** *TOML-correctness.* Two test-coverage gaps: U+007F is handled but not pinned, and
  exponent-form floats are not pinned.
- **Mi7** *readability, cross-category.* `mod.rs` calls `[inbreeding]` "the file's only cohort-sized
  axis", which spec §9's correction of 2026-08-28 contradicts.
- **Mi8** *readability.* The substitution-rate row — the other axis §9 prices at 146 bytes — writes
  at **162 characters** here, 11% over the figure §9's 62 MB estimate is built from.

### Wrong numbers in the diff's own prose (review step 8a)

Twenty claims re-derived; **three wrong and one half-wrong**, all the author's own:

| claim | truth |
|---|---|
| "nine new" tests | **ten** — nine in `to_toml.rs` plus one in `from_run_parameters.rs`; the 45 → 54 arithmetic only works with ten |
| "seventeen unit variants" | **eighteen** (4 + 4 + 3 + 4 + 3); the test covers all eighteen, only the prose was short |
| "nested two deeper" | **four braces** inside the row's own fields; max depth 5 against 1 for a plain row |
| the inbreeding rows are "142 and 158" | those are **bytes**; as characters 142 and 153, and the paragraph switched units mid-sentence |

Correct and re-derived: 269/81 lines, 982, 758/561/365, 11 lines over 120, 183, 17 curve fields,
45 → 54, and the full-suite figures.

### Nits

`an_inline_table`'s empty-slice branch is unreachable from any caller; NaN would fail the derived
`PartialEq` round trip although the TOML it writes is valid; the `toml` crate resolves to a **1.1**
parser, so its acceptance alone is not proof of 1.0 validity — which is why the TOML pass also read
every document with Python's `tomllib`.

## 7. Out of scope observations

- **Sixteen-digit floats are a real tension with no fix.** `held_out_error = 0.3333333333333333`
  sits beside a sibling `0.167`. Goal 1 requires the shortest round-tripping decimal and goal 3 pays
  for it. Worth a sentence so nobody "cleans it up".
- **36% of the file is on 4% of its lines**, and 15% of it is four copies of two curves. The
  nesting, the repetition and the seventeen curve fields are all A1's shape, not B2's format.

## 8. Missing tests to add now

`a_run_that_declared_no_batching_writes_the_flag_as_false`,
`a_shares_origin_that_fitted_nothing_writes_no_slipped_reads_key`,
`a_file_with_every_table_empty_writes_and_reads_back`,
`a_count_no_toml_integer_can_hold_still_writes_a_document_every_reader_accepts`, and a
section-anchored `every_row_of_every_table_is_one_line`.

## 9. What's good

- **`the_hand_written_words_are_serdes_words`** — the writer spells the enums itself and the test
  compares the two lists, so a rename cannot move the writer, the golden file and the reader
  together. The TOML pass confirmed no key in the writer is ever data.
- **The escape guard is exactly the specification's must-escape set** — `< 0x20 || == 0x7f`,
  checked character by character against the fetched spec.
- **Every float tested parsed as a float and never as an integer**, on all fifteen adversarial
  values including `5e-324`, `f64::MAX` and `-0.0`; and non-ASCII and C1 controls round-tripped
  unchanged.
- **Two thirds of the artefact reads like a lab record** — the reader's own words. Marking a sample
  as less inbred is a true one-line edit.
- **A duplicate row becomes a second array element, never a second key**, so the reader cannot
  silently merge two rows.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo fmt --check`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`
- `./scripts/dev.sh cargo test --lib ng::calling::parameters_file`
- `./scripts/dev.sh cargo test --lib`

Audit trail: the three per-category files in the gitignored
`tmp/review_2026-08-28_parameters-file-b2/`.
