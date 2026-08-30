# ng parameters file — C1: TOML text back into the file's shape

**Date:** 2026-08-30
**Plan:** [parameters_file.md](../../ng/impl_plan/parameters_file.md), Milestone C, step C1
**Spec:** [parameters_file.md](../../ng/spec/parameters_file.md) §4, §9
**Code:** [src/ng/calling/parameters_file/from_toml.rs](../../../../src/ng/calling/parameters_file/from_toml.rs) (new), and [mod.rs](../../../../src/ng/calling/parameters_file/mod.rs)

---

## 1. Plan

"TOML text → the file shape. Parsing, with a malformed file failing at a line number." Spec §9
promises that a malformed file "fails at read with a line number, **which is what using an
existing parser buys**", so the step is a wrapper over `toml::from_str` rather than a parser.

Milestone B had already settled that the reader can be `serde`'s: `the_written_text_reads_back_as_
the_same_file` (`to_toml.rs`) proves the hand-written writer's output round-trips through the
derived `Deserialize`, and `the_documented_inline_form_parses` (`mod.rs`) pins the exact inline
form the writer emits — a form the derived *serializer* never produces, so a serde round trip
alone cannot see it.

## 2. What the step actually builds

**The failing line, as data.** The `toml` crate renders "TOML parse error at line N, column M" in
its own prose, but the only position it exposes is a byte range (`toml::de::Error::span`); the
document text it rendered against, the key path and the line itself are all private. A caller that
wants to *use* the position — put it in a run's own error shape, point an editor at it — has to
re-derive it, and `from_toml` is the one moment the span and the text are in the same place.

So: `ParametersFile::from_toml(&str) -> Result<Self, ParametersFileError>`, where the error
carries `line: Option<usize>` beside the parser's own error as its source.

**The arithmetic is deliberately the crate's own.** `the_line_number_at` reproduces
`translate_position` (`toml-1.1.2`, `de/error.rs`), clamp included, so that this module's line and
the crate's rendered line are two derivations of one span that are meant to agree. Every failure
test asserts both. Three review agents independently proved the equivalence rather than taking the
comment's word for it: the crate counts newlines in `input[0..line_start]`, and no newline lies
between `line_start` and `safe_index`, so that count is the newline count in `input[0..safe_index]`
— which is what this module counts directly.

## 3. Measured, not assumed

- **A refusal raised after the document parsed still carries a position.** A mistyped key is not a
  syntax error — `deny_unknown_fields` fires once the document has parsed, from inside `serde`,
  whose `Error::custom` takes no span at all. The position survives because `Error::set_span` puts
  one back on the way out of every key and value that arrived without one.
- **A *missing* key has no position of its own**, so the crate reports it against whatever should
  have held it. Three grains, all measured: an inline row's key gives the row; a named section's
  key gives that section's **header** (span 4681..4702, exactly the 21 bytes of
  `[ordinary_site_prior]`); and a top-level key gives a zero-width span at byte 0, so line 1 —
  which in a file this writer produces is a comment. Spec §9's promise is kept in all three; the
  weakest one is written down rather than papered over.
- **A value spanning several lines is reported at the line it opens on.** Every failure the
  writer's own layout can produce spans one line, so this is the only fixture that says which end
  of the range the reported line comes from — and without it, reading `span.end` passes every
  other test in the file.
- **105 cuts at 97-byte intervals through the written file were every one refused with a line**,
  which is the sweep standing behind the claim that the `None` arm is unreachable in practice.
- **`cargo doc --no-deps` reports 25 unresolved links, the pre-existing baseline**, so the new
  public docs add none.

## 4. Two additions beyond the plan's line, both recorded

- **A `\r\n` test.** Spec §1.2 goal 3 is a person changing one line, and which machine they did it
  on is not the format's business. Two tests: the file reads back equal, and — the half that is
  this module's rather than the parser's — a failure in such a file names the same line. The
  second exists because the line arithmetic counts bytes: anything that normalised `\r\n` before
  parsing while leaving this counting against the original would shift every reported line, and
  the success test alone would stay green.
- **`deny_unknown_fields` on `LevelSmoothing` and `ShareSmoothing`**, the module's only two shapes
  whose variants carry fields. Raised by review as a hole in Milestone A's shape. **It is not a
  hole** — see §6.

## 5. Deliberately not done, and raised at Checkpoint C instead

**The parser's vocabulary is not translated into the file's.** For four of the five edits a person
is likely to make, the crate's diagnostics are good: a mistyped key lists the keys expected, a
warrant outside the four lists the four. Two kinds of message are not:

- a scalar of the wrong type names a Rust type — `invalid type: string "two", expected u8`;
- an `observations` table with no entry or with two reports `wanted exactly 1 element`, naming
  neither `observations` nor any of the three units.

Both are legible to a programmer and not to the geneticist goal 3 describes, and the second sits
in the path the file's own header comment invites ("delete its `observations`"). Re-wording them
turns this module from a wrapper over a parser into a translator of one, which is a decision the
owner has not been asked. Raised at Checkpoint C.

## 6. One review finding that did not reproduce

The correctness pass reported that `LevelSmoothing` and `ShareSmoothing` silently discard an
unknown key inside a `blend` or `this_periods_curve` table, contradicting the module's stated
invariant, and traced it to `serde_derive` reading `deny_unknown_fields` from the container.

**Measured, on both layouts, with the attribute and without: the key is refused either way**, with
`unexpected keys in table: note, available keys: curve_weight, curve, reach`. The `toml` crate's
own table deserialiser rejects leftover keys before `serde`'s attribute would matter. Removing
either attribute is therefore a mutation no input distinguishes.

The attributes were kept anyway, and the doc comment now says the true thing: they cost nothing
and make the module's opening claim a property of the shape rather than of one parser's
behaviour. The test added alongside them pins the **behaviour a hand-edited file meets**, not the
attribute, and its comment says so.

## 7. Changes made

- **New** `from_toml.rs`: `ParametersFile::from_toml`, `the_line_number_at`, and 15 tests.
- **`mod.rs`**: `ParametersFileError` (moved here from the reader on review, because two later
  steps in two other files add variants to it), `#[non_exhaustive]`, `line()` and
  `rendered_by_the_parser()`; `deny_unknown_fields` on the two enums; the stale "What this module
  is not, yet — no reading and no writing" section replaced; one part added to
  `a_mistyped_key_is_refused_rather_than_absorbed`.

## 8. Tests

**63 → 78 in the module.** 15 new, of which 6 assert a failing line through both derivations.

## 9. Validation

All in the dev container.

- `cargo test --lib ng::calling::parameters_file` — **78 passed, 0 failed, 2 ignored**.
- `cargo test --all-features --lib --tests` — **4,996 lib tests plus every integration binary, 0
  failed**.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo fmt --check` — clean.
- `cargo doc --no-deps` — 25 unresolved links, all pre-existing elsewhere.
- `cargo test --all-targets --all-features` still exits on the pre-existing panic in
  `benches/psp_writer_perf.rs:386`, unrelated to this feature and recorded in `PROJECT_STATUS.md`.

## 10. Mutation testing

Eight mutations; **six fail a test and two are equivalent**.

| mutation | tests failed |
|---|---|
| drop the `+ 1` from the line count | 11 |
| drop the clamp | 1 |
| `[..within]` → `[..=within]` | 2 |
| `line` always `None` | 11 |
| `span.start` → `span.end` | 1 |
| `rendered_by_the_parser` returns `""` | 10 |
| remove `deny_unknown_fields` from `LevelSmoothing` | **0 — equivalent, see §6** |
| remove `deny_unknown_fields` from `ShareSmoothing` | **0 — equivalent, see §6** |

Two of these survived a first pass and drove the two fixtures that now kill them: the multi-line
value (for `span.end`) and the newline-byte case (for the inclusive slice). Neither was reachable
by any fixture the step had before.
