# Review — ng parameters file, C1: TOML text back into the file's shape

**Date:** 2026-08-30
**Scope:** the working-tree diff of step C1 — `from_toml.rs` (new) and the `mod.rs` doc section.
**Verdict:** Request changes — **0 Blockers, 6 Majors, 15 Minors**, from four agents in isolated
worktrees at the base commit `30ab34a7` with the step applied as a patch.

The four passes: correctness; test discrimination and design fidelity; Rust API design and reuse;
and one agent whose only job was to read the error messages as the geneticist who hand-edited the
file, and to fact-check every claim in the new prose.

---

## 1. What the review confirmed rather than found

Worth recording, because two of these were claims the step's own prose made without proof:

- **`the_line_number_at` reproduces the crate's `translate_position` exactly.** Three agents proved
  the equivalence independently rather than sampling it, including at the clamp and the empty-input
  case.
- **There is no supported alternative.** `toml::de::Error` exposes `message()`, `span()` and
  `set_input()`; the input text, the key path and the rendered line are private, and `line_col()`
  has been gone since toml 0.5. The duplication is forced.
- **No panic path.** The slice is `&[u8]`, so a non-char-boundary offset cannot panic and cannot
  mis-count either — `0x0A` never appears as a UTF-8 continuation byte. `toml_parser` caps nesting
  depth, so a pathological hand-edited file is a parse error rather than a stack overflow.
- **Spans are document-relative**, re-based by `ParseError::rebase_spans`, so the text `from_toml`
  counts against is the text the span indexes.
- **The CRLF fixture is sound**: neither golden file contains a multi-line string or a backslash,
  so every `\n` in them is a terminator.

## 2. The Majors

**Ma1 — the mechanism paragraph explained the wrong mechanism.** The header said the position
survives because "the crate re-attaches the document to every error leaving the top-level
deserializer". *Every* is true, but `set_input` attaches the **text**, not the position: it is what
lets the crate *render* "at line N". What preserves the span this module reads is a different call
in different files — `Error::set_span`, in `de/deserializer/table.rs` and at six sites in
`value.rs`. The paragraph answered "serde's constructor supplies no span" with the one mechanism
that also supplies no span. Rewritten to name both and say which does what.

**Ma2 — "Every test here asserts both derivations" was false.** Four of the ten tests asserted no
line at all, and one of those four — the arithmetic test — is precisely the single-oracle case the
paragraph claimed the file did not contain. Scoped to the failure tests, with the arithmetic test
named as the deliberate exception and why it cannot consult the crate.

**Ma3 — a deleted top-level key points a caret at a comment.** Measured: the span is `0..0`, so the
line is 1, and line 1 of a written file is B3's opening note. The user is told which key is missing
but shown a caret under a `#`. Documented rather than fixed; see Ma6.

**Ma4 — Rust type names leak into the message.** `expected u8`, `expected f64`, `expected u64`.
None is a word in this file's vocabulary or in the spec.

**Ma5 — `observations` fails with a message naming nothing.** An empty table gives `wanted exactly
1 element, found 0 elements`; two entries give `more than 1 element`. Neither names `observations`
nor any of `reads` / `covered_positions` / `bases_compared`, and "element" is a word for serde's
enum representation. This is the path the file's own header comment invites — *"delete its
`observations`"* — where a user who empties the braces rather than deleting the key lands.

**Ma6 — the two above are one decision the owner has not been asked**, and it is raised at
Checkpoint C rather than settled in this step: whether this module should re-word the parser's
diagnostics into the file's vocabulary. Doing so turns a wrapper over a parser into a translator of
one. Note the contrast the reader pass drew: a *wrong unit name* already gives a good message
(`unknown variant \`read_count\`, expected one of \`reads\`, \`covered_positions\`,
\`bases_compared\``), so the gap is narrow and specific rather than general.

## 3. The Minors, grouped

**Applied — error design.** `#[non_exhaustive]` (the house pattern, five ng precedents); explicit
`#[source]` with a doc; the error moved to `mod.rs`, because two later steps in two other files add
variants and would otherwise reach sideways into the reader; `line()` kept but its doc rewritten to
say why it exists beside the public field.

**Applied — the double-print and the trailing blank line.** `#[error("…\n{source}")]` with a field
named `source` both interpolates the parser's whole diagnostic *and* exposes it as the chain's next
link, so any anyhow-style reporter prints it twice; and the crate's `Display` ends with `writeln!`,
so the wrapper's own message ended in a blank line. This crate survives the first only by accident
— `format_error_chain` skips a level whose message is a substring of the previous one. Replaced
with the house pattern for this exact type, `SampleSummaryError::ParseToml`: name the failure,
leave the detail to the source. A `rendered_by_the_parser()` accessor carries the parser's own
rendering for callers and tests that want it.

**Applied — test coverage.** Three gaps, each closed with one fixture: no failure fixture reached
end-of-input, so the clamp was pinned only against hand-written expectations; the CRLF test covered
only the success path, which is the dependency's behaviour rather than this module's; and two of the
five cases the multi-line test cited as measured had no fixture in the file. Also an unknown
*section* and an element of the wrong type, which the prose claimed and nothing exercised.

**Applied — prose.** The span of a missing top-level key is `0..0`, not "the whole document"; a
missing key in a named section spans that section's **header**, not the table; "every failure this
shape can produce spans a single line" is false for a hand-edited file and was narrowed to the
writer's own layout; goal 3 is reading and changing one line, not authoring a file from scratch;
the heading said "Three files" where the module has four; and the crate-message oracle trips on a
reword as well as on changed arithmetic, which the header now says.

**Applied — style.** `the_line_holding` renamed `the_line_number_at` (it returns a number, and the
old name read as though it returned the line's text); the helper simplified from nine lines to
three via `saturating_sub`; and the file now states its own naming rule, as both siblings do.

**Not applied — one finding that did not reproduce.** See the implementation report §6: the
reported `deny_unknown_fields` hole in `LevelSmoothing` and `ShareSmoothing` does not exist,
because the `toml` crate refuses leftover keys in a struct-variant table on its own. Measured on
both layouts, with the attribute and without. The attributes were added anyway as a property of the
shape rather than of the parser, and both the doc and the test comment say plainly that they change
no behaviour today.

## 4. Carried forward, not defects in C1

- **`ParametersFileError` will need an `Io` variant** at the first entry point that takes a path
  rather than text. Adding it later is a breaking change to any `match`, which is what
  `#[non_exhaustive]` now covers.
- **`format_version` is not checked at read.** The nearest precedent, `SampleSummary::
  from_toml_bytes`, refuses both `0` and any future version inline. Confirm C2 owns this.
- **⚑ The plan's D2 brief cannot be followed literally.** It says the three refusals must fail
  "naming the field and the two values that differ, **in the shape the census's own refusal
  uses**", and spec §9 says the same. The census's actual refusal is `Freshness::{Rebuild,
  Refused}(&'static str)` — a verdict enum carrying **only a field name**; `freshness` holds both
  differing values and discards them. So the spec describes a shape the census does not have, and
  D2 must either exceed the census (carry both values, which is what §6 and §13 test 4 want) or
  match it. **Worth settling before D2 is coded.**

## 5. Mutation testing

Eight mutations, six killed and two proved equivalent. The table is in the implementation report
§10. Two of the six survived a first pass and drove the two fixtures that now kill them.
