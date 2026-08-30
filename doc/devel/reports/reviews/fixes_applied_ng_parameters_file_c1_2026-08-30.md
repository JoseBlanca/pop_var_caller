# Fixes applied — ng parameters file, C1

**Date:** 2026-08-30
**Review:** [ng_parameters_file_c1_2026-08-30.md](ng_parameters_file_c1_2026-08-30.md)
**Outcome:** every Major and every Minor either applied or answered. **One finding was measured and
did not reproduce; one decision is escalated to Checkpoint C rather than taken here.**

---

## Applied

**Code**

- `ParametersFileError` moved from `from_toml.rs` to `mod.rs`, `#[non_exhaustive]`, and the
  parser's error carried as an explicit `#[source]` **without being interpolated into `Display`** —
  the house pattern for this exact type (`SampleSummaryError::ParseToml`). This removes both the
  double-rendering under a chain-walking reporter and the trailing blank line the crate's
  `writeln!`-terminated `Display` produced.
- `rendered_by_the_parser()` added, so a caller that wants the parser's own diagnostic without
  walking the chain has it, and the tests have an oracle independent of this module's arithmetic.
- `line()` kept, with its doc rewritten to say why it exists beside the public field.
- `the_line_holding` → `the_line_number_at`, and simplified from nine lines to three.
- `deny_unknown_fields` added to `LevelSmoothing` and `ShareSmoothing` — see "Answered" below for
  what it does and does not do.
- The file states its own naming rule, as both siblings do.

**Tests — 10 became 15 in the new file, and one existing test grew a part**

- a failure at end-of-input, which is the only fixture that puts the clamp through both derivations
  rather than through hand-written expectations;
- a failure inside a `\r\n` file, which is the half of CRLF handling that is this module's rather
  than the parser's;
- an unknown *section*, and an element of the wrong type inside an array — both were cited in the
  prose as measured and neither had a fixture;
- a sweep of 105 truncations of the written file, all refused and all naming a line, standing
  behind the claim that the `None` arm is unreachable in practice;
- the missing-key test now covers all three grains — an inline row, a named section's header, and
  the top level — with the section's header line found rather than counted back from the key's,
  because B3's note sits between them;
- `a_mistyped_key_is_refused_rather_than_absorbed` (`mod.rs`) grew a part for a key inside an
  enum's struct variant, the one shape its two existing parts did not reach.

**Prose** — six wrong claims corrected, listed in the review's §3. The one that mattered was the
mechanism (review Ma1): the header explained the position's survival by `set_input`, which
attaches the document text and no position, where the actual mechanism is `set_span` in two other
files.

## Answered, not applied

**The `deny_unknown_fields` hole does not exist.** The correctness pass reported that an unknown
key inside a `blend` or `this_periods_curve` table is silently discarded. Measured on both layouts,
with the attribute and without, the key is refused either way — `unexpected keys in table: note,
available keys: curve_weight, curve, reach` — because the `toml` crate's own table deserialiser
rejects leftover keys before serde's attribute would matter. The attributes were kept as a property
of the shape rather than of one parser, and both the doc comment and the test comment now say
plainly that they change no behaviour today. Removing either is an equivalent mutation, and the
mutation table records it as such rather than as a surviving one.

## Escalated to Checkpoint C

**Whether this module should re-word the parser's diagnostics into the file's vocabulary** (review
Ma3–Ma6). Two message shapes are legible to a programmer and not to the geneticist spec §1.2 goal 3
describes — a scalar of the wrong type names a Rust type, and an `observations` table with no entry
or two reports "wanted exactly 1 element" — and a deleted top-level key points its caret at a
comment. Fixing them turns this module from a wrapper over a parser into a translator of one, which
is a design decision rather than a coding one, and spec §4's stated reason for choosing an existing
parser is what it buys in diagnostics.

## Validation after the fixes

- `cargo test --lib ng::calling::parameters_file` — **78 passed, 0 failed, 2 ignored**.
- `cargo test --all-features --lib --tests` — **4,996 lib tests plus every integration binary, 0
  failed**.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean. `cargo fmt --check` — clean.
- `cargo doc --no-deps` — 25 unresolved links, the pre-existing baseline unchanged.
- Mutation sweep re-run on the fixed tree: six of eight fail a test, two proved equivalent.
