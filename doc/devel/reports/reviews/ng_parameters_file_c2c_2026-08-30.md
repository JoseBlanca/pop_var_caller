# Review — the two parser messages written for a programmer

*2026-08-30. Two agents in isolated worktrees at `c186097d`, each handed the step's diff as a
patch: a correctness pass and a pass by the geneticist who makes the edits. **1 Blocker, 7 Majors
and 6 Minors applied.** Module tests 120 → 126.*

## The Blocker: the file's own text could write its own error message

The translation searched `toml::de::Error`'s **rendered** output for
`invalid type: X, expected Y`. That rendering writes the offending line of the user's document
above the message — so the search read the file back. A file carrying the words
`# invalid type: string` in a trailing comment, which is exactly what a person pasting the last
error in while they edit produces, took its clause from its own comment rather than from the
failure: `reference_digest = 12 # invalid type: string` was reported as *text where text belongs*,
for a value that is a number. The same mechanism suppressed correct clauses, where a sample name
held the words `, expected `.

`toml::de::Error::message()` returns the message alone, and that is what the translation reads
now. A test whose replacement line carries the trigger text in a comment is what would notice.

## Two messages that were the wrong shape

**The unit was judged after the count.** `observations = { lanes = "four" }` reported the *type of
the number* beside a key this file has never had, because the count was read before the unit name
was matched. The unit is matched first now, so a table naming a unit the file does not count in is
told that whatever type its number is.

**The empty-`observations` message gave half the fix.** It said to delete the key rather than
empty it, and stopped there — and following it literally leaves a number somebody typed still
claiming `fitted_here`, which validates. The file's own header gives both halves; so does the
sibling refusal in `validate`. This one now does too.

## The clause named a Rust type's rule and not the reader's line

Three findings, one rewrite.

- **It had no line number**, while spec §9's promise is that a malformed file fails at read with
  one and `ParametersFileError::line()` already held it. The clause now opens with the line.
- **It said "this key"**, and a wrongly-typed *element* of a list — a sample name written as a
  number — sits on a line with no key on it, so the sentence would have been contradicted by the
  caret under it. It names the line and the two kinds, and never a key.
- **It advertised a rule it was never shown for.** "A whole number, and not a negative one" was
  the message for ten integer types — and typing a negative does not produce that message at all:
  `serde` raises `invalid value` there, not `invalid type`, so the search missed the one edit the
  clause was written for. Both shapes are translated now, and the second quotes the number back:
  *line 19 holds -1, and this key takes a whole number, zero or more*.

The word table also carried five arms nothing can reach and five for signed integers, of which the
shape has none — true today and a false statement told confidently to a geneticist the day an
`i32` field is added. It now lists only the types this file's shape uses, so a field at any other
type loses the clause rather than being described wrongly.

## Five doc claims that were false

The module header still said, sixty lines above the code that does it, that the module **does
not** translate the parser's vocabulary and that the decision had not been put to the owner. Three
comments quoted the message the derive used to give as *invalid length 0, expected exactly 1
element*, which is `serde_json`'s phrasing: the `toml` crate's own table deserialiser refuses an
empty table before `serde` sees the enum, with *wanted exactly 1 element, found 0 elements* — and
the pre-existing header two screens above quoted it correctly. One comment named a visitor
`EvidenceCountFromATable` that does not exist. And two blamed a `toml` upgrade for a message both
of whose halves are `serde`'s.

## Two the reader found in the produced file

The file's header tells you that a number you change needs `warrant = "supplied"`, and **two keys
refuse that**: the tract ladder's fallback concentration, whose warrant the file's own strata
decide, and the outlier weight, which is `defaulted` only at the caller's constant. The header now
says so. And a refusal rendered "holds 1 fitted stratum spectra".

## What the correctness pass confirmed

**The hand-written reader accepts everything the derive accepted.** Traced input by input through
`toml-1.1.2` and `serde_core-1.0.228`: both golden layouts, a negative count, a float count, a
string count, a count past `u64`, `observations` given as a number, a list, or a string. Every one
refuses in both, and the only divergence is the wording on a TOML datetime, which nobody writes
there.

## What is still the parser's, and met on three of the five likely edits

A mistyped key is *unknown field*; a warrant outside the four is *unknown variant*; a wrongly-typed
scalar still carries `expected u8` **under** the new clause, because the caret block is worth more
than the purity. "Field" and "variant" are not this file's words for a key and a value. Left alone
by the ruling's scope, and recorded here rather than fixed.
