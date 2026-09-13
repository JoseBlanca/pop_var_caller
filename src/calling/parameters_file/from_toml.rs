//! **The file's text, read back into its shape** — and spec §9's promise that "a malformed file
//! fails at read with a line number".
//!
//! # Why this is thin
//!
//! The shape derives `Deserialize` — every type but [`EvidenceCount`], whose reader is written by
//! hand for one message's sake (below) — and the writer's own text already parses through it:
//! `to_toml`'s `the_written_text_reads_back_as_the_same_file` writes the whole fixture and reads
//! it back equal, and the parent module's `the_documented_inline_form_parses` pins that the
//! one-row-a-line inline form the hand-written writer emits is a form the derived reader accepts.
//! So the parsing was settled in step B, and nothing here re-implements it.
//!
//! **What this module adds is the one thing the `toml` crate does not hand back as data: which
//! line failed.** The crate renders a line number in its own prose, but the only position it
//! exposes is a byte range ([`toml::de::Error::span`]) — the document text it rendered against,
//! the key path, and the line itself are all private. A caller that wants to *do* something with
//! the position — put it in a run's own error shape, point an editor at it — has to re-derive it
//! from that span against the text it still holds. That is
//! [`ParametersFileError::line`](super::ParametersFileError::line), and this is the one moment
//! the span and the text are in the same place.
//!
//! # Measured: a refusal raised *after* the document parsed still carries a position
//!
//! A mistyped key is not a syntax error. The document parses, and `deny_unknown_fields` refuses
//! it afterwards from inside `serde`, whose own `Error::custom` constructor takes no span at all
//! (`toml-1.1.2`, `de/error.rs`). **The position survives because the crate puts one back:
//! `Error::set_span` is called on the way out of every key and every value that arrived without
//! one** (`de/deserializer/table.rs`, in `next_key_seed` and its sibling; `de/deserializer/value.rs`,
//! at six call sites), each supplying the span of the key or the value it was reading.
//!
//! **A different call, `set_input`, attaches the document text**, and it is easy to mistake for
//! the mechanism above because it is the more visible one — it is what lets the crate *render*
//! "at line N" in its own message. It supplies no position. Without `set_span` the line this
//! module reports would be `None` for every mistyped key, and the crate's own message would lose
//! its line as well; without `set_input` only the crate's prose would suffer. The two are
//! repaired in different files, and only the first is what this module rests on.
//!
//! # What this module does not do, and where its diagnostics stop
//!
//! **It does not judge what the file says.** A document that parses and means nothing — an
//! inbreeding coefficient of 1.7, a length spectrum whose shares do not sum to one, a
//! contamination table in which no row was measured — is accepted here and refused by step C2's
//! `validate`, which runs after this and before the projection back to `RunParameters`.
//!
//! **It translates two of the parser's messages into the file's vocabulary, and only two.** Spec
//! §4 chose an existing parser partly for what its diagnostics buy, and for most of the edits a
//! person is likely to make they are good: a mistyped key lists the keys that were expected, a
//! warrant outside the four lists the four, a mistyped unit lists the three units. **Two were
//! not.** A value the key cannot hold named a Rust type — `invalid type: string "two", expected
//! u8` — where `u8`, `f64` and `u64` appear nowhere in this file, its comments or its spec; and
//! an `observations` table left empty or given two entries reported *wanted exactly 1 element*,
//! naming neither the key nor any of the three units, in the very path the file's own header
//! invites. Owner's ruling of 2026-08-30, after C1's reader pass raised both.
//!
//! **The first is a clause added to this module's own message**, with the parser's diagnostic —
//! line, column, offending line and caret — kept underneath as the error's source
//! ([`a_value_the_key_cannot_hold_in_the_files_words`]). **The second is a hand-written reader**
//! for [`EvidenceCount`] alone. Everything else the crate says is left exactly as it is.
//!
//! # One naming rule
//!
//! **A helper whose name starts with an article returns a value; a helper named for a claim
//! asserts it and returns nothing.** So `the_line_number_at` hands back a number, and
//! `the_failure_is_on_line` is a claim about an error and panics when it is false. The sibling
//! files each state their own rule and none of them is the crate's (`to_toml.rs` header).

use super::{EvidenceCount, ParametersFile, ParametersFileError};
use serde::Deserialize;
use serde::de::{self, Deserializer, MapAccess, Visitor};

impl ParametersFile {
    /// **The file, from text** — the reverse of [`ParametersFile::to_toml`].
    ///
    /// Accepts any TOML that spells this shape, not merely the layout the writer emits. Spec
    /// §1.2 goal 3 is a person reading the file their run wrote and changing one line of it, so
    /// the reader cannot be tied to the writer's own spacing.
    ///
    /// # Errors
    ///
    /// [`ParametersFileError::Malformed`] when the text is not TOML, or is TOML this shape does
    /// not accept — an unknown key, a missing one, or a value of the wrong type.
    ///
    /// **Every failure names a line, and for all but one that line is the failure's own.** A
    /// *missing* key has no position of its own, so the crate reports it against the thing that
    /// should have held it: for a key of an inline row, that row; for a key of a named section,
    /// that section's header line; and for a top-level key, a zero-width position at the first
    /// byte of the document, which is line 1. In a file this writer produced line 1 is a comment,
    /// so that one case names a line a reader cannot act on. See
    /// `a_missing_key_is_reported_at_the_thing_that_should_have_held_it`.
    pub fn from_toml(text: &str) -> Result<Self, ParametersFileError> {
        toml::from_str(text).map_err(|source| {
            let line = source
                .span()
                .map(|span| the_line_number_at(text, span.start));
            ParametersFileError::Malformed {
                line,
                in_the_files_words: a_value_the_key_cannot_hold_in_the_files_words(
                    source.message(),
                    line,
                ),
                source,
            }
        })
    }
}

/// **A value the key cannot hold, said in the file's words** — and `None` for anything else.
///
/// `serde` raises two messages for this, and the difference is not one a person editing a file
/// would draw: `invalid type: string "two", expected u8` where the *kind* is wrong, and
/// ``invalid value: integer `-1`, expected u8`` where the kind is right and the number is not one
/// the key can hold. Both end in a **Rust type** — `u8`, `f64`, `u64` — and none of those appears
/// in this file, in its comments, or in its spec. Owner's ruling of 2026-08-30 to re-word this
/// shape, and this shape only.
///
/// **It takes the message alone, never the rendering.** `toml::de::Error`'s `Display` writes the
/// offending line of the *document* above the message, so a search over the rendered text reads
/// the user's own file back — and a file with `invalid type: string` in a trailing comment, which
/// is what a person pasting the last error in while they edit produces, would have had a clause
/// fabricated out of its own text.
///
/// **It names the line and never the key.** The failing key is not in `serde`'s message, and a
/// wrongly-typed *element* of a list — a sample name written as a number — sits on a line with no
/// key on it at all, so a sentence saying *this key* would be contradicted by the caret under it.
///
/// **Every unrecognised half gives `None` for the whole sentence**, so a message this function
/// has not met reaches the reader as `serde`'s own rather than as a half-translation. Both halves
/// are `serde`'s rather than the `toml` crate's — the first from `Error::invalid_type`, the second
/// from each visitor's `expecting` — so it is a `serde` upgrade that could reword them, and one
/// that did would drop the clause rather than mistranslate it.
fn a_value_the_key_cannot_hold_in_the_files_words(
    said: &str,
    line: Option<usize>,
) -> Option<String> {
    let at = match line {
        Some(line) => format!("line {line} holds"),
        None => "one of this file's values is".to_owned(),
    };
    if let Some(after) = said.strip_prefix("invalid type: ") {
        let (found, wanted) = after.split_once(", expected ")?;
        return Some(format!(
            "{at} {}, where {} belongs",
            a_word_for_what_was_there(found.trim())?,
            a_word_for_what_a_key_holds(wanted.trim())?
        ));
    }
    // **The kind is right and the number is not** — a negative in a count, or one past what the
    // key's width holds. `serde` writes the offending value in backticks here, and quoting it
    // back is worth more than naming its kind: the reader typed it.
    let after = said.strip_prefix("invalid value: integer `")?;
    let (value, wanted) = after.split_once("`, expected ")?;
    Some(format!(
        "{at} {value}, and this key takes {}",
        a_word_for_what_a_key_holds(wanted.trim())?
    ))
}

/// The file's word for the kind of thing the parser found, from `serde`'s own phrase for it.
fn a_word_for_what_was_there(found: &str) -> Option<&'static str> {
    // `serde` writes the value itself after the kind — `string "two"`, ``integer `3` `` — so this
    // matches the opening word and drops the value, which the parser's caret already points at.
    Some(match found.split_whitespace().next()? {
        "string" => "text",
        "integer" => "a whole number",
        "floating" => "a number with a decimal point",
        "boolean" => "true or false",
        "sequence" => "a list",
        "map" => "a table",
        _ => return None,
    })
}

/// The file's word for what the key holds, from the Rust type `serde` names.
///
/// **Only the types this file's shape uses.** Every integer key in it is unsigned — the shape has
/// no signed field at all — which is what makes *zero or more* a true thing to say rather than an
/// accident. A field added at some other type falls through to `None` and keeps `serde`'s own
/// message, rather than being described by a sentence written for a type it does not have.
fn a_word_for_what_a_key_holds(wanted: &str) -> Option<&'static str> {
    Some(match wanted {
        "u8" | "u32" | "u64" => "a whole number, zero or more",
        "f64" => "a number",
        "a string" | "a borrowed string" => "text",
        "a boolean" => "true or false",
        "a sequence" => "a list, in square brackets",
        _ => return None,
    })
}

/// **Which 1-based line of `text` the byte at `offset` falls on.**
///
/// Deliberately the same arithmetic the `toml` crate renders its own message with
/// (`toml-1.1.2`, `de/error.rs`, `translate_position`), including its clamp: an offset at or past
/// the end of the text belongs to the last line rather than to a line that does not exist. The
/// two derivations are meant to agree — see this module's header — so they have to be over one
/// convention, or every test here would be comparing a line against a line-plus-one.
fn the_line_number_at(text: &str, offset: usize) -> usize {
    let within = offset.min(text.len().saturating_sub(1));
    text.as_bytes()[..within]
        .iter()
        .filter(|byte| **byte == b'\n')
        .count()
        + 1
}

/// **`observations = { reads = 812344 }`, read by hand.**
///
/// The derive reads this shape correctly and refuses it in the `toml` crate's words: *wanted
/// exactly 1 element, found 0 elements* for an empty table, and *wanted exactly 1 element, more
/// than 1 element* for one naming two units — raised by that crate's own table deserialiser
/// before `serde` sees the enum at all. Neither message contains the word `observations`, either
/// of the three unit names, or anything a reader of this file has met; "element" is the parser's
/// word for a key of a table. And the empty case sits in the path the file's own header sends a
/// person down. Owner's ruling of 2026-08-30.
///
/// **It refuses exactly what the derive refuses**, in three cases rather than one: no unit, an
/// unknown unit, and two units. The count itself is left to `serde`, whose message for a
/// non-integer there names the key and the value.
impl<'de> Deserialize<'de> for EvidenceCount {
    fn deserialize<D: Deserializer<'de>>(reader: D) -> Result<Self, D::Error> {
        /// The three units this file counts in, in the file's own spelling, for a message that
        /// has to list them.
        const THE_UNITS: &str = "`reads`, `covered_positions` or `bases_compared`";

        struct AUnitAndItsCount;

        impl<'de> Visitor<'de> for AUnitAndItsCount {
            type Value = EvidenceCount;

            fn expecting(&self, out: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(
                    out,
                    "an `observations` table naming one of {THE_UNITS} and how many of them"
                )
            }

            fn visit_map<A: MapAccess<'de>>(self, mut table: A) -> Result<Self::Value, A::Error> {
                let Some(unit) = table.next_key::<String>()? else {
                    return Err(de::Error::custom(format!(
                        "`observations` is empty. It names one of {THE_UNITS}, and how many of \
                         them stood behind the `value` beside it. If you changed that value, \
                         delete the whole `observations` key and set `warrant` to `supplied` — \
                         an empty `{{ }}` still says the number was measured, and no longer says \
                         on what"
                    )));
                };
                // **The unit is matched before the count is read**, so that a table naming a
                // unit this file does not count in is told so — rather than being told the type
                // of the number beside a key the file does not have.
                let counted: fn(u64) -> EvidenceCount = match unit.as_str() {
                    "reads" => EvidenceCount::Reads,
                    "covered_positions" => EvidenceCount::CoveredPositions,
                    "bases_compared" => EvidenceCount::BasesCompared,
                    unknown => {
                        return Err(de::Error::custom(format!(
                            "`observations` counts `{unknown}`, and the units this file counts \
                             in are {THE_UNITS}"
                        )));
                    }
                };
                let counted = counted(table.next_value::<u64>()?);
                if let Some(also) = table.next_key::<String>()? {
                    return Err(de::Error::custom(format!(
                        "`observations` names both `{unit}` and `{also}`, and it names one — \
                         delete the wrong one. The three units differ by orders of magnitude on \
                         the same cohort, so a count under the wrong one compares to nothing"
                    )));
                }
                Ok(counted)
            }
        }

        reader.deserialize_map(AUnitAndItsCount)
    }
}

#[cfg(test)]
mod tests {
    use super::super::tests::a_file_using_every_shape;
    use super::*;

    /// The line the crate's own message names, read back out of its prose.
    ///
    /// **The independent half of every failure assertion here.** It reads the number out of "TOML
    /// parse error at line N, column M", which is the crate's rendering rather than this module's
    /// arithmetic, so a test asserting both is comparing two derivations of one span.
    ///
    /// **It trips on a reworded message as well as on changed arithmetic**, and that is the
    /// intended trade: a `toml` upgrade that moves either is one somebody should look at. It
    /// cannot pass by accident — every degradation it could suffer (the phrase gone, no digits
    /// after it, an error with no span and so no context block at all) yields `None`, and `None`
    /// is compared against an expected line rather than against this module's own answer.
    fn the_line_the_crates_message_names(error: &ParametersFileError) -> Option<usize> {
        let rendered = error.rendered_by_the_parser();
        let after = rendered.split("at line ").nth(1)?;
        let digits: String = after.chars().take_while(char::is_ascii_digit).collect();
        digits.parse().ok()
    }

    /// Both derivations of the failing line, asserted to be `expected` and so to agree.
    #[track_caller]
    pub(super) fn the_failure_is_on_line(error: &ParametersFileError, expected: usize) {
        assert_eq!(
            error.line(),
            Some(expected),
            "this module's own line, from the error's byte span; the parser said:\n{}",
            error.rendered_by_the_parser()
        );
        assert_eq!(
            the_line_the_crates_message_names(error),
            Some(expected),
            "the line the toml crate's own message names; the parser said:\n{}",
            error.rendered_by_the_parser()
        );
    }

    /// The text the writer produces, with the first `find` replaced by `replace_with`.
    ///
    /// Returns the text and the 1-based line the replacement landed on, **counted from the text
    /// itself** and by a different route than [`the_line_number_at`] takes, so a fixture's
    /// expected line is a third derivation rather than an echo of the one under test.
    pub(super) fn the_written_file_with(find: &str, replace_with: &str) -> (String, usize) {
        let text = a_file_using_every_shape().to_toml();
        let at = text
            .find(find)
            .unwrap_or_else(|| panic!("the written file contains {find:?}"));
        let line = text[..at].matches('\n').count() + 1;
        (text.replacen(find, replace_with, 1), line)
    }

    /// The text the writer produces with `suffix` appended, and the line `suffix` starts on.
    fn the_written_file_followed_by(suffix: &str) -> (String, usize) {
        let text = a_file_using_every_shape().to_toml();
        let line = text.matches('\n').count() + 1;
        (format!("{text}{suffix}"), line)
    }

    /// **The file this build checked in reads back as the fixture that produced it.**
    ///
    /// It does not stand alone as evidence that the reader and the file agree — the writer's
    /// output is already pinned byte-for-byte to this file by `to_toml`'s
    /// `the_whole_shape_writes_the_documented_toml`, and that output is already read back by its
    /// `the_written_text_reads_back_as_the_same_file`, so the two together imply this one.
    ///
    /// **What it adds is the only coverage of [`ParametersFile::from_toml`] itself over a whole
    /// file.** Every pre-existing round trip calls `toml::from_str` directly and would not notice
    /// anything this entry point started doing to the text on the way — normalising line endings,
    /// stripping a byte-order mark, checking a version — which is exactly the kind of thing a
    /// later step is tempted to add here.
    #[test]
    fn the_checked_in_file_reads_back_as_the_fixture_that_wrote_it() {
        let read = ParametersFile::from_toml(include_str!("testdata/every_shape_as_written.toml"))
            .expect("the checked-in file is a parameters file");
        assert_eq!(read, a_file_using_every_shape());
    }

    /// **The same fixture in serde's own layout reads back to the same value.**
    ///
    /// The two golden files are one fixture written by two writers — array-of-table headers
    /// against one-row-a-line inline tables — so reading both is what says the reader accepts
    /// *the shape* rather than one writer's rendering of it. Like the test above, it is implied
    /// by a pre-existing pair and earns its place by going through `from_toml`.
    ///
    /// **Neither test can see a renamed key.** Both golden files are regenerated from the shape
    /// by the two ignored helpers, so a rename plus a regeneration moves the file and the reader
    /// together. `every_enum_variant_spells_as_the_file_says` is what holds the key surface.
    #[test]
    fn the_same_file_in_serdes_layout_reads_back_to_the_same_value() {
        let read = ParametersFile::from_toml(include_str!("testdata/every_shape.toml"))
            .expect("serde's own rendering of the fixture is a parameters file");
        assert_eq!(read, a_file_using_every_shape());
    }

    /// A file that is not TOML at all fails on the line that is not TOML.
    #[test]
    fn a_line_that_is_not_toml_fails_at_that_line() {
        let text = "format_version = 1\nploidy = 2\nthis line is not toml\n";
        let error = ParametersFile::from_toml(text).expect_err("that line is not TOML");
        the_failure_is_on_line(&error, 3);
    }

    /// **A key nobody knows fails at the line the typo is on, not at the top of the file.**
    ///
    /// The refusal comes from `deny_unknown_fields` after the whole document has parsed, which is
    /// where a position is easiest to lose — see this module's header.
    #[test]
    fn a_mistyped_key_fails_at_the_line_the_typo_is_on() {
        let (text, line) = the_written_file_with("reference_digest", "refrence_digest");
        let error = ParametersFile::from_toml(&text).expect_err("no such key");
        the_failure_is_on_line(&error, line);
    }

    /// A section nobody knows fails at its own header line.
    #[test]
    fn a_section_nobody_knows_fails_at_its_header() {
        let (text, line) = the_written_file_followed_by("\n[nonsense]\nanything = 1\n");
        let error = ParametersFile::from_toml(&text).expect_err("no such section");
        the_failure_is_on_line(&error, line + 1);
    }

    /// **A missing key is reported against the thing that should have held it**, which is the one
    /// failure whose line is weaker than the rest.
    ///
    /// An absence has no position, so the crate spans whatever should have contained it. Three
    /// cases, all measured: for a key of an inline row that is the row, and the line is as good
    /// as any other failure's; for a key of a named section it is the section's **header** and
    /// not the whole table; and **for a top-level key it is a zero-width position at the first
    /// byte**, which is line 1 — a comment, in a file this writer produced. Spec §9's promise is
    /// kept in all three; this is what it is worth in the weakest one, written down so that the
    /// other tests here are not read as saying more than they do.
    #[test]
    fn a_missing_key_is_reported_at_the_thing_that_should_have_held_it() {
        let (top_level, _) = the_written_file_with("ploidy = 2\n", "");
        let error = ParametersFile::from_toml(&top_level).expect_err("ploidy is not optional");
        the_failure_is_on_line(&error, 1);

        let (a_row, line) = the_written_file_with("warrant = \"borrowed\", ", "");
        let error = ParametersFile::from_toml(&a_row).expect_err("a warrant is not optional");
        the_failure_is_on_line(&error, line);

        // The section's **header**, which is several lines above the key that was removed —
        // B3's note sits between them — so the header's line is found rather than counted back
        // from the key's.
        let (a_section, key_line) = the_written_file_with("rung = \"fitted_curve\"\n", "");
        let header = a_section
            .find("[ordinary_site_prior]")
            .expect("the written file has that section");
        let header_line = a_section[..header].matches('\n').count() + 1;
        assert!(
            header_line < key_line - 1,
            "the point of this case is that the header is not the line above the key: header \
             {header_line}, key {key_line}"
        );
        let error = ParametersFile::from_toml(&a_section).expect_err("a rung is not optional");
        the_failure_is_on_line(&error, header_line);
    }

    /// A scalar of the wrong type fails at the line the value is on.
    #[test]
    fn a_value_of_the_wrong_type_fails_at_the_line_it_is_on() {
        let (text, line) = the_written_file_with("ploidy = 2", "ploidy = \"two\"");
        let error = ParametersFile::from_toml(&text).expect_err("a ploidy is not a string");
        the_failure_is_on_line(&error, line);
    }

    /// An element of the wrong type inside an array fails at that element, not at the array.
    #[test]
    fn an_element_of_the_wrong_type_fails_at_that_element() {
        let (text, line) = the_written_file_with("\"TS-1\",", "1,");
        let error = ParametersFile::from_toml(&text).expect_err("a sample name is not a number");
        the_failure_is_on_line(&error, line);
    }

    /// **A warrant the file does not define is refused**, rather than read as one that is.
    ///
    /// The four warrants are the file's compatibility surface, and an unknown one is the case
    /// where a reader falling back to a default would report a number the run guessed as one it
    /// had fitted.
    #[test]
    fn a_warrant_that_is_not_one_of_the_four_is_refused() {
        let (text, line) = the_written_file_with("\"fitted_here\"", "\"fitted_over_there\"");
        let error = ParametersFile::from_toml(&text).expect_err("there are four warrants");
        the_failure_is_on_line(&error, line);
    }

    /// **A value that spans several lines is reported at the line it opens on**, not the line it
    /// closes on.
    ///
    /// Every failure the *writer's own layout* can produce spans one line — measured, over a
    /// mistyped key, a missing key at each of its three grains, an unknown section, an element of
    /// the wrong type and a scalar of the wrong type. A hand-edited file is not so limited: a
    /// value written across lines, or a row a person split in two, spans several. This is the
    /// fixture that says *which end* of that range the reported line comes from. Without it,
    /// reading the span's end instead of its start passes every other test in this file.
    #[test]
    fn a_value_that_spans_several_lines_is_reported_at_the_line_it_opens_on() {
        let (text, line) =
            the_written_file_with("ploidy = 2", "ploidy = \"\"\"\ntwo\nplants\n\"\"\"");
        let error = ParametersFile::from_toml(&text).expect_err("a ploidy is not a string");
        the_failure_is_on_line(&error, line);
    }

    /// **A failure at the very end of the text is reported at the last line**, which is the clamp
    /// the arithmetic inherits from the crate.
    ///
    /// The only fixture here whose span reaches the end of the document, so it is the only one
    /// that puts the clamp through both derivations rather than through this module's own
    /// hand-written expectations.
    #[test]
    fn a_failure_at_the_end_of_the_text_is_reported_at_the_last_line() {
        let text = "format_version = 1\nploidy = 2\nsamples = [";
        let error = ParametersFile::from_toml(text).expect_err("that array never closes");
        the_failure_is_on_line(&error, 3);
    }

    /// **A file whose lines end `\r\n` reads back the same**, which is a file edited on Windows.
    ///
    /// Goal 3 is a person changing one line, and which machine they did it on is not something
    /// the format gets to care about. TOML's grammar allows either ending; this pins that nothing
    /// in this shape re-introduces the restriction. The fixture is safe because neither golden
    /// file contains a multi-line string or an escape, so every `\n` in it is a terminator.
    #[test]
    fn a_file_edited_on_a_machine_that_ends_lines_differently_reads_the_same() {
        let text = a_file_using_every_shape().to_toml().replace('\n', "\r\n");
        let read = ParametersFile::from_toml(&text).expect("CRLF is a TOML line ending");
        assert_eq!(read, a_file_using_every_shape());
    }

    /// **And a failure in such a file names the same line it would have named otherwise.**
    ///
    /// The success path above is the parser's behaviour rather than this module's. This is the
    /// half that is this module's: the line arithmetic counts bytes, so anything that normalised
    /// `\r\n` to `\n` before parsing while leaving the text this counts against unnormalised
    /// would shift every reported line by the number of carriage returns above it, and the
    /// success test would stay green.
    #[test]
    fn a_failure_in_a_file_that_ends_lines_differently_names_the_same_line() {
        let (text, line) = the_written_file_with("reference_digest", "refrence_digest");
        let error = ParametersFile::from_toml(&text.replace('\n', "\r\n"))
            .expect_err("no such key, whatever ends the lines");
        the_failure_is_on_line(&error, line);
    }

    /// **However the file is truncated, the failure still names a line.**
    ///
    /// Spec §9 promises a line number without qualification, and the line is an `Option` because
    /// the crate's span is. Nothing in this suite reaches the `None` arm and nothing measured on
    /// this reader ever has; this sweeps a whole family of broken files rather than arguing from
    /// the crate's source that no such file exists.
    ///
    /// **105 cuts at 97-byte intervals through the written file, and all 105 were refused with a
    /// line.** The count is asserted loosely rather than exactly, so that a change to the golden
    /// file's length does not fail a test about something else.
    #[test]
    fn however_the_file_is_truncated_the_failure_still_names_a_line() {
        let whole = a_file_using_every_shape().to_toml();
        let mut refused = 0;
        for cut in (0..whole.len()).step_by(97) {
            let Some(text) = whole.get(..cut) else {
                continue; // not a character boundary
            };
            if let Err(error) = ParametersFile::from_toml(text) {
                refused += 1;
                assert!(
                    error.line().is_some(),
                    "a truncation at byte {cut} named no line; the parser said:\n{}",
                    error.rendered_by_the_parser()
                );
            }
        }
        assert!(
            refused > 50,
            "the sweep is only evidence if most truncations are refused, and only {refused} were"
        );
    }

    /// The line arithmetic, at the edges the crate's clamp and its exclusive count exist for.
    ///
    /// **These expectations are hand-derived from the crate's `translate_position` rather than
    /// taken from a live error**, because the two branches that matter — an empty text and an
    /// offset past the end — are not reachable through a document. What holds the arithmetic
    /// against the crate on reachable inputs is every test above that calls
    /// [`the_failure_is_on_line`].
    #[test]
    fn the_line_number_at_an_offset_is_what_the_crate_counts() {
        assert_eq!(the_line_number_at("", 0), 1, "an empty text has one line");
        assert_eq!(the_line_number_at("a\nb\nc", 0), 1);
        assert_eq!(
            the_line_number_at("a\nb", 1),
            1,
            "the newline ending a line belongs to the line it ends, as the crate counts it"
        );
        assert_eq!(
            the_line_number_at("a\nb\nc", 2),
            2,
            "the first byte of line 2"
        );
        assert_eq!(the_line_number_at("a\nb\nc", 4), 3);
        assert_eq!(
            the_line_number_at("a\nb\nc", 99),
            3,
            "an offset past the end belongs to the last line, as the crate clamps it"
        );
    }
}

#[cfg(test)]
mod the_two_messages_written_for_a_programmer {
    use super::super::tests::a_file_using_every_shape;
    use super::tests::{the_failure_is_on_line, the_written_file_with};
    use super::*;

    fn refused(text: &str) -> ParametersFileError {
        ParametersFile::from_toml(text).expect_err("this edit is not a parameters file")
    }

    /// **A value the key cannot hold is reported in the file's words, and the parser's own
    /// diagnostic is still underneath it.**
    ///
    /// `serde` says `invalid type: string "two", expected u8`, and `u8` appears nowhere in this
    /// file, its comments or its spec. Owner's ruling of 2026-08-30: re-word this shape, and this
    /// one only.
    ///
    /// **The line comes out of the fixture rather than being typed here**, and is asserted
    /// through `the_failure_is_on_line`, which compares this module's arithmetic against the
    /// number the crate renders in its own prose — so a line added to the writer's header moves
    /// the expectation with it rather than failing the test.
    #[test]
    fn a_wrong_type_is_reported_in_the_files_own_words() {
        let (text, line) = the_written_file_with("ploidy = 2", "ploidy = \"two\"");
        let error = refused(&text);
        assert_eq!(
            error.to_string(),
            format!(
                "the parameters file could not be read: line {line} holds text, where a whole \
                 number, zero or more belongs"
            )
        );
        the_failure_is_on_line(&error, line);
        assert!(
            error.rendered_by_the_parser().contains("expected u8"),
            "the parser's own diagnostic is kept as the source, caret and all: {}",
            error.rendered_by_the_parser()
        );

        // A number key given text, and a text key given a number: the same shape both ways round.
        let (text, line) = the_written_file_with(
            "reference_concentration = 1.0",
            "reference_concentration = \"one\"",
        );
        let error = refused(&text);
        assert_eq!(
            error.to_string(),
            format!(
                "the parameters file could not be read: line {line} holds text, where a number \
                 belongs"
            )
        );
        let (text, line) = the_written_file_with("reference_digest = ", "reference_digest = 12 #");
        let error = refused(&text);
        assert_eq!(
            error.to_string(),
            format!(
                "the parameters file could not be read: line {line} holds a whole number, where \
                 text belongs"
            )
        );
    }

    /// **A number of the right kind that the key still cannot hold says which number it was.**
    ///
    /// `serde` raises this one as `invalid value` rather than `invalid type` — a different message
    /// for what a person editing a file would call the same mistake — so a translation written for
    /// the first alone would miss the edit most likely to make it: typing a negative, or a count
    /// past the width of its key.
    #[test]
    fn a_number_the_key_cannot_hold_names_the_number() {
        let (text, line) = the_written_file_with("ploidy = 2", "ploidy = -1");
        let error = refused(&text);
        assert_eq!(
            error.to_string(),
            format!(
                "the parameters file could not be read: line {line} holds -1, and this key takes \
                 a whole number, zero or more"
            )
        );
        the_failure_is_on_line(&error, line);

        let (text, line) = the_written_file_with("ploidy = 2", "ploidy = 300");
        assert_eq!(
            refused(&text).to_string(),
            format!(
                "the parameters file could not be read: line {line} holds 300, and this key takes \
                 a whole number, zero or more"
            )
        );
    }

    /// **A message this module has no better word for reaches the reader as the parser's own.**
    ///
    /// The clause is added where it helps and omitted everywhere else, so a `serde` upgrade that
    /// reworded its diagnostics loses the clause rather than mistranslating it — and the shapes
    /// spec §4 chose this parser *for* are left exactly as they are.
    ///
    /// **The last case is the one that matters.** The translation reads the parser's *message*
    /// and never its rendering, because the rendering echoes the offending line of the user's own
    /// file above the message — so a file carrying `invalid type: string` in a comment, which is
    /// what a person pasting the last error back in while they edit produces, would otherwise have
    /// had a clause fabricated out of its own text.
    #[test]
    fn every_other_failure_keeps_the_parsers_own_message() {
        for (find, replace_with) in [
            // A mistyped key, which the parser answers by listing the keys it expected.
            ("ploidy = 2", "ploidyy = 2"),
            // A warrant outside the four, which it answers by listing the four.
            (
                "warrant = \"fitted_here\"",
                "warrant = \"fitted_over_there\"",
            ),
            // A rung outside its own four, which is the same shape one vocabulary along.
            (
                "rung = \"fitted_curve\"",
                "rung = \"a curve somebody fitted\"",
            ),
            // Not TOML at all.
            ("ploidy = 2", "ploidy ="),
            // A file whose own text carries the words the translation looks for. The failure is
            // the mistyped key; the comment is what a person pasting the last error back in
            // produces, and a search over the parser's *rendering* would have read it.
            (
                "ploidy = 2",
                "ploidyy = 2 # invalid type: string \"two\", expected u8",
            ),
        ] {
            let (text, _) = the_written_file_with(find, replace_with);
            let error = refused(&text);
            assert_eq!(
                error.to_string(),
                "the parameters file could not be read",
                "nothing is added to a message this module has no better word for; the parser \
                 said:\n{}",
                error.rendered_by_the_parser()
            );
        }
    }

    /// **An `observations` table left empty names the key, its three units, and both halves of
    /// the fix.**
    ///
    /// The derive answers it with *wanted exactly 1 element, found 0 elements* — no key, no unit,
    /// and "element" is the parser's word for a key of a table. **And this is the path the file's
    /// own header sends a reader down**: it says that a number you change needs its warrant set
    /// to `supplied` and its `observations` deleted, and a reader who empties the braces rather
    /// than deleting the key lands here. So the message says *both* halves — deleting the
    /// observations alone leaves a number somebody typed claiming it was fitted here, which
    /// validates.
    #[test]
    fn an_emptied_observations_table_names_the_key_and_its_units() {
        let (text, _) =
            the_written_file_with("observations = { reads = 812344 }", "observations = {}");
        let said = refused(&text).rendered_by_the_parser();
        assert!(said.contains("`observations` is empty"), "{said}");
        assert!(
            said.contains("`reads`, `covered_positions` or `bases_compared`"),
            "it names the three units the file counts in: {said}"
        );
        assert!(
            said.contains("delete the whole `observations` key")
                && said.contains("set `warrant` to `supplied`"),
            "and both halves of the edit the header's own instruction meant: {said}"
        );
        assert!(
            !said.contains("element"),
            "and it does not speak of elements: {said}"
        );
    }

    /// **A table naming two units, or a unit this file does not count in, says so and says what
    /// to do.**
    ///
    /// The unknown unit is checked before the count is read, so a table naming a unit the file
    /// does not have is told that — rather than being told the type of the number beside a key
    /// the file has never had.
    #[test]
    fn an_observations_table_that_names_the_wrong_number_of_units_says_which() {
        let (text, _) = the_written_file_with(
            "observations = { reads = 812344 }",
            "observations = { reads = 812344, bases_compared = 7 }",
        );
        let said = refused(&text).rendered_by_the_parser();
        assert!(
            said.contains("`reads`") && said.contains("`bases_compared`"),
            "{said}"
        );
        assert!(said.contains("delete the wrong one"), "{said}");

        for count in ["4", "\"four\""] {
            let (text, _) = the_written_file_with(
                "observations = { reads = 812344 }",
                &format!("observations = {{ lanes = {count} }}"),
            );
            let said = refused(&text).rendered_by_the_parser();
            assert!(
                said.contains("counts `lanes`"),
                "the unit is judged before the count, so a unit this file does not have is what \
                 the message is about whatever type its number is: {said}"
            );
            assert!(
                said.contains("`reads`, `covered_positions` or `bases_compared`"),
                "{said}"
            );
        }
    }

    /// **The hand-written reader accepts every shape the derive did**, which is what stops the
    /// re-wording from narrowing what a file may say.
    #[test]
    fn the_hand_written_reader_takes_every_unit_the_file_writes() {
        let file = a_file_using_every_shape();
        let read = ParametersFile::from_toml(&file.to_toml()).expect("the writer's own text");
        assert_eq!(read, file);

        // All three units, each in the row the file writes it in.
        let counts: Vec<EvidenceCount> = [
            read.base_quality_calibration.by_read_group[0]
                .error_probability_multiplier
                .observations,
            read.inbreeding.by_sample[0]
                .inbreeding_coefficient
                .observations,
            read.repeat_tracts.substitution_rate_by_stratum[0]
                .rate
                .observations,
        ]
        .into_iter()
        .flatten()
        .collect();
        assert_eq!(
            counts,
            vec![
                EvidenceCount::Reads(812_344),
                EvidenceCount::CoveredPositions(180_600_412),
                EvidenceCount::BasesCompared(40_122),
            ]
        );

        // **And the layout `serde`'s own serialiser writes**, which is the other of the two golden
        // files: a `[[…]]` header table rather than an inline one, on the same three units.
        let serdes = toml::to_string(&file).expect("the shape serialises");
        assert_eq!(
            ParametersFile::from_toml(&serdes).expect("and reads back"),
            file
        );
    }
}
