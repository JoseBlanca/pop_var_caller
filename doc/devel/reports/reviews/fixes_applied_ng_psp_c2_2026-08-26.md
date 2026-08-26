# Fixes applied — ng psp C2, the record head and the skip

*2026-08-26. Answers [the review](ng_psp_c2_2026-08-26.md) of commit `b3177cdd`, finding by finding.
Branch `ng-psp-encoding`.*

---

## What changed, in one paragraph

**Both Blockers are fixed and so are all thirteen Majors**, and three of them changed the shape of
the code rather than adding a test: the offset base is now a type the encoder owns and advances
itself; the head's byte-count width is written once and derived everywhere; and a fault inside a
body the head already bounded is damage rather than a request for more bytes. The record test count
goes from 53 to **65** and the library's from 4,634 to **4,646**. Three claims of this project's own
were measured wrong and are corrected in place, one of them a number carried out of C1's fixes after
it had stopped being true.

## The number I restated instead of re-measuring

**C1's golden fixture is 77 bytes, not 75.** The fixes that closed C1's Blocker widened its
reference bases by one and its first observation's bases by one; I then wrote "the fixture's body is
still 75 bytes — the two cancel" into the C1 report, the C1 fixes report and `PROJECT_STATUS.md`.
They do not cancel, they add. Corrected in all three, with the C1 impl report's own note saying what
happened. **This is the failure mode this project has recorded before and it happened again the same
way**: a number about my own fixture, restated from memory rather than re-run, in prose rather than
in an assertion. The assertion that would have caught it is the golden byte list, which was in the
same commit.

## Finding by finding

### Blockers

**B1 — a fault inside a bounded body reported as the retry signal.** *Fixed.*
`RecordDecodeError::inside_a_bounded_body` re-expresses every fault the body reported as a fault in
the record containing it: a `Truncated` becomes `Malformed` with the reason *"it runs past the N
bytes the head declared for the body"*, and every offset shifts to record-relative. `decode_record`
is the only caller, because it is the only place that knows the body was bounded —
`decode_record_body` on its own is legitimately handed short buffers, which C1's own cut sweep
tests. `a_fault_inside_a_bounded_body_is_damage_and_more_bytes_do_not_help` runs the same damaged
record twice, once with the buffer holding exactly it and once with 4,096 further bytes, and asserts
the two errors are **equal** — which is the property that makes a retry pointless, stated as an
assertion rather than as prose.

**B2 — the coordinate ceiling.** *Fixed on both sides.* `record_span` derives the width without
`GenomeRegion::len`, so the only overflow left is the whole-axis region; a region whose last base is
the last coordinate there is gets its own refusal, `EndsAtTheCoordinateCeiling`. The reader refuses
the same regions, so **no file this module writes can contain one and no other writer's file can
smuggle one in** — which matters because the value is not the reader's to keep: it is handed out in
a public `RecordHead` and detonates in whoever asks its width.

### Majors

**M1 — the offset base.** *Fixed, and further than the review asked.* `OffsetBase` has two
constructors — `at_block_start(Position)` and `after(&RecordHead)` — so the slip that reads most
naturally, handing a record its own start, no longer compiles. And the encoder **holds its own base
and advances it**: `RecordEncoder::for_block(first_position)`, `start_block`, and an
`encode_record` that takes no base at all. Two of the three measured failure modes are now
unrepresentable rather than merely tested. The third — a *reader* holding a stale base — is still
possible, because holding an old `OffsetBase` is sometimes deliberate, so
`a_base_stale_by_one_record_moves_every_record_after_it` prices it: the walk consumes the run
exactly, reports no error, and lands the third record early by exactly the gap the base never
absorbed. Milestone D1 moves the reader's half the same way.

**M2 — the declared encodings were pinned by nothing.** *Fixed.*
`record_fields_declares_the_names_a_written_file_carries` now asserts `(name, encoding)` pairs for
all twenty-three fields, so C1's Blocker is closed on the axis its own fix missed.

**M3 — the head the encoder hands back.** *Fixed.*
`the_head_the_encoder_hands_back_is_the_head_it_wrote` compares it with the head read from the
bytes, over all three fixture records — which differ in region, body length and non-reference count,
so no single mutation of the returned struct survives.

**M4 — the head's `u32` guards.** *Fixed.*
`a_head_count_too_large_for_its_field_is_refused_rather_than_narrowed`, the head's copy of the body's
own test.

**M5 — two origins for one offset.** *Fixed.* Everything out of `read_record_head` and
`decode_record` is now record-relative; `decode_record_body` called on its own stays body-relative
because it has no head in front of it, and `bytes_in`'s doc says both. The head-against-body
mismatch asserts its offset.

**M6 — half of C1's classification pair.** *Fixed.*
`a_head_declaring_a_body_no_buffer_holds_is_truncated_not_malformed` pins the decision **and its
reason**: spec §8 refuses a fixed maximum record size, so there is no length at which a declared
body becomes damage — what bounds it is the head's own `u32` and Milestone D's block.

**M7 — a fourth thing the head cannot describe.** *Fixed* with B2.

**M8 — the body's width in three places.** *Fixed.* `BodyByteCount` is the type; the head's field,
the decoder's ceiling and the encoder's guard all read the width from it, and
`a_body_is_at_most_four_gibibytes_wherever_that_is_written` spells the number out. `declared_body_bytes`
is lifted out of the encoder so `BodyTooLong` has a test at all — reaching it through the encoder
needs a four-gibibyte body no test may allocate.

**M9 — the skip test compared a walk against itself.** *Fixed, and it changes what C3 will be.* The
test now carries a ⚠ saying its evidence is the comparison against the fixture, not the agreement
between the walks. Beside it,
`a_walk_that_skips_every_other_record_matches_a_full_decode_on_the_ones_it_keeps` builds only every
other record and reads the rest by head alone — **C3's oracle in miniature**, so that C3's is
written in a shape that can fail.

**M10 — no arbitrary-input property over the new entry points.** *Fixed.*
`arbitrary_bytes_read_as_a_record_or_are_refused_but_never_panic`, over eight base values including
the ceiling, and it **asks the decoded region for its width on purpose** — that call is the
assertion, which is why it is not written as `is_empty`.

**M11 — the head's read count unchecked against the body.** *Fixed.* `decode_record` derives it from
the body it just built and refuses a disagreement. It costs one more pass over observations already
in hand, which is small beside building them, and the disagreement that matters reads low: a varying
position whose head says zero is one the cohort's first pass skips and nothing looks at again.

**M12 — the wire name.** *Fixed:* `body-bytes` → `record-body-byte-count`. It says it is a count,
and it does not put a field called `body-bytes` inside the TOML body whose own length `header.rs`
calls `body_bytes`. **On the spec:** §4.3 calls this field `record_length`, which is a different
inaccuracy — the head is part of the record and this field measures only the body, which §4.3's own
prose says. Both specs put the byte layout and the names in the implementation's hands
(psp_record_encoding §1.3, psp_file_format §1.2), so this is not a departure needing a ruling; it is
recorded here so the next reader of §4.3 knows which name the files carry.

**M13 — what the count counts.** *Fixed* on the public field's own doc, which is what a filter
author reads.

### Minor and Nits — what was taken

Applied: the module diagram now spells the head's fields as the manifest does and says the chain-id
changes are not there yet; the duplicate manifest test is gone (the golden name list subsumes it);
`three_records_in_order`'s first gap is now 13 against a span of 7, so
`positions_are_rebuilt_…` exercises the case its doc claims — **it did not before, and the review
measured that**; the head's field machinery moved out from under the body's banner and became
`RECORD_HEAD_FIELDS` / `record_head_field`, with the body's helper `body_field`; `BodyReader` is
`FieldReader` and its doc covers both halves; `bytes_read` on the two whole-record types is
`record_bytes`, so the body's consumed count keeps the name that means it; `RecordInBuffer` is
`LocatedRecord` and carries `#[must_use]`; `push` is `encode_record`, so `decode_record` has a
visible twin; `StartsBeforeThePreviousRecord`'s fields are `previous_position` and `offered_start`
and its message no longer names a record that may not exist; `DECLARED_FIELD_COUNT` and
`declared_fields` say whose set they mean; the docs the `RecordLayout` rename left behind are
rewritten; "nothing is written before refusing" is stated on `encode_record` and asserted at both
refusals; a zero offset is pinned as legal with the reason; and the three measurement-quoting slips
are fixed — the "99 positions in 100" now names the corner it was measured on, the duplicated 2.06×
sentence points at the module doc, and the 20 M records/s is attributed to decoding with the
writer's own rate stated as unmeasured.

**Also corrected: "changing the head's encoding later is a manifest change, not a format change"**,
which appeared in the doc comment, the C2 report and `PROJECT_STATUS.md`. It is wrong — this reader
refuses a file declaring anything but what it knows, so a build that switched would reject every psp
written before it and vice versa. The deferral to D2 still holds; its cost was understated. It is
free only while no psp exists, which is until Milestone F.

**Not taken:** `RecordEncoder::for_block` deliberately has no `Default` and no `new` — a writer with
no block to write is not a value — which also settles the review's Minor about a derived `Default`
reappearing. And the unknown-field error still cannot name *which* later field ran out: the declared
names are `String` and this module's errors name a `&'static str`, so the byte offset the error
carries, plus the declared order, is what identifies it.

## Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::psp::record` | **65 passed** (53 before) |
| `cargo test --lib ng::psp` | **109 passed** |
| `cargo test --lib --bins --tests --examples` | library **4,646 passed**, 14 ignored; every other target green except `examples/ng_generic_loci_dump` (11 failures, all `NotFound { ".../ref.fa.repeats.parquet" }`, pre-existing) |
