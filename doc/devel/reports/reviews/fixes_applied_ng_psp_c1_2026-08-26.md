# Fixes applied — ng psp C1, the record body codec

*2026-08-26. Answers [the review](ng_psp_c1_2026-08-26.md) of commit `74563cf7`, finding by
finding. Branch `ng-psp-encoding`.*

---

## What changed, in one paragraph

**The Blocker is closed by two assertions that do not come from the code they test** — the
fixture's exact 75 bytes, and the nineteen declared names written out longhand — and by a doc
comment that stops claiming the array drives the codec. **Twelve of the thirteen Majors are fixed;
one is deferred to Milestone D with its hazard named.** Three claims the commit message led with
were measured wrong and are corrected here rather than quietly dropped. The record test count goes
from 22 to **39**, the library's from 4,603 to **4,620**, and
`cargo clippy --lib --tests --all-features -- -D warnings` — red on this branch before, which is
why C1's test code had never been linted — is now green and joins the gate.

## The one wrong byte the Blocker's own fix caught

The golden-bytes test was written with the fixture's summed log-error worked out by hand: I put
`255, 170, 21` where the encoder writes `255, 159, 21`. The slip was one step in the LEB128
division (2,720 where the remainder makes it 2,719) and the test failed on its first run. **That is
the test doing exactly the job it was added for**, one commit early — and it is why the byte string
in it must never be regenerated to make a test pass.

## Finding by finding

### Blocker

**B1 — `BODY_FIELDS` and the codec are two lists.** *Fixed.*
Two golden tests: `the_fixture_encodes_to_these_exact_bytes` (75 bytes, annotated field by field)
and `record_body_fields_declares_the_names_a_written_file_carries` (the nineteen names longhand).
Between them, reordering the codec fails the first and reordering the array fails the second.
The array's doc comment no longer says it drives the codec; it says what the array is — **a file's
fingerprint of its own layout, which nothing in the types can hold to the bytes** — and names the
test that does hold it, and says plainly that changing both together is a format change.

The `macro_rules!` that would emit both the array entry and the codec call from one line was
considered and not built: it needs the cardinality the manifest does not carry, so it is a header
change at the earliest.

### Major

**M1 — a declared length no body could hold reported as the retry signal.** *Fixed.* A new
`MOST_BYTES_A_BODY_CAN_DECLARE` (the `u32` the head's `body_bytes` is) bounds every declared byte
length and every declared count, and past it the fault is `Malformed`. `read_count` divides that
ceiling by the entry's byte floor, so a count is bounded by what a body could actually hold rather
than by the raw `u32`.

**M2 — nothing pins the `Truncated` / `Malformed` split.** *Fixed.* The 75-cut sweep now asserts
`Truncated` at every cut, and `a_length_no_body_could_hold_is_malformed_and_never_truncated` covers
the three doors into the other class — an over-long reference, an over-long observation count and an
over-long run count. `RecordDecodeError`'s doc names both tests as what holds the line.

**M3 — the `u16` narrowing guard has no test.** *Fixed.*
`a_witness_coordinate_past_a_locus_position_is_refused_rather_than_narrowed` supplies 65,536 as a
run start and again as a run length, and asserts the field and the number in the message.

**M4 — `encode_body` reads fields one at a time.** *Fixed.* Both record types are destructured
exhaustively with **no `..`**, so a field added to either is a compile error on the write side.
`chain_ids` is named rather than swept up, so a rename shows. And because the compiler's own repair
invites a neutral value, `the_fixture_leaves_no_field_at_a_value_a_decoder_could_invent` asserts
every count in the fixture is non-zero — the fixture needed four values nudged to satisfy it
(`num_fwd` on the second observation, `q_sum` and `placed_left` on the third, and `read_group` on
the first), which is the point rather than a side effect.

**M5 — witness runs rewritten into a different witness.** *Fixed, and the decision went the other
way from the C1 report.* Runs are now required to ascend and never touch, checked as they are read;
out-of-order, overlapping and touching pairs are all `Malformed`. The report had called the
normalisation deliberate. It is not defensible: read depth is raised over a witness, so a merged run
credits a read with positions it never saw, the byte count stays right so C2's length check cannot
see it, and it was the one place in the decoder where corruption produced a different valid record.
Three cases tested.

**M6 — `take`'s overflow check untested.** *Fixed.*
`a_length_prefix_past_the_address_space_is_refused_not_indexed` supplies 2⁶⁴−1 as a length prefix at
offset 0 and again at a non-zero offset — the case where a cursor that added without checking lands
behind itself.

**M7 — a locus kind from a later minor version read as corruption.** *Fixed.*
`RecordDecodeError` gains a third variant, `Unsupported { field, bytes_in, tag }`, whose instruction
is *upgrade the reader* — the same one `header.rs` gives for a newer format, and the reason the
error type's doc now opens by saying there are three things a caller might do rather than two.

**M8 — the doc promised a refusal no code performs.** *Fixed, by making the doc true.* It now says
plainly that this reader cannot tell a per-record field from a per-observation one, that such a file
is **accepted and decoded into plausible nonsense from the second observation onwards**, and that a
later writer adding a per-observation field must raise the format version. The same correction is
made in the C1 implementation report and is recorded in `PROJECT_STATUS.md`, because the original
claim reached the commit message and the status file.

**M9 — `#[derive(Default)]` on `RecordBodyLayout`.** *Fixed.* Dropped. `current()` is renamed
`as_this_build_writes_it()` and carries a `# When this is right` paragraph naming the one case
(bytes this process encoded itself) and saying which checks it skips; `decode_record_body`'s doc
points at both constructors.

**M10 — the declared step can be severed from the type's.** *Fixed.*
`the_declared_step_is_the_types_step_and_is_four_thousand_and_ninety_six` pins both the link and the
number, and the refusal test now derives its wrong value from the right one
(`STEPS_PER_NAT / 4`) so it cannot be satisfied by a coincidence of literals.

**M11 — a kind can be added write-only.** *Fixed.* `every_locus_kind_round_trips` walks a list whose
membership an exhaustive `match` enforces, so a kind added to `LocusKind` is a compile error in the
test as well as in `put_kind`. `the_locus_kind_tags_are_the_numbers_the_files_carry` pins all three
tag bytes, which no test did — two of the three were pinned only by accidents of other fixtures.

**M12 — the walked-past test uses one unknown field.** *Fixed.* A two-field case with different
encodings, so the order matters, and a cut case asserting `Truncated` at every cut inside the
trailing field. The fixed-point fixture is now two bytes rather than one, which is what makes that
encoding's skip load-bearing — a one-byte value cannot tell a varint skip from a one-byte skip.

**M13 — neither new error can fold into `PspReadError` without loss.** **Deferred to Milestone D,
deliberately.** The two variants the review drafts need a block index, a record index and a path,
and D and F are what settle what context those carry; adding a public enum variant now and changing
it then is worse than adding it once. **The hazard is recorded where it will be met**: a
`RecordDecodeError` must not be folded through `CorruptBlock`, whose `#[source]` is a
`std::io::Error` and whose documented meaning is *the file is damaged* — the wrong instruction for a
short buffer, and the one thing D4's retry has to branch on.

### Minor and Nits — what was taken

Applied: the truncation sweep also runs under a layout with a trailing unknown field;
`a_squared_mapping_quality_sum_past_u32_round_trips`; `entries_to_reserve` gains an absolute
ceiling (`MOST_ENTRIES_RESERVED = 4,096`) beside its relative one, with a direct test — a 1 MiB
body reserved 13.0 MB before, and the two byte floors are measured by a test rather than asserted in
a comment; a long, non-`ACGT` sequence through a multi-byte length prefix;
`the_fields_a_writer_declares_survive_a_header_round_trip`; the fixture's reference bases now cover
its region (seven over seven, where it carried six); the two duplicate tests are one; the
three-run witness test drops its redundant whole-record assertion and says what it adds; the decode
errors now name their fields **with the manifest's own names**, so a message and a header entry are
the same string and the third vocabulary is gone; `decode_record_body` returns a named
`DecodedRecordBody { record, bytes_read }` rather than a tuple whose second element `_` discards;
`BodyBytes` is `BodyReader` with a `bytes_read` field and every advancing method prefixed `read_` or
`skip_`; `room_for` is `entries_to_reserve`; `RecordLayoutError` is `RecordBodyLayoutError`;
`encode_body`/`decode_body` are `encode_record_body`/`decode_record_body`, because `body` in this
module already means the header's TOML body; the module's summary line covers both halves;
`bytes_left` saturates and says why; both error enums carry the line saying why they are
`#[non_exhaustive]`; the design argument moved off a test doc into a pointer; `put_signed_varint`
joins `put_varint` so the sink is on the same side on adjacent lines; the errors carry `bytes_in`,
the offset into the body, which was already held and discarded.

**Four wire field names changed**, which is free now and a format break later:
`reads-placed-left` → `reads-starting-left-of-the-locus` (the old name did not say left of what),
`reads-discarded-by-cap` → `reads-discarded-by-the-depth-cap` (the header never defined the cap),
and `ssr-motif` / `ssr-left-flank` / `ssr-right-flank` → `repeat-*`, which is the word the error
messages and the tests already used. `mapq-sum` and `mapq-sum-of-squares` are **kept**: MAPQ is the
standard term for the quantity and the review's alternative spelled it out at the cost of the word a
reader would search for. No psp has ever been written.

**Not taken:** indexing the error's field by position in `BODY_FIELDS` behind a `BodyField` newtype
— passing the manifest's own `&'static str` name achieves the same single vocabulary in four lines
rather than forty. And the unknown-field error still cannot name *which* later field ran out: the
declared names are `String` and this module's errors name a `&'static str`; the byte offset the
error now carries, plus the declared order, identifies it.

## Validation

Run in the container on the tree being committed:

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --lib --all-features -- -D warnings` | clean |
| `cargo clippy --lib --tests --all-features -- -D warnings` | **clean — new to the gate**, and red on this branch before these fixes |
| `cargo test --lib ng::psp::record` | **39 passed** (2 head, 37 body; 22 before) |
| `cargo test --lib --bins --tests --examples` | library **4,620 passed**, 14 ignored; every other target green except `examples/ng_generic_loci_dump` (11 failures, all `NotFound { ".../ref.fa.repeats.parquet" }`, pre-existing) |

Three of the clippy fixes are in test code from steps A1–A3 rather than from C1 — a hex literal's
digit grouping, a complex type behind a named alias, and a `vec!` where an array does. They are here
because the gate they were blocking is one C2 and C3 want.
