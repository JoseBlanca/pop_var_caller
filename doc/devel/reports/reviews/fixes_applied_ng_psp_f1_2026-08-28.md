# ng psp store — F1 fixes applied

*2026-08-28. Applies [the F1 review](ng_psp_f1_2026-08-28.md) to the tree at `34047878`,
branch `ng-psp-encoding`. Every finding is addressed; none is deferred.*

---

## 1. Finding by finding

| | finding | what was done |
|---|---|---|
| **B1** | the ordering scan only ever ran on its first pair | every ordering fixture is **eight entries with the break walked across all seven pairs**. The `.take(1)` mutation now fails 3 tests. ⚠ The first attempt still could not fail: starting the run at position 1 made the broken entry *equal* to its predecessor at the first break, and equal positions are legal — the fixture starts at 1,001 now, and says why. |
| **B2** | `Display` panics at entry 0 | both variants carry `previous_entry` beside `entry_number`; **no arithmetic in any message template**. `every_refusal_renders_at_entry_zero_without_panicking` renders all six variants at entry 0 and asserts none is empty and none prints a struct dump. |
| **M1** | `LARGEST_ENTRY_BYTES` is 18, the widest entry 23 | corrected to `5 + 10 + 8`, with the doc saying why the 18 was there and what it was carried from. **And made checkable**: `the_width_constants_are_the_widths_the_encoder_produces` encodes the widest and narrowest entries and asserts both constants against the bytes — without it, restoring the 18 was caught by nothing, because nothing observable reads a reservation. |
| **M2** | no fixture resets position at a contig boundary | `a_position_that_resets_at_a_new_contig_is_accepted` — five entries across three contigs, each starting at base 1 again. The raw-position mutation now fails 2 tests. |
| **M3** | the checksum was pinned only by avalanche | a golden value, `683_841_834`, measured from the fixture. Swapping in FNV-1a now fails. |
| **M4** | an over-long varint reported as a truncation | a new `Overlong` variant (the enum is `#[non_exhaustive]`, so adding one is free), and `take_varint` matches on `VarintError` rather than discarding it. Tested with ten continuation bytes in a buffer that does not end early. |
| **M5** | no message was pinned | `each_refusal_says_what_is_wrong_and_where` asserts all six rendered messages as whole sentences. Cutting one down, or mislabelling a field, now fails. |
| **M6** | a field added to the entry is silently dropped by the encoder | `encode_index` destructures with no `..`, so a new field is a compile error there too — the same guard `BlockBuilder::from_manifest` already uses on the manifest. |
| Minor | `offset` collides with the base-pair sense | renamed `block_offset`, matching production and the footer's `index_offset`/`trailer_offset`, with the collision recorded on the field. |
| Minor | stringly-typed field names | `IndexEntryField { Contig, FirstPosition, BlockOffset }` with one `Display`, so a message cannot name a field the format does not have. `FieldTooLarge` becomes `ContigNumberTooLarge` — the only field that can raise it. |
| Minor | `entry`/`at` overloaded | `entry_number` for ordinals, `byte_cursor` for the cursor, `index_entry` for the test constructor. |
| Minor | the offset's width stated three times, one an `expect` | `first_chunk::<8>()` makes it a type: no `checked_add`, no `expect`, and the mutation that turned eight tests into panics no longer compiles. |
| Minor | the checksum doc's claim was false of ng's layout | rewritten: all four regions outside the blocks are uncompressed, and what distinguishes the index is that the other three carry framing damage shows up in — a length and sentinel, a magic — while the index is a bare run of entries. |
| Minor | the decode's doc did not say what it leaves unchecked | it now says a contig number is checked against the width of a `ContigId` and **not** against the header's contig list, and an offset not against the file's length, and that both belong to `open`. |
| Minor | `TrailingBytes`' half-names | `trailing_bytes` / `entries_read`. |
| Nit | the plan's F1 wording describes F3 and F4 | the entry now says the codec is F1's and the building and reading are F3's and F4's. |

## 2. Not done, and why

- **A writer-side `debug_assert` that entries are ordered.** Suggested by one checklist, tried,
  and **reverted**: the codec's own tests must be able to encode a disordered index in order to
  prove the decoder refuses one, and the assert fired in three of them. The obligation is real
  but belongs to **F3's `finish`**, which is what assembles an index from `BlockBuilder`'s
  output; it is recorded there rather than dropped.
- **Refusing padded LEB128.** An index therefore has more than one byte form, and its checksum
  moves between them. Bounded in practice: ng always re-encodes the index from decoded entries
  and rewrites the footer with it, so a foreign padded index becomes canonical on append with a
  checksum that matches. Recorded on `decode_index`; refusing it outright is a format-level rule
  and would be its own change.
- **Splitting `decode_index`'s three jobs into three functions.** Judged and left: the
  allocation bound and the ordering pass both read better where the doc argues for them.

## 3. Verification

Seven defects injected into the fixed tree, **seven caught**:

| defect | caught by |
|---|---|
| only the first pair is order-checked (`.take(1)`) | 3 |
| raw positions compared, contigs ignored | 2 |
| the checksum becomes FNV-1a | 1 |
| an over-long varint called a truncation again | 1 |
| the truncation message stops naming the field | 1 |
| a field spelled with the wrong word | 2 |
| the entry width shrinks back to 18 | 1 — **and by nothing at all before the constants test was added** |

Gate, in the container:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo test --lib --bins --tests --examples --no-fail-fast` — **library 4,811 passed, 0
  failed**, against 4,804 at `34047878`. `ng::psp::index` is **20 tests**, against 13. The 21
  example failures are the known pre-existing `ref.fa.repeats.parquet` breakage.
