# ng psp store — G4 fixes applied

*2026-08-28, on top of `395cfa26`. Answers [the G4 review](ng_psp_g4_2026-08-28.md).*

---

## 1. The Blocker, and where it belonged

**B1 — a footer putting the block index inside the header let `append` truncate the header
away.** Reproduced here before anything was changed: `PspReader::open` accepted the file, `append`
accepted it, and a 3,742-byte psp became 109 bytes with `finish` returning `Ok`.

**The fix is in `reader.rs`, not in `append`.** `open` now checks `index_offset >= header_bytes`
as a rule about the *file*, before the per-entry loop — because **the per-entry rule cannot carry
it: on an empty index there are no entries to check.** That is the whole defect in one sentence,
and it is the sentence now on the check. Every operation that starts from `open` gets the bound,
which is where G3's `replace_trailer` fix should have put its sibling.

`a_footer_that_puts_the_block_index_inside_the_header_is_refused` pins it from `append`'s side,
and asserts the file is byte-identical afterwards.

## 2. Finding by finding

| finding | what was done |
|---|---|
| **M1** the seam test | `a_record_inside_the_last_blocks_span_is_refused` — a record between the last block's *first* and *last* records, which is the case the old fixture could not reach. It asserts the fixture's own precondition (that the last block spans more than one position) before proving anything |
| **M2** the silently-ignored level | `ManifestRefusal::UnreadableLevel` — **absent** is the one shape that falls back, because it is what a file written before the parameter existed looks like; every other shape is refused. `a_level_recorded_in_a_shape_this_writer_cannot_read_is_refused` walks a string, a float and a boolean, and asserts that absent still falls back |
| **M3** the invented level | `ManifestRefusal::LevelPastAnyLevel { recorded: i64 }` carries the file's own number. `an_append_is_refused_on_a_recorded_level_zstd_will_not_take` covers 99 and ±5,000,000,000, and asserts the refusal's text **contains the recorded value** — which is the half that was wrong |
| **M4** the blockless psp | `a_psp_that_holds_no_records_is_appendable`, asserting its own precondition (zero blocks) first |
| **M5** the `WriteStats` docs | rewritten to say which population each field counts, and `an_appended_file_reads_back_as_one_file` now pins all three: 6 records written, 9 blocks in the file, and `bytes` equal to the file's length |
| **Minor** the `Reopen { .. }` assertion | matches `Incomplete` specifically, with a ⚠ naming G3's M2 |
| **Minor** the seam's cost | `the_seam_is_found_by_reading_only_the_last_block` — the walk sees 5 records, not 40. A seam walk that quietly became whole-file would fail it |
| **Minor** the interruption window | `an_append_stopped_part_way_leaves_a_file_no_reader_accepts` — every truncation from the block index onwards, and the complete file opens |
| **Minor** `push`'s destructure | every field named, no `..`, and `spent` and `path` reached through the destructured bindings rather than through `self` |
| **Minor** the names | `find_the_last_record_in` and `build_the_compressor_the_header_records` — verbs, and the second no longer claims a file it does not take |
| **Minor** the arch tree | `writer.rs` names `append` |
| **the unasserted claims** | `a_record_appended_into_the_old_last_cell_opens_a_second_block_for_it` pins both halves of `continuing_after`'s doc |
| **`ManifestRefusal::Compressor`'s ⚠** | rewritten: G4 made it reachable — an append on a file whose recorded level zstd will not take — and the doc still said nothing reaches it. **G3's review asked for that note and G4 is what made it false**, which is the shape worth naming |
| **the three wrong numbers** | corrected in the G4 implementation report, with what was measured beside what was claimed |

## 3. Not done, and why

- **Two copies of "write a block and give it its index entry"** now exist: `push`'s inline one and
  `put_block`, which `finish` still calls for the last block. Merging them means handing
  `put_block` the borrow split that removed the copy, which is the type split G4 argued against.
  Recorded as the smallest live piece of that argument.
- **`UnsupportedManifest`'s `Display` claims the manifest** for a failure in
  `header.writer.parameters`, which is provenance. A separate variant is the honest fix; it is a
  message, not a behaviour, and the milestone ends here.

## 4. Verification

In the container, on the tree being committed:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo test --lib ng::psp` — **391 passed, 0 failed**, against 381 before the fixes and 373
  before the step.
- `cargo test --lib` — **4,928 passed, 0 failed, 14 ignored**.

**Ten tests added net**, eight of them the review's own bodies. The measurement the review
corrected was re-taken here: a fifth section between the trailer and the footer fails **82 tests
with the byte-accounting assertion and 66 without it** — it accounts for sixteen, where the report
had claimed none.
