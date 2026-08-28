# ng psp store — F4 fixes applied

*2026-08-28. Applies [the F4 review](ng_psp_f4_2026-08-28.md) to the tree at `6b7bfc02`,
branch `ng-psp-encoding`. This is the last step of Milestone F.*

---

## 1. Finding by finding

| | finding | what was done |
|---|---|---|
| **B1** | the sections rule was tested on one side | `sections_that_stop_short_or_run_past_the_footer_are_both_refused` — the fixture nudges the trailer's declared length **both ways**. Relaxing `!=` to `<` now fails it. |
| **M1** | the budget's derivation was wrong twice | rewritten on the constant: the buffers are §4.4's, not §5.3's, and **the 190 kB already contains a 32 kB window**, so the room for one is `500 − 190 − 16 − 16 + 32 = 310 kB` rather than 278. 2^18 is the largest power of two under both — **a wrong derivation for a right number**, and the doc now says so rather than quietly correcting itself. The §4.2 citation moves to §7's table, which says *raise the budget, or rewrite the file*. |
| **M2** | the test's name overclaimed | renamed `opening_decompresses_no_block`, which is what it shows, with a note that an `open` reading the blocks without decompressing passes it and that showing *that* needs `open` drivable over an arbitrary source. **The stronger statement is now a test**: `the_opener_cannot_reach_any_block_decoding_code` reads the module's own imports and fails if either the block or the record module is named. |
| **M3** | a block offset into the header opened | the range is `header_bytes..index_offset` — bounded at **both** ends. `read_header_and_its_length` surfaces the length `Header::decode` already returned and `read_header` was throwing away. Two tests: one for an entry at byte 0, one for an entry **exactly at** `index_offset`, which the old fixture missed by one byte. |
| **M4** | the allocation bound is not where it looks | the comment now states where the bound actually comes from — `decode_footer`'s abut rule plus `read_footer`'s sections rule give `index_offset + index_bytes <= file_bytes` before this function is entered — and says plainly that **the check below it cannot fire**, that it was written believing it was the bound, and why it stays. |
| Minor | `trailer()` read twice only by accident | `the_trailer_reads_the_same_twice`; dropping the seek now fails three tests. |
| Minor | the budget is not stored | `PspReader::look_back_window_budget_bytes()`. At several thousand open samples the budget is the number that multiplies, and a caller tuning it should be able to read back what it got. |
| Minor | a second `open(2)` per sample | the header is read from the handle `open` already holds, through `read_header_from`, which rewinds first so it does not depend on where the caller left the cursor. |
| Minor | the vocabulary trap fired a third time | the module doc's list of redefined words now carries it, and names the sharpest case: **production has a `PspReader` too, and its `trailer()` returns the fixed tail — which here is `footer()`.** Kept rather than renamed, because spec §6.2 fixes both words; what was missing was the warning. |

## 2. The three wrong numbers

- **The library was 4,841 before, not 4,836.**
- **`ng::psp::` is 319 after and 304 before, not 314 and 299** — and the report contradicted its
  own commit message, which is the tell that neither was measured.
- **The `tests_support` lift was described as part of this diff**; it landed in the parent commit.

Re-measured by checking out the parent and running it. **The ten-defect table was right** — the
first in this milestone to survive re-scoring intact, after F2's had a duplicate row and F3's had
totals two too high.

## 3. Verification

Six defects injected into the fixed tree, every anchor asserted to match exactly once:

| defect | caught by |
|---|---|
| the abut rule only refuses sections that stop short | 1 |
| a block exactly at the index's first byte is allowed | 1 |
| a block inside the header is allowed | 1 |
| the trailer's seek is dropped | 3 |
| the opener reaches for the block module | 1 |
| the index length is not bounded by the file | **0 — and that is the point of M4**: the bound is real and lives elsewhere, so this check cannot fire. Recorded, not hidden. |

Gate, in the container:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo test --lib --bins --tests --examples --no-fail-fast` — **library 4,861 passed, 0
  failed**, against 4,856 at `6b7bfc02` and 4,841 at its parent. `ng::psp::` is **324 tests**. The
  21 example failures are the known pre-existing `ref.fa.repeats.parquet` breakage.

## 4. Left for later, with a home

- **`Damaged { reason: String }`** joins the six other `reason: String` variants F3's review
  routed to Milestone G, where `append` and `replace_trailer` build on the same error surface.
- **The 500 kB budget may not fit what an open sample holds** — the block index alone is about
  336 kB for a whole genome. **H4**, and it is written on the constant.
- **`blocks()` returns what the spec calls the block index**, and `head_magic` is a noun for a
  function that reads. Both are renames inside a surface G extends; taken there rather than at the
  last step of a milestone.
