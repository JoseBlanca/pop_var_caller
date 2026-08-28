# ng psp store — F3 fixes applied

*2026-08-28. Applies [the F3 review](ng_psp_f3_2026-08-28.md) to the tree at `238c0b49`,
branch `ng-psp-encoding`.*

---

## 1. Finding by finding

| | finding | what was done |
|---|---|---|
| **B1** | a lost block gives a valid-looking file | `PspWriter` carries **`spent: Option<&'static str>`**, set on each of the three failure paths after the builder has handed a block over — the head not reading back, the compressor refusing, the write failing. `finish` refuses a spent writer before it writes anything. ⚠ **The first attempt set it on one path of three**: the replacement for the other two never matched and applied silently, which the anchor check caught. |
| **B2** | the index test never checked the index | `every_index_entry_points_at_the_block_it_names` now **inflates each block and decodes its head**, comparing contig and first position against the entry. The reviewer's mutation — every entry's position shifted 100 bases — now fails it. |
| **M1** | the flush gap was closable here | `a_failed_flush_is_surfaced_rather_than_swallowed`, on `/dev/full`, Linux-only and skipped where the device is absent. **F3 recorded this as needing H2; it did not.** Swallowing the flush error now fails it. |
| **M2** | the compression level never reached the file | `create` records `zstd-compression-level` in the header's writer parameters, which is what `block.rs`'s own doc assigns to F3 by name — and what `append` will need to match bytes already written. |
| **M3** | `create` truncates an existing psp silently | **Documented, not changed** — the fix the reviewer proposes adds a parameter to the signature arch §4.2 fixes, so it is the owner's call. The doc now says plainly that creating over a finished psp destroys it, that this is what `File::create` means and what a re-run wants, that the spec warns about the milder case in `append` and not about this, and that a caller who must not destroy one checks first. Raised in the review's open questions. |
| **M4** | the readable-check could move after the writes | `a_finish_that_refuses_leaves_no_finished_file_behind` — a refused `finish` must leave no file ending in the footer magic. Moving the guard now fails it. |
| Minor | the readable-check omitted the index checksum | it verifies the checksum too; computing it over the trailer now fails **16 tests** |
| Minor | `head_of` is a noun phrase for a fallible decode, with the path first | → `decode_the_head_of(payload, path)` |
| Minor | `check_it_is_readable`'s "it" has no antecedent | → `check_the_index_and_footer_read_back` |
| Minor | `push`'s doc claimed a refusal costs the file nothing | corrected: that holds for a coordinate refusal, and the doc now says every other failure is unrecoverable and why — including that the old "an I/O failure is terminal" reasoning holds only while the failure persists |
| Minor | `finish`'s doc said three durability steps while the code does two | the explanation moved into the doc comment, where a docs reader meets it |

## 2. Not done, and where it went

- **Splitting `PspWriter` to remove the `Option<BlockBuilder>` and the per-block `to_vec()`.** The
  reviewer compiled it and it is a genuine improvement — both `expect`s and both copies gone. It
  is a structural change to the type at the last step of a milestone, and the allocation it
  removes is one per block against a walk that compressed the block. **Recorded for G**, where
  `append` reopens the same machinery and the split will pay for itself twice.
- **Turning three `reason: String` variants into `#[source]` chains.** Right, and it touches the
  error surface that F4 and G both build on. Same home.
- **A `debug_assert` that `written` accounts for every section.** The reviewer measured it as
  taking a fifth-section defect from 1 failing test to 7. Also G, with the split.

## 3. Verification

Seven defects injected into the fixed tree, every anchor asserted to match exactly once:

| defect | caught by |
|---|---|
| a spent writer is allowed to finish | 1 |
| the guard runs after the writes instead of before | 1 |
| the footer's checksum is computed over the trailer | 16 |
| the flush error is swallowed | 1 |
| an index entry's coordinate is shifted past its block | 1 |
| an out-of-order record is accepted (re-run from F3) | 1 |
| the readable-check is skipped (re-run from F3) | 1 |

Gate, in the container:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo test --lib --bins --tests --examples --no-fail-fast` — **library 4,856 passed, 0
  failed**. `ng::psp::` is **319 tests**. The 21 example failures are the known pre-existing
  `ref.fa.repeats.parquet` breakage.

## 4. The wrong number

F3's report and commit message said "twelve defects injected, eleven caught". **The table names
ten and nine were caught.** Every individual kill count was right; both totals were two too high,
and the survivor was called the twelfth when it was the tenth. Corrected in the F3 report.

**Third step running with a miscounted defect table.** F1's was right, F2's had ten rows for nine
defects, and this one had ten rows counted as twelve. The counts come from a script; the totals
were typed. **They are now read off the table rather than written beside it.**
