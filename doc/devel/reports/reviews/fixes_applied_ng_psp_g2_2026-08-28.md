# ng psp store — G2 fixes applied

*2026-08-28, on top of `d2b29285`. Answers [the G2 review](ng_psp_g2_2026-08-28.md), whose
findings are numbered there.*

---

## 1. Finding by finding

| finding | what was done |
|---|---|
| **M1** the never-decoded claim | `a_declined_records_body_is_never_decoded` — the agent's test, adopted. It is the only shape that can hold the claim: a body a full walk *cannot* read, so a walk that declines it must walk past without looking. The fixture it needs already existed inside another test and is now a helper, `a_psp_whose_record_this_build_cannot_read` |
| **M2** the weaker damage reader | a fourth ⚠ paragraph on the type, **with a number measured here rather than quoted**: on a three-record block of 102 payload bytes, every byte flipped in turn, a full walk refuses 93 and a walk declining every body accepts 72 of those 93. `a_declining_walk_accepts_damage_a_full_walk_refuses` is that measurement, and it asserts the inequality rather than the counts |
| **M3** the constant-false predicate | replaced by `the_cohorts_first_pass_predicate_keeps_the_records_where_something_varied`, whose fixture makes every third record carry a read that disagrees — and which asserts that the fixture separates the run before it proves anything |
| **M4** the untested `current_block` | `a_selective_walk_from_a_later_block_names_that_blocks_own_ordinal` — the agent's test, adopted |
| **M5** the sweep that misses the branch | `a_psp_with_damaged_blocks_walks_selectively_without_panicking`: 400 psps with one to four bytes flipped, walked with a predicate that declines every body |
| **Minor** the inert `Debug` | written out by hand, naming the predicate and skipping it. **It repaid itself at once**: the agent had to avoid `expect_err` in one of its own supplied tests because the derive did not apply, and that test now uses it |
| **Minor** the `RecordHead` re-export | removed; `reader.rs` names the type through the module root, where its public path already is |
| **Minor** the forwarded accessors | `SelectiveRecordIter::walk()` added, so a question added to `RecordIter` later is reachable without also being added here; the three spellings stay as the ones a selective walk is asked often |
| **Minor** the names | `only_where` → `building_only_where`, and `SelectiveIter` → `SelectiveRecordIter`. The old pair read as the filter the doc spends two paragraphs denying |
| **Minor** arch §4.1 and §1 | `records_where` returns a `Result` with the reason, `RecordIter::building_only_where` is written out with the decision behind it, and §1's tree names `SelectiveRecordIter` |
| **Minor** untested paths | `records_where_refuses_a_manifest_this_build_cannot_read` (and that the predicate was never called) and `records_where_on_a_psp_with_no_records_is_empty` |
| **Nit** the wrecked-block preamble | `wreck_the_block` joins `a_finished_psp`, `rewrite` and `footer_of` in `tests_support`; three tests use it |
| **Nit** the typo | `the_first_passs_…` is gone with the test it named |

## 2. One suggested fix did not compile as written

The `RecordHead` fix came with a replacement comment for `reader.rs` containing the literal
`` `record::` ``, and the agent reported the suite green with it applied. **It is not**: the
guard test scans the whole non-test half of the file, comments included, and `record::` is on its
forbidden list — the reason G1's review widened that list to the braced import forms. Applied
verbatim it fails `the_opener_cannot_reach_any_block_decoding_code`. The comment is reworded to
say the same thing without spelling a path, and the note itself now says the scan covers comments.

## 3. Not done, and why

- **A predicate that panics and is resumed through `catch_unwind`** re-applies the record's
  chain-id changes and can then report a sound file as corrupt. A panicking predicate is a caller
  bug, the recovery needs `catch_unwind` around an iterator the caller does not own, and the
  behaviour lives in `block.rs`'s `next_record_where`, which G2 did not write. Recorded here.
- **The chain-id live-set test keeps all three records in one block**, so the reset at a block
  boundary is not seen under a predicate. `a_walk_from_any_block_is_the_tail_of_a_walk_from_the_first`
  holds the reset for a full walk; the selective case is worth a fixture and is not urgent.
- **`RecordIter::next` always yields `record: Some(..)` yet its item type carries the `Option`.**
  Pre-existing from Milestone D3 and shared with `BlockStream`; changing it means two item types.

## 4. Verification

In the container, on the tree being committed:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo test --lib ng::psp` — **361 passed, 0 failed**, against 355 before the fixes and 348
  before the step.
- `cargo test --lib` — **4,898 passed, 0 failed, 14 ignored**, against 4,892 and 4,885.
- `cargo test --doc ng::psp` — 1 passed.

**Six tests added net**, four of them the agents' own bodies. The four defects of the step were
re-injected against the fixed tree and all four are still caught; **the two the review found
surviving are now caught too** — the decode-then-discard mutant by two tests, and the dropped
block addend by three.
