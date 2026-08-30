# ng psp store — G1 fixes applied

*2026-08-28, on top of `ffb0e2a3`. Answers [the G1 review](ng_psp_g1_2026-08-28.md), whose
findings are numbered there.*

---

## 1. The Blocker, and the rule that replaces it

**B1 — `records_from` skipped records when two blocks share a first position.**

The search is now `partition_point(|entry| entry.first_position < at)` and then one back: **the
block before the first block that starts at or after `at`**.

**Neither agent's proposed fix was taken, and the reason is a second hole both left open.** Both
proposed entering a run of equal positions at its head, which repairs the reproduced case. It does
not repair the case with no run at all: a block's *last* record may begin on the base the next
block begins on — the same byte ceiling produces it — so a block whose first position is strictly
below `at` can still hold a record starting at `at`.

The rule that misses nothing is the one above, and it is short to justify. Block *k*'s records
start in `[e_k, e_{k+1}]`, both ends inclusive, so a record starting at `at` can be in block *k*
whenever `e_k ≤ at ≤ e_{k+1}`; the earliest such *k* is one less than the first entry at or after
`at`. **Its cost is one extra block, and only when `at` falls exactly on a block's first
position** — which is exactly when the previous block might hold a record there. For a coordinate
strictly inside a block, which is the everyday case, it picks the same block the old rule did.

That is a change to what `records_from` returns at an exact block boundary, and it is inside the
contract spec §6.2 and arch §4.1 already fix: the walk starts at or before the coordinate asked
for. Both documents now say which block, and why.

`records_from_a_blocks_first_position_starts_that_block` asserted the old rule and is replaced by
`records_from_a_blocks_first_position_enters_the_block_before_it`, which pins the new one for all
eight of the fixture's blocks. The reproduction is
`a_walk_from_a_position_two_blocks_share_starts_at_the_first_of_them`.

## 2. Finding by finding

| finding | what was done |
|---|---|
| **B1** | the rule above; two tests, one replacing the one that asserted the old rule |
| **M1** the catch-all | `refuse` is one exhaustive match with no `_` arm, hoisted to a free function so G2's selective walk can share it. `BufferCeilingUnderTheBuffer` now maps to the caller's-mistake class. `walk_from`'s own `map_err` names the two variants `BlockStream::new` can return |
| **M2** the ceiling with no knob | `PspReader::with_a_record_buffer_ceiling` exists, and the message names it. The rule lives in `walk.rs` beside the buffer it is about, because `reader.rs` may not name `block` |
| **M3** the unpinned ordinal | `a_walk_from_a_later_block_names_the_failing_block_by_its_own_ordinal` — the agent's test, adopted verbatim |
| **M4** the out-of-range ordinal | the field doc says which case is not an index and to use `block_index().get(block)`; `a_fault_in_the_bytes_introducing_a_block_names_one_past_the_last` pins it |
| **M5** the transposition | a `WalkStart` struct with named fields; the call site is a struct literal |
| **M6** the record-level upgrade arm | `a_record_naming_a_locus_kind_this_build_does_not_know_asks_for_a_newer_reader` — adopted, and it locates the kind byte by diffing two payloads rather than by a magic offset |
| **M7** `live_reads` untested | `live_reads_names_the_reads_live_at_the_record_just_yielded`, on a fixture whose records carry ids that arrive and depart. **The first version of this test was mine and was the failure the finding describes**: written against `a_finished_psp`, where every `chain_ids` is empty, it passes for every wrong implementation. It is replaced, and the doc says why |
| **M8** the introducing-bytes arm | covered by M4's test, which asserts the source variant as well as the ordinal |
| **Minor** message and naming | `UnsupportedRecordEncoding` carries *upgrade the reader*; `damaged_by`'s `reason` names the section and the cause says what was wrong with it; `blocks_read` → `blocks_begun`; `current_block()` added; `NoSuchBlock`'s fields are `ordinal_asked_for` and `blocks_in_the_file`; the `blocks` field is `block_index` |
| **Minor** `Box<dyn Error>` | replaced by a typed `DamageFound` — `Footer` or `BlockIndex`, transparent — so a caller matches instead of downcasting. The test that had to downcast now matches |
| **Minor** the guard test | the doc says what it checks; the forbidden list gains `block::` and `record::`, the braced forms that walked past it |
| **Minor** the arch doc | `walk.rs` is in §1's tree with the reason it is separate; §4.1 has `block_index()`, a fallible `records()` with why, `records_from(at: GenomePosition)` and `records_from_block` |
| **Minor** untested corners | a coordinate past every block; a contig the file does not carry; `records_from_block(0)` on a psp with no blocks |
| **Minor** no damaged walk | both hostile sweeps adopted: 400 psps with flipped block bytes, 300 with a hostile block index |
| **Nit** the footer preamble | `footer_of` joins `a_finished_psp` and `rewrite` in `tests_support`; three tests use it |

## 3. Not done, and why

- **The `BlockOrdinal` newtype** one agent proposed. `WalkStart` removes the transposition at the
  one call site that exists, and the newtype would reach three public error variants and the
  ordinal's `usize`/`u64` split. It is the right change if a fourth caller appears; it is not
  worth a public-API change in a step that has none.
- **The doubled `: {source}` tails on four `BlockReadError` variants.** Pre-existing (Milestone D),
  and G1 is only what first puts them in front of a caller. Whoever next touches `block.rs`'s
  error text should drop them.
- **`refuse`'s `Io` arm remains untested.** The agent that found it could not build a fixture: it
  needs a `read(2)` failure on a `Take<&mut File>` mid-block, and truncating the file under an
  open walk gives an end of file rather than an error. Recorded rather than papered over, and it
  is H2's — where a writer is killed for real.
- **`tracing` for applied defaults.** ng has no `tracing` anywhere; that is a module-wide decision
  and not G1's.

## 4. Verification

In the container, on the tree being committed:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo test --lib ng::psp` — **348 passed, 0 failed**, against 337 before the fixes and 324
  before the step.
- `cargo test --lib` — **4,885 passed, 0 failed, 14 ignored**, against 4,874 and 4,861.
- `cargo test --lib --bins --tests --examples --no-fail-fast` — the only failures are
  `examples/ng_generic_loci_dump.rs` (11) and `examples/ng_ssr_loci_dump.rs` (10), both
  pre-existing and on a missing `ref.fa.repeats.parquet` fixture. Neither imports `psp`.

**Eleven tests added net**, four of them the reviewing agents' own bodies adopted as supplied.
