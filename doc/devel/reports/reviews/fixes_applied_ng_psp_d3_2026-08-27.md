# Fixes applied — ng psp D3 (the review of `91f2d15f`)

*2026-08-27. Answers [the review](ng_psp_d3_2026-08-27.md). Step D3 of
[`../../ng/impl_plan/psp_file_format.md`](../../ng/impl_plan/psp_file_format.md), branch
`ng-psp-encoding`.*

---

## 1. The Blocker, and it is the defect D3 claimed to have fixed

**A refused stream did not end.** `fail` emptied the reader's own compressed buffer and set no
terminal state, so what actually stopped a refused reader was that its next read hit the end of
the source — which lasts exactly as long as the whole file fits in one 16 kB fill. Past that, the
reader read on from wherever the source stood, took four arbitrary bytes for a block length, and
carried on.

**Four of the five agents found it independently**, and the measurements are worse than a leak:

- on a 69,769-byte, 841-block file with one bit flipped, over a plain `&[u8]` — the exact shape
  every committed test used — it refused after 129 records and then **handed back 3,585 more**;
- on a 2,960-byte file with one corrupt byte it refused after 5,980 records, **handed over 5,681
  more, and ended cleanly**;
- through a source that returns short reads, a four-block file gave **nine further real records
  and no second refusal at all**.

That is a damaged file reported as damaged and then read as though it were not. It contradicts
arch §4.1 — *"An iterator that fails yields `Err` once and then `None`"* — and D3's commit
message, report §5 and `PROJECT_STATUS.md` all list it as a **fixed survivor**. It was not fixed;
the mutant died on the fixtures' size.

**Fixed** with a terminal `refused` flag that every entry point checks and nothing clears, and
with a fixture that can see it: a file of 2,999 incompressible records over more than four read
chunks, damaged in its first block, tried through both a whole-file source and one that dribbles.
**Removing the flag, or removing the guard that reads it, each fails that test now.**

**And making the flag load-bearing took a second pass.** The first version also freed the
compressed buffer, and reading into a zero-length buffer returns zero — which looks exactly like
the end of the source. So the flag was still redundant and its removal still passed. `fail` now
frees the *rolling* buffer, which may have grown past what the reader budgeted for, and keeps the
read chunk, which is the budget. The flag is the only thing that stops a refused reader, and the
tests say so.

## 2. The other Majors

| what | what was done |
|---|---|
| **A block whose length prefix claims more than its frame used was unguarded** — dropping that condition passed all 66 tests while a two-block file read back as one, with records from a block nobody asked for | the three ways a block can fail to end where it said are three variants now, not one reporting `bytes_left: 0` for two of them; and the untested direction has a test |
| **`buffered_bytes` — the number Milestone H4 measures against the 500 kB budget — was pinned by nothing.** A stub returning `0` passed every check the module made, and so did one that dropped the read chunk | it is asserted exactly, on a fresh reader and on one mid-file |
| **The retry arm had one witness.** Putting spec §8's exact defect into it — advance the coordinate before pumping — failed **one** of 66 tests | the walk that crosses the arm most often now compares the records it built instead of counting them, and a new test reads a whole file through a source that returns **one byte at a time**, so every record is retried, most of them many times. The mutant now fails two |
| **No test used a source that returns short reads or fails with `Interrupted`** — every fixture was a `&[u8]`, which returns everything asked for and never fails, so `read_exactly`'s refill loop and the `Interrupted` arm never ran | a source that dribbles and interrupts, at four settings, reading the same records |
| **`Zstd` was the only variant spanning two instructions** — a corrupt frame and a too-narrow window differ *only* in a 20-digit code | it carries zstd's own resolved account: "Restored data doesn't match checksum" against "Frame requires too much memory for decoding" |
| **`rolling_at` was per-block state outside the group the compiler enforces**, and the doc said the opposite | it is in `BlockCursor`, and `BlockCursor::between_blocks` is written out longhand rather than `..Self::opening(0)`, so a field added to it is a compile error at **both** constructors rather than silently inheriting one's choice |
| The buffer's shrink-back was asserted by a **substring** of a `Debug` string — `"163840"` contains `"16384"`, so a ten-times-larger buffer passed | asserted as a number |

## 3. ⚠ A tension between two spec sections that only bites on corrupt input — for Checkpoint D

**Measured by a review agent: a 4,132-byte block drove the reader to hold 67,125,248 bytes**, 2,048
times the 32 kB budget. A psp block's *inflated* size is not bounded by its size on disk, and when
no record can be parsed out of one the rolling buffer doubles until the frame runs out.

**Nothing is sized from a declared length** — 1,473,500 fuzzed inputs found no allocation driven by
one, no panic and no hang — so this is bounded by the data rather than by an attacker's
arithmetic, and `fail` now releases it the moment the block is refused. But:

- **spec §8** says a record larger than the buffer must make it grow, and that *"a fixed maximum
  record size is not a safe assumption to bake in"*;
- **spec §1.1** puts an open sample at 500 kB, and says so is what makes a thousand-sample run fit.

**Both cannot hold on a corrupt file, and which gives is the owner's call.** It is written at the
line that grows and raised at Checkpoint D; nothing here decides it.

## 4. ⚠ Three more numbers of mine were wrong

The review re-derived every figure in D3's report, commit message and status entry. Correct: the
line counts, 66, 182, 17 new tests, 49→66, 4,702→4,719, the `E0063` text verbatim, and nine of the
thirteen falsification rows. **Wrong:**

1. **"a refused stream is not ended → 1"** — it is **2**; the narrow-window test asserts the same
   property. And more to the point, §1 above: neither of them was catching what the row claimed.
2. **"the block-end check is dropped → 2"** — **3** on the tree D3 committed. (On the tree this
   commit ships it is 2 again, because the check is now three separate guards.)
3. **"the decoder is not reset → no behaviour to change"** — **wrong when it was written**.
   Removing that line changed the outcome on **11,664 of 60,000** damaged inputs — and **none of
   them before the first refusal**, which is to say all of them in exactly the states §6 called
   unreachable *because they already end the stream*. They did not end the stream; that was the
   Blocker. **With §1 fixed the claim is true, and I re-measured it on this tree: removing the
   reset leaves all 187 `ng::psp` tests green.** A claim that was false has become true by fixing
   something else, which is worth saying rather than quietly keeping.

## 5. What the review confirmed rather than found

- **No input makes this reader loop.** 1,473,500 fuzzed inputs across three shapes — damaged
  files, damaged payloads recompressed so the reader meets a well-formed frame with damaged
  contents, and lying block heads — with 0 hangs, 0 panics and no allocation from a declared
  length. The largest payload tried was 255,940 bytes, 15.6× the rolling buffer.
- **No block ever returned a record count other than the true one** across 13,500 lying heads.
- **The `block/` directory split is still right to defer**, and this time the checklist said so on
  its own terms: arch §1 assigns all three concerns to this file verbatim, the checklist has no
  line-count rule, and its one directory rule points the other way.

## 6. Two documentation debts paid, and one misquote corrected

- The module doc still said *"Two halves"* and *"nothing reads one back through a rolling buffer
  (Milestone D3)"* — **the exact debt the D2 review attached to the previous deferral, re-incurred
  in the commit that discharged it.** It now describes three.
- D3's report §2.6 said `module_structure` had *"judged the deferral right rather than merely
  tolerable"*. What the D2 review says is that it **did not object**. Corrected.
- The timings in `READ_CHUNK_BYTES` and `next_record_where` read as measurements of this reader.
  They are the measuring prototype's, and nothing has timed `BlockStream`; both now say so, and
  Milestone H is named as where that happens.

## 7. Validation

| command | result |
|---|---|
| `cargo fmt --check` | exit 0, no output |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::psp::block` | 71 passed; 0 failed |
| `cargo test --lib` | 4,724 passed; 0 failed; 14 ignored |
| `cargo test --lib --bins --tests --examples` | every target green except `examples/ng_generic_loci_dump` (11 failures, pre-existing) |

Every defect the review found surviving was re-injected against the strengthened tests:

| defect | tests failed of 71 |
|---|---|
| a refused stream is not marked | 1 |
| the guard that reads that mark is removed | 1 |
| a block over-declaring its compressed length is not checked | 1 |
| `buffered_bytes` is stubbed to zero | 2 |
| the retry arm advances the coordinate before pumping | 2 |
| the block-end check is dropped | 2 |

Six for six.
