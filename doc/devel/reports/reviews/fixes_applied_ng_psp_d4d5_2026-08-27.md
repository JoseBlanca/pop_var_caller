# Fixes applied — ng psp D4 and D5 (the review of `4258c4b6..3a278a84`)

*2026-08-27. Answers [the review](ng_psp_d4d5_2026-08-27.md). Steps D4 and D5 of
[`../../ng/impl_plan/psp_file_format.md`](../../ng/impl_plan/psp_file_format.md), branch
`ng-psp-encoding`.*

---

## 1. The surviving mutation, and the gap it sat in

**The reader's new ceiling was measured against the wrong thing, and no test could tell.** `pump`
refuses once the rolling buffer has grown as far as this reader allows. Written against the
buffer's **capacity** rather than what it **holds**, all 81 tests passed — while a well-formed
2,100,075-byte block of three 700,000-base records was refused with *"a record needs more than the
1048576 bytes this reader allows one to hold"*. A valid file, reported as a record too large for
the reader.

**The gap was that no reader test decoded a block past the ceiling at all.** The largest block any
of them met was 210,894 bytes, and at the shipped 100 kb block size a fully covered block is about
1.76 MB decompressed — so the untested case is the ordinary one.

**Fixed** with two tests that come at it from both sides, and the second is the killer:

- `a_well_formed_block_larger_than_the_ceiling_reads_back_every_record` — a block twice the
  ceiling, read record for record, with the buffers ending at the 32 kB budget. The ceiling bounds
  one *record*; a block may be any size.
- `a_block_whose_records_each_grow_the_buffer_to_the_ceiling_still_reads` — three records each
  needing more than half the ceiling, so the buffer's capacity sits **at** the ceiling while what
  it holds is far below. Measured: this is the only test of the 85 that separates the two
  spellings, and under the mutation it is the only one that fails.

## 2. The other Majors

| what | what was done |
|---|---|
| **`BlockCursor` and `PerBlockState` are `Copy`**, so Milestone E's chain-id live set — a collection — gives `error[E0204]` before the `E0063` that is the point, and the cheapest way past it is `BlockStream`, where the field compiles, passes 197 tests, and is never reset | `Copy` dropped from both, with the measurement written at each type. What is left is the error that sends the coder to the right place |
| **The refusal reported its own limit twice** — `held_bytes` can only ever equal `allowed_bytes`, so the message read *"needs more than 1048576 … it had 1048576"* — and never located the record | the variant is `RecordLargerThanTheReaderAllows { block, records_read, allowed_bytes }`: which block, how many of its records were already read, and the one number that is a knob. The redundant field is gone |
| **The ceiling was documented as a knob and there was no knob** — a `pub const` read straight out of `pump`, the third reader budget in this type with no way in | `BlockStream::with_a_buffer_ceiling(ceiling)` sets it and `buffer_ceiling()` reads it back, refusing a ceiling at or under the rolling buffer (`BufferCeilingUnderTheBuffer`). `a_record_past_the_ceiling_is_read_by_raising_it` reads a record the default refuses |
| **The measurement justifying 1 MiB did not reproduce, and priced the bottom of the range.** "Roughly 30 kB" measures **18,292 bytes**; a thousand samples is the bottom, and at three thousand 1 MiB each is **3.07 GB** against spec §1.1's 1.5 GB | **the value moved to 512 KiB**, which is 10× the largest record this caller's depth cap can produce (48,693 bytes at a 150-base span) and gives 1,572,864,000 bytes at three thousand samples — §1.1's budget to within 5 %. The arithmetic is now a test, so it moves when the constant does |
| **`a_refill_inside_a_block_head_is_retried` tested something the design makes unreachable.** A reader clears its buffer at a block start and zstd emits an internal block whole, so a refill lands *before* a head, never inside one — over 108,746 head restarts the largest partial head was **zero bytes** | it is `a_block_head_is_read_after_however_many_restarts_it_takes`, and the reader counts head restarts so the test can assert one happened. The docstring says what is real and what the old one claimed |
| **The every-boundary sweep's docstring claimed the property its fixture cannot reach** — 0 retries over all 837 schedules, while the sibling's docstring said so correctly | kept, because it killed three mutations, with the claim replaced and `assert_eq!(retries, 0)` pinning the fact the sibling rests on |
| **The ceiling's error was tested as a `Display` impl** — swapping the two fields in the message failed nothing, and the only reader that emitted it met a corrupt block | a well-formed oversized record now produces the refusal, and raising the ceiling then reads it. That is the whole distinction the variant exists to draw: the file is damaged, or the reader is budgeted too small |
| **"The measured price of not merging is about 7 %"** cited two rows of the same non-merging cut rule at different block sizes | **what merging would have saved was never measured**, and the comment says that. The ruling stands on the alignment argument, which needs no number |

## 3. ⚠ Seven numbers of mine were wrong

The review re-derived every figure in the three commit messages, both reports and the status
entry. **Everything in D4's retry table reproduced exactly**, all eight rows, along with the record
and block counts, both test counts, D5's eleven-failure injection and every other falsification
row. **Wrong:**

1. **"roughly 30 kB" for a record at the depth cap** — the stated construction measures **18,292
   bytes**; a 150-base span at the same depth is 48,693. The ratio drawn from it was out by nearly
   a factor of two.
2. **"blocks of up to 210,894 bytes"** — it was **one** block. All 1,999 records fell in cell 0 of
   the grid the fixture named, so the sweep retried 58,778 times and never crossed a block
   boundary while doing it. The fixture now cuts: **three blocks, the largest 73,849 bytes, and
   14 retries even when the whole file arrives at once** — the table is re-measured and in the
   test's own docstring.
3. **"about 7 % is the price of not merging"** — those two figures compare block sizes on the same
   non-merging writer. §2 above.
4. **"twelve thousand of them are 1,823 bytes on disk"**, justifying the incompressible fixture —
   it is **4,543 bytes** at the grid its caller uses and **165** at the 100 kb default. Neither is
   1,823.
5. **"thousands of bases from where it belongs"** for the far-apart fixture — cells 0, 7, 41 and
   900 of a 100 kb grid put it about **86 million** bases out. The code comment understated what
   the commit message got right.
6. **"2,048 times what an open sample budgets for"** — 2,048× is against the 32 kB the two buffers
   hold between them; against spec §1.1's 500 kB per open sample it is **131×**. Both are now in
   the sentence.
7. **D3 §2.6's "`record.rs` … is still the smaller of the two"** — `block.rs` is now the larger by
   601 lines. The deferral of the `block/` split still stands, on the two structural premises that
   have not moved.

Corrected forward in the code and in [D4](../implementations/ng_psp_d4_2026-08-27.md),
[D5](../implementations/ng_psp_d5_2026-08-27.md) and `PROJECT_STATUS.md` — not by rewriting
history.

## 4. What the review confirmed rather than found

- **No input makes this reader panic, hang or grow.** 1,523,400 fuzzed inputs across three
  corpora — 1,520,000 byte-level damages, 3,000 lying counts and lengths recompressed, and 400
  inflation bombs. 1,514,438 refused, 8,962 ended cleanly, none did anything else.
- **The ceiling does what it was added for.** Over bombs decompressing to as much as 67,108,864
  bytes, the reader added **1,200,128 bytes** of resident size — sampled every 2 ms, because a
  reading taken after a case is worthless once `fail` has released the buffer.
- **And the harness that says so can fail**, shown twice: with the ceiling disabled it fails on
  memory, and with the retry arm's refill removed the watchdog aborts the process at 20 s.

## 5. Two things left as they are, and said out loud

- **`>=` rather than `>` in the ceiling check is an unpinned one-byte convention.** Measured:
  swapping it fails none of the 201 `ng::psp` tests. What it decides is whether the buffer may
  ever *hold* exactly the ceiling or must stop one below; an oracle would have to search for the
  exact rolling length at the deciding pump, which pins the test to zstd's emission sizes to prove
  one byte on a 512 kB budget. Recorded at the line rather than given a brittle test.
- **The `block/` directory split is still deferred.** The two structural premises are unchanged —
  arch §1 assigns all three concerns to this file verbatim, and the checklist's one directory rule
  points the other way — and `module_structure` filed **no findings**. The quantitative premise
  D3 offered has flipped and is dropped rather than repeated.

## 6. Two spec debts, now scheduled rather than conditional

Both are Checkpoint D items, listed in `PROJECT_STATUS.md` rather than left to *"when that document
is next touched"*:

- **Spec §4.1 and §12 question 3 say the merge-across-empty-spans rule ships**, and the owner
  ruled against it on 2026-08-27. Until they are corrected the design authority instructs the
  opposite of the code, and the next person implementing from §4.1 will implement merging.
- **The reader's record ceiling has no row in spec §6.7 or §7**, which are where each error class
  is paired with what the user has to do — and where Milestone F4's `PspReadError` mapping reads
  them from. The row it wants is *a record needs more of the buffer than the reader allows → raise
  the reader's ceiling, or the block is corrupt*.

## 7. Validation

| command | result |
|---|---|
| `cargo fmt --check` | exit 0, no output |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::psp::block` | 85 passed; 0 failed |
| `cargo test --lib ng::psp` | 201 passed; 0 failed |
| `cargo test --lib` | 4,738 passed; 0 failed; 14 ignored |
| `cargo test --lib --bins --tests --examples --no-fail-fast` | 4,846 passed; two example targets fail on the missing `ref.fa.repeats.parquet` fixture (`ng_generic_loci_dump` 11, `ng_ssr_loci_dump` 10), neither importing `psp` |

Every defect the review found surviving was re-injected against the strengthened tests, one at a
time, on a clean copy:

| defect | tests failed of 85 |
|---|---|
| the ceiling is measured on the buffer's capacity, not what it holds | 1 |
| the ceiling is raised out of reach | **compile error** — the memory fixture's `const` guard |
| a block-head restart is not counted | 1 |
| a record restart is not counted | 1 |
| the rolling buffer is never drained | 13 |
| the buffer never shrinks back at a block | 2 |
| `with_a_buffer_ceiling` accepts a ceiling under the buffer | 1 |
| `>=` → `>` in the ceiling check | **0 — and left that way on purpose**, §5 |

Seven for eight, with the eighth reported as an unpinned convention rather than as a survivor.
