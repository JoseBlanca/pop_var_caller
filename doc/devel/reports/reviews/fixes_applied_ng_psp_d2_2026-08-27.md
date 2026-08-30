# Fixes applied — ng psp D2 (the review of `2200001f`)

*2026-08-27. Answers [the review](ng_psp_d2_2026-08-27.md). Step D2 of
[`../../ng/impl_plan/psp_file_format.md`](../../ng/impl_plan/psp_file_format.md), branch
`ng-psp-encoding`.*

---

## 1. What the review found

Eight checklists across five worktrees. Between them **35 further mutations**, **1.4 million
fuzzed inputs** to `compressed_block_bytes`, and an independent re-run of the head-encoding
measurement on both corpora.

**One Blocker, eight Majors, and five numbers of mine measured wrong.** The Blocker and the
largest Majors are all the same shape, and it is the shape D1's review found too: **the
compressor's four settings each had one test, and none of those tests exercised the
configuration a writer actually ships with.**

## 2. The Blocker — the window cap was never tested at the shipped settings

The cap is the format's central design decision. Both tests that proved it reached zstd ran at
`with_level(MIN_LOOK_BACK_WINDOW_LOG, 1)` — a 1 kB window at level 1 — and every test that used
the shipped path compressed about a kilobyte, which is *under* the 32 kB window, so zstd narrows
the frame's own declaration to fit and the cap is inert. Measured: guarding the `WindowLog` call
with `if level < ZSTD_COMPRESSION_LEVEL` compiles, raises no clippy warning, and passes all 42
tests — while every file a writer then produces needs a wider window than its own manifest
declares, which is the file D3's decoder refuses block by block.

**And half of the older test could not fail at all.** At the minimum window, "one exponent less"
asks zstd for a window below its own floor, so the *parameter* is refused and the frame is never
looked at — an assertion satisfied by any bytes whatever. Proved by replacing the frame with
`b"not a frame at all"`: still green.

**Fixed** by `the_shipped_settings_cap_a_payload_larger_than_the_window`, which runs at the
default window and level over a payload measured to exceed the window — and by giving the test
helper a two-state error so it can assert *the frame* was refused rather than the cap. The mutant
now fails four tests.

## 3. The Majors

| what | what was done |
|---|---|
| **The level reached zstd untested** — four agents, independently. Clamping every block to level 3, or forcing every caller to 9, passed all 42 tests, because the only level-mentioning test compared two compressors at the *same* level | `the_chosen_level_reaches_the_bytes` asserts level 19 < level 9 < level 1 in block bytes, and that the writer's constructor writes at the shipped level. Both mutants now fail it |
| **`compressed_block_bytes` returned an unvalidated length** — four agents. Under fuzzing it returned a length longer than the slice it was given 656,785 times, the largest 4,294,967,299 from eight bytes in, and the module's own walk sliced on it | it is `compressed_block_at`, returning `Whole` / `NotAllHere` / `NoLengthYet`. The distinction is the one this module types everywhere else, and it is what a D3 reader pulling bytes from a file actually branches on |
| **`BlockCompressor` never saw a `Manifest`**, so a window could come from anywhere and a manifest field it must honour had no compiler-flagged home. Adding a field to `Manifest` flagged three sites and none was the compressor | `BlockCompressor::from_manifest`, destructuring with no `..`, is what a writer uses. A file whose frames need more window than its own manifest promises is one every reader refuses with a zstd error code — the unactionable failure spec §4.2 says the declaration exists to prevent |
| **Nothing compressed a payload zstd cannot shrink.** Replacing `compress_bound` with the payload's own length passed all 42, and fails on a block holding one record — which is ordinary at sparse coverage | `a_payload_zstd_cannot_shrink_still_compresses`, over a one-record block and an empty one. The oversized fixtures are what hid it: every compression test built forty records or more |
| **The framing had no golden anchor.** Two coordinated changes — writing *and* reading the length big-endian — were killed only by one test's incidental byte literals; rewritten as a round trip, they went green | `an_on_disk_block_is_a_little_endian_length_then_a_zstd_frame` pins the endianness and the frame's magic without pinning zstd's own output, and asserts the fixture's frame is under 256 bytes so a big-endian length could not agree with a little-endian one |
| **`with_level` did not check the argument it exists for.** zstd clamps silently: a level of 100 gives exactly the level-22 bytes, and −131,073 gives a block 169× the shipped level's | refused, against `zstd::compression_level_range()`, with both bounds in the message |
| **The level was unrecoverable** — not in the file, not from the writer's version, not from a running compressor | the compressor keeps it and `compression_level()` hands it back. **Recording it in the header's writer-provenance parameters is F3's**, and the constant now says so |
| **`BlockCompressError`'s message and source were untested** for two of three variants | `a_compressors_refusals_say_what_broke_and_what_the_bounds_are` pins all three, including that zstd's own account survives as the cause |

## 4. A trap the review found that is now written down

**A frame's own declared window is never *wider* than the file's, and is routinely narrower.**
zstd narrows it to fit a payload smaller than the cap, so an ordinary block under a 32 kB
declaration writes a frame needing 1 kB. Spec §4.2 calls a mismatch between our declaration and
the frame's "a corruption worth detecting rather than tolerating" — **and the check that
implements that has to be `≤`, not `=`**, or Milestone D3 rejects almost every block of every
file this writer produces. `a_frames_own_window_is_never_wider_than_the_file_declares` is the
test, and it is there for D3 to find.

## 5. The measurement, and the harness that produces it

**Every published figure reproduced exactly** when an agent rebuilt the harness and re-ran it on
both corpora — all eight bytes-a-record figures, all four percentages, both record counts, both
block counts, and the refused field. **And the conclusion is stronger than D2 claimed**: twelve
further fixed widths were tried, sixteen arms in all, plus four grids and six windows, and no
fixed width beats varint on either sample. Varint's head is already 4.000 bytes a record on
tomato — one byte per field, the floor — which is why nothing can undercut it.

**Two things about the harness were fixed.**

- **The fixed-width arm had no oracle.** The varint arm checks itself byte-for-byte against
  `BlockBuilder`; the fixed arm — whose bytes the conclusion turns on — was checked by nothing.
  It has a test module now: a fixed head reads back field for field, the width bound is exact at
  every boundary, an eight-byte width holds `u64::MAX`, and the varint arm writes LEB128 against
  an independent expectation. A `[[example]]` entry puts them in a plain `cargo test`.
- **The corpus now describes itself.** Two of the labels in §5 below were mine, from the specs,
  and one was wrong.

## 6. ⚠ Five numbers of mine were wrong

1. **The library runs 4,695 tests at `2200001f`, not 4,694.** Two agents measured it
   independently and I confirmed it. D2's report and commit message both say 4,694; the report's
   own preceding figure (4,685 after D1, plus ten new tests) already implied 4,695.
2. **The tomato corpus is 10.25 reads a position, not 3.** The label is inherited from spec §4.1
   and §4.3 — arch §7 already says 11.4 and disagrees with them — so D2 did not invent it, but it
   is stated in D2's report as a fact about the corpus and it is not one. HG002 chr21 measures
   279.99, which matches. **So the measurement's depth range is 10× to 280×, and neither corpus
   reaches the 3× end of the range this caller is committed to.** The harness prints
   `mean-reads-a-position` now, so the label cannot go stale again.
3. **"644 regions" is not re-derivable from the corpus.** It is quoted from spec §4.1 and the
   file gives 1,217 maximal runs of consecutive covered positions in 281 occupied 100 kb cells.
   The sentence it was supporting — why the narrow fixed arm cannot encode the human sample — now
   uses what the harness measures: **15 records have a within-block position offset over 65,535,
   the largest 90,467**. On tomato the widest is 21,115 and none passes it, which is why the same
   arm fits there.
4. **"one allowed a single exponent less refuses it"** was not demonstrated — §2 above.
5. **"every one of the frame's bytes is refused"** is true for whole-byte inversion and not for
   single-bit damage: over one 69-byte frame, 543 of 552 bit flips were refused and **9 were
   accepted and inflated to the payload unchanged**, none to anything different. The guarantee is
   **damage is refused or harmless, never silently different**, and the test now asserts that —
   over every bit rather than every byte.

All five are corrected in [the D2 report](../implementations/ng_psp_d2_2026-08-27.md) and in
`PROJECT_STATUS.md`. Fixed forward, not by rewriting history.

## 7. What was not done

- **`block.rs` is not yet split into a directory.** module_structure judged the deferral right
  and did not object: arch §1 assigns all three concerns to this file by design, the re-exports
  follow the documented convention, and nothing D2 added belongs in `header.rs`. What the
  deferral owed is done — the module doc no longer says "nothing here compresses".
- **The `psp` → `locus_generation` back-reference** — `psp` calls itself a peer and imports
  `SampleLocusObservations` from a pipeline stage. Pre-existing, public API, and 30 files under
  `src/` plus 18 examples consume that type across four stages. Noted for Milestone F, which is
  the last cheap moment.
- **The window's range predicate is spelled in both `header.rs` and `with_level`.** Re-checking is
  justified — `Manifest`'s fields are public — but the *predicate* could be a function in
  `header.rs`, which owns the bounds. Left as it is; both spellings use the same two constants,
  so they cannot drift in value, only in wording.

## 8. Validation

| command | result |
|---|---|
| `cargo fmt --check` | exit 0, no output |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo clippy --example ng_psp_head_encoding --all-features -- -D warnings` | clean |
| `cargo test --lib ng::psp::block` | 49 passed; 0 failed |
| `cargo test --example ng_psp_head_encoding` | 4 passed; 0 failed |
| `cargo test --lib` | 4,702 passed; 0 failed; 14 ignored |
| `cargo test --lib --bins --tests --examples` | every target green except `examples/ng_generic_loci_dump` (11 failures, pre-existing) |

*`cargo clippy --examples` is red on `examples/psp_row_stream_roundtrip.rs`, which this branch
does not touch — verified by `git diff HEAD --stat` on that path being empty. The project's gate
is `--lib --tests`.*

**Every defect the review found surviving was re-injected against the strengthened tests**, one
at a time, on a clean copy:

| defect | tests failed of 49 |
|---|---|
| the window is capped only below level 9 | 4 |
| every caller's level is forced to 9 | 1 |
| the frame buffer reserves the payload's length, not zstd's bound | 1 |
| the level is clamped to 3 | 1 |
| the framing is big-endian on both sides | 10 |
| a declared length is believed unconditionally | 1 |
| the level's range check is dropped | 2 |

Seven for seven.
