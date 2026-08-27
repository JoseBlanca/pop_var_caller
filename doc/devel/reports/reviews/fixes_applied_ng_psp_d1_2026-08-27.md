# Fixes applied — ng psp D1 (the review of `a5eccee9`)

*2026-08-27. Answers [the review](ng_psp_d1_2026-08-27.md) finding by finding. Step D1 of
[`../../ng/impl_plan/psp_file_format.md`](../../ng/impl_plan/psp_file_format.md), branch
`ng-psp-encoding`.*

---

## 1. What the review was, and what it found

Eight checklists, each in its own worktree. Between them: **72 mutations applied and run**,
**1.4 million fuzzed inputs** to the two parsers and **501,434 builder pushes**, and the first
independent re-derivation of every number the commit claimed.

**Two Blockers, nine Majors.** Both Blockers were found by more than one agent independently, and
**neither was a defect in the code** — both were properties the code has and no test held:

- **the cut is by a record's *start*, not its end**, and no fixture had a record whose span
  crossed a grid multiple, so changing one identifier passed all 21 tests;
- **`BlockRecords::split` handed a whole-block caller `Truncated`**, the class whose instruction
  is *fetch more bytes* — the C2 review's Blocker one level up, and this one *was* a defect.

**The parsers themselves came out sound.** Over 700,000 arbitrary-byte decodes of
`BlockHead::decode`, 700,000 of `BlockRecords::split`, 400,005 round trips and 40,000 random
record streams: no panic, no arithmetic overflow, and nothing ever sized from a declared length —
a head declaring `record_count = u64::MAX` in front of three bytes allocates nothing. Two
properties Milestone D3's retry loop will rest on held on every input: **every strictly shorter
prefix of a decodable head came back `Truncated`**, never `Ok` and never `Malformed`; and
**`Malformed` is stable under adding bytes** — 21,286 damaged inputs, twelve further random bytes
each, none of which turned the refusal into anything else. And the refusal rollback held
byte-for-byte across 91,247 interleaved refusals.

## 2. The two Blockers

### B1 — nothing held the cut to a record's start

`cell_of` measured the grid cell from `region.start`, which its own doc called a decision.
Changing it to `region.end` passed all 21 tests, because every fixture in the file was a record
of span 1 or a record whose span stayed well inside its 100,000-base cell.

**Why it is not cosmetic.** A record's span is *sample-dependent* — a deletion widens a locus in
one sample and not in another — so a cut taken from the end makes a block boundary depend on
which sample is being written. That destroys the one property the grid exists for, and it
destroys it silently: every file still reads back self-consistently, and only a cohort read
across samples would ever notice.

**Fixed** by giving the fixtures the case: a new test over a record starting at 99,998 and ending
at 100,002; a record widened across each cell's own upper boundary in the whole-cut oracle; a
widened record in the sparse arm of the cross-sample test, which now also asserts that the two
samples cut the *same number* of blocks; and a property test asserting that every record in a
block falls in that block's own grid cell, over random grids and spans. **The mutant now fails
five tests.**

### B2 — `split` handed a whole-block caller the retry class

`BlockRecords::split` forwarded `BlockHead::decode`'s error unchanged. `Truncated` is right for a
growing buffer, and that is what Milestone D3 will hold — but `split`'s caller holds the block
entire, so no quantity of further bytes changes the answer and the retry is a fixed point.

**Fixed** by converting at the boundary where the length becomes known, mirroring what
`record.rs` already does for a record inside a bounded body, with a test at every cut of a
nine-byte head. Damage passes through unchanged, which is its own test.

## 3. The Majors

| what | what was done |
|---|---|
| **A varint fault misfiled as `Truncated` passed all 21 tests**, and the arm was a wildcard over a `#[non_exhaustive]` enum where `record.rs` matches it exhaustively | `VarintError` is matched exhaustively, the overflow's wording is this format's rather than the codec's, and a test covers all three head fields — and asserts the refusal survives 4,096 further bytes |
| **`BlockHead::record_count` was a `u64` guarded by a doc comment**, so `encode` wrote heads `decode` refused | it is a `NonZeroU64`; the zero has no representation, and `OpenBlock` opens with the record that opened it already counted |
| **The same-cell encode had no rollback and no test reached it** — the existing rollback test's codec refusal always landed on the cut path | the live buffer is truncated back on refusal, and a test offers a refusal inside the open block's own grid cell |
| **`from_manifest` was tested only with no ceiling**, so dropping the ceiling entirely passed | the manifest test now exercises both halves of the declared cut rule |
| **`from_manifest` ignored the manifest's declared field layout**, though an append writes records into a file that already declares one | it refuses a layout this build cannot write, through `RecordLayout::from_manifest` |
| **A zero byte ceiling was accepted** where `header.rs` refuses it in a manifest | refused, with the header's own reason corrected: a zero ceiling does not close blocks empty, it gives every record a block of its own |
| **Nothing pinned "no ceiling" as the default** — a hidden 1 MiB fallback passed all 21 tests | `DEFAULT_BLOCK_BYTE_CEILING` is a named constant with the open question beside it, and a test over 89,999 records in one grid cell asserts one block |
| **No property test**, on two round-trip laws with `proptest` already a dev-dependency | both laws are property tests now |
| **The per-block reset had no compiler-flagged home** — a field added to `RecordEncoder` and initialised in `for_block` alone builds clean, is never reset, and is the exact shape Milestone E's chain-id live set will take | `RecordEncoder`'s per-block state is a `PerBlockState` struct that `start_block` rebuilds whole, and `BlockBuilder`'s fields are sorted into what lives for the file and what lives for one block |

**The last of those buys less than it looks, and the doc now says so.** Measured both ways: a
field added *inside* `PerBlockState` is `error[E0063]` at its only constructor, so it cannot be
left out of the reset; a field added to `RecordEncoder` *beside* it and initialised in `for_block`
still compiles and still is never reset — 148 tests green. The split does not force the choice.
It makes the reset automatic once the choice is made, and names the two lifetimes so a wrong
choice is visible.

## 4. Minors and Nits taken

Names, mostly, and all of them applied and compiled: `closed`/`open`/`last_accepted` were three
of the bare participles the checklist names outright, and are now `closed_block_payload`,
`open_block`, `last_accepted_region`; `BlockCell` named its containment backwards — a grid cell
contains a block, not the other way round — and is `GridCell { contig, cell_index }`;
`BlockHead::decode` returned an unnamed `(Self, usize)` against a decision `record.rs` had
written down in the opposite direction, and returns a `DecodedBlockHead`;
`byte_ceiling_reached` is `has_reached_byte_ceiling`; `BlockHeadError` is
`BlockHeadDecodeError`, since only decoding can fail; `BlockWriteError::Record` is
`RecordRefused`.

Also: `BlockCutRuleError` is now its own type, because building a builder and pushing a record
fail in disjoint ways and F3 would otherwise be handed a variant that cannot occur at its call
site; `BlockHead::encode` and `from_manifest` destructure exhaustively, so a field added to
either type is a compile error at the place that has to decide what to do with it;
`ContigOutOfOrder`'s message names its arguments rather than passing them positionally, and a
test renders it; the truncation test asserts the field and the offset it was reporting, both of
which could be hard-coded to constants with the suite green; `Manifest::as_this_build_writes_it`
is the one place the four defaults are assembled, where four sites had been building the literal
by hand; and `PspReadError::CorruptBlock`'s hazard note now names `BlockHeadDecodeError` beside
`RecordDecodeError`.

**One finding was refused, and its claim withdrawn instead.** `a_record`'s doc said records of
the same span carry different reference bases, and the review showed the oracle's records aliased
in pairs. The suggested fix was to key the bases on something that cannot alias — but **no
keying can deliver that**: over a four-letter alphabet, two one-base records must sometimes carry
the same base, and the fixture is full of one-base records. The attempted fix failed on its first
run, which is how this was found. The claim is gone, replaced by what is true: the builder cuts
on coordinates and never reads inside a record, so no fixture here could tell one record's
payload from another's, and `record.rs`'s own property test over arbitrary bytes is what owns
that.

## 5. ⚠ A number in D1's commit message was wrong

The commit message reads: *"about 10 % of the file on a patchy sample at a 5 kb grid, 17.557
bytes a record against 16.444"*. **The 5 kb figure is 18.242**; 17.557 is the figure at 100 kb.
Spec §4.1's table: 5 kb 18.242, 20 kb 18.084, 100 kb 17.557, 1,000 kb 16.444. The 10 % is
18.242 against 16.444, and what 100 kb recovers is 17.557.

[The implementation report](../implementations/ng_psp_d1_2026-08-27.md) §2.6 has all three
right, so this is the commit message alone. **Fixed forward here rather than by rewriting
history**, which is this project's rule.

The review re-derived every other number the commit and the report claimed and reproduced them:
21 tests, 137 in `ng::psp`, 4,653 → 4,674, 63 records over three contigs and three cells, and
nine of the ten per-mutant kill counts. The tenth — *"the cut is a running total, not a grid"*,
reported as 5 — the agent's own implementation of that mutation killed 3; the description admits
more than one implementation and neither of us can claim the other's number.

## 6. What the review found and this commit did *not* act on

- **`block.rs` will want splitting into a directory when D2 and D3 land** (module_structure).
  Noted for D2 rather than done blind: the split should be named by what the files own, and D2's
  compressor is the first sibling.
- **`FieldReader` is now written twice**, once in `record.rs` and once as `BlockHeadReader`
  here. Deliberate, and stated at the type: `FieldReader` is hard-wired to `RecordDecodeError`,
  and making it generic over the fault would put a trait call or a callback on the path that
  decodes about twenty million records a second, to serve a parser that runs about a hundred and
  sixty times per sample. What the duplication carries over is the invariant
  (`bytes_read <= bytes.len()`), which is what makes the slicing total.
- **⛦ Spec §12 question 3 says the accumulate-across-empty-spans rule *ships*.** D1's report and
  commit message read it as leaving the rule itself open; it leaves only its *threshold* open.
  `Manifest` has no field for the rule and the plan's D1 does not list it, so it is still not
  built — but it is a **decision for the owner at Checkpoint D**, not a thing this fix commit
  should have quietly settled. Raised there.

## 7. Validation

| command | result |
|---|---|
| `cargo fmt --check` | exit 0, no output |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::psp::block` | 32 passed; 0 failed |
| `cargo test --lib ng::psp` | 148 passed; 0 failed |
| `cargo test --lib --bins --tests --examples` | library 4,685 passed; 0 failed; 14 ignored. Every other target green except `examples/ng_generic_loci_dump` (11 failures, pre-existing) |

**Every defect the review found surviving was re-injected against the strengthened tests**, one
at a time, on a clean copy:

| defect | tests failed of 32 |
|---|---|
| the cut takes a record's end, not its start | 5 |
| every varint fault reported as `Truncated` | 1 |
| `split` forwards `Truncated` to a whole-block caller | 1 |
| the ceiling fires on `>` rather than `>=` | 2 |
| the ceiling counts twenty bytes of block head against itself | 1 |
| `from_manifest` drops the declared ceiling | 2 |
| a `None` ceiling silently becomes 1 MiB | 1 |
| the same-cell refusal leaves a stray byte behind | 1 |
| the same-cell path counts a record before encoding it | 1 |
| `Truncated` always reports offset zero | 1 |

Ten for ten. The eleventh — `start_block` assigning one field of `PerBlockState` rather than the
whole struct — is a **no-behaviour-change** mutation while the struct has one field, and is
listed as such in §3 rather than as a survivor.
