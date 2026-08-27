# Fixes applied — ng psp E3 (the review of `4f138292`)

*2026-08-27. Answers [the review](ng_psp_e3_2026-08-27.md). Step E3 of
[`../../ng/impl_plan/psp_file_format.md`](../../ng/impl_plan/psp_file_format.md), branch
`ng-psp-encoding`.*

---

## 1. The first Blocker is E1's Blocker, one level up

**`read_record_head` applied the chain-id changes and only then bounded the record's body.** A body
that stops early is a `Truncated`, and that class means *fetch more bytes and re-parse this record
from its first byte*. `BlockStream` does exactly that, under a comment saying the restart resumes
*"against state this arm has not touched"* — which this step falsified.

**All three review agents found it independently.** Two measurements, and the second is the worse
one:

- **A well-formed file is rejected as damaged.** Over 1,999 records each naming twenty live reads,
  in blocks larger than the 16 kB rolling buffer, the reader refuses at record 149 with *"id 150,
  which is already live"*: the retry meets an arrival the first pass already applied.
- **A record that only departs reads retries to `Ok` with the wrong set.** A departure is a
  *position* in the live set. The first pass shrinks it; the second resolves the same positions
  against the shrunken set and takes out different reads. Measured, `[1, 2]` where the truth is
  `[1, 2, 4]` — no error at all, and every later record in the block composed against a set short
  by one read.

**And it is the same defect E1's review found, one level down.** There, `read_changes` applied the
departures between reading them and reading the arrivals. Here, `read_record_head` applied the
whole thing before it bounded the body. Both times the fix is the same shape and both times the
suite was silent.

**Fixed** by splitting the parse from the apply — `LiveSetReader::parse_changes` reads without
moving the set, `apply_the_changes_just_parsed` moves it — with `reader.take` on the body in
between. The comment on the line now carries both failures, so the third occurrence has to be
written past them.

## 2. The second Blocker, and the reason the suite was silent about both

**A writer that named only the first observation's reads passed all 4,770 library tests.** A
record's observations are split by allele, by witness and by read group, so a locus with two
alleles from two lanes is four observations with the reads spread across them — and no psp fixture
had ever put chain ids on more than one.

That is the same blindness that hid the first Blocker: **every fixture on the restartable-parse
path had `chain_ids: Vec::new()`**, so each record's changes were the two bytes `0, 0`, and
applying them twice is applying them never. The one test that provably retries — 37,209 restarts
at a byte a read — walked records that named no reads at all.

**Fixed by making the fixtures name reads**, in the three places it matters:

| test | what changed |
|---|---|
| `records_that_straddle_the_buffer_name_their_reads_exactly_once` | **new** — 1,399 incompressible records each naming twenty live reads, in blocks past the rolling buffer, read whole and then a byte at a time. Asserts the blocks really are larger than the buffer and that a byte a read really does retry more often than there are records |
| `a_record_cut_short_is_truncated_at_every_cut` | its record names reads now, and it asserts the live set has not moved at any of the cuts |
| `the_reads_a_record_names_are_the_union_over_its_observations` | **new** — ids on three observations, checked on both sides of the codec |

## 3. The Majors

| what | what was done |
|---|---|
| **`read_record_head`'s `&mut LiveSetReader` carried an undocumented "exactly once a record" contract**, written only in an inline comment rustdoc does not render — and `decode_the_body_of`, the pairing that satisfies it, was not exported, so the only route for an outside caller holding a head was `decode_record`, which parses it again | the contract and an `# Errors` section are in the rendered doc, and `decode_the_body_of` is in `mod.rs`'s re-export list |
| **Nothing tested `further_in`'s re-basing**, and nothing tested that it keeps the fault's class where `inside_a_bounded_body` converts it. Both mutations passed everything | `a_sub_parsers_fault_is_re_based_and_keeps_its_class`, over all three classes and then through a real record head. The two mutations fail 1 and 3 tests now |
| **Nothing tested `skip_unknown_field`'s `ChainIdChanges` arm** — replacing its arrival-count read with a constant passed all 4,770 | `a_later_writers_chain_id_changes_field_is_measured_and_stepped_over` |
| **The new encoding's doc said it "cannot be walked past by anybody"** while the same commit walks past it in `skip_unknown_field` | corrected: it *can* be measured and stepped over, and what must never be stepped over is the copy in a record's **head**, which is a field this reader knows by name and refuses a file that moves |

## 4. The Minors worth naming

- **`FieldReader::skip` was the one advancing method with no bound**, on a value that came out of
  parsing untrusted bytes, in a type whose doc rests on `bytes_read <= bytes.len()` being total —
  and `read_varint` slices on that unconditionally, so an overshoot is a panic in a decoder whose
  contract is that corrupt input gives an error. Today's one caller cannot overshoot; E4 adds
  another. It is bounded now, with a `debug_assert` naming the invariant.
- **`RecordEncoder::for_block` is still a second way to open a block**, called with a placeholder
  coordinate by `BlockBuilder`, while `encode_record_starting_a_block`'s doc claimed to be "the
  only way". The claim is now "the only way to *re*start one", which is what is true, and
  `for_block` says what it is and what reaching `encode_record` through it would cost.
- **Doc drift the commit created and did not reconcile**: five places still counted four head
  fields and a twenty-three-field manifest; `PerBlockState`'s doc still argued the chain-id live
  set belongs *inside* it, which is not where this step put it (and the doc now says why —
  `LiveSetWriter` owns scratch that must survive a boundary, and one entry point resets both);
  `decode_record`'s doc still argued against being two calls after the commit made it two; and a
  header test still said "that is not one of the six".

## 5. What the review confirmed rather than found

- **12,000,000 fuzzed inputs across eight seeds** through the record decoder, each under a deadline
  and a 2 GiB address-space limit and with the harness shown able to abort: no panic, no hang, and
  no allocation driven by a declared length.
- **The step's two design claims hold, traced rather than taken.** All three of
  `encode_record_starting_a_block`'s refusals fire before either reset, `write_a_record` is
  infallible after them, and `BlockBuilder` has no other path needing a rollback — the mutation
  that reverses the order is caught by `a_refusal_at_a_cut_leaves_the_open_blocks_live_reads_alone`.
  And the changes really do sit in front of `body_bytes`'s reach.
- **The header side is complete**: the spelling, `ALL_ENCODINGS`, the parse match, the parameter
  check and the round trip over every scheme.
- **`further_in` is right to preserve the class**, and the reasoning was checked: the chain-id
  parser is handed the rest of the record, not a bounded slice, so more bytes is exactly what would
  help.
- **Every figure in the implementation report re-derived and holds.**

## 6. One guard that is defensive, and says so

`apply_the_changes_just_parsed` is a no-op unless a parse is waiting. **Removing that guard passes
every test in the module**, because no caller applies twice today. It stays because `read_changes`
parses *and* applies while `read_record_head` calls the two separately — the two spellings sit
beside each other and the mistake is one line away. That is the same reason `begin_next_block`
resets zstd's decoder, and it is recorded as defensive rather than left implying the line is
load-bearing.

## 7. Validation

| command | result |
|---|---|
| `cargo fmt --check` | exit 0, no output |
| `cargo clippy --lib --tests --all-features -- -D warnings` | clean |
| `cargo test --lib ng::psp` | 237 passed; 0 failed |
| `cargo test --lib` | 4,774 passed; 0 failed; 14 ignored |

**Six defects re-injected against the strengthened tests:**

| defect | tests failed |
|---|---|
| **the changes are applied before the body is bounded** (Blocker 1) | 2 |
| **only the first observation's reads are named** (Blocker 2) | 1 |
| the skip arm ignores the arrival count | 1 |
| a sub-parser's fault is not re-based | 1 |
| a sub-parser's `Truncated` is converted to damage | 3 |
| applying the changes twice moves the set twice | **0 — reported as defensive**, §6 |

Five for six, with the sixth named rather than counted.
