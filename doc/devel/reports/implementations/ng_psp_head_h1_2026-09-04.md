# psp record head — H1: the head carries the keep rule's denominator, and the locus kind's tag

**Date:** 2026-09-04
**Plan step:** [psp_head_compared_reads.md](../../ng/impl_plan/psp_head_compared_reads.md) Milestone H, step H1
**Spec:** [psp_head_compared_reads.md](../../ng/spec/psp_head_compared_reads.md) §3, §3.1, §6
**Branch:** `ng-psp-mode`

## Plan

The cohort merge keeps a locus when some single sample shows at least
`max(floor, share × its compared reads)` non-reference reads. The head carried the numerator and
not the denominator, so at three reads a position the floor decided everything and the head
sufficed, and at three hundred the share decided and the head could not answer. This step adds
`reads-compared-with-reference` to the head, and moves the `locus-kind` tag forward from the body
so a reader that has not built a record can still say what kind of locus it is.

No consumer: nothing in this branch reads either field. The change lands now because a head
layout change costs nothing while no psp exists outside this crate's tests and costs a format
version from the run driver's stored-file milestone onwards.

## Assumptions and recorded choices

- **The kind's tag is in the head and its detail stays in the body.** `RecordHead` gains a new
  `LocusKindTag` — the three kinds without a repeat tract's motif and flanks — and the body keeps
  the motif and the two flanks, present exactly when the head's tag says repeat tract. It is a
  move and not a copy, so there is no second answer to check the first against; the body decoder
  is handed the kind it is decoding under.
- **Placement inside the head (the spec marks it soft).** The denominator goes where the spec
  puts it, between `non-reference-reads` and `record-body-byte-count`. The tag goes immediately
  after `reference-span`, because those are the two fields the cohort's width bound reads
  together — it governs generic loci only, and a repeat tract's span may lawfully be wider.
- **`put_kind` splits into `put_kind_tag` and `put_kind_detail`**, and `read_locus_kind` into
  `read_locus_kind_tag` (the head) and `read_kind_detail(tag)` (the body). **A kind added to
  `LocusKind` is still a compile error, in two places** — `LocusKindTag::of` and
  `put_kind_detail`, both exhaustive with no wildcard. It is no longer a compile error in the
  tag's writer, which matches on `LocusKindTag`; what holds that side in step is
  `every_locus_kind_round_trips`, and the tag's reader can never be exhaustive because it matches
  a number read from a file.
- **⛦ An unknown tag is now refused earlier and by more readers, and that is the owner's to
  rule at Checkpoint H.** It is `Unsupported` as before — *upgrade the reader*, not *rebuild the
  file* — but a walk that would have **skipped** the record no longer can, because it meets the
  tag while deciding whether it wants the record at all. **It is a choice, not a consequence of
  the move**: the tag is a self-delimiting integer, so an `Unknown(u64)` arm would let a skipping
  walk carry on and leave the refusal to whoever built the body, which is what a psp from a later
  writer met before this step. What argues for refusing is what the tag is in the head for — two
  of the cohort merge's pre-assembly decisions read it, and the width bound is applied differently
  to a repeat tract than to a generic locus, so a reader that cannot classify a record cannot
  correctly decide it does not want it either. Recorded at `read_locus_kind_tag`; the spec does
  not settle it, and it cost one test its mechanism (below).

## Changes made

- `RECORD_HEAD_FIELDS` 5 → 7: `locus-kind` after `reference-span`, and
  `reads-compared-with-reference` after `non-reference-reads`. `BODY_FIELDS` 22 → 21: the
  `locus-kind` entry is gone, the three repeat-tract fields stay last.
- `RecordHead` gains `reads_compared_with_reference: u32` and `locus_kind: LocusKindTag`.
  `LocusKindTag` is new and public, with `LocusKindTag::of(&LocusKind)`.
- `write_a_record` keeps both halves of `non_reference_and_compared_reads()` — the derivation
  already computed the denominator and the writing line dropped it — and writes both, plus the
  tag.
- `read_record_head` reads the tag and both counts, and **refuses a head whose varying count
  exceeds its compared count** without touching a body. The two are summed over the same
  observations, so the first counts a subset of the second.
- `decode_the_body_of` compares the denominator against the body it just built, beside the
  numerator's existing comparison. One derivation call returns both, so the extra check costs no
  extra pass.
- `decode_record_body` takes the locus kind as an argument; it is what says whether the body's
  last three fields are there at all.

## Tests

`cargo test --lib 'ng::psp'` goes from **417 to 421**; the library target from **6,153 to 6,157**.

Four new:

- **The head carries the compared count**, over two shapes: the rich fixture, where only one of
  three observations spanned the whole locus, so the head reads 0 varying out of 137 compared and
  the partial observations' 3 and 1 reads are in neither count; and the same record with its
  complete observation removed, where both counts are zero and only the rule's floor can decide
  the locus.
- **A head reading more varying reads than it compared is refused without its body** — 5 out of
  3, 5 out of 0 and 1 out of 0 refused; 5 out of 5 accepted, which is what stops the check being
  written as `>=`, and 0 out of 0 accepted, which is the record whose reads all stopped inside
  the locus.
- **A head whose compared count disagrees with its body is refused either way** — 138 and 136
  against a body of 137, by moving one byte of the head's LEB128.
- **A manifest from before this change is refused by name**, in both of its shapes: a whole old
  manifest disagrees at position 2, where this build expects `locus-kind` and the file declares
  `non-reference-reads`; one that stops after the varying-read count names
  `reads-compared-with-reference` as the field it is missing. **The first is the one an old
  writer would actually produce**, and the field it names is not the one spec §6 predicted.

Changed, and each for a reason:

- **`every_locus_kind_round_trips` now goes through the head.** It used the body-only round trip,
  which after the move would hand the decoder the kind it was about to check — proving nothing.
  It now writes and reads the whole record and asserts the head's tag as well as the record.
- **`a_locus_kind_from_a_later_writer_says_upgrade_the_reader` forges the head's third byte**
  rather than the body's last.
- **`the_fixture_encodes_to_these_exact_bytes`** loses the body's trailing kind byte, and
  **`the_head_this_version_writes_is_these_exact_bytes`** is 8 bytes rather than 6 — over a new
  record, because its old one made all three of the changed bytes a literal zero (see the review
  below). Both are the pinned byte strings, moved deliberately.
- **`a_head_whose_non_reference_count_disagrees_with_its_body_is_refused_either_way`** gained its
  second direction, and its name says so.
- **`a_declined_records_body_is_never_decoded` needed a new fault.** Its whole point is that a
  declined record's body is never decoded, and it proved it with a body a full walk cannot
  read — a locus-kind tag no build knows. That fault is in the head now, where a declining walk
  meets it too. The fixture is now a repeat-tract record whose stored motif has been emptied to
  no bases: one byte's value, not a length, so the body is exactly as long as its head declares
  and everything before the motif reads normally. The test also now asserts the half that gives
  the other half its meaning — that a full walk over that file *is* refused, naming
  `repeat-motif`.
- **`a_decode_refilled_at_every_source_byte_boundary_gives_the_same_records` asserted that the
  sweep never retries a parse, and that stopped being true.** With one more byte a record, 281
  parses restart over the sweep's 1,188 refill schedules — 34 at one byte a read, 17 at two, 2 at
  seventeen, still 1 at 1,183 bytes of the 1,188-byte file, and none when the whole file arrives
  at once. The mechanism is not established and the docstring now says so: thirteen of the file's
  fourteen block payloads are 93 to 100 bytes and the last is 38, and raising the writer's block ceiling until they are 127
  to 128 bytes takes every schedule back to zero, which is the opposite direction from *larger
  payloads are cut more often*. The sweep's real property — every schedule returns the same
  records — is untouched, and the sibling test remains the only one where a block is larger than
  the 16 kB rolling buffer. The assertion is now that the sweep does retry, so the count cannot
  quietly stop being taken.
- Byte-index constants in **seven** places moved with the head: the synthetic block payload in
  `block.rs`, four forged or asserted head positions in `record.rs`, the forged head in the
  selective-walk fixture, and the head re-encoder in the unknown-trailing-field walk.

### Mutation pass after the fixes

Re-run on the fixed tree (`tmp/h1_mutations.sh`, four mutations, all killed):

| mutation | result |
|---|---|
| the head-only check exempts a zero denominator | 1 failure — was 0 before the fix |
| the numerator's head-vs-body check refuses only when the head reads low | 1 failure — was 0 before the fix |
| the denominator's head-vs-body check refuses only when the head reads low | 1 failure |
| the kind tag and the varying count swap places in **both** the writer and the reader | 5 failures, and the head's pinned byte string is now one of them — it was not before the fixture changed |

## Validation results

Run in the container from the worktree:

- `cargo fmt` clean; `cargo clippy --all-targets --all-features -- -D warnings` exit 0.
- `cargo test --lib` — **6,157 passed**, 0 failed, 14 ignored.
- `cargo test --tests` — every integration target green.
- `cargo build --all-targets` — no errors.

## Review

Two reviewers, each in its own worktree over this step's patch. No Blockers.

### What the correctness reviewer found — three tests that could pass with the feature broken

It ran six mutations. Four were killed; **two survived the whole psp suite**, and both were real
gaps rather than false alarms:

- **The head-only check exempted a zero denominator and nothing noticed.** Guarding it with
  `compared > 0` left 421 tests green. What that lets through is the worst shape there is: a head
  saying five reads varied out of none compared, whose share the keep rule computes as zero, so
  the locus is kept on a denominator saying nothing was ever looked at. The refusal test now runs
  three shapes — 5 out of 3, 5 out of 0, 1 out of 0 — and asserts that 0 out of 0 is accepted,
  which is the record whose reads all stopped inside the locus.
- **The numerator's head-against-body check was provoked in one direction only.** Refusing solely
  when the head reads *low* also left 421 green, and a head declaring 100 varying reads over a
  body showing 19 was accepted. That test now runs both directions, matching the denominator's,
  which had both from the start.

It also found that **the head's pinned byte string was not pinning the fields this step added**.
Its fixture had no observations and a generic kind, so the kind tag, the varying count and the
compared count were all a literal `0` — indistinguishable from each other and from a reordering.
Measured: swapping the writer's and the reader's order of `locus-kind` and `non-reference-reads`
*together*, which is exactly the silent format drift that test exists to catch, left it green. The
fixture is now a repeat tract with 2 reads varying out of 5 compared, so the three bytes are 1, 2
and 5.

One clause of the spec came back **weakened rather than implemented**: §6 expects a manifest from
before this change to be refused naming `reads-compared-with-reference` as missing, and a genuine
pre-change file is refused one field earlier, at `locus-kind` — a consequence of putting the tag
before the counts, which the spec marks soft. The test says so now and is renamed accordingly.

Categories 1 to 3 — the encode/decode pair, the two refusals' placement and boundary, and the
kind's move — came back clean, clause by clause.

### What the fallout reviewer found

**Fixed here.** The module's own record diagram still drew a four-field head; `put_kind_tag`'s doc
argued a compile-time guarantee that moved to `LocusKindTag::of` when the match changed type, and
`put_kind_detail` had the guarantee and not the note; the new head-only check called itself *the
one validity check a reader can make without touching a body* three lines below three others that
also touch no body; the module doc quoted the head's 9.2 % and 5.8 % without saying they were
taken on a five-field head; `benches/ng_psp_perf.rs` enumerated the head's decoded fields and was
two short; the retry-sweep note said fourteen payloads of 93 to 100 bytes when the last is 38; and
this report undercounted the moved byte-index sites. `LocusKindTag` also gained
`#[non_exhaustive]`, which `LocusKind` has and which it mirrors one-for-one.

**Left for H2, which is the next step.** The owning specs still describe a four-field head:
[`psp_file_format.md`](../../ng/spec/psp_file_format.md) §2's definition list, §4.3's diagram, its
field-by-field list and its 0.077-bytes-a-record line;
[`psp_record_encoding.md`](../../ng/spec/psp_record_encoding.md) §2.3;
[`cohort_merge.md`](../../ng/spec/cohort_merge.md)'s position-summary entry; and
[`run_streaming.md`](../../ng/spec/run_streaming.md) §3.3's ⚑, which is now met and still reads as
open. One new code comment cites §4.3 for the kind tag, which §4.3 will only say once H2 lands.

**Recorded, not fixed — pre-existing at the branch point.** This module's doc says a record's
chain ids are dropped and join at Milestone E, and `BODY_FIELDS`' doc says the same; both were
already false before this step, since the chain ids landed at E4. `chain_ids.rs`'s
*"Nothing calls this yet"* is false the same way. `PROJECT_STATUS.md`'s C2 entry spells the
four-field head, and is a dated log entry rather than a description of the format today.

## Tradeoffs and follow-ups

- **`examples/ng_psp_head_encoding.rs` still measures a four-scalar head** and so under-reports
  it. Extending it to the sixth field and re-measuring on both corpora is step H3, which is the
  next step of this milestone.
- **The head grew from 32 to 40 bytes in memory** (`size_of::<RecordHead>()`), which its own test
  bounds. It is read once per record and holds no allocation; nothing measured this step's cost
  on the wire, which is H3's.
- **No consumer**, by design: the cheap-numbers read that folds heads across samples belongs to
  the successor plan.
