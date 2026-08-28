# ng — building the psp store: implementation plan

*Draft, 2026-08-26. Turns the settled design — [`../spec/psp_file_format.md`](../spec/psp_file_format.md)
(the container), [`../spec/psp_record_encoding.md`](../spec/psp_record_encoding.md) (the record),
[`../spec/psp_chain_id_encoding.md`](../spec/psp_chain_id_encoding.md) (the chain ids) and
[`../arch/psp_file_format.md`](../arch/psp_file_format.md) (the types) — into build order.
**Not a place for new design:** a question this plan cannot answer from those four goes back to
them.*

*The experiments that produced the design are
[`psp_encoding_experiments.md`](psp_encoding_experiments.md); its Milestone Z is done and its results
are in [the memory review](../../reports/reviews/psp_memory_milestone_z_2026-08-25.md).*

---

## Scope

**In:**

- `src/ng/psp/` — the whole module: header, block, record, index, footer, chain ids
  ([arch §1](../arch/psp_file_format.md)).
- **The three approximated quantities becoming integers at the type**, which reaches outside
  `psp/`. The owner put this work here on 2026-08-25 rather than in a plan of its own.
- The oracles: a parity check against the measuring prototype, restart-equals-sequential, refusal of
  an interrupted file, worker-count invariance.

**Out — with a home each:**

- **Wiring the store into the locus generator or the cohort merge.** This plan builds a module and
  proves it against a fixture; making the pipeline write and read one is
  [`../spec/run_streaming.md`](../spec/run_streaming.md)'s, which owns the run objects.
- **Deciding the trailer's contents.** The container stores opaque bytes
  ([spec §3.4](../spec/psp_file_format.md)); what goes in them follows the statistical work.
- **Any change to production's `.psp`.** ng's store is new code beside it; production stays as it is.
- **A newtype for `ChainId`** — arch §7 raises it and routes it to whoever owns ng's chain-id
  minting.

## Principles (how the order was chosen)

- **Types first, then implementation**, inside every milestone (project rule).
- **The algorithmic heart before the plumbing.** A record round-trips in memory (C) before anything
  compresses (D) or writes a file (F). The hardest correctness work is then testable with no I/O.
- **Rounding moves upstream before the encoder depends on it.** B changes types outside `psp/`, and
  it comes early because C's `FixedPoint` encoding stores what those types produce. Building C first
  would mean writing an encoder that rounds, then removing the rounding.
- **Reuse over rewrite.** `encode_u64_leb128` and friends are called as-is
  ([arch §6](../arch/psp_file_format.md)); this plan writes no second varint.
- **Verify against something outside the code.** The prototype
  (`examples/psp_row_stream_roundtrip.rs`) already reads and writes this shape and produced every
  measurement the specs quote; it is the oracle, not self-consistency.
- **Isolate the steps that fail silently.** E is the chain ids, where a wrong answer is a genotype
  rather than a crash. Its steps are their own commits with the oracle green either side.
- **Incremental, with pauses.** One milestone, then stop.

## Preconditions (already in place)

- `SampleLocusObservations` and `SequenceObservation`
  ([`src/ng/locus_generation/mod.rs:40,295`](../../../../src/ng/locus_generation/mod.rs)) and
  `ReadWitness` ([`witness.rs:259`](../../../../src/ng/locus_generation/witness.rs)) — what is
  written and what must come back.
- The varint codec ([`src/psp/varint.rs:46,83,119`](../../../../src/psp/varint.rs)).
- ng's newtypes: `GenomeRegion`, `ContigId`, `Bp`, `ReadGroupId`
  ([`src/ng/types.rs:79,13,185,210`](../../../../src/ng/types.rs)).
- The measuring prototype, which writes and reads the settled shape and whose `verify` walks it
  against a source in lockstep.
- A tomato accession and a 279-reads-a-position human sample, both used throughout the specs.

---

## Milestone A — the vocabulary and the header

✅ **A1 — the types, no logic.** `RecordHead`, `Header`, `Manifest`, `FieldSpec`, `FieldEncoding`,
`BlockIndexEntry`, `Footer`, `PspReadError`, `PspWriteError`. Doc comments carry the invariant;
no function bodies.
***Source:*** [arch §3, §4.3](../arch/psp_file_format.md).

✅ **A2 — the header: build, encode, parse, validate.** Magic, `u64` length, TOML body, sentinel.
Round-trips, and `head` on the file shows the body.
***Depends:*** A1. ***Source:*** [arch §3.2](../arch/psp_file_format.md), [spec §3.1](../spec/psp_file_format.md).

✅ **A3 — `read_header`, and refusal of a bad one.** Reads a header without a footer, which is what it
is for. Unknown major version is `UnsupportedVersion`, not a parse failure.
***Depends:*** A2. ***Source:*** [spec §6.6, §6.7](../spec/psp_file_format.md).

> **Checkpoint A: the vocabulary is fixed and a header round-trips. Pause for review.**

## Milestone B — the approximated quantities become integers upstream

**This milestone changes code outside `psp/`**, and it is here because the owner put it here
(2026-08-25). It is what keeps direct mode and psp mode bit-identical: rounding where the value is
computed means both routes see the same number, and the store then stores an integer it was handed
rather than approximating anything.

⛦ **B1 — the window's GC fraction becomes an integer type.** 0–100, rounded at construction.
Terminal per-window statistic; its consumer bins its input.
***Source:*** [spec psp_record_encoding §5.1.1](../spec/psp_record_encoding.md).
**BLOCKED, 2026-08-26 — ng does not compute this quantity anywhere.** It exists only in the frozen
production tree (`src/sample_summary/coverage.rs:283`). It is the per-window statistic
[`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md) §2.2 *proposes* and calls
"the first new accumulator step 4 would add", and that accumulator is not built. Nothing can round
where nothing computes.

⛦ **B2 — the window's mean coverage becomes an integer type**, in quarter-reads. Same argument.
***Depends:*** B1 (the pattern). ***Source:*** as B1.
**BLOCKED for the same reason as B1.**

✅ **B3 — the summed log-error becomes an integer type**, in units of 1/4,096 of a natural log.
**Own commit, do not bundle**: it enters a likelihood, so a wrong scale is a wrong genotype and not a
crash. The oracle is the calling comparison — the same records called before and after, with the
count of changed genotypes and the movement in quality scores reported.
***Depends:*** B2. ***Source:*** [spec psp_record_encoding §5.1.1](../spec/psp_record_encoding.md).

> **Checkpoint B: the three quantities are integers where they are computed, and calls are
> unchanged. Pause for review.**

## Milestone C — one record to bytes and back

No compression, no file — a record encoded into a `Vec<u8>` and decoded from it. **This is where the
correctness work is, and it is all testable in memory.**

✅ **C1 — encode and decode a record body**, every field except the chain ids: the observations, their
sequences, `read_witness`, `read_group`, the support moments, and the two counts of reads that showed
nothing. Round-trips exactly.
***Depends:*** A1, B3. ***Source:*** [arch §2](../arch/psp_file_format.md), [spec psp_record_encoding §5](../spec/psp_record_encoding.md).

✅ **C2 — the record head, and the skip.** `position_offset | reference_span | non_reference_reads |
body_bytes | body`, and a decoder that takes the head and either builds the body or advances past it.
***Depends:*** C1. ***Source:*** [spec §4.3](../spec/psp_file_format.md).

✅ **C3 — the body stands on its own.** Every difference the body carries restarts at the record, so a
skipped body strands nothing. **Own commit, do not bundle**: getting this wrong is the failure the
head exists to avoid, and it is silent — the records after a skipped one decode plausibly and wrong.
The oracle is a decode that skips every other record and still matches a full decode on those it
keeps.
***Depends:*** C2. ***Source:*** [spec §4.3](../spec/psp_file_format.md), [spec psp_record_encoding §2.3](../spec/psp_record_encoding.md).

> **Checkpoint C: a record round-trips, and a skipped record costs a pointer advance. Pause for
> review.**

## Milestone D — the psp block

✅ **D1 — cut blocks on the genomic grid.** A block ends when a position crosses into the next
multiple of `genomic_block_size_bp`, and never crosses a contig. The optional byte ceiling closes one
early.
***Depends:*** C3. ***Source:*** [spec §4.1](../spec/psp_file_format.md).

✅ **D2 — compress a block with a capped look-back window**, and record the window in the header.
***Depends:*** D1. ***Source:*** [spec §4.2](../spec/psp_file_format.md).

✅ **D3 — the streaming reader.** Pull decompressed bytes into a rolling buffer, parse a record, hand
it over, keep nothing. 16 kB buffers.
***Depends:*** D2. ***Source:*** [spec §5.1, §4.4](../spec/psp_file_format.md).

✅ **D4 — the restartable parse.** A record that straddles the rolling buffer is retried from its
start with the running state restored, and the buffer grows for a record larger than it.
**Own commit, do not bundle**: a parse that half-advances the state before failing corrupts every
record after it, plausibly. The oracle is a decode forced to refill at every possible boundary.
***Depends:*** D3. ***Source:*** [spec §8](../spec/psp_file_format.md).

✅ **D5 — every running difference resets at a block boundary.** The oracle is
restart-equals-sequential: reading from an arbitrary block gives what a full read gives from there.
**Own commit** — the failure is silent and plausible, because coverage is smooth.
***Depends:*** D4. ***Source:*** [spec §3.2](../spec/psp_file_format.md).

> **Checkpoint D: blocks cut on the grid, stream back, and a reader can start at any of them. Pause
> for review.**
>
> ✅ **Reviewed** — eight checklists across three worktrees, 21 mutations, 1,523,400 fuzzed inputs:
> [the review](../../reports/reviews/ng_psp_d4d5_2026-08-27.md),
> [the fixes](../../reports/reviews/fixes_applied_ng_psp_d4d5_2026-08-27.md). **The two spec
> corrections it owed are made**: §4.1 and §12 question 3 carry the owner's 2026-08-27 ruling
> against merging near-empty blocks instead of the withdrawn rule, and §6.7/§7 have the row for the
> reader's record ceiling that Milestone F4's `PspReadError` mapping reads them for.

## Milestone E — the chain ids

**The whole milestone fails silently.** A wrong live set gives a record an id nobody folded, and the
merge composes an allele for a read that was never there. Every step is its own commit with its
oracle green either side.

✅ **E1 — the live set and its changes.** Per position, who arrived and who left; the set restated at
each block's start.
***Depends:*** D5. ***Source:*** [spec psp_chain_id_encoding §4](../spec/psp_chain_id_encoding.md).

✅ **E2 — re-entry.** An id may go live, stop, and go live again, because a pair's mates rarely
overlap — **83 % of ids on the human sample and 91 % on tomato**. A stream that assumes one stretch
per id silently loses the second mate of nine reads in ten. **Own commit**, and the oracle counts
two-stretch ids on the fixture and asserts they all come back.
***Depends:*** E1. ***Source:*** [spec psp_record_encoding §6](../spec/psp_record_encoding.md).

✅ **E3 — the changes ride in the record head; the exception lists stay in the body.** The changes
carry state, so a skipped body must not strand them.
***Depends:*** E2, C3. ***Source:*** [spec psp_record_encoding §6](../spec/psp_record_encoding.md).

✅ **E4 — the derived residual, and the check that guards it.** The residual observation's ids are the
live set minus the others; reads the depth cap discarded and reads that produced no observation must
not leak in. **Own commit.** The oracle is the inequality the walk already asserts: a derived list is
at most the observation's read count and at least half of it, since at most two mates share an id.
***Depends:*** E3. ***Source:*** [spec psp_chain_id_encoding §5](../spec/psp_chain_id_encoding.md).

> **Checkpoint E: chain ids round-trip exactly, including ids that go live twice. Pause for review.**
>
> ✅ **Reviewed** — four rounds of eight checklists (E1–E4), each in isolated worktrees:
> [E1](../../reports/reviews/ng_psp_e1_2026-08-27.md) ·
> [E2](../../reports/reviews/ng_psp_e2_2026-08-27.md) ·
> [E3](../../reports/reviews/ng_psp_e3_2026-08-27.md) ·
> [E4](../../reports/reviews/ng_psp_e4_2026-08-28.md), with their fix reports beside them.
> **What Milestone E leaves unmeasured is what the column costs**: spec
> [`psp_chain_id_encoding.md`](../spec/psp_chain_id_encoding.md) §7 wants the compressed bytes, the
> reader's peak resident live set, and a full scan against a seek, on both corpora. None of that
> can be taken until Milestone F opens a file, and every figure quoted through E is the spec's own,
> from a prototype over alignments rather than from this writer.

## Milestone F — the container

✅ **F0 — a field's cardinality joins the manifest.** Added by the owner on 2026-08-28, when the
question A1 left open was raised at the milestone that closes it. Spec §4.5 asks the manifest to
carry each field's name, cardinality and encoding; it carried two of the three. **This is the last
step where adding the third is free** — from F1 a manifest change costs a format version.
The key is `shape`, not `cardinality`: production's store already uses that word for how often
a field appears, and ng's manifest deliberately does not carry that — **the owner ruled against
carrying it, 2026-08-28**, because the fields queued to be added appear once per record and a
format-version bump is free until a psp file exists that someone keeps.
***Depends:*** nothing in F. ***Source:*** [spec §4.5](../spec/psp_file_format.md).

✅ **F1 — the block index's codec.** Its bytes, its checksum, and the refusals that keep a seek
honest; *building* one while writing is F3's and *reading* one at open is F4's.
***Depends:*** D5. ***Source:*** [arch §3.4](../arch/psp_file_format.md).

✅ **F2 — the footer**, magic last, with the index checksum.
***Depends:*** F1. ***Source:*** [arch §3.4](../arch/psp_file_format.md), [spec §3.3](../spec/psp_file_format.md).

✅ **F3 — `create`, `push`, `finish`.** `push` rejects an out-of-order record; `finish` writes index,
trailer, footer, then flushes, surfaces the buffered writer's errors, and syncs.
***Depends:*** F2, E4. ***Source:*** [spec §6.3](../spec/psp_file_format.md).

✅ **F4 — `open`, and refusal.** Footer, index, header, no block touched. A file with no footer is
`Incomplete`; a window larger than the reader's budget is `WindowTooLarge`.
***Depends:*** F3. ***Source:*** [spec §6.2, §6.7](../spec/psp_file_format.md).

> **Checkpoint F: a file can be written, closed, reopened and read end to end. Pause for review.**
>
> ✅ **Reviewed** — five rounds of eight checklists (F0–F4), each in isolated worktrees:
> [F0](../../reports/reviews/ng_psp_f0_2026-08-28.md) ·
> [F1](../../reports/reviews/ng_psp_f1_2026-08-28.md) ·
> [F2](../../reports/reviews/ng_psp_f2_2026-08-28.md) ·
> [F3](../../reports/reviews/ng_psp_f3_2026-08-28.md) ·
> [F4](../../reports/reviews/ng_psp_f4_2026-08-28.md), with their fix reports beside them.
> **Four Blockers, and every one was a test that could not fail rather than code that was
> wrong** — a suite that proved the encoding table consistent with itself, an ordering scan
> exercised on one pair, a writer whose recovered failure produced a valid file with records
> missing, and a sections rule tested on one side. **What Milestone F leaves unmeasured is what an
> open sample costs.** F4 ships a refusal threshold carved out of the 500 kB budget by arithmetic
> on the spec's figures, and the block index alone is about 336 kB of that budget for a whole
> genome — H4 is where both are measured.

## Milestone G — the rest of the surface

✅ **G1 — `records_from`**, turning a coordinate into a block with one index lookup.
***Depends:*** F4. ***Source:*** [spec §6.2](../spec/psp_file_format.md).
> **Reviewed** — nine checklists across three agents:
> [the review](../../reports/reviews/ng_psp_g1_2026-08-28.md),
> [the fixes](../../reports/reviews/fixes_applied_ng_psp_g1_2026-08-28.md). **The design question
> the step was handed is settled by the documents**: spec §6.2 matches the coordinate against
> where records *start* and arch §4.1 starts reading at a block's first record, so `records_from`
> is block selection and not an overlap query. **The Blocker was the index search**: two blocks
> may share a first position, and *the last block starting at or before the coordinate* entered
> that run at its end and lost the earlier block's records with no error.

✅ **G2 — `records_where`**, the head-driven skip as a public iterator.
***Depends:*** F4, C2. ***Source:*** [spec §6.2](../spec/psp_file_format.md).
> **Reviewed** — nine checklists across two agents:
> [the review](../../reports/reviews/ng_psp_g2_2026-08-28.md),
> [the fixes](../../reports/reviews/fixes_applied_ng_psp_g2_2026-08-28.md). No Blocker; **every
> Major was about what the tests did not hold**. Nothing held the claim the type exists for — a
> decode-then-discard implementation passed all 355 tests — and the property that a declining walk
> is a weaker reader of damage than a full one was correct and unwritten: **of 93 single-byte
> corruptions a full walk refuses on one small file, a walk declining every body accepts 72**.

☐ **G3 — `replace_trailer`**, rewriting from the trailer's offset and leaving blocks and index alone.
***Depends:*** F3. ***Source:*** [spec §6.5](../spec/psp_file_format.md).

☐ **G4 — `append`**, truncating at the index offset and keeping the manifest. Enforces coordinate
order across the seam and **fails on a manifest it cannot honour**.
***Depends:*** F3. ***Source:*** [spec §6.4](../spec/psp_file_format.md).

> **Checkpoint G: all five operations exist. Pause for review.**

## Milestone H — the oracles, and the numbers

☐ **H1 — parity with the prototype.** A store written by `psp/` and one written by
`examples/psp_row_stream_roundtrip.rs` yield the same records, compared with the strictness each
field requires: integers, sequences, witnesses and chain-id lists **exactly**; the fixed-point fields
inside their own step. **A blanket tolerance would pass while a chain-id list was being corrupted.**
***Depends:*** G4. ***Source:*** [arch §8](../arch/psp_file_format.md).

☐ **H2 — an interrupted write is refused.** Kill a writer before `finish`; the reader rejects the
file rather than reading the blocks that reached disk.
***Depends:*** F4. ***Source:*** [spec §6.3](../spec/psp_file_format.md).

☐ **H3 — worker-count invariance.** One sample gathered at 1, 2, 4, 8 and 16 workers gives
byte-identical files apart from the header's timestamp.
***Depends:*** F3. ***Source:*** [spec §7](../spec/psp_file_format.md).

☐ **H4 — the memory number, measured not asserted.** N samples open and walked in lockstep at 1, 8,
62 and 5,000, peak resident against the 500 kB per-open-sample budget. **Report the per-sample slope,
not one point** — the budget is about the slope.
***Depends:*** G2. ***Source:*** [spec §5.2](../spec/psp_file_format.md).

☐ **H5 — what the head is worth at depth, which nobody knows.** A walk keeping one record in a
hundred against one keeping all, on tomato at 3 reads a position **and on the 279-reads-a-position
sample**. The chain-id changes ride in the head and grow with depth — 0.432 bytes a position at 11.4
reads, 6.42 at 293 — so the head grows while the skip's value shrinks, and **how much of the 2.06×
survives is an open question the arch doc records** (arch §7).
***Depends:*** H4. ***Source:*** [spec psp_record_encoding §6](../spec/psp_record_encoding.md).

> **Checkpoint H: the store is proven against the prototype and its costs are measured. Pause for
> review.**

---

## Verification summary

| milestone | proven by |
|---|---|
| A — vocabulary and header | header round-trip; `head` on the file shows the body; a bad version is refused |
| B — quantities as integer types | the calling comparison: same records called before and after, changed genotypes and quality movement reported |
| C — the record codec | in-memory round-trip, exact per field; a decode skipping every other record matches a full one |
| D — the psp block | restart equals sequential from an arbitrary block; a decode forced to refill at every boundary |
| E — the chain ids | exact chain-id lists; two-stretch ids counted on the fixture and all recovered; the derived residual's count inequality |
| F — the container | write, close, reopen, read end to end; a file with no footer is refused |
| G — the surface | each operation on a fixture; `append` refuses a manifest it cannot honour |
| H — the whole store | parity with the prototype record for record; peak resident against the budget at 5,000 open samples |

## Out of scope — next plans

- **Wiring the store into the run** — [`../spec/run_streaming.md`](../spec/run_streaming.md) owns the
  run objects that would write and read one.
- **The trailer's contents** — follows the statistical work
  ([spec §3.4](../spec/psp_file_format.md)).
- **A `ChainId` newtype** — arch §7, to whoever owns ng's chain-id minting.
- **Correcting [`../spec/run_streaming.md`](../spec/run_streaming.md) §7.2**, whose "tens of
  kilobytes" the 500 kB budget supersedes, and
  **[`../arch/module_layout.md`](../arch/module_layout.md)'s** note that ng has no `.psp` yet.
