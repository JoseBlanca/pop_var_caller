# ng — storing the chain ids in a psp file, and the experiment that decides how

*Status: design spec, 2026-08-18; **the experiment it defined is no longer needed — the owner
settled it on 2026-08-25: the changes-only form (§4) ships and the intermediate distances form is
not built.** What survives here is §4's design, §5's traps and §8's rejected alternatives. §3's three
arms and §7's experiment are kept as the record of what the decision was taken against, not as work
to do. **One thing §4 did not anticipate**: the records are skippable now
([`psp_record_encoding.md`](psp_record_encoding.md) §2.3), so the live-set changes go in a record's
head and only the exception lists stay in its body — that document's §6 has it. It is a
piece of the psp encoding, whose spec is still owed
([`run_streaming.md`](run_streaming.md) §10); that document owns the byte layout, blocks,
compression, the index and versioning, and this one owns a single column inside it. Reads on
[`cohort_merge.md`](cohort_merge.md) §14 Q2 (the ruling that made the column big) and
[`run_streaming.md`](run_streaming.md) §12.1 (the file must not depend on worker count).*

*Downstream: [`../arch/psp_file_format.md`](../arch/psp_file_format.md) §3 and
[`../impl_plan/psp_file_format.md`](../impl_plan/psp_file_format.md) Milestone E, which is where this
document's silent-failure traps are built and guarded.*

---

## 1. Why this column suddenly matters

**A chain id names the read that produced a piece of evidence.** One `u64` per read, or per read
*pair* with its mates collapsed onto one id, allocated monotonically from zero and never reused
(`src/ng/locus_generation/pileup/chain_id_allocator.rs:1-10`).

Until 2026-08-17 both callers stored the ids of the reads that **disagreed** with the reference and
dropped the rest. Production drops them positionally — its REF bucket, `allele_index == 0` — and its
own comment says what that saves: **the dropped ids are about 96.6% of all chain ids on real
cohorts** (`src/pileup/walker/open_record.rs:150-160`).

**The owner's ruling of 2026-08-17 reverses it, and the cohort merge is why.** A cohort locus can
span several of one sample's records, and a read's allele over that locus is what it showed at each
of them — so a read that covered a position and agreed with the reference has to be distinguishable
from one that never reached it. Unnamed, those are the same absence. The rule the merge now works
to admits no third case: *either we know the read covered the whole locus and its allele is
elongated with what it showed, or we know it did not and it is removed as evidence; being unable to
decide is an error that must never happen* ([`cohort_merge.md`](cohort_merge.md), Preconditions).
ng's mint therefore names every read it folds
(`src/ng/locus_generation/pileup/open_record.rs`, `KeyedObservation::chain_ids`).

So a column that carried about 3.4% of the ids now carries all of them, and it does so at **every
position of the genome rather than at the positions where something varies.** The owner accepted
that cost against the alternative — discarding every read that cannot prove it spanned a locus,
which at a few reads a position is most of the evidence. This document is about paying less for it
in the file without changing what the merge is given.

### 1.1 The one property that makes a cheaper encoding possible

**The merge never uses an id's value. It only asks whether two records name the same read.**
Identity within one cohort locus is the whole requirement, and a locus never crosses a segment
boundary ([`run_streaming.md`](run_streaming.md) §4.3). So the stored form is free to be anything
that preserves equality over a window of one locus — the compression is a storage concern that
touches no calling code.

---

## 2. What this document is for

**Goals**

1. Define a **differential encoding** of the chain-id column precisely enough to build.
2. Define the **experiment** that decides whether it ships: what to implement, on what data, what
   to measure, and against which baselines.
3. Record what was already tried and rejected, so the branch does not re-explore it.

**Non-goals**

- **The psp encoding as a whole.** Byte layout, block sizing, compression settings, the index, the
  trailer and versioning are [`run_streaming.md`](run_streaming.md) §10's spec, not this one.
- **Changing what the merge consumes.** The in-memory `SequenceObservation::chain_ids` stays as it
  is; this is about the file. If the experiment ships, the reader reconstructs the same lists.
- ~~**Deciding the outcome.**~~ **Decided 2026-08-25 by the owner, without the experiment: the
  changes-only form ships.** This document's job is now to say what it is and what will bite the
  coder, not how to choose it.

**It does not** re-open the ruling, touch the in-memory shape, or propose a read-major file.

---

## 3. The current encoding, and the honest baseline

**Today's column is a list per allele row: a LEB128 length, then one raw little-endian `u64` per
id** (`encode_list_column`, `src/psp/block.rs:327-334`; the CSR variant the writer actually calls,
`:352-361`). The payload sits inside a **zstd-compressed block** (`src/psp/block.rs:5-31`), so the
bytes on disk are already smaller than the bytes encoded.

ng has no psp writer yet; production's is the reuse candidate and the shape a straightforward port
inherits.

**The experiment needs three arms, not two, or it will overstate the result.** Monotone `u64`s
compress; a differential scheme compared only against raw-and-zstd would be credited with a saving
that a much simpler change also delivers.

| arm | what it stores | why it is in the experiment |
|---|---|---|
| **A — as today** | LEB128 length + raw LE `u64` per id, zstd on the block | the port's default; the number everything else is measured against |
| **A′ — delta + varint** | ids sorted ascending, first as a varint, then deltas as varints | the cheap change: no reader state, a few lines, and it exploits monotone allocation |
| **B — differential** | §4 | the proposal |

---

## 4. The differential encoding

**The idea.** A read of length L is named at every one of the L positions it
covers, so the column stores each read about L times. What changes between adjacent positions is
only which reads arrived and which left. Store *that* stream instead, and a read costs about two
entries for its whole length rather than L. Everything else about the column then follows from a
second observation: at nearly every position nearly every read agrees with the reference, so the
reference observation's id list is **the reads that are live minus the reads the other observations
name** — and need not be stored at all.

So the column becomes two things:

- **The live-read stream, per position.** The ids that begin covering this position and the ids
  that stop. Arrivals are ascending and near-dense, so delta-varint applies; departures reference
  ids already live and can be spelled as positions in the live set rather than as ids.
- **The exception lists, per record.** Every observation *except* the designated residual one keeps
  an explicit id list, as today. That is the ~3.4% of ids production used to store, so this half is
  no worse than what shipped before the ruling.

**The residual observation is named, not inferred.** The writer marks which observation of a record
carries the residual; the reader computes its ids as the live set minus every other observation's.
Naming it costs one small integer per record and removes any guessing about which "the reference
one" is when observations split by witness and by read group.

**Random access is preserved by restating, not by chaining across blocks.** Each block begins with
the full live set; the differential stream runs only inside the block. A reader that seeks to a
block reconstructs from its start, which is what it already does for every other column. The block
cut rule is unchanged and stays a function of the observation stream alone
([`run_streaming.md`](run_streaming.md) §12.1) — the encoding must not make a cut depend on how the
live set happens to look, or the file stops being worker-count invariant.

### 4.1 What it should cost, as arithmetic rather than measurement

At depth `d` and read length `L`, per position per sample, **before compression**:

- **arm A** carries `d` ids: `8d` bytes, plus framing. At `d = 30`, 240 bytes.
- **arm B** carries the arrivals and departures, about `2d/L` entries — at `d = 30`, `L = 100`,
  fewer than one entry a position — plus the exception lists at roughly `0.034 × 8d` bytes.

**These are estimates from the shape of the data, not measurements**, and two things will cut them
down: zstd already exploits much of the redundancy in arm A, and real reads do not all have one
length. The experiment exists because the gap between this arithmetic and the compressed bytes is
exactly what is unknown.

---

## 5. What will bite the coder

- **A chain id is a fragment, not a read, so its coverage is not one interval.** Mates collapse onto
  one id (`chain_id_allocator.rs:1-10`), and the two mates need not be adjacent — the unsequenced
  gap between them means an id can go live, stop, and go live again. An arrivals-and-departures
  stream that assumes one interval per id will silently lose the second mate.
- **Records are not positions, and one sample's records can overlap.** A record's footprint spans
  the event it holds, and the table keeps heterogeneous spans open at once — *"a wide deletion
  record may stay open while shorter records open and close around it"*
  (`src/ng/locus_generation/pileup/open_record.rs:1463-1466`). The live-set stream is indexed by
  position; the id lists belong to records. The mapping between them is the part to get right
  before writing any bytes.
- **Two classes of read are counted and never named**, and they must not leak into the residual:
  `reads_without_observation` and `reads_discarded_by_cap`
  (`src/ng/locus_generation/mod.rs`, `SampleLocusObservations`). A read the depth cap discarded is
  not in any observation, so if it is in the live set the residual gains an id nobody folded.
- **The residual arithmetic fails silently.** If the reader derives one id too many, the reference
  observation gains a read that does not exist, and the merge will happily compose an allele for
  it. Give the reader something to check against: an observation's derived id count must be at most
  its `num_obs` and at least half of it, since at most two mates share an id — the same inequality
  the walk's own differential against production already asserts
  (`src/ng/locus_generation/pileup/parity.rs:2245`).
- **Compare after compression, never before.** The column's cost on disk is post-zstd. A raw-byte
  comparison will make arm B look better than it is, and it is the disk and decode numbers the
  decision turns on.

---

## 6. Cross-cutting concerns

**Memory.** Arm B gives the reader state it does not have today: the live set, which is depth-sized
per sample, held while a block is scanned. That is small beside the observations themselves, but it
is new, and at the top of the committed cohort range it is `samples × depth` ids resident. The
experiment must report the reader's peak, not only the file size.

**Errors.** A differential stream can be internally inconsistent in ways a list cannot — a departure
for an id that is not live, a residual that comes out negative. Both are corrupt-input failures and
belong with the psp reader's error type, not with the merge; the merge must never see a
half-reconstructed list.

**Concurrency.** None of this changes: blocks stay independently decodable, and the live set is
per-reader state inside one block.

---

## 7. The experiment

**When.** On a branch, after [`run_streaming.md`](run_streaming.md) §10's psp encoding spec lands
and ng has a writer to modify. Before that there is nothing to measure.

**What to build.** One reader/writer interface, three implementations behind it (A, A′, B from §3),
so the same walk can write all three and the same merge can read all three. The merge's output must
be **identical across the three** — that is the correctness gate, and it is the same
byte-identity oracle the module already uses.

**On what data.** Both committed corners, one chromosome each: the **63-sample tomato panel at
about three reads a position**, and **HG002 at its own depth**. They differ in the two variables the
arithmetic in §4.1 turns on — depth and cohort size — so a scheme that wins on one and loses on the
other will show it.

**What to report**, per arm and per corner:

1. the chain-id column's **compressed** bytes, and the whole file's;
2. **write** wall time and peak RSS;
3. **read** wall time for a full sequential scan, and for a seek to a random block — the second is
   where arm B pays;
4. the **merge's end-to-end** wall time on the same input, which is what the operator feels;
5. the reader's **peak resident live set**.

**What decides it.** The owner, on those numbers. This spec deliberately sets no threshold: the
trade is a file-size saving against reader state and decode complexity, and how much complexity is
worth how many bytes is not a question a document written before the measurement should pretend to
answer.

---

## 8. Alternatives, and why they are not in the experiment

- **A local index instead of the identifier** — number the reads live in a block 0, 1, 2… and store
  that, one byte where there is one `u64` now. **Rejected on production's experience** (the owner,
  2026-08-18): it was tried there, the code came out considerably more complex, and it cost
  performance. Recorded rather than re-derived; nothing in ng makes it cheaper than it was.
- **Dropping the reference reads' ids**, as production does. This is the thing the 2026-08-17
  ruling reversed, and reversing it back would put the merge into the state where it cannot tell a
  read that agreed from a read that was never there.
- **A read-major file** — each read stored once with its events, which removes the redundancy
  completely. It is the wrong shape for everything else the psp does: the streaming reader, the
  two-phase decode and the column access all want positions.

---

## 9. Deferred, with a recommended home

- **The byte layout of the two streams** — how the live set and the exception lists sit inside a
  block, and which column each is — to the psp encoding spec
  ([`run_streaming.md`](run_streaming.md) §10), which owns block structure. This document fixes what
  is stored and why, not where the bytes go.
- **Whether the in-memory shape should change too** — a narrower id, an inline buffer at low depth —
  to the measurement the cohort merge spec already owes for what one observation costs
  ([`cohort_merge.md`](cohort_merge.md) §8). Nothing here depends on it.

---

## 10. Open questions

1. **Does the saving survive zstd?** — OPEN, and it is the question the experiment exists for.
   Arm A's ids are monotone and highly compressible already. *Leaning:* a real saving remains,
   because zstd cannot know that consecutive positions repeat the same set, but the size of it is a
   guess until measured. **Settled by:** §7's arm A against arm B, compressed.
2. **Is arm A′ enough?** — OPEN. If delta-varint alone gets most of the way, it wins on complexity
   by a wide margin: no reader state, no residual arithmetic, no new error class. **Settled by:**
   the same run, which is why A′ is in it.
3. **Should the residual observation be derived at all, or should every observation keep its list?**
   — OPEN. Deriving it is where most of the saving is and where the silent failure is. *Leaning:*
   derive, with the count check of §5. **Settled by:** the size difference between the two variants
   of arm B, which the experiment can report from one implementation by writing both.
4. **How often does an id go live twice** — the mate gap of §5? — OPEN, and it decides whether the
   stream needs a re-entry form or can treat every arrival alike. **Settled by:** counting it on the
   two corners' alignments, which needs no encoding work and could be done before the branch.
