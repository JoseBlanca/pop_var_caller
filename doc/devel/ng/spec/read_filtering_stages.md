# ng — read filtering in stages: two filters and a converter

**Date:** 2026-08-03 · **Status:** design draft, **no code yet**. Four open questions for the
owner in §8.
**Supersedes** the single-type shape in [`read_filtering.md`](read_filtering.md) §5. That
document's filter *policy* — which nine filters run, their thresholds, their order — is
unchanged and stays there; this one is only about how the work is divided.
**Companions:** [`alignment_cursor.md`](alignment_cursor.md) (the reader above this),
[`ref_seq.md`](ref_seq.md) (the reference the last filter consults).
**Code-facing companion:** [`../arch/read_filtering_stages.md`](../arch/read_filtering_stages.md).

---

## 1. What this is

Today one type, `ReadFilter`, does three jobs in one loop
([`filtering.rs:895-951`](../../../../src/ng/read/filtering.rs#L895)): it rejects records on
their flag and mapping quality, decodes the survivors into ng's own read type, then rejects
some of *those* on their length, their CIGAR and how well they match the reference.

Splitting it gives three named things:

```
record readers → region narrowing        raw records, one reused buffer   ← unchanged
      ↓  admission on the raw record     flag, mapping quality — no reference
      ↓  decode                          raw record → AlignedRead
      ↓  admission on the read           length, CIGAR, mismatch fraction
      ↓  the cursor keeps what survives
```

**Goals.**

- A converter that only converts, reachable on its own.
- A first filter that needs **no reference**: today nothing can filter on flag and mapping
  quality without supplying a `RawRefSeq` and paying a probe fetch for every contig in the
  header ([`filtering.rs:688-702`](../../../../src/ng/read/filtering.rs#L688)).
- The pre-decode/post-decode boundary expressed in types rather than in a comment.
- **Byte-identical output.** Same reads kept, same drops charged to the same reasons.

**Non-goals.**

- **No change to the filter policy.** The nine filters, their thresholds, their evaluation
  order and the `DropReason` set stay exactly as `read_filtering.md` settled them.
- **No change to what the cursor keeps.** It keeps decoded, filtered reads, and
  [`alignment_cursor.md`](alignment_cursor.md) §5 gives the reason — a read served by about a
  dozen consecutive regions is transformed once, not a dozen times. This split does not reopen
  that.
- **No performance change sought.** Nothing here should get faster or slower. If a measurement
  moves, something is wrong.
- **Not a new module.** This is a reorganisation inside `read/filtering.rs`; see §6.

**It does not** introduce a trait with competing implementations, change the error taxonomy's
*meaning*, add a configuration knob, or touch `AlignedRead`.

## 2. Why the work is divided the way it is

**The division is not a matter of taste — it falls out of what each filter reads.** Six of the
nine filters read two integers off the undecoded record. The other three read fields that the
decode *produces*.

| filter | reads | where that comes from |
|---|---|---|
| #1 duplicate, #3 supplementary, #4 secondary, #5 unmapped, #6 QC-fail | the SAM flag | on the raw record |
| #2 low mapping quality | the SAM mapping quality | on the raw record |
| #7 too short | the read's sequence length | decode uppercases the sequence into a `Vec<u8>` |
| #9 bad CIGAR | the read's CIGAR | decode converts noodles' CIGAR into ng's `CigarOp` |
| #8 high mismatch fraction | CIGAR, sequence, base qualities, position, contig — plus a reference fetch | same |

Cited: `verdict_pre_decode` ([`filtering.rs:210-238`](../../../../src/ng/read/filtering.rs#L210))
and `verdict_post_decode` ([`filtering.rs:269-322`](../../../../src/ng/read/filtering.rs#L269));
the CIGAR conversion and the uppercase are in `decode_record`
([`aligned_read.rs:102`](../../../../src/ng/read/aligned_read.rs#L102) and
[`:107`](../../../../src/ng/read/aligned_read.rs#L107)).

**So "all nine filters run on raw records" is not available.** Reaching it would mean either
re-implementing the mismatch rule and the CIGAR scan against noodles' types — a second copy of
two rules, which is the failure this module guards against hardest — or running the uppercase
and the CIGAR conversion before the filter, which is the decode under another name. Two filter
stages with the converter between them is the shape the code already has; this makes it visible.

**The decode is not a formality.** Besides the two conversions above it computes the adaptor
boundary from the flag, the mate's placement and the template length
([`aligned_read.rs:112-120`](../../../../src/ng/read/aligned_read.rs#L112)). That is real work,
and it is work a filter has no business doing.

## 3. What will bite you

**The raw record is one reused buffer, so a stage cannot return one.** `RecordSource::read_next`
fills a caller-owned buffer in place and the whole pass allocates exactly one record
([`filtering.rs:366-382`](../../../../src/ng/read/filtering.rs#L366)). There is one live raw
record at any moment. So the first stage cannot have the shape `fn(RawRecord) -> Option<RawRecord>`
without cloning a `RecordBuf` per record — measured elsewhere in this codebase at ~680 bytes
across seven allocations
([`record_reader/container.rs:13-17`](../../../../src/ng/read/input/record_reader/container.rs#L13)).
It has to be a **verdict about the buffer**, not a value returned from it. The existing trait
doc calls this the lending-iterator problem and says the buffer shape was chosen for exactly
this reason ([`filtering.rs:359-365`](../../../../src/ng/read/filtering.rs#L359)).

**A drop must be charged to a read group, and the first stage runs before any read exists.**
This is already solved and easy to miss: the raw record carries its read group, stamped by the
region narrowing ([`region_records.rs:222`](../../../../src/ng/read/input/region_records.rs#L222)),
and `RawRecord::read_group` exists for the tally rather than for any filter
([`filtering.rs:350-357`](../../../../src/ng/read/filtering.rs#L350)). The first stage can charge
its own drops correctly without decoding anything.

**Filter order inside stage two is deliberate and is not the spec's table order.** It runs #7,
then #9, then #8, so that a read failing the cheap checks never pays the reference fetch — and
so a read failing both #9 and #8 is charged to the root cause rather than the symptom
([`filtering.rs:240-268`](../../../../src/ng/read/filtering.rs#L240)). One test pins it, using a
fixture record that fails both
(`a_walk_charges_every_drop_reason_by_hand_count`, `input/cursor.rs`). Preserve the order or
that test fails, which is what it is for.

**Three things stop the same way and one of them must not restart.** The filter's state is
`Running` / `EndOfInput` / `Failed`
([`filtering.rs:646-658`](../../../../src/ng/read/filtering.rs#L646)). A long-lived filter reaches
`EndOfInput` at the end of *every* region, because the region narrowing reports a region boundary
the only way a record source can, and repositioning undoes it. `Failed` is permanent. Whatever
owns the stages has to keep that distinction; collapsing it into one flag is how a cursor ends up
reading on from a file already known to be broken.

**The `Option<ReadGroupId>` on the buffer is a guard, not a convenience.** It is cleared when the
buffer is handed out and set again after each record is filled, so a source that forgets to stamp
produces a refusal rather than attributing the read to the previous record's library
([`filtering.rs:479-497`](../../../../src/ng/read/filtering.rs#L479)). Any restructure that hands
the buffer between stages must keep that order.

## 4. What this buys

Three things, and only the first is about tidiness.

**The reference stops being a precondition for filtering at all.** `ReadFilter::new` takes a
`RawRefSeq` and fetches a zero-length window on every contig in the header before it will build
([`filtering.rs:688-702`](../../../../src/ng/read/filtering.rs#L688)). With the stages split,
anything that only wants flag and mapping-quality filtering — a coverage histogram, an
insert-size pass, a read-group pre-pass — can have it without a reference. **Nothing in ng needs
this today**; it is a capability the current shape forecloses, not a request anyone has made.

**Two vestiges dissolve rather than needing their own cleanup.** `with_validated_contigs` (the
probe-free constructor, [`filtering.rs:746`](../../../../src/ng/read/filtering.rs#L746)) and
`ReadFilterBuffers` (a lend-and-reclaim seam) both exist because a filter used to be built per
region by a pooled caller. That caller is gone; `ReadFilter::new` is the only route into
`with_validated_contigs`, and the only value anything passes for the buffers is
`Default::default()`. If the reference enters only at stage two, both stop having a reason to
exist.

**The ordering becomes a type instead of a comment.** Today the only thing stopping someone
hoisting `decode()` above the flag checks is prose. That change would keep every test green and
pay a full decode — name, sequence, qualities, CIGAR, adaptor boundary — for every duplicate,
secondary and low-mapping-quality read.

## 5. What could go wrong with it

**You may end up with the same type, renamed, plus three more.** Someone still has to own the
three-way stop, the error taxonomy, and the loop that pulls until something survives. A driver
owning all three stages will look a great deal like today's `ReadFilter`. **The split is only
worth doing if the parts are separately nameable and separately testable** — if the design ends
up with a driver that no one can use the pieces of, it has added three types and bought nothing.
That is the question §8's first two items exist to settle before code.

**This is the module every dump's byte-identity runs through.** Step 1 decides which reads reach
every locus. The evidence that a restructure changed nothing is not the unit suite — it is the
four acceptance dumps (§7).

## 6. Where it lives

`src/ng/read/filtering.rs`, unchanged. The module-layout rule is that a step with no competing
implementations is a file rather than a folder
([`../arch/module_layout.md`](../arch/module_layout.md) principle 1a), and splitting one type
into three inside a file does not create alternatives to sit side by side. If the file grows past
readability once the three are in it, that is a reason to make `read/filtering/` a folder — a
mechanical move, decided at code time, not here.

## 7. How we know it works

**The unit suite is necessary and not sufficient.** 45 tests live in `filtering.rs` and most
drive `ReadFilter` as one thing; they will need re-pointing at whichever type keeps their
subject.

**The evidence is output identity on real data:**

| anchor | what it proves |
|---|---|
| `ng_generic_loci_dump` and `ng_ssr_loci_dump`, HG002 chromosome 21 (BAM) | 251,792 and 4,406 lines, byte-identical |
| the same two on tomato SL4.0ch01 (CRAM) | 1,718,914 and 11,945 lines, byte-identical |
| `ng_generic_walk_probe`, chromosome 21 | `loci=236081 observations=251786 reads_admitted=54709` |

A read admitted or dropped differently shows up in all of them. `reads_admitted` is the direct
one: it is step 1's keep count for the run.

**Two properties need a test each, because output identity cannot see them:**

- **The first stage builds without a reference.** The capability §4 claims is not exercised by
  any existing caller, so it needs a test that constructs it with no `RawRefSeq` at all —
  otherwise the claim is untested and will quietly stop being true.
- **The decode still happens only for survivors.** This is a *work* property, not an answer
  property: hoisting the decode changes no output. Assert it directly — a fixture of records that
  all fail stage one, and a converter that has been asked for nothing.

## 8. Open questions — for the owner, before any code

**Q1. Is the first stage a type, or does the free function stay a free function?**
`verdict_pre_decode(flag, mapq, &config)` is already a pure function
([`filtering.rs:210`](../../../../src/ng/read/filtering.rs#L210)). It needs to become a type only
if it must hold something — the configuration, or its own tally. *Leaning: a type, because it has
to charge drops to read groups and that means owning a tally.* Settles §5's "same type renamed"
risk: if stage one is only a function, there is no second type and less to justify.

**Q2. Who owns the per-read-group tally when two stages both drop reads?** One shared value the
driver owns and lends, or one per stage, merged on read. *Leaning: one, owned by the driver* —
the tally is read as a whole (`SampleCursor::read_group_counts`) and merging two on every read
would be work for nothing. But this is the decision that most shapes whether the stages are
usable apart, so it is the owner's.

**Q3. What is the driver called, if `ReadFilter` becomes the post-decode stage?** The repo's
naming rule is that a name says what the value *is*. Candidates: keep `ReadFilter` for the
composite and name the stages for what they admit; or make the stages `RecordAdmission` /
`ReadAdmission` and the driver `FilteredReads`. *No leaning — this is a taste call and it is
yours.*

**Q4. Does the cursor consume a driver, or compose the three stages itself?** It holds a
`ReadFilter<RegionRecords, R>` today
([`input/cursor.rs:241`](../../../../src/ng/read/input/cursor.rs#L241)) and drives it through
four methods: `next`, `has_failed`, `restart_after_end_of_input`, `source_mut`. If the cursor
composes the stages, that surface grows and the cursor gains knowledge of the split. *Leaning:
a driver, so the cursor's surface is unchanged* — which also means the cursor's tests do not move.

## 9. Deferred, with a home

- **Removing `with_validated_contigs` and `ReadFilterBuffers`.** They dissolve as a *consequence*
  of this split (§4), so they belong to whichever step lands stage two — not to a separate
  cleanup, and not to this document to schedule.
- **A whole-file read path.** The capability §4 opens (filtering without a reference) is not the
  same as a way to *read* a whole file, which ng does not have: `AlignmentFile::open` requires an
  index. If one is ever wanted it arrives as a `RecordReader` arm, beside the others
  ([`../arch/alignment_cursor.md`](../arch/alignment_cursor.md) §1.3).
- **`read/filtering/` as a folder.** See §6; decided at code time on the file's real size.

## 10. Reuse map

Nothing here is new code in the sense of new logic. Every rule already exists and is being
re-homed.

| what | existing code | how it is reused |
|---|---|---|
| the six flag/mapping-quality filters | `verdict_pre_decode` [`filtering.rs:210`](../../../../src/ng/read/filtering.rs#L210) | **move**, body unchanged |
| the three decode-dependent filters | `verdict_post_decode` [`filtering.rs:269`](../../../../src/ng/read/filtering.rs#L269) | **move**, body unchanged — including the #7/#9/#8 order |
| the converter | `RawRecord::decode` [`filtering.rs:349`](../../../../src/ng/read/filtering.rs#L349) → `decode_record` [`aligned_read.rs:67`](../../../../src/ng/read/aligned_read.rs#L67) | **reuse as-is**; it is already a standalone function |
| the raw record and its buffer contract | `RawRecord` / `RecordSource` [`filtering.rs:334`](../../../../src/ng/read/filtering.rs#L334), [`:366`](../../../../src/ng/read/filtering.rs#L366) | **reuse as-is** — `read_group()` already exists for the tally |
| the tally | `ReadFilterCounts` / `ReadGroupCounts` [`filtering.rs:122`](../../../../src/ng/read/filtering.rs#L122), [`:661`](../../../../src/ng/read/filtering.rs#L661) | **reuse as-is**; Q2 decides who holds it |
| the three-way stop | `FilterState` [`filtering.rs:646`](../../../../src/ng/read/filtering.rs#L646) | **reuse as-is**, on the driver |
| the error taxonomy | `ReadFilterError` [`filtering.rs:578`](../../../../src/ng/read/filtering.rs#L578) | **reuse as-is** — its three variants already name the three stages: `Source`, `Decode`, `Reference` |

**The parity oracle is the four dumps** (§7). There is no second implementation to differ from —
the oracle is this code before the change.
