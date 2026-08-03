# ng — read filtering in stages: two filters and a conversion

**Date:** 2026-08-03 · **Status:** design settled, **no code yet**. Both questions §9 raised are
resolved (owner, 2026-08-03).
**Revises** [`read_filtering.md`](read_filtering.md) §5, which gave step 1 a single type. That
document says which nine filters run, what their thresholds are, in what order they are
evaluated, and what each drop is called. All of that is unchanged and stays there. This one is
only about how the work is divided and where each piece lives.
**Companions:** [`alignment_cursor.md`](alignment_cursor.md) (the reader below this),
[`ref_seq.md`](ref_seq.md) (the reference the last filter consults).
**Code-facing companion:** [`../arch/read_filtering_stages.md`](../arch/read_filtering_stages.md).
**Build order:** [`../impl_plan/read_filtering_stages.md`](../impl_plan/read_filtering_stages.md).

---

## 1. What this is

Today one type, `ReadFilter`, does three jobs in one loop
([`filtering.rs:895-951`](../../../../src/ng/read/filtering.rs#L895)): it rejects reads on their
flag and mapping quality, converts the survivors into ng's own read type, then rejects some of
*those* on their length, their CIGAR and how well they match the reference. Only two of the
three are filtering.

**Two words, and the design is about the boundary between them.** After the rename in §6 the
types carry the distinction themselves:

- a **raw aligned read** (`RawAlignedRead`) is what comes off the file undecoded — a SAM flag
  and a mapping quality readable without unpacking anything else;
- an **aligned read** (`AlignedRead`) is the converted one: sequence uppercased, CIGAR in ng's
  own operations, adaptor boundary worked out
  ([`aligned_read.rs:67`](../../../../src/ng/read/aligned_read.rs#L67)).

The conversion turns one into the other, and the question this document answers is **which
filters run on which side of it, and who owns the loop.**

```
AlignedReadsReader         finds and unpacks raw aligned reads, one reused buffer
      ↓
RegionRawAlignedReads      this region's raw aligned reads only
      ↓
AlignmentCursor            ┌ filter the raw aligned read   (flag, mapping quality)
                           ├ convert it                    (raw aligned read → aligned read)
                           ├ filter the aligned read       (length, CIGAR, mismatch fraction)
                           └ keep it, and serve it to every region that overlaps it
```

**Goals.**

- Filtering stops being something a conversion happens inside. The two filters and the
  conversion each get a name, and each can be tested on its own.
- `read/filtering.rs` holds **the keep-or-drop rules and the thresholds they use, and nothing
  else**: it opens no files, converts nothing, and drives no loop over reads.
- The first filter needs **no reference**: today nothing can filter on flag and mapping quality
  without supplying a `RawRefSeq` and paying a probe fetch for every contig in the header
  ([`filtering.rs:688-702`](../../../../src/ng/read/filtering.rs#L688)).
- **Byte-identical output.** Same reads kept, same drops charged to the same reasons.

**Non-goals.**

- **No filter changes.** Which nine filters run, their thresholds, the order they run in and
  what each drop is called are `read_filtering.md`'s and stay exactly as they are.
- **No change to what the cursor keeps** (§3).
- **Speed is a constraint here, not a goal.** The change is for a simpler shape, and it is
  worth doing on that alone. A *significant* slowdown would be a reason to reconsider the
  design; an improvement is welcome and needs no defending. Neither is what this is for.
  Either way, measure it — §8.

**It does not** introduce a trait, add a type, change the meaning of any error, add a
configuration knob, or touch `AlignedRead`.

## 2. Why the division falls where it does

Step 1 runs nine filters; [`read_filtering.md`](read_filtering.md) §3 is where they are defined.
**Six of them read two values that are on the raw aligned read already** — the SAM flag and the
mapping quality. **The other three read fields that only the conversion produces.**

| filter | what it reads | which side of the conversion |
|---|---|---|
| #1 duplicate, #3 supplementary, #4 secondary, #5 unmapped, #6 QC-fail | the SAM flag | the raw aligned read |
| #2 low mapping quality | the SAM mapping quality | the raw aligned read |
| #7 too short | the sequence length | the aligned read — the conversion uppercases the sequence into a `Vec<u8>` |
| #9 bad CIGAR | the CIGAR | the aligned read — the conversion turns noodles' CIGAR into ng's `CigarOp` |
| #8 high mismatch fraction | CIGAR, sequence, base qualities, position, contig, plus a reference fetch | the aligned read, for the same two reasons |

Cited: `verdict_pre_decode` ([`filtering.rs:210-238`](../../../../src/ng/read/filtering.rs#L210)),
`verdict_post_decode` ([`filtering.rs:269-322`](../../../../src/ng/read/filtering.rs#L269)), and
the two conversions ([`aligned_read.rs:102`](../../../../src/ng/read/aligned_read.rs#L102),
[`:107`](../../../../src/ng/read/aligned_read.rs#L107)).

**Converting is not free, and that is what fixes the order.** Building an aligned read copies the
name, uppercases the sequence, rebuilds the CIGAR as ng's own operations and works out the
adaptor boundary ([`aligned_read.rs:67-120`](../../../../src/ng/read/aligned_read.rs#L67)). So:
reject on what the raw aligned read already carries, convert what survives, then reject on what
the conversion produced. **A read dropped by the first six never pays for a conversion.**

The last three cannot move earlier — they read fields the conversion makes. Taking them off the
raw read instead would mean writing the mismatch rule and the CIGAR scan a second time against
noodles' types, and two copies of one rule is the thing this module guards against hardest.

## 3. Why the filters and the conversion sit below what the cursor keeps

Once the three pieces are separable, there are three places the cursor's edge can fall. All
three produce the same reads. They differ in how much work they spend doing it, and the two we
did not choose spend it in opposite directions.

**(a) The cursor keeps raw aligned reads and the filters sit above it.** Then a read is converted
once per region that serves it. Not twice: consecutive regions overlap by about 93 % because the
caller widens each one by 5,000 bases against regions averaging 390, so a read is returned by
**about a dozen** consecutive regions
([`input/cursor.rs:12`](../../../../src/ng/read/input/cursor.rs#L12),
[`:245`](../../../../src/ng/read/input/cursor.rs#L245)). Twelve conversions per surviving read.
On CRAM it is worse again — keeping raw reads means keeping a whole decoded container, which
[`alignment_cursor.md`](alignment_cursor.md) §5 already rules out.

**(b) The cursor converts everything and the filters sit above it.** Then every raw aligned read
is converted before the flag and mapping-quality filters can speak. On a BAM with 30 % duplicates
— ordinary for real cohort input after duplicate marking — 30 % of conversions are thrown away,
and the cursor's kept set holds reads that are about to be discarded, so the memory goes too.

**(c) The filters sit below what the cursor keeps** — today's shape, and this design's. Converted
**once**, and only for reads that already cleared flag and mapping quality.

**So we choose (c): the two filters and the conversion all sit below the kept set, and the cursor
owns the loop.** Given what this change is for — a simpler shape, and no significant slowdown
(§1) — (c) is the one that costs nothing to get. (a) and (b) are not wrong; they would produce
the same output, and if the goals were different they might be the better trade.

It is written down because separating three things that were fused is exactly when someone
rearranges them, and moving any of them above the kept set brings back (a) or (b) without
anything failing.

**A rationale considered and rejected, recorded so nobody re-proposes it.** *"Splitting the
pieces lets us convert fewer reads."* It does not. The conversion already sits after the cheap
filters and before the expensive ones
([`filtering.rs:920`](../../../../src/ng/read/filtering.rs#L920), *"Decode only the pre-decode
survivors"*), and this design does not move it. Reads dropped by #1–#6 already pay no
conversion; reads dropped by #7/#8/#9 pay one and always will, because those filters read what
the conversion produces. Not one conversion is saved. The reasons to do this are §5's.

## 4. What will bite you

**The raw aligned read is one reused buffer, so a filter cannot return one.** The reader fills a
caller-owned buffer in place and the whole pass allocates exactly one
([`filtering.rs:366-382`](../../../../src/ng/read/filtering.rs#L366)). A filter shaped
`fn(RawAlignedRead) -> Option<RawAlignedRead>` would clone a noodles `RecordBuf` per read —
measured elsewhere in this codebase at ~680 bytes across seven allocations
([`aligned_reads_reader/container.rs:13-17`](../../../../src/ng/read/input/aligned_reads_reader/container.rs#L13)).
A filter takes a **borrow** and returns a **verdict**, and the verdict carries the drop reason
because the tally is keyed on it.

**A drop must be charged to a read group, and the first filter runs before any aligned read
exists.** Already solved, and easy to undo by accident: the raw aligned read carries its read
group, stamped by the region narrowing
([`region_records.rs:222`](../../../../src/ng/read/input/region_records.rs#L222)), and the
`read_group` accessor exists for the tally rather than for any filter
([`filtering.rs:350-357`](../../../../src/ng/read/filtering.rs#L350)).

**The second filter's order is #7 → #9 → #8, not the numbering order.** The cheap checks run
first so a doomed read never pays the reference fetch, and a read failing both #9 and #8 is
charged to the root cause rather than the symptom
([`filtering.rs:240-268`](../../../../src/ng/read/filtering.rs#L240)). Pinned by
`a_walk_charges_every_drop_reason_by_hand_count` (`input/cursor.rs`), whose fixture's last
record fails both.

**An unmapped read travels through `RawAlignedRead`.** Filter #5 drops it, so it exists as one
before being rejected — and an unmapped read is not aligned. `AlignedRead` has no such case: the
conversion refuses a record with no reference id or no position. The name is still right — SAM
calls every line an alignment record, unmapped ones included, and this type only ever comes from
an alignment file — but say so in the type's doc rather than leave it to be discovered and
"fixed".

**Deleting the source trait takes three fatal-error tests with it if nobody notices** — §6.

## 5. What this buys

**The three-way stop collapses to a boolean.** `FilterState::EndOfInput` exists **only** because
the filter is a separate type that cannot tell why its source stopped. Its own doc says so
([`filtering.rs:635-645`](../../../../src/ng/read/filtering.rs#L635)): a long-lived filter reaches
`EndOfInput` at the end of *every* region, because the region narrowing reports a region boundary
the only way it can, and repositioning has to undo that. The cursor is the thing that *causes*
region ends, so once it owns the loop it never has to ask. `FilterState`,
`restart_after_end_of_input`, `has_failed` and `source_mut` all go, replaced by one `failed` flag.

**No new types — this is a subtraction.** With the tally on the cursor (§7) the first filter
holds nothing, and is already a plain function today. The second needs the reference bases and
the buffer they are read into, and the cursor can hold both. So `ReadFilter` goes, the source trait goes, and what
remains is the two filters, the conversion, and a loop in the type that was already driving all
three.

**The reference stops being a precondition for filtering at all.** Nothing in ng needs that today
— it is a capability the current shape forecloses, not a request anyone has made — but a coverage
histogram, an insert-size pass or a read-group pre-pass would each want flag and mapping-quality
filtering without a reference, and cannot have it.

**Two vestiges dissolve rather than needing their own cleanup.** `with_validated_contigs` and
`ReadFilterBuffers` exist because a filter used to be built per region by a pooled caller. That
caller is gone; both stop having a reason to exist once the cursor owns the loop.

**`read/filtering.rs` ends up holding only the rules and their thresholds.** It stopped
reading files on 2026-08-03, when the two whole-file readers it owned were deleted; this takes
the conversion and the loop out too. What is left is: what counts as a read worth keeping.

## 6. Naming, and where each piece lives

**The rename is part of this design, not a tidy-up.** `RawRecord` names its topic rather than its
value — a record of *what*? A reference sequence, a read, an observation? Renaming it pairs it
with the type it converts into, so the relationship the whole design turns on becomes visible in
the names instead of needing a paragraph of prose.

| now | becomes |
|---|---|
| `RawRecord` (trait) | `RawAlignedRead` |
| `NoodlesRawRecord` | `NoodlesRawAlignedRead` |
| `RecordReader` (enum) | `AlignedReadsReader` |
| `BamRecordReader` / `CramRecordReader` / `InMemoryRecordReader` | `BamAlignedReadsReader` / `CramAlignedReadsReader` / `InMemoryAlignedReadsReader` |
| `record_reader/` (module) | `aligned_reads_reader/` |
| `RegionRecords` | `RegionRawAlignedReads` |
| `RecordSource` (trait) | **deleted** — below |

**The readers leave "raw" out of their names** (owner's call, 2026-08-03). They yield undecoded
reads, so each reader's doc comment has to say so where a caller meets it; the type name does not
carry it.

**The source trait is deleted, not renamed.** It has exactly one production implementation
(`RegionRecords`) and, once the cursor owns the loop, no generic consumer at all — it would
survive only so two test doubles could feed the filter. **The consequence is real and must not be
missed:** those doubles exist to raise the fatal errors a real file cannot conveniently raise
(`read_filter_source_read_error_is_fatal`, `read_filter_decode_error_is_fatal`,
`read_filter_reference_error_mid_stream_is_fatal`). With the trait gone the in-memory reader must
be able to yield a scripted error — which is the better arrangement anyway, because the error
path then runs through the real chain instead of through a fake that bypasses two layers.

**Module homes.**

- `read/aligned_read.rs` — `RawAlignedRead`, `NoodlesRawAlignedRead`, `AlignedRead`, and the
  conversion between them. One thing in two states, and now named so; the conversion already
  lives here.
- `read/filtering.rs` — the rules and their thresholds: `ReadFilterConfig`, `DropReason`,
  `FilterVerdict`, `ReadFilterCounts`, and the two filters themselves.
- `read/input/aligned_reads_reader/`, `read/input/region_raw_aligned_reads.rs` — unchanged in
  substance, renamed.
- `read/input/cursor.rs` — gains the loop.

## 7. The tally

**One tally per read group met, held by the cursor, cumulative until the caller says otherwise.**
A drop rate is a read group's property: one bad library shows up as an anomalous mapping-quality
or mismatch rate, and summing across read groups erases exactly that signal.

**`AlignmentCursor::reset_counts()` lets the caller choose the window.** On the cursor, not on
`AlignmentFile`, and deliberately: that file held a `Mutex<Vec<ReadGroupCounts>>` until three
commits ago, fed by each per-region stream on `Drop`, and it went with the rest of that path.
Putting it back means either a lock on the hottest loop or a fold at cursor drop — which makes
the number unreadable while the cursor is alive and loses it entirely if the cursor leaks. A
run-wide total is aggregation, not shared state: read each cursor's counts before dropping it,
which is what `SampleCursor::read_group_counts` already does across a sample's files.

**Not per region.** Regions overlap by ~93 % and a read is filtered once, when first read off the
file — never again when replayed. A per-region tally would record *where the reader happened to
be when it met a bad read*, not how bad the region is; the numbers would not sum to the
chromosome's total and would not be comparable between regions.

## 8. How we know it works

**The unit suite is necessary and not sufficient.** 45 tests live in `filtering.rs` and most
drive `ReadFilter` as one thing; each moves to whichever piece keeps its subject.

**The evidence is output identity on real data:**

| anchor | what it proves |
|---|---|
| `ng_generic_loci_dump`, `ng_ssr_loci_dump` — HG002 chr21 (BAM) | 251,792 and 4,406 lines, byte-identical |
| the same two — tomato SL4.0ch01 (CRAM) | 1,718,914 and 11,945 lines, byte-identical |
| `ng_generic_walk_probe`, chr21 | `loci=236081 observations=251786 reads_admitted=54709` |

`reads_admitted` is the direct one: step 1's keep count for the run.

**And time the walk, because the constraint in §1 needs a number to be checked against.** The
same probe prints `seconds` and `loci_per_second`; on the development machine chromosome 21 runs
in **≈1.9 s** (one run, `seconds=1.876`, `loci_per_second=125834`, 2026-08-03 — a single
measurement on one machine, not a benchmark). Compare before and after on the same machine in
the same session. A few per cent either way is noise. A large regression is a reason to
reconsider the design rather than to accept it; a large improvement is a pleasant surprise and
changes nothing about why the change was made.

**Three tests do not exist yet**, and each covers something no output comparison can see:

- **The first filter runs with no reference at all.** Untested, §5's capability quietly stops
  being true.
- **The conversion is asked for nothing when every read fails the first filter.** A *work*
  property: hoisting the conversion changes no output.
- **A scripted read error still surfaces as fatal, through the real chain.** Replaces the three
  test-double tests the deleted trait takes with it (§6).

## 9. Resolved decisions & open questions

**Q1 — resolved (owner, 2026-08-03): the cursor holds them, and both filters stay plain
functions.** The reasoning is below, kept because it is what the answer was chosen against.

**Q1. Where do the reference bases and the buffer they are read into live?**

The **first filter** — the one that looks at the SAM flag and the mapping quality — needs nothing
but the read in front of it and the thresholds. It is a plain function today, one that stands on
its own rather than belonging to an object, and it can stay one.

The **second filter** needs two more things: the reference bases, to compare the read against for
the mismatch check, and a buffer to read those bases into. The buffer is reused from read to read
so the check costs no allocation. Something has to hold both. Either the second filter becomes a
small object that holds them, or the cursor holds them and hands them over on each call — which
is what happens today, one layer up.

*Chosen: the cursor holds them, and both filters stay plain functions.* It is the smaller change,
and it leaves `read/filtering.rs` holding only rules and thresholds, with nothing to construct and
no state to get wrong.

**Q2 — resolved. The check compares two contig tables instead of opening 2,580 contigs.**

**Three different checks are involved, and it helps to keep them apart.**

1. **The reference itself is sound.** `read_reference_info`
   ([`reference_info.rs:239`](../../../../src/ng/reference_info.rs#L239)) builds the contig
   table, checks the FASTA against its `.fai`, and rejects duplicate names. Everything that opens
   a reference goes through it, and this design does not touch that.
2. **The alignment file agrees with the reference.** `AlignmentFile::open` compares the file's
   `@SQ` list against that contig table with `ContigList::first_disagreement` — names, lengths
   and, where both sides carry one, digests
   ([`open_bam.rs:206-211`](../../../../src/ng/read/input/open_bam.rs#L206)).
3. **The accessor handed to `cursor()` is over that same reference.** Nothing checks this
   directly. What stands in for it is a loop that asks the accessor for a zero-length window on
   every contig in the header ([`filtering.rs:688`](../../../../src/ng/read/filtering.rs#L688)) —
   which proves each contig *resolves* and nothing about its length.

**Check 3 is the weak one and the expensive one.** It opens every contig in turn and discards
each: ~2,580 on GRCh38 at roughly 52 µs per open with a shared index
([`ref_seq.rs:622-655`](../../../../src/ng/ref_seq.rs#L622)), an estimated **~130 ms per cursor**,
paid once per file per chromosome. (Arithmetic on a documented micro-measurement, not a
measurement of this check.) And for all that it never notices an accessor whose contigs have the
right names and the wrong lengths.

**The fix: ask the accessor which table it is over, and compare.** Every accessor already answers
that — `ContigTable::contigs()`, implemented by all three of them
([`ref_seq.rs:257`](../../../../src/ng/ref_seq.rs#L257)) — and `ContigList` is comparable, with
`first_disagreement` for the message. So the cursor's constructor does what the open gate does,
against the accessor instead of the reference:

```rust
pub fn cursor<R: RawRefSeq + ContigTable>(
    self: &Arc<Self>, contig: ContigId, reference: R,
) -> Result<AlignmentCursor<R>, AlignmentFileError> {
    self.contigs
        .first_disagreement(reference.contigs())
        .map_err(|detail| AlignmentFileError::ContigReconcile { … })?;
    …
}
```

**This is strictly better on both counts.** It costs one comparison of ~2,580 entries — integer
and string compares in memory, microseconds — against ~2,580 file opens. And it proves more:
names *and* lengths agree, which is the same thing the gate proves for the file, so the file, the
reference and the accessor are all held to one table instead of two-and-a-half.

**Why not check it once, in whatever builds the cursors?** Because at that moment the accessor
does not exist. `AlignmentFile::open` and `SampleReads::open` run before any accessor is made:
the caller supplies one per file, through a factory, and that is not incidental —
`WindowedRefSeq` holds an open per-contig reader behind a `RefCell` and is `Send` but not `Sync`,
so a sample's k files must each have their own or they share one file position and one sliding
window. Moving the factory into `open` would let the check run once, but it would not fully close
the hole — validating one accessor a factory returns does not prove the next one is over the same
FASTA — and it changes `open`'s signature and every caller. With the comparison above costing
microseconds there is nothing left to save.

**The one thing this does not do** is stop a caller passing an unrelated accessor *that happens
to carry a matching table*. Closing that means the file owning its reference and handing out
accessors itself, which is a larger change — §10.

**Impl note:** the added `+ ContigTable` bound is satisfied by every accessor in the tree
(`InMemoryRefSeq`, `ResidentRefSeq`, `WindowedRefSeq`, and the three test spies) and propagates to
`SampleReads::cursor` and to both generators. Mechanical, but it does touch their signatures.

## 10. Deferred, with a home

- **Removing `with_validated_contigs` and `ReadFilterBuffers`.** They dissolve as a consequence
  of this change (§5), so they belong to whichever step lands it — not to a separate cleanup.
- **A whole-file read path.** Filtering without a reference (§5) is not the same as a way to
  *read* a whole file, which ng does not have: `AlignmentFile::open` requires an index. If one is
  wanted it arrives as an `AlignedReadsReader` arm, beside the others.
- **The file owning its reference, and handing out accessors itself.** Q2 closes the hole that
  matters — an accessor over a different reference — by comparing contig tables. What remains is
  a caller passing an unrelated accessor whose table happens to match. Removing the possibility
  means `AlignmentFile` holding an `OpenReference` for BAM as well as CRAM (today it is `None`
  for BAM) and building accessors itself, which would also drop the factory from
  `SampleReads::cursor` and the type parameter from both generators. That is its own change, with
  its own reasons, and it is not this one's to make.
- **`read/filtering/` as a folder.** A step with no competing implementations is a file
  ([`../arch/module_layout.md`](../arch/module_layout.md) principle 1a), and this change makes
  `filtering.rs` *smaller*. Revisit only if that stops being true.

## 11. Reuse map

No new logic. Every rule exists and is being re-homed or renamed.

| what | existing code | action |
|---|---|---|
| the six flag/mapping-quality filters | `verdict_pre_decode` [`filtering.rs:210`](../../../../src/ng/read/filtering.rs#L210) | **keep**, renamed for the vocabulary |
| the three conversion-dependent filters | `verdict_post_decode` [`filtering.rs:269`](../../../../src/ng/read/filtering.rs#L269) | **keep**, renamed; the #7/#9/#8 order unchanged |
| the conversion | `RawRecord::decode` [`filtering.rs:349`](../../../../src/ng/read/filtering.rs#L349) → `decode_record` [`aligned_read.rs:67`](../../../../src/ng/read/aligned_read.rs#L67) | **reuse as-is**, called by the cursor |
| the loop | `ReadFilter::next` [`filtering.rs:895`](../../../../src/ng/read/filtering.rs#L895) | **move** into `AlignmentCursor::next_read` |
| the tally | `ReadFilterCounts` / `ReadGroupCounts` [`filtering.rs:122`](../../../../src/ng/read/filtering.rs#L122), [`:661`](../../../../src/ng/read/filtering.rs#L661) | **reuse as-is**, owned by the cursor |
| the errors | `ReadFilterError` [`filtering.rs:578`](../../../../src/ng/read/filtering.rs#L578) | **reuse as-is** — its three variants already name the three pieces |
| the raw read and its buffer contract | `RawRecord` [`filtering.rs:334`](../../../../src/ng/read/filtering.rs#L334) | **rename and move** to `aligned_read.rs` |
| the region narrowing | `RegionRecords` [`region_records.rs`](../../../../src/ng/read/input/region_records.rs) | **rename**; its `RecordSource` impl becomes inherent methods |
| the three-way stop | `FilterState` [`filtering.rs:646`](../../../../src/ng/read/filtering.rs#L646) | **delete** — a `failed` flag on the cursor replaces it (§5) |
| the source trait and its doubles | `RecordSource` [`filtering.rs:366`](../../../../src/ng/read/filtering.rs#L366), `FakeSource`, `ErroringSource` | **delete**; the in-memory reader gains a scripted error (§6) |

**The parity oracle is the four dumps** (§8). There is no second implementation to differ from —
the oracle is this code before the change.
