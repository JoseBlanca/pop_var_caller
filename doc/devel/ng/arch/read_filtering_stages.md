# ng — read filtering in stages: types & interfaces

*Architecture draft, 2026-08-03. Code-facing companion to
[`../spec/read_filtering_stages.md`](../spec/read_filtering_stages.md) — every **why** points
there and is not re-argued here. Revises the single-type shape in
[`read_filtering.md`](read_filtering.md) §1; the filter policy that doc pins down is unchanged.
Under [`module_layout.md`](module_layout.md). Naming per
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md). Signatures are
illustrative; the **contract** is the deliverable.*

**Two words, used exactly as the spec defines them** (spec §1): a **raw record** is the
undecoded reused buffer (`RawRecord` / `NoodlesRawRecord`), a **read** is ng's `AlignedRead`
after the conversion. The two stage types are named for which of the two they judge.

**Four of the shapes below are not settled** — spec §8 holds them as questions for the owner.
They are marked `OPEN:` where they bite. Do not build past one without an answer.

---

## Module home

`src/ng/read/filtering.rs`, unchanged — spec §6. Nothing moves file.

## 1. The types

### 1.1 What a stage answers

Both stages answer the same question about a candidate — keep it, or drop it for a named
reason — and the answer type already exists.

```rust
// Unchanged, in place: filtering.rs:96-118.
pub enum FilterVerdict { Keep, Drop(DropReason) }
pub enum DropReason { Duplicate, LowMapq, Supplementary, Secondary, Unmapped,
                      QcFail, TooShort, HighMismatchFraction, BadCigar }
```

**No new verdict type, and no `Option<T>` in a stage's return.** `Option` would lose the reason,
and the reason is what the tally is keyed on.

### 1.2 The first stage — admission on the raw record

Filters #1–#6, on the raw record's SAM flag and mapping quality. It needs no reference and no contig
probe; that is the capability the split exists to open (spec §4).

```rust
/// Decides whether an undecoded record is worth decoding, on its flag and its
/// mapping quality alone.
///
/// **Takes a borrow and returns a verdict, never a record.** The record lives in a
/// single buffer the source refills in place, so returning one would mean cloning it
/// (spec §3).
pub struct RecordAdmission {
    config: ReadFilterConfig,
    /// One tally per read group met, in first-seen order.
    /// OPEN: spec §8 Q2 — this may belong to the driver instead.
    counts: Vec<ReadGroupCounts>,
}

impl RecordAdmission {
    pub fn new(config: ReadFilterConfig) -> Self;

    /// Judge one record and charge the drop, if any, to the read group the record
    /// carries. Infallible: nothing here reads a file or the reference.
    pub fn admit(&mut self, record: &impl RawRecord) -> FilterVerdict;

    pub fn counts(&self) -> &[ReadGroupCounts];
}
```

**Contract.** Infallible. Reads `flag()`, `mapq()` and `read_group()` off the record and nothing
else. Charges every drop to a reason and a read group; a record whose source never stamped a
group is charged to `None`, never to an arbitrary one. Holds no reference, no buffer and no file.

`OPEN:` spec §8 Q1 — whether this is a type at all, or whether `verdict_pre_decode` stays a free
function and the driver keeps the tally. The block above assumes a type; if Q1 lands the other
way, delete it and keep `verdict_pre_decode(flag, mapq, &config)` as it is today.

### 1.3 The converter

Already exists and already stands alone: `RawRecord::decode`
([`filtering.rs:349`](../../../../src/ng/read/filtering.rs#L349)) calling `decode_record`
([`aligned_read.rs:67`](../../../../src/ng/read/aligned_read.rs#L67)). **No new type.** What
changes is that it is called by the driver between two stages rather than from inside one.

Its failure stays fatal to the run, not a drop — a record that cleared stage one is mapped, so a
decode failure means a corrupt record (spec, `read_filtering.md` §7).

### 1.4 The second stage — admission on the read

Filters #7, #9 and #8, on the decoded read, **in that order** — the cheap checks first so a doomed read never pays the
reference fetch, and a read failing both #9 and #8 charged to the root cause (spec §3).

```rust
/// Decides whether a decoded read survives, on its length, its CIGAR, and how well
/// it matches the reference.
///
/// The only stage that touches a reference, and the only one that can fail.
pub struct ReadAdmission<R: RawRefSeq> {
    reference: R,
    config: ReadFilterConfig,
    /// Reused scratch for the mismatch check's reference fetch; touched only when
    /// that filter runs.
    ref_buf: Vec<u8>,
    /// OPEN: spec §8 Q2, as above.
    counts: Vec<ReadGroupCounts>,
}

impl<R: RawRefSeq> ReadAdmission<R> {
    /// Fallible, because it probes every contig of `header` against the reference so
    /// an in-loop fetch failure means corrupt input rather than a mismatched
    /// reference. This is `ReadFilter::new`'s probe, moved.
    pub fn new(header: &sam::Header, reference: R, config: ReadFilterConfig)
        -> Result<Self, RefSeqError>;

    /// `Err` is a reference-fetch failure and is **fatal to the run**, never a drop.
    pub fn admit(&mut self, read: &AlignedRead) -> Result<FilterVerdict, RefSeqError>;

    pub fn reference(&self) -> &R;
    pub fn counts(&self) -> &[ReadGroupCounts];
}
```

**Contract.** The probe runs once, at construction, over the header's `@SQ` list. `admit` reads
only the decoded read and the reference. `reference()` exists so a long-lived caller can release
the bases a walk has gone past — the cursor uses it today
([`input/cursor.rs:553`](../../../../src/ng/read/input/cursor.rs#L553)).

**`with_validated_contigs` and `ReadFilterBuffers` are deleted here**, not in a separate cleanup
— spec §9.

### 1.5 The driver

Someone has to pull records until one survives, own the three-way stop, and turn three failure
kinds into one error. That is this, and it keeps the name callers already use.

```rust
/// One sample's reads, filtered — the composition of the two stages and the
/// converter, as an `Iterator<Item = Result<AlignedRead, ReadFilterError>>`.
pub struct ReadFilter<S: RecordSource, R: RawRefSeq> {
    source: S,
    /// The single record buffer reused across the whole pass.
    record_buf: S::Record,
    records: RecordAdmission,
    reads: ReadAdmission<R>,
    state: FilterState,
}
```

**Contract, unchanged from today.** Lazy; one record resident. Fused — a fatal error is yielded
once and then `None`. Only `EndOfInput` is undone by `restart_after_end_of_input`; `Failed` is
permanent ([`filtering.rs:646-658`](../../../../src/ng/read/filtering.rs#L646)).

The surface the alignment cursor drives is **unchanged** — `next`, `has_failed`,
`restart_after_end_of_input`, `source_mut`, `counts`, `reference`
([`input/cursor.rs:241`](../../../../src/ng/read/input/cursor.rs#L241) and its uses). That is
what keeps this change inside one file.

`OPEN:` spec §8 Q3 (is `ReadFilter` still the right name for the composite, now that two of its
three jobs have names of their own) and Q4 (does the cursor keep consuming a driver, or compose
the stages itself — the block above assumes a driver).

## 2. Errors

**No new error type, and no change of meaning.** `ReadFilterError`'s three variants already name
the three stages ([`filtering.rs:578`](../../../../src/ng/read/filtering.rs#L578)):

| variant | raised by |
|---|---|
| `Source` | the record source, under stage one |
| `Decode` | the converter |
| `Reference` | stage two's mismatch check |

`RecordAdmission::admit` cannot fail at all, which is why it returns a bare verdict.

## 3. Design decisions — decided

- **Two filter stages, not one** — spec §2.
- **A stage takes a borrow and returns a verdict, never a record** — spec §3.
- **The drop reason travels in the verdict**, so an `Option` return is ruled out: the tally is
  keyed on reason and read group.
- **The read group is read off the raw record, not the read** — it is already stamped there for
  exactly this ([`filtering.rs:350-357`](../../../../src/ng/read/filtering.rs#L350)).
- **Stage two owns the reference and the contig probe; stage one has neither** — spec §4.
- **Stage two keeps the #7 → #9 → #8 order**, not the numbering order — spec §3, pinned by
  `a_walk_charges_every_drop_reason_by_hand_count`.
- **The driver keeps the three-way `FilterState`** — spec §3.
- **No trait, no bake-off:** no competing implementations, so concrete types in one file —
  `module_layout.md` principle 1a.

## 4. Reconciliation with existing code

Every row read at the cited line, 2026-08-03. This is a re-homing: no row is new logic.

| what | existing code | action |
|---|---|---|
| stage one's body | `verdict_pre_decode` [`filtering.rs:210`](../../../../src/ng/read/filtering.rs#L210) | **move** into `RecordAdmission::admit`, body unchanged |
| stage two's body | `verdict_post_decode` [`filtering.rs:269`](../../../../src/ng/read/filtering.rs#L269) | **move** into `ReadAdmission::admit`, body unchanged |
| the converter | `RawRecord::decode` [`filtering.rs:349`](../../../../src/ng/read/filtering.rs#L349) → `decode_record` [`aligned_read.rs:67`](../../../../src/ng/read/aligned_read.rs#L67) | **reuse as-is**; called by the driver |
| the contig probe | `ReadFilter::new` [`filtering.rs:688`](../../../../src/ng/read/filtering.rs#L688) | **move** to `ReadAdmission::new` |
| the raw record and its buffer contract | `RawRecord` / `RecordSource` [`filtering.rs:334`](../../../../src/ng/read/filtering.rs#L334), [`:366`](../../../../src/ng/read/filtering.rs#L366) | **reuse as-is**, unchanged |
| the verdict and the reasons | `FilterVerdict` / `DropReason` [`filtering.rs:96`](../../../../src/ng/read/filtering.rs#L96), [`:104`](../../../../src/ng/read/filtering.rs#L104) | **reuse as-is** |
| the tally | `ReadFilterCounts` / `ReadGroupCounts` [`filtering.rs:122`](../../../../src/ng/read/filtering.rs#L122), [`:661`](../../../../src/ng/read/filtering.rs#L661); `tally_for_current_record` [`:846`](../../../../src/ng/read/filtering.rs#L846) | **reuse as-is**; `OPEN:` Q2 decides the owner |
| the three-way stop | `FilterState` [`filtering.rs:646`](../../../../src/ng/read/filtering.rs#L646) | **reuse as-is**, on the driver |
| the errors | `ReadFilterError` [`filtering.rs:578`](../../../../src/ng/read/filtering.rs#L578) | **reuse as-is** |
| the probe-free constructor + lent buffers | `with_validated_contigs` [`filtering.rs:746`](../../../../src/ng/read/filtering.rs#L746), `ReadFilterBuffers` | **delete** — their only caller was the deleted reader pool |
| the caller | `AlignmentCursor` holds `ReadFilter<RegionRecords, R>` [`input/cursor.rs:241`](../../../../src/ng/read/input/cursor.rs#L241) | **unchanged**, if Q4 lands on a driver |

## 5. Open items

**Genuine open design questions** — spec §8 holds the reasoning and a leaning for each:

- `OPEN: Q1` — is stage one a type, or does `verdict_pre_decode` stay a free function?
- `OPEN: Q2` — one tally on the driver, or one per stage merged on read?
- `OPEN: Q3` — what the composite is called if `ReadFilter` names stage two instead.
- `OPEN: Q4` — driver, or the cursor composing the three stages?

**Impl-time confirmations, not decisions:**

- Whether `RecordAdmission::admit` takes `&impl RawRecord` or a concrete `&S::Record`. The
  generic form is stated above; the concrete one may inline better and changes no contract.
- Whether `ref_buf` stays a field of stage two or moves to the driver. It is scratch for one
  filter; keeping it beside that filter is the assumption.
- Whether `filtering.rs` needs to become `read/filtering/` once three types live in it. Decide on
  the file's real size (spec §6).

## 6. Test & bench shape

Tests stay beside the code, in `filtering.rs`'s own module. 45 tests live there now and most
drive `ReadFilter` as one thing; each moves to whichever type keeps its subject.

**Two tests do not exist yet, and spec §7 says why they are the point of the change:** stage one
building with no reference, and the converter being asked for nothing when every record fails
stage one.

**The regression anchors are output identity on real data, not the unit suite** — spec §7 has the
four dumps and the walk probe's figures.
