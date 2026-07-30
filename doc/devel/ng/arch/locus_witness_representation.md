# ng — what a locus observation says a read witnessed: types & interfaces

*Architecture draft (2026-07-30), the code-facing companion to
[`../spec/locus_witness_representation.md`](../spec/locus_witness_representation.md). It gives the
types and signatures as they appear in code; **every *why* points back to the spec**, which is the
authority and is not re-argued here. It **amends** [`locus_generation.md`](locus_generation.md) —
the arch doc that owns `SampleLocusObservations`, `SequenceObservation` and `ReadWitness` — rather
than restating it: §1 below names the three blocks of that doc this change supersedes, and
everything else there stands. The generic generator's own interfaces are in
[`locus_generation_pileup.md`](locus_generation_pileup.md). Signatures are illustrative; the
**contract** is the deliverable.*

**No code yet.** The spec settles the design; this settles the shapes; the implementation plan
orders the work.

## Module home

A new file, `src/ng/locus_generation/witness.rs`, holding **every type that answers "what did this
read see"** — `ReadWitness`, `LocusLen`, `WitnessedLocusPositions`, `UnwitnessedBases` — and
re-exported from `locus_generation` so no consumer's import path changes.

Why a file rather than more of `mod.rs`: the canonical-representation invariant (spec §3.3) is the
one thing here that fails silently, it spans four types, and private fields need one module
boundary to be private *to*. In a 1,700-line `mod.rs` it has no home.

```
src/ng/locus_generation/
  mod.rs      – SampleLocusObservations, SequenceObservation, LocusGenerator, dispatcher, NoLoci
  witness.rs  – ReadWitness, LocusLen, WitnessedLocusPositions, UnwitnessedBases   ← new
  ssr.rs      – the STR generator
  pileup/     – the generic generator
```

---

## 1. The types

### 1.1 What a read witnessed, on the locus axis

A read's witness is a **set of locus positions**, not one run (spec §3.1). It is held as canonical
half-open runs — sorted, non-empty, non-adjacent, non-overlapping — because that is the only form
in which two equal sets compare equal, and observation identity depends on that comparison
(spec §3.3).

```rust
/// The locus positions one read witnessed — **a set, in locus coordinates**, and never
/// empty (a read that witnessed nothing inside a footprint does not fold into it).
///
/// *Invariant, enforced at construction and the reason the field is private:* the runs
/// are sorted, half-open, non-empty, and separated by at least one gap position. One set,
/// one representation — so `Eq`/`Hash` may be derived and two reads with the same witness
/// share one observation (spec §3.3).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct WitnessedLocusPositions(/* private; encoding OPEN, §4 */);

impl WitnessedLocusPositions {
    /// From runs in any order, normalised: sorted, then adjacent and overlapping runs
    /// merged. `None` if `runs` is empty or any run is.
    pub fn new(runs: impl IntoIterator<Item = (u16, u16)>) -> Option<Self>;

    /// The one-run case, which is what the STR path mints and the common generic one.
    pub fn one_run(offset_in_locus: u16, positions_covered: u16) -> Option<Self>;

    /// The runs, canonical order, half-open `[start, end)` in locus coordinates.
    pub fn runs(&self) -> impl ExactSizeIterator<Item = (u16, u16)> + '_;

    /// How many locus positions in total — what `num_obs_along_locus` sums over.
    pub fn positions_covered(&self) -> u32;

    /// Prefix / suffix constraints, unchanged in meaning: the first run starts at 0, the
    /// last run ends at the locus length (spec §3.1).
    pub fn is_flush_left(&self) -> bool;
    pub fn is_flush_right(&self, locus_len: LocusLen) -> bool;
}
```

**`ReadWitness` keeps both variants and its constructors' signatures.** `Complete` stays so
`complete_observations` remains an equality test and the STR path's call sites do not move
(spec §1 goal 4, §3.1).

```rust
pub enum ReadWitness {
    /// The read reached both borders and witnessed every position between them.
    Complete,
    /// The positions it did witness — one run, or several.
    Observed { positions: WitnessedLocusPositions },
}

impl ReadWitness {
    /// Unchanged signatures, now building a one-run set. Still clamp-then-derive.
    pub fn from_left(positions_covered: u16, locus_len: LocusLen) -> Self;
    pub fn from_right(positions_covered: u16, locus_len: LocusLen) -> Self;
    /// **New** — the interior run neither constructor could express, which is what the
    /// deferred note on the variant asked for (spec §1).
    pub fn from_run(offset_in_locus: u16, positions_covered: u16, locus_len: LocusLen) -> Self;

    pub fn is_flush_left(&self) -> bool;
    pub fn is_flush_right(&self, locus_len: LocusLen) -> bool;
}
```

*Contract.* All three constructors clamp into `locus_len` and return `Complete` when the result
covers the whole locus — so a run reaching both borders can never masquerade as partial. The
clamp remains a convention rather than a type invariant, for the reason already recorded on the
variant: a run clamped against *some* `LocusLen` proves nothing about the locus it is finally
attached to, so the real check lives in `num_obs_along_locus`, where the region is in hand.

### 1.2 Which bases the read did not supply, on the read axis

When an indel's anchor position has no `Match` event the fold takes that base from the reference so
the indel has an anchor, and `bases` cannot currently say which base that was. It fires **once per
indel event**, not once per observation, so this is a set and not a flag (spec §2).

```rust
/// Which entries of `bases` came from the reference rather than from the read —
/// **indices into `bases`, the read axis**, not locus positions.
///
/// Empty is the common case and costs no allocation. *Invariant:* strictly ascending, so
/// one set has one representation and `Eq`/`Hash` are derivable (spec §3.3).
#[derive(Debug, Clone, Default, PartialEq, Eq, Hash)]
pub struct UnwitnessedBases(/* private; encoding OPEN, §4 */);

impl UnwitnessedBases {
    /// From indices in any order; deduplicated and sorted.
    pub fn new(indices: impl IntoIterator<Item = u16>) -> Self;
    /// The common case, and what `Default` gives.
    pub fn none() -> Self;
    pub fn is_empty(&self) -> bool;
    pub fn contains(&self, index_into_bases: u16) -> bool;
    pub fn indices(&self) -> impl ExactSizeIterator<Item = u16> + '_;
}
```

**The field on the observation is not an `Option`.**

```rust
pub struct SequenceObservation {
    pub bases: Box<[u8]>,
    /// Which of `bases` the read did not supply. Empty on an observation the read
    /// witnessed whole — the common case (spec §3.2).
    pub unwitnessed_bases: UnwitnessedBases,
    // … read_witness (renamed, §3), then unchanged: read_group, num_obs, num_fwd,
    //   q_sum, mapq_sum, mapq_sum_sq, placed_left, chain_ids
}
```

*Why the two live on different axes and are not one type:* `read_witness` is in locus positions,
`bases` is allele content in read coordinates. An insertion adds bases without positions and a
deletion positions without bases, so one index cannot address both (spec §3.2).

---

## 2. What the fold hands over

The fold works in **reference** coordinates and `finalise` resolves into **locus** coordinates
against the record's final footprint. Those are two axes, and the type system separates them for
the same reason `LocusLen` exists: two same-shaped quantities that can be transposed get two types
so the compiler refuses the mix.

```rust
/// The reference positions one read witnessed inside a record, canonical runs, never
/// empty. The fold-time counterpart of `WitnessedLocusPositions`; `witness_of` is the
/// only conversion. Replaces `RefSpan`, which could hold one run only.
pub(super) struct WitnessedRefPositions(/* private */);

/// The read's witness, or `None` only when it witnessed nothing inside the record at all.
/// **A hole no longer yields `None`** — that is the second gap this change closes
/// (spec §1, §4).
pub(super) fn apply_events_into(
    allele_seq: &mut Vec<u8>,
    record_pos: u32,
    ref_seq: &[u8],
    events: &[ReadEvent],
    unwitnessed: &mut UnwitnessedBases,   // the two borrow sites record into this
) -> Option<WitnessedRefPositions>;

/// Resolve a witness against the **final** footprint — once, at `finalise`, never during
/// the fold. Clamps into the footprint at both ends before measuring.
pub(super) fn witness_of(
    witnessed: &WitnessedRefPositions,
    record_pos: u32,
    record_end_exclusive: u32,
) -> ReadWitness;

/// The sort key `finalise` orders observations by. **No longer `(u8, u16, u16)`** — it
/// has to be a total order over sets, so it borrows. `ReadWitness` still gets no `Ord`
/// impl of its own: that would export this file's sorting convention to every consumer.
pub(super) fn witness_order(witness: &ReadWitness) -> impl Ord + '_;
```

*Contract.* `apply_events_into` returns `None` only for a read with no witnessed position inside
the record; every other read now yields a witness, and the drop path narrows to that one case
(spec §4). `unwitnessed` is cleared by the callee, like `allele_seq`, so the caller reuses one
buffer per fold rather than allocating per read.

**Errors: none new, and the out-of-range case is already unreachable.** A footprint wider than a
`u16` run can describe is rejected at construction by
`PileupGeneratorConfig::check` ([generator.rs:123-131](../../../../src/ng/locus_generation/pileup/generator.rs#L123)),
which returns `PileupGeneratorConfigError::RecordSpanExceedsCoverageRun` when `max_record_span`
exceeds `MAX_RECORD_SPAN_CEILING` ([:43](../../../../src/ng/locus_generation/pileup/generator.rs#L43)).
So the spec's "must fail loudly rather than clamp" (§3.4, §5) is met by a gate that already exists
and a `debug_assert` restating the envelope in `witness_of` — not by a new error variant. The
step still mints no error of its own ([`locus_generation.md`](locus_generation.md) §4).

---

## 3. Decisions — decided (why in the spec)

- **`ReadCoverage` is renamed `ReadWitness` — decided (owner, 2026-07-30).** "Coverage" reads as
  *depth* to a geneticist, which is why the type's own doc has to correct it — "one read's span,
  not depth" — and after this change the payload is a witness outright. The names that follow are
  mechanical: the field `read_coverage` → `read_witness`, and `coverage_of` / `coverage_order` →
  `witness_of` / `witness_order`. **193 uses of the type, 91 of the field, across 12 files**, none
  of it behaviour. Not raised by the spec; it is this doc's call because it is a name.
- **`ReadWitness::Observed` carries a set, not a run** — one shape for one run, two runs, N runs;
  the flush predicates and `positions_covered` become derivations from it. *Rejected:* one
  observation per run (breaks `num_obs` as a read count and needs a read identity `chain_ids`
  cannot supply), and a bare `SmallVec` of runs on the variant (leaves `bases` a concatenation the
  read never showed as one sequence) — spec §3.1.
- **`Complete` stays a variant** — keeps `complete_observations` an equality test and the STR
  path's call sites unmoved (spec §3.1). *Note:* the spec supports this with an observation count
  that is not recorded in the Milestone D report; the decision stands on the call-site and
  equality-test grounds, which are independent of it (§4).
- **Canonicality is enforced at construction, so `Eq`/`Hash` are derived** — the spec asks for
  hand-written impls "not derived over a representation that admits duplicates" (§3.3); a private
  field with a normalising constructor removes the ambiguity at the source, after which derived
  impls are correct *and* cheaper than hand-written ones. Same end, one fewer thing to get wrong.
- **`witness_order` borrows instead of returning `(u8, u16, u16)`** — a total order over a set
  cannot be a fixed-width tuple. `ReadWitness` still gets no `Ord`, unchanged reasoning
  ([open_record.rs:255-259](../../../../src/ng/locus_generation/pileup/open_record.rs#L255)).
- **`BaseProvenance` is named `UnwitnessedBases`** — "provenance" names the topic; the value is the
  set of bases the read did not witness, which is what a consumer decides on. It also keeps one
  term for one concept: positions witnessed on the locus axis, bases unwitnessed on the read axis.
  *Live alternative:* `BorrowedBases`, the spec's own verb (§2), equally accurate; it lost on the
  term-family point alone.
- **The field is `UnwitnessedBases`, not `Option<UnwitnessedBases>`.** The spec's `Option` buys "8
  bytes and no allocation" (§3.2) only for a heap-backed type; inline, an empty set already costs
  no allocation and the `Option` adds a discriminant and a second spelling of "nothing". No tension
  with `Option`-means-absent: an empty set is a set, not a sentinel. **Reopen if the encoding turns
  out heap-backed** (§4).
- **Two witness types, one per coordinate axis** — public for the locus, `pub(super)` for the
  fold's reference axis. The reasoning that minted `LocusLen`
  ([mod.rs:241-250](../../../../src/ng/locus_generation/mod.rs#L241)): same-shaped quantities that
  can be transposed get two types so the compiler refuses the mix.

*Boy-scout note.* `RecordWitness` ([open_record.rs:121](../../../../src/ng/locus_generation/pileup/open_record.rs#L121))
is a tally of witness outcomes, not a witness, and after this change there are two real `*Witness`
types to confuse it with. `RecordWitnessCounts` when the file is next touched — `pub(super)`, so
the rename is local.

---

## 4. Open items

- `OPEN:` **The encoding of both sets, and the inline bound.** Interface above is written so the
  encoding is private and swappable: canonical runs behind `runs()`, ascending indices behind
  `indices()`. Candidates are a `SmallVec` of runs (leaning — the common witness is one run and the
  STR path never mints more) against a bitmask (`u128` covers 128 positions but pays 16 bytes to
  say "positions 0..40"). **Settled by the distribution of runs per witness and of unwitnessed
  bases per observation on real data — neither measured** (spec §8). Do not pin the field types
  before those two counts exist.
- *Impl-time confirmations, not design items:* whether `runs()` returns `(u16, u16)` pairs or a
  `Run` newtype (a newtype if any call site ever holds one loose); whether `UnwitnessedBases`
  indexes `bases` as `u16` or `u32` (bounded by `bases.len()`, which an insertion can grow past a
  footprint — check against the widest observed allele before pinning `u16`).

---

## 5. Reconciliation with existing code

Every row read at the cited line (2026-07-30), in the `ng-generic` worktree. Convergence on
existing types, not new types beside them.

| this doc's name | existing code | action |
|---|---|---|
| `ReadWitness` | `ReadCoverage` [mod.rs:213](../../../../src/ng/locus_generation/mod.rs#L213) | **rename** — 193 uses of the type and 91 of the field `read_coverage` across 12 files (§3) |
| `WitnessedRefPositions` | `RefSpan` [open_record.rs:65](../../../../src/ng/locus_generation/pileup/open_record.rs#L65) | **replace** — same role (the fold's witnessed extent, reference coordinates, never empty), widened from one run to a set. Its "no empty span" invariant carries over verbatim |
| `WitnessedLocusPositions` | `ReadCoverage::Observed { offset_in_locus, positions_covered }` [mod.rs:216-224](../../../../src/ng/locus_generation/mod.rs#L216) | **replace the payload** — the two `u16`s are the one-run case of the set |
| `is_flush_left` / `is_flush_right` | [mod.rs:327-348](../../../../src/ng/locus_generation/mod.rs#L327) | keep signatures; body derives from the first and last run |
| `num_obs_along_locus`'s clamp | [mod.rs:69-104](../../../../src/ng/locus_generation/mod.rs#L69) | **keep** — the comment there explains why the bound is not expressible on the type, and a set does not change that; the loop iterates runs instead of one range |
| `complete_observations` | [mod.rs:123](../../../../src/ng/locus_generation/mod.rs#L123) | unchanged — `Complete` is still a variant, still an equality test |
| the two borrow sites | Insertion [open_record.rs:1217-1224](../../../../src/ng/locus_generation/pileup/open_record.rs#L1217), Deletion [:1228-1238](../../../../src/ng/locus_generation/pileup/open_record.rs#L1228) | each records the index it pushes into `UnwitnessedBases`; the `offset >= consumed_until` guard is already exactly "no `Match` emitted this position" |
| the drop path | [open_record.rs:1335-1345](../../../../src/ng/locus_generation/pileup/open_record.rs#L1335) | **narrows** to "witnessed nothing"; the set-of-read-ids mechanism and the subtract-prior-contribution step survive untouched — a read still becomes non-contiguous when the window widens |
| `witness_of` | `coverage_of` [open_record.rs:184-212](../../../../src/ng/locus_generation/pileup/open_record.rs#L184) | renamed with its type; resolves a set against the footprint. The both-ends clamp and the `Complete` short-circuit are kept |
| `witness_order` | `coverage_order` [open_record.rs:261](../../../../src/ng/locus_generation/pileup/open_record.rs#L261) | renamed with its type; signature changes to borrow, and `pub(super)` for the differential stays |
| the observation identity | `ObservationKey` [open_record.rs:228](../../../../src/ng/locus_generation/pileup/open_record.rs#L228) | same three fields; the witness one grows a set. Renamed with its sibling — spec §4 |
| the out-of-range gate | `PileupGeneratorConfig::check` [generator.rs:123-131](../../../../src/ng/locus_generation/pileup/generator.rs#L123), `MAX_RECORD_SPAN_CEILING` [:43](../../../../src/ng/locus_generation/pileup/generator.rs#L43) | **reuse as-is** — already rejects a footprint wider than a run can describe, so §2 needs no new error |
| the STR generator's mints | [ssr.rs:770](../../../../src/ng/locus_generation/ssr.rs#L770), [:821-822](../../../../src/ng/locus_generation/ssr.rs#L821), [:889](../../../../src/ng/locus_generation/ssr.rs#L889), [:989](../../../../src/ng/locus_generation/ssr.rs#L989) | **call sites unchanged** — `from_left`/`from_right` keep their signatures, which is why they keep them |

---

## 6. Test & bench shape

**Where they live.** The invariant tests in `witness.rs`, the fold's in `pileup/open_record.rs` and
`pileup/tests.rs`, the differential and the census in `pileup/parity.rs`. **No new bench** — no
bake-off here, and the cost question is an allocation count, not a wall time.

**Spec §7 lists the five acceptance checks; only one of them is this doc's to shape.** The
canonicality property is the failure that is silent, so it needs a property test rather than a
fixture: two `WitnessedLocusPositions` built from different orderings of the same runs compare
equal and hash equal, for any generated multiset. Its mutation is the real check — remove the
normalising step in the constructor and observation counts must inflate. Nothing else in the change
fails quietly; the rest are fixtures the spec already enumerates.
