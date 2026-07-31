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
read see"** — `ReadWitness`, `LocusLen`, `WitnessedLocusPositions` — and
re-exported from `locus_generation` so no consumer's import path changes.

Why a file rather than more of `mod.rs`: the canonical-representation invariant (spec §3.3) is the
one thing here that fails silently, it spans four types, and private fields need one module
boundary to be private *to*. In a 1,700-line `mod.rs` it has no home.

```
src/ng/locus_generation/
  mod.rs      – SampleLocusObservations, SequenceObservation, LocusGenerator, dispatcher, NoLoci
  witness.rs  – ReadWitness, LocusLen, WitnessedLocusPositions   ← new
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
    Partial { positions: WitnessedLocusPositions },
}

impl ReadWitness {
    /// Unchanged signatures, now building a one-run set. Still clamp-then-derive.
    pub fn from_left(positions_covered: u16, locus_len: LocusLen) -> Option<Self>;
    pub fn from_right(positions_covered: u16, locus_len: LocusLen) -> Option<Self>;
    /// **New** — the positions witnessed, which subsumes the interior run neither of the
    /// two could express (spec §1), and the only constructor that may answer `Complete`.
    pub fn from_witnessed_runs(
        runs: impl IntoIterator<Item = (u16, u16)>,
        locus_len: LocusLen,
    ) -> Option<Self>;

    pub fn is_flush_left(&self) -> bool;
    pub fn is_flush_right(&self, locus_len: LocusLen) -> bool;
}
```

*Contract — **revised at D3 (owner, 2026-07-31)**. This paragraph asked for a third constructor
`from_run(offset, covered, locus_len)` and for **all** constructors to return `Complete` when the
clamped run covers the whole locus. Both were replaced, by a split on **what the caller claims**:*

- ***A reach*** *—* `from_left` */* `from_right`*: "the read got at least this far from this
  border". A lower bound, and on the STR path it is counted in* **read** *bases against a locus
  measured in* **reference** *positions (`ssr.rs`'s `reach` against `locus.segment.tract_len()`),
  two rulers that diverge under stutter. A reach at or past the locus length therefore says the
  read* ran out of read*, not that it reached the far border, so these never answer `Complete` —
  which matters because `Complete` gates `complete_observations`, i.e. what a likelihood may score
  as an* exact *length. Implementing the original contract would have scored a lower bound as a
  measurement, and moved two columns of the STR byte-identity oracle while doing it.*
- ***A witnessed set*** *—* `from_witnessed_runs`*: "these are the positions the read witnessed",
  on the locus's own ruler. Completeness is then arithmetic rather than inference, decided on the*
  total *positions covered and never on the outer edges (a set flush at both borders can still
  have a hole). The precondition is that ruler and it is the caller's to keep; the shape is what
  protects it, since a producer holding a reach and a border cannot build locus-coordinate runs by
  accident.*
- `from_run` *is* not *built: an interior run is* `from_witnessed_runs([(3, 7)], len)`*, and two
  spellings of one run differing only in whether they decide completeness is a coin-flip for the
  caller.*
- `Complete` *stays a* bare *variant a caller writes when it knows structurally — the STR
  delimiter reporting both borders of the tract anchored in this read. It gains no payload: a
  stored span would be a claim about a locus the type cannot see (below), and 1,646,289 of
  1,647,161 observations on the chr1 run are `Complete`, so every one of them would build, compare
  and hash a run that says nothing new.*

The clamp remains a convention rather than a type invariant, for the reason already recorded on the
variant: a run clamped against *some* `LocusLen` proves nothing about the locus it is finally
attached to, so the real check lives in `num_obs_along_locus`, where the region is in hand.

### 1.2 `SequenceObservation` gains no field

Its only change is the rename of `read_coverage` to `read_witness` (§3). A borrowed anchor base
stays indistinguishable from a sequenced one: the spec deferred that half on the measurement — **8
borrowed bases in 225 million event-folds, never two in one observation** (spec §8) — and the
design for it, when a measurement asks, is recorded in spec §6, not here. An arch doc that typed a
field nobody is building would be the drift this one exists to prevent.

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
/// **A hole no longer yields `None`** — that is the gap this change closes (spec §1, §3.1).
/// The signature is otherwise today's: the two borrow sites are untouched, since that half
/// is deferred (§1.2).
pub(super) fn apply_events_into(
    allele_seq: &mut Vec<u8>,
    record_pos: u32,
    ref_seq: &[u8],
    events: &[ReadEvent],
) -> Option<WitnessedRefPositions>;

/// Resolve a witness against the **final** footprint — once, at `finalise`, never during
/// the fold. Clamps into the footprint at both ends before measuring, rebases the runs from
/// reference onto locus coordinates, and hands them to `ReadWitness::from_witnessed_runs`,
/// which owns the `Complete` decision (D3): this is the caller that states positions rather
/// than a reach, and the footprint's width is the emitted locus's length exactly.
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
(spec §4). The runs are accumulated into a buffer the caller owns and the callee clears, like
`allele_seq`, so a fold allocates nothing per read.

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
- **`ReadWitness::Partial` carries a set, not a run** — one shape for one run, two runs, N runs;
  the flush predicates and `positions_covered` become derivations from it. *Rejected:* one
  observation per run (breaks `num_obs` as a read count and needs a read identity `chain_ids`
  cannot supply), and a bare `SmallVec` of runs on the variant (leaves `bases` a concatenation the
  read never showed as one sequence) — spec §3.1.
- **The variant is `Partial`, not `Observed`** — "observed" is not a contrast with `Complete`,
  since a complete witness was observed too, and once the enum says *witness* the word carries
  nothing. `PartialLeft` / `PartialRight` were removed for being side-tagged; this payload is a set
  of positions with no side, so the word is free again (spec §3.1).
- **`Complete` stays a variant** — keeps `complete_observations` an equality test and the STR
  path's call sites unmoved (spec §3.1).
- **Canonicality is enforced at construction, so `Eq`/`Hash` are derived** — the spec asks for
  hand-written impls "not derived over a representation that admits duplicates" (§3.3); a private
  field with a normalising constructor removes the ambiguity at the source, after which derived
  impls are correct *and* cheaper than hand-written ones. Same end, one fewer thing to get wrong.
- **`witness_order` borrows instead of returning `(u8, u16, u16)`** — a total order over a set
  cannot be a fixed-width tuple. `ReadWitness` still gets no `Ord`, unchanged reasoning
  ([open_record.rs:255-259](../../../../src/ng/locus_generation/pileup/open_record.rs#L255)).
- **No field for the borrowed base — deferred on the measurement, not on the design** (spec §3.2,
  §6). When it is built: named for what it holds (`UnwitnessedBases`, since "provenance" names a
  topic; `BorrowedBases` was the equally accurate alternative), and **not** an `Option` — the
  spec's "8 bytes and no allocation" argument holds only for a heap-backed type, and inline an
  empty set already allocates nothing.
- **`ObservationRow` becomes `KeyedObservation`; `ObservationKey` keeps its name.** The accumulator
  holds the same information as the public type in the fold's layout — its key grouped rather than
  flattened — so once the public one is an observation, "row" is the last place in the crate where
  the word names something that is not a table. `KeyedObservation { key: ObservationKey, support,
  chain_ids }` says what it is and what it is keyed by. The key was already named for what it is and
  does not move. **26 uses across 4 files** (spec §4, which left this decision here).
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

**The encoding is no longer open: runs, two inline.** Every witness in 225 million DNA-seq
event-folds is one run, and the RNA-seq case that motivates the change is two (spec §8). A bitmask
loses on that evidence — it would pay 16 bytes to say "positions 0..40" for a case that is one run
in every observation measured. The encoding stays private behind `runs()` so it can still move.

- `OPEN:` **How often the hole fires on real RNA-seq.** The spec's fixture proves reachability, not
  frequency, and the aligner's junction behaviour decides it — some emit a short deletion where the
  truth is a short intron, which is exactly the shape that widens a record across a junction. **Not
  a design question**: the failure is demonstrated and the shape does not change with the rate. It
  is the implementation plan's input, for ordering.
- *Impl-time confirmations, not design items:* whether `runs()` returns `(u16, u16)` pairs or a
  `Run` newtype (a newtype if any call site ever holds one loose).

---

## 5. Reconciliation with existing code

Every row read at the cited line (2026-07-30), in the `ng-generic` worktree. Convergence on
existing types, not new types beside them.

| this doc's name | existing code | action |
|---|---|---|
| `ReadWitness` | `ReadCoverage` [mod.rs:213](../../../../src/ng/locus_generation/mod.rs#L213) | **rename** — 193 uses of the type and 91 of the field `read_coverage` across 12 files (§3) |
| `WitnessedRefPositions` | `RefSpan` [open_record.rs:65](../../../../src/ng/locus_generation/pileup/open_record.rs#L65) | **replace** — same role (the fold's witnessed extent, reference coordinates, never empty), widened from one run to a set. Its "no empty span" invariant carries over verbatim |
| `WitnessedLocusPositions` | `ReadCoverage::Observed { offset_in_locus, positions_covered }` [mod.rs:216-224](../../../../src/ng/locus_generation/mod.rs#L216) | **replace the payload**, and rename the variant to `Partial` — the two `u16`s are the one-run case of the set |
| `is_flush_left` / `is_flush_right` | [mod.rs:327-348](../../../../src/ng/locus_generation/mod.rs#L327) | keep signatures; body derives from the first and last run |
| `num_obs_along_locus`'s clamp | [mod.rs:69-104](../../../../src/ng/locus_generation/mod.rs#L69) | **keep** — the comment there explains why the bound is not expressible on the type, and a set does not change that; the loop iterates runs instead of one range |
| `complete_observations` | [mod.rs:123](../../../../src/ng/locus_generation/mod.rs#L123) | unchanged — `Complete` is still a variant, still an equality test |
| the two borrow sites | Insertion [open_record.rs:1217-1224](../../../../src/ng/locus_generation/pileup/open_record.rs#L1217), Deletion [:1228-1238](../../../../src/ng/locus_generation/pileup/open_record.rs#L1228) | **unchanged** — deferred (§1.2). Recorded here because the guard `offset >= consumed_until` is already exactly "no `Match` emitted this position", so the day it is built, this is where |
| the drop path | [open_record.rs:1335-1345](../../../../src/ng/locus_generation/pileup/open_record.rs#L1335) | **narrows** to "witnessed nothing"; the set-of-read-ids mechanism and the subtract-prior-contribution step survive untouched — a read still becomes non-contiguous when the window widens |
| `witness_of` | `coverage_of` [open_record.rs:184-212](../../../../src/ng/locus_generation/pileup/open_record.rs#L184) | renamed with its type; resolves a set against the footprint. The both-ends clamp and the `Complete` short-circuit are kept |
| `witness_order` | `coverage_order` [open_record.rs:261](../../../../src/ng/locus_generation/pileup/open_record.rs#L261) | renamed with its type; signature changes to borrow, and `pub(super)` for the differential stays |
| the observation identity | `ObservationKey` [open_record.rs:228](../../../../src/ng/locus_generation/pileup/open_record.rs#L228) | **name unchanged**; the witness field inside it grows a set |
| `KeyedObservation` | `ObservationRow` [open_record.rs:237](../../../../src/ng/locus_generation/pileup/open_record.rs#L237) | **rename** — the fold's accumulator, same fields (§3) |
| the out-of-range gate | `PileupGeneratorConfig::check` [generator.rs:123-131](../../../../src/ng/locus_generation/pileup/generator.rs#L123), `MAX_RECORD_SPAN_CEILING` [:43](../../../../src/ng/locus_generation/pileup/generator.rs#L43) | **reuse as-is** — already rejects a footprint wider than a run can describe, so §2 needs no new error |
| the STR generator's mints | [ssr.rs:770](../../../../src/ng/locus_generation/ssr.rs#L770), [:821-822](../../../../src/ng/locus_generation/ssr.rs#L821), [:889](../../../../src/ng/locus_generation/ssr.rs#L889), [:989](../../../../src/ng/locus_generation/ssr.rs#L989) | **call sites unchanged** — `from_left`/`from_right` keep their signatures, which is why they keep them |

---

## 6. Test & bench shape

**Where they live.** The invariant tests in `witness.rs`, the fold's in `pileup/open_record.rs` and
`pileup/tests.rs`, the differential and the census in `pileup/parity.rs`. **No new bench** — no
bake-off here, and the cost question is an allocation count, not a wall time.

**Spec §7 lists the acceptance checks; two of them are this doc's to shape.**

The **canonicality property** is the failure that is silent, so it needs a property test rather
than a fixture: two `WitnessedLocusPositions` built from different orderings of the same runs
compare equal and hash equal, for any generated multiset. Its mutation is the real check — remove
the normalising step in the constructor and observation counts must inflate.

The **spliced fixture** is the regression anchor for the change's whole purpose, and it exists
already: spec §8's 15 bp intron with a 20 bp deletion widening the record across it. It belongs in
`pileup/tests.rs` rather than in a dump example, because it asserts a fold outcome — the read
appears with a two-run witness where today it is absent from the record. Keep the one-position
sensitivity in the test's comment: at a 16 bp deletion the footprint stops short of exon 2 and the
read is recorded normally, which is what makes the fixture a knife-edge rather than a decoration.
