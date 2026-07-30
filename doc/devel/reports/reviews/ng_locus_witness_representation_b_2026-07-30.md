# Code Review: ng locus witness representation — Milestone B (the witness types)
**Date:** 2026-07-30
**Reviewer:** rust-code-review skill (orchestrator) + 4 category sub-agents, one isolated worktree each
**Scope:** commits `89030aa..11de107` on `ng-pileup-generator`
**Status:** Approve-with-changes — all applied in `19aa2a2`

---

## 1. Scope

- **What:** the four commits that build the witness types and wire none of them in — `d82a0e8`
  (move `ReadWitness`/`LocusLen` into `witness.rs`), `0058c2e` (one `sort_key`, crate-absolute
  imports), `53c9af5` (`WitnessedLocusPositions`), `11de107` (`WitnessedRefPositions` + the
  shared `canonicalise_runs`).
- **Reviewed against:** `11de107`.
- **In scope:** `witness.rs` (new), `locus_generation/mod.rs`, `pileup/{open_record, parity,
  genome_walk}.rs`, `ssr.rs`, `proptest-regressions/ng/locus_generation/witness.txt`.
- **Out of scope:** `src/pileup/` (production, frozen); Milestone C's fold change (not
  written); Milestone A's renames (reviewed separately).
- **Categories dispatched:** `reliability` (the invariant is the deliverable and the types are
  unused, so "can it be violated?" is the only real question), `refactor_safety` (a claimed
  *pure move* and a comparator consolidation), `module_structure` + API-fitness (placement,
  visibility, and whether the API can support C1), `naming` + `smells` (two constructors, two
  conventions; long docs that must earn their length).
- **Skipped:** `errors` (no error type, no `Result`, no `?` in the diff), `defaults`,
  `unsafe_concurrency`, `tooling`, `idiomatic` (folded into naming/smells for a 4-agent
  fan-out proportionate to a 2,052-line diff).

## 2. Verdict

**Approve-with-changes.** Both of the milestone's mechanical claims survived independent
checking: B1 *is* a pure move (two lines changed, both intra-doc links gaining `super::`), and
the `sort_key` consolidation changed no emitted order (both deleted bodies byte-identical to
the survivor, all four call sites keyed and oriented the same). The canonical form itself
survived a deliberate attack — reversed runs, duplicates, containment, touching,
`(0, u16::MAX)`, 500 descending runs, `(u32::MAX-1, u32::MAX)`, `from_positions(u64::MAX)` —
with no panic, no overflow, and no second representation. Derived `Eq`/`Hash` over `SmallVec`
is genuinely spill-invariant, verified by asserting `merged.spilled() && !direct.spilled()` in
the existing equality test.

What the review found instead: **the tests were thinner than the invariant**, and **the API
could not do what the arch says it must**.

## 3. Execution status

| command | result |
|---|---|
| `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings` | clean, at review time and after fixes |
| `cargo test --release --lib ng::locus_generation` | 287 → **292** passed, 0 failed |
| `cargo test --lib --bins --tests --examples --all-features` | 2,818 → **2,822** passed, 0 failed |
| STR dump on the tomato CRAM | byte-identical throughout |
| `cargo doc --no-deps --all-features` | unresolved links unchanged against the pre-B tree |

Findings labelled "Needs verification": **0** — every finding below was demonstrated by a
mutation or a probe in an isolated worktree, with its output quoted.

## 4. Open questions

1. **Should `WitnessedRefPositions` live with the fold at all?** Its field is private to a
   ~3,100-line module, and the C1 fold will be written *inside* that module — so direct
   construction, bypassing `canonicalise_runs`, is available exactly where the risk is. That
   contradicts the rationale `witness.rs` was created under ("a private field needs one module
   boundary to be private to"), while the arch's *Module home* list puts only the three public
   types in `witness.rs`. **Not resolved here** — it is the arch's call, and it is the one
   item this review leaves open. Raised at Checkpoint B.
2. **Is `pub(super)` still the right visibility for the walk's internals?** Milestone A's
   review proposed demoting three of them and the resolution cited `open_record::witness_order`
   as the load-bearing contrast — which B1 has since **deleted**. Recommended: settle it once
   as policy rather than per milestone.

## 5. Top 3 priorities (all applied)

1. **The property test could not see a normaliser that loses positions** — the guard the plan
   names for B2, blind to the one class of defect it exists to catch.
2. **The API could not satisfy arch §2's caller-owned buffer**, so C1 would have allocated per
   (read × widen) on the multi-junction RNA-seq case the milestone exists for.
3. **`start()` / `end_exclusive()` invited the exact C1 mistake** the milestone is meant to
   remove — convergent across two categories.

## 6. Findings

### Major

**M1: `witness.rs` — the property test asserts the canonical form's shape, never its content**
- **Categories:** reliability
- **Confidence:** High. Mutating the merge's `open_end.max(end)` to `end` — truncating a
  containing run, losing witnessed positions — left `a_set_is_the_same_set_whatever_order_its_runs_arrive_in`
  **passing**, and both reference-axis tests with it. Only the single `contained` fixture line
  caught it. Sorting precedes merging, so both orderings lose the same positions and
  order-independence still holds.
- **Applied:** the test now compares the position sets of input and output (`BTreeSet` over the
  flattened runs) and the `positions_covered` total. Verified: fails under the mutation, passes
  clean.

**M2: `open_record.rs` — the constructor set cannot satisfy arch §2's "a buffer the caller owns and the callee clears"**
- **Categories:** module_structure
- **Confidence:** High. The reviewer wrote the fill-in-place signature against the type and got
  four `E0599`s (no `clear`, `push`, `canonicalise`, `default`), then measured the
  value-returning fallback with `spilled()`: `two-run witness len=2 spilled=false, three-run
  witness len=3 spilled=true`. The cost is per (read × widen) because `refold_live_reads`
  rewrites `FoldedReadState.witnessed` on every widen.
- **Applied:** `take_from(&mut WitnessedRefRuns)` and `refill_from(&mut self, &mut
  WitnessedRefRuns)`, the latter swapping so the buffer inherits the old set's storage. A test
  asserts the buffer comes back empty and that a failed refill leaves both sides unchanged.

**M3: `open_record.rs` — `start()` / `end_exclusive()` name the enclosing span on the type whose point is that the span lies**
- **Categories:** refactor_safety, naming — *convergent*
- **Confidence:** High. `witness_of` computes `past_last` as `witnessed.end.saturating_add(1)`
  against the **inclusive** `RefSpan`; substituting `end_exclusive()` type-checks and is wrong
  in both that expression and the neighbouring `debug_assert`. The type's own test asserted
  `end_exclusive() - start() == 21` where six positions were witnessed — the file demonstrated
  the wrong number and exposed the two accessors to compute it.
- **Applied:** both deleted (they were dead code); `positions_covered()` added in their place.
  C1 walks `runs()`. **Recorded in the plan:** when C2 takes a set, `witness_of`'s `+ 1`
  disappears and its `debug_assert` comparison becomes `>`.

**M4: `witness.rs` — `one_run`'s overflow rejection is not discriminated by any test**
- **Categories:** reliability
- **Confidence:** High. `checked_add` → `saturating_add` left all 287 tests green. The only
  boundary case, `one_run(u16::MAX, 1)`, returns `None` under checked, wrapping *and*
  saturating add, because each sum lands `<= start` and `start >= end` rejects it.
- **Applied:** `(60_000, 10_000)` added, where saturation would return a 5,535-position set.
  Verified: fails under the mutation.

**M5: `witness.rs` — `LocusLen::from_positions` has no test at all**
- **Categories:** reliability
- **Confidence:** High. Replacing the saturation with `positions as u16` left 287 green. A
  65,536-position region becomes `LocusLen(0)`, on which every witness is flush-right and
  zero-long — a wrong depth everywhere, no panic.
- **Applied:** saturation and `of_region` both tested. Verified against the mutation.

**M6: `witness.rs` — two constructors, two conventions, one primitive shape**
- **Categories:** naming
- **Confidence:** High. `new` took half-open pairs, `one_run` took offset+length, both `u16`
  pairs of locus positions. Feeding the first an offset/length pair is **not rejected** —
  `(4, 5)` is a valid canonical set covering one position where the caller meant five. Only the
  reversed spelling errors. This is the hazard `LocusLen`'s own doc says it was minted to
  remove, reintroduced on the neighbouring constructor's parameter list.
- **Applied:** `from_half_open_runs` and `one_run_from_offset_and_length`. Newtypes for the two
  components would close it completely and are recorded for whenever `ReadWitness::Partial`'s
  fields are revisited — they are the same two `u16`s, and sealing that variant is deferred
  (spec §6).

### Minor (all applied)

- **`canonicalise_runs` reported success as an ignorable `bool`** while mutating in place:
  `canonicalise_runs(&mut runs); Some(Self(runs))` compiles, builds a degenerate set, and no
  lint objects (verified). Now takes and returns the buffer.
- **The `expect(dead_code)` tripwire was all-or-nothing.** A block-level `expect` is fulfilled
  if *any* item inside still triggers the lint, so wiring one accessor at C1 would have been
  silent — proven by wiring `new` + `start` and getting a clean build. Now per-method.
- **A new unresolved intra-doc link** at `open_record.rs:75`, invisible to the default doc build
  because `pileup`'s items are private (`--document-private-items` shows it). Now a full path.
- **The module doc described the file as it was two commits earlier** — naming two of its four
  residents, and justifying it by an invariant belonging to a type it never mentioned.
- **"No consumer's import path names this module" was falsified by B3 itself**, which imports
  `witness::canonicalise_runs` by path. Both claims qualified.
- **`is_flush_right`'s `saturating_add` → `wrapping_add` survived** — reachable because
  `Partial`'s fields are public and its clamp is a convention. Now tested.
- **`start > end` was never tested**, only `start == end`. Added.
- **The reference axis had no containment fixture**, so it would not have noticed the very
  drift `canonicalise_runs` is shared to prevent. Added.

### Nits (recorded, not applied)

- `witnessed` names a half-open pair in one scope and an inclusive `RefSpan` in the next, inside
  one file — C1 rewrites that function and should collapse the two.
- `parity.rs`'s `ObservationIdentityWithoutGroup` hardcodes `sort_key`'s `(u8, u16, u16)`, which
  C4 replaces with a borrowing order.
- `witness_of` returns `ReadWitness`, not `Option<ReadWitness>`, while every set constructor
  returns `Option` — at C2, clamping a set against the footprint can in principle empty it.
- No test produces a canonical set of more than two runs; the property test reaches only the
  `u16` instantiation of `canonicalise_runs`.

## 7. What's good

- **The canonical form is genuinely robust.** A reviewer attacked it with adversarial input
  across both instantiations and could not produce a second representation, a panic, or an
  overflow.
- **`SmallVec` was the right encoding for a reason beyond size:** it compares and hashes as a
  slice, which removes the inline-vs-spilled equality hazard a hand-rolled array would have
  carried — and the milestone's own test exercises exactly that path (a merge down to two
  runs against two built directly), which a reviewer confirmed is not vacuous.
- **Sharing one normaliser between the two axes** pre-empted the drift that two copies of
  `witness_order` had just produced, one milestone earlier.
- **B1 being a genuinely pure move** is what let two reviewers verify it mechanically rather
  than by reading.

## 8. Commands to re-verify

```
/Users/jose/devel/pop_var_caller-ng-pileup/scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings
/Users/jose/devel/pop_var_caller-ng-pileup/scripts/dev.sh cargo test --release --lib ng::locus_generation
/Users/jose/devel/pop_var_caller-ng-pileup/scripts/dev.sh cargo test --lib --bins --tests --examples --all-features
```

Each new test must fail under its own mutation, not merely pass:

| test | mutation |
|---|---|
| `a_set_is_the_same_set_whatever_order_its_runs_arrive_in` | `open_end.max(end)` → `end` |
| `an_empty_input_or_an_empty_run_is_rejected_rather_than_dropped` | `checked_add` → `saturating_add` |
| `locus_len_saturates_rather_than_truncating` | `positions.min(u16::MAX as u64) as u16` → `positions as u16` |
| `a_run_whose_end_would_overflow_is_still_flush_right` | `saturating_add` → `wrapping_add` |
| `a_run_contained_in_another_neither_survives_nor_shortens_it` | as the first |
