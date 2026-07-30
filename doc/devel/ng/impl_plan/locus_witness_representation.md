# ng — what a locus observation says a read witnessed: implementation plan

**Status:** draft, 2026-07-30. The build order for the witness representation: the vocabulary
renames, the `witness.rs` module, the two witnessed-set types, the fold change that stops
discarding a holed read, and the consumers that follow. Design is settled in
[`../spec/locus_witness_representation.md`](../spec/locus_witness_representation.md) (spec) and
[`../arch/locus_witness_representation.md`](../arch/locus_witness_representation.md) (types &
interfaces). **This plan turns that design into order; it is not a place for new design** — every
open question is resolved in the spec §8, and the one item still `OPEN:` in the arch doc (how often
the hole fires on real RNA-seq) is an input to *ordering*, not to any step's content.

Follows the generic locus generator's plan
([`locus_generation_pileup_generator.md`](locus_generation_pileup_generator.md)), complete through
its Milestone D, whose measurements are what made this design decidable.

---

## Scope

**In:** the crate-wide rename of the observation and witness vocabulary; a new
`src/ng/locus_generation/witness.rs`; `WitnessedLocusPositions` and `WitnessedRefPositions`;
`apply_events_into` returning a set instead of a span and no longer discarding a holed read; the
narrowed drop path; `witness_of` and `witness_order`; the `from_run` constructor; the consumers
(`num_obs_along_locus`, the flush predicates, both dump tools, the divergence census); and the two
regression anchors — the canonicality property and the spliced fixture.

**Out (with a home, nothing dropped):**

- **Saying which bases the read supplied** (`UnwitnessedBases`) — deferred on the measurement,
  design recorded in **spec §6**. 8 borrowed bases in 225 million event-folds.
- **Consuming partial evidence** — step 7's censored likelihood, owned by
  [`locus_generation_pileup.md`](../spec/locus_generation_pileup.md) §10.
- **Sealing `ReadWitness`'s fields** — spec §6, revisit with a later arch pass.
- **Whether `reads_without_observation` survives at all** — spec §6; it needs comparing against
  `reads_silent_over_footprint` first.
- **Parallelism** — deferred whole (`locus_generation.md` §9).

## Principles (how the order was chosen)

- **The renames land first, alone, and behaviour-free.** They are by far the largest diff — 193 uses
  of one type, 102 of one field — and entangling them with a behaviour change would make every
  later diff unreadable and `git bisect` useless. Milestone A changes no bytes of output.
- **Types first, then implementation**, within every milestone (project rule).
- **Compute before you act.** C1 builds the witnessed set and *still* returns `None` on a hole, so
  it is provably byte-identical; only C3 changes what the walk emits. This split is not
  hypothetical: the §8 measurement probe did exactly C1's work and all 275 locus-generation tests
  stayed green, which is the evidence that the shape is behaviour-neutral.
- **Isolate the silent failures.** Two steps here fail quietly rather than loudly — B2
  (a non-canonical set inflates observation counts instead of erroring) and C3 (a read appears or
  vanishes from a record). Each is **its own commit with its oracle green before and after**.
- **Verify against ground truth, not self-consistency.** The STR dump's byte-identity and the
  generic parity anchor are external oracles; the spliced fixture is drawn from a real failure.
- **Incremental, with pauses.** One milestone, then stop.

## Preconditions (already in place — the executor confirms these before A1)

1. **Worktree `/Users/jose/devel/pop_var_caller-ng-pileup`, branch `ng-pileup-generator`.** It
   carries both the docs and the generator code since the `ng-generic` merge, so no branch is
   created and nothing is read out of the sibling worktree. Confirmed by the counts each step
   cites, which were taken in this tree.
2. `cargo test --release --lib ng::locus_generation` green — 275 tests, including
   `ng_agrees_with_production_where_production_fabricated_nothing` and
   `ng_emits_the_same_bytes_in_a_second_process`.
3. The STR dump byte-identity oracle runs on the committed fixture and on a tomato CRAM.
4. **No measurement probe in this tree.** The §8 instrumentation was throwaway and lives
   uncommitted in the *other* worktree (`pop_var_caller-ng-generic`); its numbers are in spec §8
   and its spliced fixture is rebuilt properly at D6. **Do not port it across** — a probe-free tree
   is what keeps A1's diff to the rename.
5. The four specs and the arch doc are committed and use the new vocabulary, so the code is what is
   out of step — not the docs.

---

## Milestone A — the vocabulary (no behaviour change)

Every step here is mechanical and must leave output byte-identical. Any test that changes
expectations is a step that did more than rename.

- ✅ **A1.** `ObservedSequence` → `SequenceObservation`, and the field
  `SampleLocusObservations::observed_sequences` → `observations`. 39 uses of the type across 7
  files, 102 of the field. *Depends:* —. *Source:* spec §1, §4.
- ✅ **A2.** `ReadCoverage` → `ReadWitness`, and the field `read_coverage` → `read_witness`. 193 and
  91 uses across 12 files, including the `read_coverage` column in both dump tools' TSV output.
  *Depends:* A1. *Source:* arch §3.
- ✅ **A3.** `coverage_of` → `witness_of`, `coverage_order` → `witness_order`. Signatures unchanged
  at this point. *Depends:* A2. *Source:* arch §3.
- ✅ **A4.** The variant `ReadWitness::Observed` → `Partial`. 55 match sites. *Depends:* A2.
  *Source:* spec §3.1, arch §3.
- ✅ **A5.** `ObservationRow` → `KeyedObservation`; `ObservationKey` keeps its name. 26 uses across
  4 files. *Depends:* A1. *Source:* arch §3.
- ✅ **A6.** The code's own doc comments: "row" and "cell" → "observation", where they name this
  type. **Leave the aligner's ~340 matrix rows and the dump tools' TSV rows alone** — those are
  real tables. *Depends:* A1–A5. *Source:* spec §6.

> **Checkpoint A: the vocabulary is one word, and nothing moved.** The STR dump is byte-identical
> apart from its renamed column header; the generic parity anchor and the second-process
> byte-identity test are green. Pause for review.

## Milestone B — the witness types (types first, nothing wired)

- ☐ **B1.** Create `src/ng/locus_generation/witness.rs` and **move** `ReadWitness` and `LocusLen`
  into it, re-exported from `locus_generation` so no import path changes. Pure move, no edits.
  *Depends:* A6. *Source:* arch *Module home*.
  **Two items the Milestone A review deferred here, because the move is where they cost
  nothing** ([review](../../reports/reviews/ng_locus_witness_representation_a_2026-07-30.md)
  Mi11, Mi13): (a) `witness_order` exists **twice**, byte-identical, in `pileup/open_record.rs`
  and in `ssr.rs`'s tally — `open_record`'s comment justifies withholding an `Ord` impl because
  it "would export *this file's* sorting convention", which the STR copy refutes, so the move
  should absorb the comparator rather than carry two copies across (a reviewer verified that
  deriving `Ord` and reducing both call sites leaves 275 tests passing); (b) three files import
  the witness vocabulary through `super::super::` while two use the crate-absolute path, and in
  `open_record.rs` the same spelling resolves to two different modules — converting them here
  makes `grep crate::ng::locus_generation::witness` answer "who depends on this".
- ☐ **B2.** `WitnessedLocusPositions`: private field, `new` normalising (sort, merge adjacent and
  overlapping), `one_run`, `runs`, `positions_covered`, `is_flush_left`, `is_flush_right`. Encoding
  is runs with two inline. Derived `Eq`/`Hash`, sound because construction canonicalises. **Own
  commit, do not bundle** — a non-canonical representation inflates observation counts silently
  rather than failing. Its guard is the property test in the same commit: sets built from different
  orderings of the same runs compare and hash equal, over generated multisets. *Depends:* B1.
  *Source:* arch §1.1, spec §3.3, §8 (encoding).
- ☐ **B3.** `WitnessedRefPositions`, `pub(super)`, the fold-axis counterpart: canonical runs in
  reference coordinates, never empty. Replaces `RefSpan`'s role but is not yet used by the fold.
  Carries `RefSpan`'s "no empty span" invariant verbatim. *Depends:* B2. *Source:* arch §2.

> **Checkpoint B: both set types exist, canonical, unused.** Their own tests pass; the walk is
> untouched and every oracle from Checkpoint A is still green. Pause for review.

## Milestone C — the fold (the behaviour change)

- ☐ **C1.** `apply_events_into` accumulates runs into a caller-owned buffer and returns
  `WitnessedRefPositions` — **and still returns `None` on a hole.** Byte-identical by construction;
  the §8 probe demonstrated this exact shape against the full suite. *Depends:* B3. *Source:*
  arch §2.
- ☐ **C2.** `witness_of` resolves a `WitnessedRefPositions` against the final footprint into a
  `ReadWitness`, keeping the both-ends clamp and the `Complete` short-circuit. Still one run in
  practice, so still byte-identical. *Depends:* C1. *Source:* arch §2.
- ☐ **C3.** **Stop discarding a holed read.** `apply_events_into` returns `None` only when the read
  witnessed nothing inside the record; the drop path narrows to that case, keeping the
  set-of-read-ids mechanism and the subtract-prior-contribution step. **Own commit, do not
  bundle** — this is the step where a read appears in a record it was absent from, and the failure
  mode is a wrong number rather than a crash. Its oracles: the spliced fixture (D6) flips from
  "read absent" to "read present with a two-run witness", the STR dump stays byte-identical, and
  the generic census gains its own class rather than absorbing the change into an existing one.
  *Depends:* C2. *Source:* spec §1 goal 1, §3.1, §4; arch §2.
- ☐ **C4.** `witness_order` borrows instead of returning `(u8, u16, u16)`, and `finalise` sorts
  with it. `ReadWitness` still gets no `Ord` of its own. *Depends:* C3. *Source:* arch §2, §3.

> **Checkpoint C: a read with a hole is recorded instead of discarded.** The spliced fixture
> passes, output is deterministic across two processes, and the STR path has not moved. Pause for
> review.

## Milestone D — consumers and surfaces

- ☐ **D1.** `num_obs_along_locus` iterates the runs instead of one range. **Its clamp stays** — the
  comment there explains why the bound is not expressible on the type, and a set does not change
  that. *Depends:* C4. *Source:* spec §4, arch §5.
- ☐ **D2.** `is_flush_left` / `is_flush_right` derive from the first and last run; signatures
  unchanged. *Depends:* D1. *Source:* arch §1.1.
- ☐ **D3.** `ReadWitness::from_run`, the interior-run constructor the deferred note on the variant
  asked for. `from_left` / `from_right` keep their signatures, so the STR generator's four call
  sites do not move. *Depends:* D2. *Source:* arch §1.1, spec §1 goal 3.
- ☐ **D4.** Both dump tools print the set rather than one run, and the generic dump's invariant
  check (`offset_in_locus + positions_covered <= footprint`) becomes a check per run. *Depends:*
  D3. *Source:* spec §4, §8. **Deferred here by the Milestone A review** (Mi15): `witness_label`
  exists in **three** example dumps, identical down to a shared seven-line comment, and the three
  have already drifted — `ng_ssr_loci_dump` emits `partial:left`/`partial:right` where the other
  two emit `partial_left`/`partial_right` but keep `partial:interior`. The two research dumps'
  labels are their own output and must not move by accident, so D4 should share the *derivation*
  and let each tool spell its own strings — deciding the drift rather than inheriting it. D4 also
  owns the generic dump's `observed:<offset>+<positions>` value label and its
  `rows_observed`/`reads_observed` counter keys, the last user-visible uses of A4's retired
  variant name.
- ☐ **D5.** The divergence census counts holed witnesses and the positions inside them — the
  counters the §8 measurement used, kept this time rather than thrown away. *Depends:* D4.
  *Source:* spec §4.
- ☐ **D6.** The spliced fixture as a permanent test in `pileup/tests.rs`: a 15 bp intron plus a
  20 bp deletion widening the record across it, asserting the read appears with a two-run witness.
  Its comment records the knife-edge — at 16 bp the footprint stops one position short of exon 2
  and the read is recorded normally either way. *Depends:* C3. *Source:* spec §7, §8; arch §6.

> **Checkpoint D: every consumer reads the set, and the surfaces show it.** Pause for review.

## Milestone E — verification at scale

- ☐ **E1.** STR dump byte-identity on the committed fixture and a tomato CRAM, across the whole
  change. *Depends:* D6. *Source:* spec §7.1.
- ☐ **E2.** The generic parity anchor green on real data, and the census's new class counted and
  floored on the HG002 and tomato runs from spec §8. *Depends:* E1. *Source:* spec §7.2.
- ☐ **E3.** Allocations per observation on the chr1 run, measured against Milestone D's baseline of
  1,647,161 observations at 461 MB peak RSS. The requirement: no allocation on an observation that
  witnessed one run. *Depends:* E2. *Source:* spec §5, §7.5.
- ☐ **E4.** **When a spliced BAM is available:** run the same probe over it and record the hole
  rate. This does not gate any step — it is the number that tells us how much the change bought.
  *Depends:* E3. *Source:* spec §8 (open), arch §4 (`OPEN:`).

> **Checkpoint E: the change is proven at scale.** Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A — vocabulary | STR dump byte-identical apart from its renamed column; generic parity anchor green; second-process byte-identity green. A test whose expectations change means the step did more than rename |
| B — the set types | the canonicality property test over generated multisets, plus its mutation: remove the normalising step and observation counts inflate |
| C — the fold | C1/C2 byte-identical against the full suite; C3 against the spliced fixture (read absent → present with two runs), the STR dump, and a census class that is counted rather than absorbed |
| D — consumers | `num_obs_along_locus` on a fixture mixing complete and multi-run witnesses; the dumps' per-run invariant asserted on every locus of every run |
| E — at scale | STR byte-identity on a tomato CRAM; the generic anchor on HG002 chr1 and two tomato samples; allocations per observation against the 461 MB baseline |

## Out of scope (next plans)

- **`UnwitnessedBases`** — spec §6 holds the design and the measurement that deferred it. It gets a
  plan when a count asks for one.
- **Step 7's censored likelihood** — the first real consumer of a partial witness, and the reason
  this evidence is being kept. Home: the step 7 spec.
- **Sealing `ReadWitness`'s fields** and **the fate of `reads_without_observation`** — spec §6.
