# ng — the alignment cursor, Milestone D1+D2: review

*Three agents, each in its own worktree over the same staged diff (base `ee0c94b`), on
**correctness of the per-region reset**, **design fit and API**, and **test reliability**.
Implementation report: [ng_alignment_cursor_d_2026-08-02.md](../implementations/ng_alignment_cursor_d_2026-08-02.md).
Fixes: [ng_alignment_cursor_d_fixes_2026-08-02.md](../implementations/ng_alignment_cursor_d_fixes_2026-08-02.md).*

**Two Blockers, both real in the sense that mattered, and only one of them a defect in the
code.** The reliability agent proved the feature could be switched off with the whole suite
green. The design agent's Blocker was a correct reading of the evidence that measurement then
refuted — and the measurement it forced found something larger that predates this change.

---

## Blockers

### The whole of D2 could be switched off and nothing failed

`open_walk`'s chromosome comparison replaced by `if true`, so that **every region mints a fresh
cursor** — precisely the per-region query D2 exists to remove — left **1,557 tests passing**.

They had to. A generator that keeps a cursor and one that rebuilds it emit exactly the same
loci; that is the correctness requirement. The only thing that moves is what the reader
*avoided*, and `CursorCounts` — documented in `cursor.rs` as "the only way to tell whether it is
working" — was read nowhere outside its own module's tests. Spec §11.5 asks for exactly this
("the saving is asserted, not assumed"), and at C there was no caller to assert it through.

**Fixed** by `PileupGenerator::cursor_counts()`, which sums retired chromosomes' cursors and
asks the live one, plus two tests. Deliberately *not* on `PileupGeneratorCounts`, which the dump
tools print verbatim against byte-identical baselines.

### The reference accessor's lifetime, and the trap perf-review L2 names

The design agent found that this change moves the mismatch filter's accessor from *per region*
to *per chromosome*, that `RawChromReader::fetch` extends its window forever unless
`evict_before` is called, that nothing in `src/ng/read/` calls it, and that L2 warns by name
that a shared accessor walking a contig grows to ~250 MB on chromosome 1. It also observed
correctly that the tandem-repeat-targeted fixture — 0.64 % of positions covered — reposition
rather than extends, so it **cannot show the term**.

That is sound reasoning from the evidence to hand, so it was measured rather than argued: a
synthetic contiguously-covered 20 Mb contig at 30×, with the fixture build in a separate process
so its gigabytes could not mask the walk.

**Refuted on magnitude.** The delta is +3.4 MB on a 25.6 MB baseline at 20 Mb, it does not scale
with contig length (it is non-monotone across four spans and smallest at the largest), and the
walk is 1.32× faster there too. The change adds a roughly fixed ~6 MB, not a contig-scaled
buffer. Full table in the implementation report.

**⚠ And the measurement found something bigger, which is not D's.** *Both* sides grow at ~1 byte
per base. `PileupGenerator`'s own `reference: Arc<R>` and the preparer's accessor are
**run**-lifetime and can never be evicted — `Arc` cannot give the `&mut self` `evict_before`
needs — so L2's ~250 MB is already there at `ee0c94b`, and ng's whole-genome memory claim
(30.1 MB for 18.5 M loci) was measured on the fixture that hides it. **Raised for the owner at
Checkpoint D**; it belongs with the deferred per-chromosome reference registry (spec §12).

---

## Major, fixed

- **A region abandoned half-walked leaked a record into the next region.** Deleting
  `open_records.drain_all()` / `reset()` from `WalkerState::begin_region` passed all 1,555
  tests. Two agents found it independently, and both proved the same failure: the leftover
  record is *finalised by the next region's walk*, so that region leads its output with a locus
  at a coordinate nobody asked about — with reads folded under `read_id`s since restarted at
  zero. `begin_segment` on a half-drained walk reaches it; none of the six new tests abandoned
  a walk. Fixed by a test.
- **A failure after the first region on a chromosome was blamed on the first.**
  `PreparedSampleReads::region` is what every shed error names, and it is set when the
  chromosome's stream is built; the reassignment in `move_to_region` could be deleted with the
  suite green. The one shipped test for a preparation failure walks a single region, where the
  two are the same region. Fixed by a test.
- **`RegionReadSource`'s contract omitted the property that makes the walker's look-ahead
  discard lossless.** The walker throws its peeked read away at every region boundary, and that
  is safe only because the source *replays*. The trait said only "point the source at `region`".
  An implementor that consumed would drop one read per boundary, silently — the failure class
  that cost 3,830 loci with 1,471 tests green. D3 adds the second implementor. Fixed in the
  contract.
- **`fold_region_walk` justified its plain sum by a property this change removed.** "each
  region's walk owns its active set" is no longer true; the sum is safe because
  `ActiveReads::begin_region` zeroes `silent_exits` where `reset` preserves it. The comment a
  maintainer reads while deciding whether the fold is safe asserted the old reason. Fixed.

## Minor, fixed

`IngestError::Cursor`'s claim that "every `CursorError` names the path" (one variant does not,
and it duplicates a condition the enum already has); the box made the generator unconditionally
`!Send` for no gain (`+ Send` added); `OpenReadQuery` now reachable by a condition the split
exists to separate (`AfterFailure`), documented; `enter_chromosome`'s `chain_ids.reset()`
claimed a purpose `begin_region` already serves (it earns its place on the failure path
instead); the ordering rationale in `PileupWalker::move_to_region` was untrue as written; a
chromosome test discarded its middle walk; eleven stale doc sites across `src/ng/`, the bench
and both example tools; `genome_walk.rs`'s "changed in one respect" header.

## Known and accepted

Mutations that survive, each defensive code with no reachable failure the agents could
construct — recorded rather than papered over with a contrived fixture:

- `pending.clear()` in `PileupWalker::move_to_region`. Instrumented across two fixtures,
  `pending` is always empty at a region boundary; it becomes reachable only when one walker tick
  closes several records at once (a widened record over a long deletion) *and* the caller
  abandons mid-drain.
- `*chrom_id`, `*walker_pos`, `*last_admitted_chrom_id` in `begin_region`. `enter_chrom`
  overwrites all three whenever the new region has a read, so they matter only for a region with
  no reads at all, where nothing is observable.
- `next_read_id = 0` and `by_read_id.clear()` in `ActiveReads::begin_region`. The first is
  fidelity with a fresh walker, not correctness. The second is unreachable because the
  allocator's `pending_mates` is reset, so no mate cross-link can name a stale id.
- `PreparedSampleReads::move_to_region`'s error arm is dead through this generator:
  `WrongChromosome` is guarded by `open_walk`'s comparison and `AfterFailure` by the `failed`
  latch.

## Deferred to the owner

- **`arch/locus_generation_pileup.md` and `spec/locus_generation_pileup.md` now misdescribe this
  module** — including an "**As built**" block reproducing the deleted three-parameter
  signature, "One `PileupWalker` per region", "one read query per segment", and arch §5's
  explicit "No spec fold-in is owed here". A fold-in is now owed. Not done here: this skill does
  not edit the design docs, and Milestone D's plan does not include a fold-in. **Raised at
  Checkpoint D.**
- **The generic path's per-read-group drop tallies are unreachable.** `RegionReads::Drop` folded
  each query's `ReadFilter` tally back into the `AlignmentFile`, which is what made
  `SampleReads::counts()` report step-1 drop rates; a cursor has no `Drop` and `SampleCursor`
  exposes `CursorCounts` instead. Spec §3 promises "a caller reads `cursor.counts()` whenever it
  likes" for exactly this tally. Nothing consumes it outside its own module today. **Belongs on
  Milestone F's inventory.**
