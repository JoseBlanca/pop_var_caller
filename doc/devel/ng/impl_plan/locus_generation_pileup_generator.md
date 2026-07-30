# ng generic locus generator — the changes and the generator (plan 3 of 3)

**Status:** draft, 2026-07-28. Turn the proven-identical copy into ng's generic `LocusGenerator`:
stop the fabrication, add what a `PileupRecord` cannot say, wrap it in a region walk, and measure
what changed. Design settled in [`locus_generation_pileup.md`](../spec/locus_generation_pileup.md)
(spec) and [`../arch/locus_generation_pileup.md`](../arch/locus_generation_pileup.md) (types &
interfaces). This turns that design into build order; it is **not** a place for new design.

Follows [`locus_generation_pileup_prerequisites.md`](locus_generation_pileup_prerequisites.md) and
[`locus_generation_pileup_port.md`](locus_generation_pileup_port.md).

---

## Scope

**In:** taking ownership of the copies outright (A0 — the reference adaptor deleted, `copy_fidelity`
narrowed), the behaviour changes inside them (the rest of Milestone A), the ng-shaped output
(Milestone B), the generator and its region walk (Milestone C), and the measurements that decide
whether the change was worth making (Milestone D).

**Out (deferred, with homes in spec §10):** the cohort merge; the windowed depth/GC statistics; the
`SsrBundle` generator; parallelism; BAQ; cross-region phasing; consuming partial observations
(step 7's freebayes `1/k` scheme); an aggregating accessor on the shared type; moving
`PreparedRead`/`CigarOp` out of `pileup/walker/`.

## Principles (how the order was chosen)

- **The algorithmic heart before the plumbing.** The no-fabrication rule (A) is the reason this step
  exists; it is built and measured before the generator, the region walk, or the dispatch that feed
  it.
- **Types first, then implementation** — `RefSpan` and the `FoldedReadState` fields (A1) before the
  builder that produces them (A2).
- **Isolate the silent failures.** A2, A3 and C2 produce *quietly wrong numbers* rather than
  crashes — a mis-derived extent is a wrong depth, a mis-clamped walk is missing evidence. Each is
  **its own commit with its oracle green before and after**, so `git bisect` can find it.
- **Verify against ground truth, then against a definition.** A's oracle is still production, on the
  loci ng does not change. **That class is narrower than "reads that witnessed the whole footprint",
  which this line used to name:** production can be wrong about a read that witnessed everything
  (spec §3, class 6), so the oracle holds on a fixture where production fabricates nothing rather
  than on a filter over one where it does. Beyond it there is no oracle, so the new behaviour is
  pinned by fixtures written to fail the *wrong* implementation.
- **Ungated / container builds.** All `cargo` via `./scripts/dev.sh`; a native host build at
  completion.

## Preconditions (already in place)

- **Plans 1 and 2 are complete.** The shared locus type is final; the region stream is owned; ng's
  walker is proven to compute what production computes — 1.5 M synthetic records and 437 k real ones,
  zero divergences — and the differential has been shown to fail six ways. **That baseline cannot be
  rebuilt after A2**, which is the first commit that makes the two walkers differ on purpose. A0 sits
  before it precisely because it is the last step the full differential can still verify.
- **`copy_fidelity.rs` guards all eight copies today.** This plan retires it file by file, starting
  at A0; that narrowing is part of each step, not a cleanup at the end.
- **The production defect plan 2 found is fixed** (`5f32a62`), in production and in ng's copy. It sat
  in `apply_events_to_ref_into` — the function A2 replaces — so the fix is inherited rather than
  redone, and the three-read case that pins it is a ready-made fixture for the witnessed-extent rule.
- The dump-tool precedent [`examples/ng_ssr_loci_dump.rs`](../../../../examples/ng_ssr_loci_dump.rs)
  exists and D2 follows it.

---

## Milestone A — take ownership, then stop the fabrication

A0 is a pure refactor with no behaviour change; A1–A5 are the reason this step exists:
**nothing is written into an observation that its read did not witness** (spec §4).

- ✅ **A0 — delete the reference adaptor. Pure refactor, own commit, and it goes FIRST.**
  ng's copies speak production's `MultiChromRefFetcher` for one reason only: they were
  transcribed verbatim, so their signatures are production's. That was never a design
  choice, and from this plan on the two walkers diverge — so `open_record.rs` takes a
  `RefSeq` directly and `RefSeqFetcher` is **deleted**, not renamed (owner, 2026-07-29).

  What goes with it: `to_chrom_ref_fetch_error` and its two lossy spots (a contig *name*
  rendered as an id, the `u64 → u32` narrowing); `WalkerError::Fasta`'s source becomes
  `RefSeqError`; and ng stops importing `MultiChromRefFetcher` / `ChromRefFetchError`
  entirely — one fewer coupling to frozen production. Use **`fetch_into`**, not `fetch`:
  arch §4 already noted the latter allocates a `Vec<u8>` per `open_new`/`widen` where the
  former writes into a caller buffer, and the buffer comes free with the removal.

  **Why first, and not later: it is the last moment the refactor can be proven free.** A0
  changes no behaviour, so the stage-1 differential must stay green across it — and every
  step after A2 deliberately breaks that differential. Done later, this refactor would land
  with no oracle at all.

  *Two things this step must handle, not discover:*
  - **The differential's two sides stop sharing one accessor.** Production's walker takes
    `MockFasta`, ng's takes an `InMemoryRefSeq` over the same bytes. They agree only because
    the generator draws its reference from `ACGT`, where `MockFasta`'s raw passthrough and
    `InMemoryRefSeq`'s canonicalisation coincide. **Pin that**, in `parity.rs`, with a test
    asserting the two serve identical bytes over the fixture contig — otherwise a later
    generator change that introduced a lower-case or ambiguity-coded base would show up as
    "the walkers disagree" and be chased into the walk.
  - **`copy_fidelity.rs` must be narrowed, not retired.** A0 is the first commit that makes
    a copied file genuinely diverge. Drop `open_record.rs` (and `errors.rs`) from its
    checked set **in this commit**, leaving `cigar_cursor.rs`, `decompose.rs`,
    `active_read_set.rs`, `chain_id_allocator.rs` and `tests.rs` still guarded as the
    verbatim copies they still are. Each later step drops the file it changes, so the guard
    keeps protecting what is still a copy instead of being switched off wholesale at the
    first change.

  *Depends:* —. *Source:* arch §1.3, §4 (the `fetch_into` note); owner decision 2026-07-29.

- ✅ **A1 — the state the rule needs.** `RefSpan { start, end }` (1-based inclusive, reference
  coordinates); `FoldedReadState` gains `witnessed: RefSpan` and `read_group: ReadGroupId`;
  `AlleleSupportStats` becomes ng's copy, production's **minus `placed_start`** (`placed_left` is
  kept — it feeds production's QUAL penalty). Types only. *Depends:* —. *Source:* arch §1.2, spec §6.
- ✅ **A2 — the builder stops filling. Own commit, do not bundle.** `apply_events_into` emits only
  what the events cover and returns the extent; `None` when the witnessed positions are
  non-contiguous. Three traps the spec names: it still needs `ref_seq` for an **indel's anchor base**
  when no `Match` emitted it (the one recorded residual); `events_overlapping` does **not** clip a
  deletion to the window, so the extent must be intersected with `[record_pos, record_end)`; and
  `bases.len()` is **not** `positions_covered` (an insertion adds bases without positions, a
  deletion the reverse). *Oracle:* reads whose events tile the footprint must stay byte-identical to
  production — there are no gaps to fill — so the stage-1 differential still passes **on that class**
  and fails only where it should. *Depends:* A1. *Source:* spec §4, §8.
- ✅ **A3 — `widen` extends the REF bucket only. Own commit, do not bundle.** `alleles[0]` grows;
  no other bucket does. Evict `num_obs == 0` buckets at widen — without production's
  append-to-every-bucket the re-fold no longer lands in its old bucket, so stranded empties
  accumulate against a `find_allele_index` that is a linear byte-compare scan. *Depends:* A2.
  *Source:* spec §4, §7.
- ✅ **A4 — coverage at `finalise`.** `coverage_of(witnessed, record_pos, record_end)` →
  `Complete` | `Observed { offset_in_locus, positions_covered }`, resolved once, against the
  **final** footprint. A read `Complete` when it folded becomes `Observed` after a widen with nothing
  about the read having changed. *Depends:* A3. *Source:* spec §4, §6.
- ✅ **A5 — the no-observation path.** A non-contiguous witness yields no observation and the read is
  recorded — as a **per-record set of read ids, not a counter**: the path is reached at every
  position the record is affected at, so a counter multiplies by the footprint length. And a read
  that folded contiguously can *become* non-contiguous when the window widens, so the path must
  **subtract its prior contribution** before dropping it, or a live contribution is stranded in a
  bucket for a read with no row. *Depends:* A4. *Source:* spec §4, §6.

> **Checkpoint A: the fabrication is gone, and the class ng does not change still matches production
> byte for byte.** Pause for review.

---

## Milestone B — the ng-shaped output

- ✅ **B1 — the bucket key.** `ReadContribution` carries the read group; the key becomes
  `(bases, read_coverage, read_group)`. With one read group the row count must be unchanged — that
  is what "free at one read group" has to mean. *Depends:* A5. *Source:* spec §6.

  *Implemented as the row identity realised at `finalise`, not as a fold-time bucket key* —
  `read_coverage` is only knowable against the final footprint (A4), and arch §1.2 puts the
  bucketing at `finalise` for exactly that reason. So the fold keys buckets on bases alone, and
  `observation_rows` re-derives rows **per read**. The step keeps emitting `PileupRecord`, projecting
  the rows back onto the positional allele list, so the stage-1 differential still proves the
  re-derivation faithful — the same "keep the oracle alive across the risky change" that put A0
  first. B2 deletes the projection.
- ✅ **B2 — `finalise` returns `SampleLocusObservations`.** Rows are re-derived from `folded_reads`
  (per-read), not from per-bucket totals — coverage and group are per-read facts. **Sort by
  `(bases, read_coverage, read_group)` before emitting**: `folded_reads` is an `AHashMap` with a
  per-process seed, so fold order would make the output non-deterministic run to run. Chain ids are
  dropped **per read** — none for a read that agreed with the reference across everything it
  witnessed — because row splitting means "the REF row" is no longer a unique row. *Depends:* B1.
  *Source:* spec §3, §6.

  **B2 owns the determinism test, and it must run the walk in *separate processes* — owner,
  2026-07-29.** Milestone A's review found two mechanisms whose only observable effect is the order
  rows are created in: `refold_live_reads`' `ids.sort_unstable()` (delete it and all 151 tests pass,
  in four separate processes too) and the contributor skip in the same function. Neither is pinnable
  from inside one process, because `ahash`'s seed is fixed *within* a process and both walkers run in
  the same one — the tests compare ng against production, never one ng run against another.

  Rather than write a cross-process test for a mechanism this step makes redundant, the decision is
  that **B2's sort is the guarantee and B2's test is where determinism is proven** (owner: "enough
  that B2 makes it so"). Once rows are ordered by `(bases, read_coverage, read_group)`, creation
  order is invisible in the output and the sort at `refold_live_reads` is belt-and-braces, kept
  because it costs nothing on a rare path. The test to write here is therefore: **the same input
  walked in two separate processes emits byte-identical output** — which is the property spec §7 and
  §13 actually claim, and which no test in Milestone A could have made.
- ✅ **B3 — the per-record counters.** `reads_without_observation` (A5's set) and
  `reads_discarded_by_cap`. The cap truncates in the walk, **before any record exists**, so the
  truncated ids must be plumbed into the fold and registered per affected record; the count is those
  absent from `folded_reads` at `finalise`. Two cases have no clean answer and are recorded as such:
  a read truncated where no record is open, and a truncated read whose deletion would have widened a
  record. *Depends:* B2. *Source:* spec §6.

> **Checkpoint B: ng emits its own locus type, deterministically.** Pause for review.

---

## Milestone C — the generator

- ✅ **C1 — `PileupGenerator` + config + counts.** The struct, `PileupGeneratorConfig` (five knobs,
  defaulting to production's `pub const`s **by name**), `PileupGeneratorCounts` (production's
  `RunSummary` fields plus `reads_silent_over_footprint` and `records_outside_region`). *Depends:*
  B3. *Source:* arch §1.1.

  **`max_record_span` is capped at `u16::MAX` (65,535), and the config rejects more — owner,
  2026-07-29.** A `ReadCoverage` run is expressed as two `u16`s, narrowed through
  `LocusLen::from_positions`, which **saturates**. A footprint wider than 65,535 therefore makes a
  partial witness report a truncated `positions_covered` — a wrong number, no error. Production's
  own `--max-record-span` is an unbounded `u32` with a default of 5,000, so inheriting the knob "by
  name" would inherit the hazard; **this is the one knob where ng's constant is not simply
  production's.**

  *Why capping is the right answer rather than widening the run to `u32`:* the cap is not a
  constraint in practice. A locus is at most ~100 bp of reference, and a 5,000 bp record is already
  unreachable with Illumina reads — the existing default is generous by a factor of fifty, and
  65,535 by a factor of six hundred. Widening the run would touch the shared locus type and the STR
  generator that also mints coverage, to buy a range no data can occupy (owner, 2026-07-29).

  A4 left a `debug_assert` in `coverage_of` stating the envelope; **this step is what makes it
  provable** rather than hopeful, so the assert stays as documentation of the invariant and the
  rejection is the enforcement.
- ✅ **C2 — the region walk. Own commit, do not bundle.** `begin_segment` records the region and
  opens nothing (it cannot fail; opening a query can, so the first `next_locus` is where an
  `IngestError` surfaces). The query is `[region.start, region.end + max_record_span]` — the halo,
  without which a record whose footprint crosses the boundary silently loses the support lying
  beyond it. The walk **stops** once `walker_pos > region.end` and no open record is anchored at or
  before `region.end`, or the halo is walked in full at every boundary. Records are dropped on their
  **anchor**. *Depends:* C1. *Source:* spec §2.
- ✅ **C3 — the allocator across segments.** One `ChainIdAllocator` on the generator, `reset()` at
  each region end — it preserves `next_id` and clears `pending_mates`/`active_count`, without which
  a pending mate cross-pairs between contigs and `active_count` leaks toward
  `ActiveReadsExhausted`. **Snapshot its counters at `begin_segment` and fold the delta**: `reset()`
  preserves them and `RunSummary` takes them by *assignment*, so summing per-region summaries would
  triangular-sum `chain_allocations` and `mate_lookup_evictions`. *Depends:* C2. *Source:* spec §8.
- ✅ **C4 — wire it in.** `impl LocusGenerator<()> for PileupGenerator`; fill `GeneratorSet`'s
  `generic` slot, replacing `Unfilled(NotImplemented)`. `LocusGenerationError` gains a `Walker`
  variant. *Depends:* C3. *Source:* arch §2.1, spec §7.

> **Checkpoint C: ng mints generic loci end to end.** Pause for review.

---

## Milestone D — prove it, then measure it

- ✅ **D1 — stage-2 parity.** `project(PileupRecord) -> SampleLocusObservations` in the test module;
  every divergence falls in one of the **six** named classes and is counted, not excused. Build the
  **permanent** anchor here too — *every* locus of a fixture on which production fabricates nothing;
  that is what replaces the stage-1 differential this plan retires. *Depends:* C4. *Source:* spec §3.
  **As written this step said "five" classes and defined the anchor as "loci where every folded read
  witnessed the whole footprint"; D1 disproved both and spec §3 now carries the corrected form.**

  **Two findings changed the step, both recorded in the
  [D report](../../reports/implementations/ng_locus_generation_pileup_generator_d_2026-07-29.md):**

  - **There are six classes, not five.** ✅ **Spec §3 now carries the sixth row and the
    mechanism (2026-07-30).** Production's `widen` appends reference bases to every bucket, and a
    read that was **not a contributor at the widening step** never has that tail revised — so
    production holds that read's haplotype against a stale footprint. The read can be a
    **complete** witness and production still wrong about it, so this is not class 1, and §13.2
    wants the two counts separately.

    **The mechanism as first written here was too strong and is corrected in §3:** production does
    *not* "re-fold nothing" — `process_position` re-folds every contributor into an affected
    record, and production's own comment argues the append is equivalent to a re-fold *"modulo the
    new bytes never being event-modified by this read"*. `refold_live_reads` (ng's, A3) is what
    covers the reads that modulo excludes.
  - **The anchor's predicate is insufficient**, for the same reason: "every folded read witnessed
    the whole footprint" does not select an agreeing class. The anchor is therefore built on a new
    fixture — every read on a contig sharing one event set, so every read is re-folded by its own
    copy of every widening event and none is left stale — and the fixture's property is asserted
    rather than argued.

  `to_pileup_record` is off the differential and off `mock_reference.rs`; its 44 inherited-test
  callers stay, which is a decision carried to Checkpoint D (report §8).
- ✅ **D2 — the dump tool and the fixtures that must fail a wrong implementation.**
  `examples/ng_generic_loci_dump.rs`, following `ng_ssr_loci_dump.rs`. Six new fixtures, because
  **no inherited test exercises the defect**: a read adaptor-masked over part of a multi-base record
  (must be `Observed`, not `Complete`); an interior `N` and a ref-skip (no observation, counted); a
  record widened past an already-expired read (bases must not have grown); an indel whose anchor base
  was masked (the residual, pinned at one base); two read groups on one allele (two rows summing to
  the single-group total); and a deletion at a region boundary (support must match a single-region
  walk). *Depends:* D1. *Source:* spec §12, §13.

  **It forced the two counters that were defined and never incremented**:
  `reads_silent_over_footprint` (a per-read `ever_contributed` flag in the active set — which
  released `active_read_set.rs` from `copy_fidelity.rs`, the fourth of eight) and
  `reads_declined_by_preparer`, which reads zero on every real run because no v1 preparer
  declines anything, and exists so that *is* a statement. **Two of the six fixtures could not
  fail** and were rewritten: the chain-id one had only whole-footprint witnesses, where ng's
  per-read rule and production's positional one coincide; the boundary one had no read starting
  past the boundary, so removing the halo entirely changed nothing.
- ✅ **D3 — the measurements.** Two numbers, both deliverables rather than by-products. **The size of
  production's defect:** how many loci, reads and reference bases production credits to reads that
  never sequenced them — the number that turns the indel-deficit hypothesis into a result or kills
  it. **Throughput:** the dump over one human chromosome against production's `pileup` on the same
  CRAM. There is no candidate-only fallback, so a bad number is a performance problem to solve, not a
  design to reconsider — and the allele-list growth (spec §7) is the first place to look.
  *Depends:* D2. *Source:* spec §7, §13.

  **The defect:** on HG002 at 300×, chr1:1–6 Mb, production credited **871 reads over 162 loci
  (0.33 %) with 1,550 reference bases they never sequenced**; a 10× window over the same
  coordinates gives **zero**, and a tomato CRAM 6 reads over 4 loci. Class 6 (production's stale
  widen) is **zero on real data**. **Throughput:** chr1, same 30× BAM, single-threaded both
  sides — ng **34.7 s** against production's **10.2 s**, so **3.4× slower**. **The peak-RSS half
  of this line is withdrawn** (review 2026-07-30): ng's 461 MB is dominated by the *dump tool's*
  whole-run row buffer — `ObservationRow` is 152 B and the run emitted 1,541,788 loci, so the
  vector spine alone is ≥ 234 MB — and that buffer is region-length-shaped, the shape §7 forbids
  the generator to be. The number sizes the tool, not the generator; re-measuring needs the tool
  to stream first. See the D report §2. Region typing (a pass production does not make) is 5.9 s
  of it; the
  regions average **391 bp** and 36 % of the records the walk finalised are discarded at the
  clamp, so **the first lever is the region grain, not the fold** — `column_depth_truncations` is
  0 and the allele list never grows here. **And the measurement found the differential's own
  `q_sum` tolerance wrong**: a fixed 1e-9 *absolute* grain, which 300× depth (`q_sum ≈ −3,360`)
  reported as an unlisted divergence one grain wide. Now a *relative* tolerance, and a tolerance
  rather than a rounding.

> **Checkpoint D: the generic path mints loci, the defect is measured, and the cost is known.**
> Pause for review — and for the decision whether the read-group grain stays per-`@RG`.

---

## Verification summary

| milestone | proven by |
|---|---|
| A0 | the stage-1 differential still green **in full** — A0 changes no behaviour, and this is the last step at which that is true |
| A1–A5 | the stage-1 differential still green **on reads whose events tile the footprint**; failing elsewhere by design, with each failure enumerated |
| B | one-read-group row counts unchanged; output byte-identical across repeated runs (the sort) |
| C | loci from two adjacent regions concatenate coordinate-sorted with no gap at the join; `records_outside_region` zero-sum across neighbours |
| D1 | every divergence in one of **six** named classes (D1 found the sixth); the fabrication-free fixture green as a permanent anchor |
| D2 | six fixtures, each written so a span-derived or fill-preserving implementation fails it — **mutation-verified**, and two of them rewritten after passing the mutation they were written for |
| D3 | the defect size reported on two organisms and two depths; throughput measured against production's `pileup` at 3.4×, and the cost attributed. **Peak RSS is *not* measured** — the number taken sized the dump tool's row buffer, and is withdrawn (D report §2) |

## Out of scope (next work)

- **Step 7's use of partial observations** — freebayes' `1/k` prefix/suffix scheme, decided but not
  built here (spec §10).
- **The cohort merge**, which is what consumes these loci; **the windowed statistics**, which slide
  over their per-position depth; **parallelism**, deferred whole.
- **A second generic generator** — the active-region or haplotype-window definition of a locus, which
  the `LocusGenerator` trait exists to let sit beside this one.
