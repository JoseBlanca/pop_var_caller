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

**In:** the behaviour changes inside the copied walker (Milestone A), the ng-shaped output
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
  class ng does not change (reads that witnessed the whole footprint). Beyond that class there is no
  oracle, so the new behaviour is pinned by fixtures written to fail the *wrong* implementation.
- **Ungated / container builds.** All `cargo` via `./scripts/dev.sh`; a native host build at
  completion.

## Preconditions (already in place)

- **Plans 1 and 2 are complete.** The shared locus type is final; the region stream is owned; ng's
  walker is proven byte-identical to production's and the differential has been shown to
  discriminate. **That baseline cannot be rebuilt after A2** — the first commit of this plan makes
  the two walkers differ on purpose.
- The dump-tool precedent [`examples/ng_ssr_loci_dump.rs`](../../../../examples/ng_ssr_loci_dump.rs)
  exists and D2 follows it.

---

## Milestone A — stop the fabrication

The one behaviour change, and the reason for the whole step: **nothing is written into an
observation that its read did not witness** (spec §4).

- ☐ **A1 — the state the rule needs.** `RefSpan { start, end }` (1-based inclusive, reference
  coordinates); `FoldedReadState` gains `witnessed: RefSpan` and `read_group: ReadGroupId`;
  `AlleleSupportStats` becomes ng's copy, production's **minus `placed_start`** (`placed_left` is
  kept — it feeds production's QUAL penalty). Types only. *Depends:* —. *Source:* arch §1.2, spec §6.
- ☐ **A2 — the builder stops filling. Own commit, do not bundle.** `apply_events_into` emits only
  what the events cover and returns the extent; `None` when the witnessed positions are
  non-contiguous. Three traps the spec names: it still needs `ref_seq` for an **indel's anchor base**
  when no `Match` emitted it (the one recorded residual); `events_overlapping` does **not** clip a
  deletion to the window, so the extent must be intersected with `[record_pos, record_end)`; and
  `bases.len()` is **not** `positions_covered` (an insertion adds bases without positions, a
  deletion the reverse). *Oracle:* reads whose events tile the footprint must stay byte-identical to
  production — there are no gaps to fill — so the stage-1 differential still passes **on that class**
  and fails only where it should. *Depends:* A1. *Source:* spec §4, §8.
- ☐ **A3 — `widen` extends the REF bucket only. Own commit, do not bundle.** `alleles[0]` grows;
  no other bucket does. Evict `num_obs == 0` buckets at widen — without production's
  append-to-every-bucket the re-fold no longer lands in its old bucket, so stranded empties
  accumulate against a `find_allele_index` that is a linear byte-compare scan. *Depends:* A2.
  *Source:* spec §4, §7.
- ☐ **A4 — coverage at `finalise`.** `coverage_of(witnessed, record_pos, record_end)` →
  `Complete` | `Observed { offset_in_locus, positions_covered }`, resolved once, against the
  **final** footprint. A read `Complete` when it folded becomes `Observed` after a widen with nothing
  about the read having changed. *Depends:* A3. *Source:* spec §4, §6.
- ☐ **A5 — the no-observation path.** A non-contiguous witness yields no observation and the read is
  recorded — as a **per-record set of read ids, not a counter**: the path is reached at every
  position the record is affected at, so a counter multiplies by the footprint length. And a read
  that folded contiguously can *become* non-contiguous when the window widens, so the path must
  **subtract its prior contribution** before dropping it, or a live contribution is stranded in a
  bucket for a read with no row. *Depends:* A4. *Source:* spec §4, §6.

> **Checkpoint A: the fabrication is gone, and the class ng does not change still matches production
> byte for byte.** Pause for review.

---

## Milestone B — the ng-shaped output

- ☐ **B1 — the bucket key.** `ReadContribution` carries the read group; the key becomes
  `(bases, read_coverage, read_group)`. With one read group the row count must be unchanged — that
  is what "free at one read group" has to mean. *Depends:* A5. *Source:* spec §6.
- ☐ **B2 — `finalise` returns `SampleLocusObservations`.** Rows are re-derived from `folded_reads`
  (per-read), not from per-bucket totals — coverage and group are per-read facts. **Sort by
  `(bases, read_coverage, read_group)` before emitting**: `folded_reads` is an `AHashMap` with a
  per-process seed, so fold order would make the output non-deterministic run to run. Chain ids are
  dropped **per read** — none for a read that agreed with the reference across everything it
  witnessed — because row splitting means "the REF row" is no longer a unique row. *Depends:* B1.
  *Source:* spec §3, §6.
- ☐ **B3 — the per-record counters.** `reads_without_observation` (A5's set) and
  `reads_discarded_by_cap`. The cap truncates in the walk, **before any record exists**, so the
  truncated ids must be plumbed into the fold and registered per affected record; the count is those
  absent from `folded_reads` at `finalise`. Two cases have no clean answer and are recorded as such:
  a read truncated where no record is open, and a truncated read whose deletion would have widened a
  record. *Depends:* B2. *Source:* spec §6.

> **Checkpoint B: ng emits its own locus type, deterministically.** Pause for review.

---

## Milestone C — the generator

- ☐ **C1 — `PileupGenerator` + config + counts.** The struct, `PileupGeneratorConfig` (five knobs,
  defaulting to production's `pub const`s **by name**), `PileupGeneratorCounts` (production's
  `RunSummary` fields plus `reads_silent_over_footprint` and `records_outside_region`). *Depends:*
  B3. *Source:* arch §1.1.
- ☐ **C2 — the region walk. Own commit, do not bundle.** `begin_segment` records the region and
  opens nothing (it cannot fail; opening a query can, so the first `next_locus` is where an
  `IngestError` surfaces). The query is `[region.start, region.end + max_record_span]` — the halo,
  without which a record whose footprint crosses the boundary silently loses the support lying
  beyond it. The walk **stops** once `walker_pos > region.end` and no open record is anchored at or
  before `region.end`, or the halo is walked in full at every boundary. Records are dropped on their
  **anchor**. *Depends:* C1. *Source:* spec §2.
- ☐ **C3 — the allocator across segments.** One `ChainIdAllocator` on the generator, `reset()` at
  each region end — it preserves `next_id` and clears `pending_mates`/`active_count`, without which
  a pending mate cross-pairs between contigs and `active_count` leaks toward
  `ActiveReadsExhausted`. **Snapshot its counters at `begin_segment` and fold the delta**: `reset()`
  preserves them and `RunSummary` takes them by *assignment*, so summing per-region summaries would
  triangular-sum `chain_allocations` and `mate_lookup_evictions`. *Depends:* C2. *Source:* spec §8.
- ☐ **C4 — wire it in.** `impl LocusGenerator<()> for PileupGenerator`; fill `GeneratorSet`'s
  `generic` slot, replacing `Unfilled(NotImplemented)`. `LocusGenerationError` gains a `Walker`
  variant. *Depends:* C3. *Source:* arch §2.1, spec §7.

> **Checkpoint C: ng mints generic loci end to end.** Pause for review.

---

## Milestone D — prove it, then measure it

- ☐ **D1 — stage-2 parity.** `project(PileupRecord) -> SampleLocusObservations` in the test module;
  every divergence falls in one of the **five** named classes and is counted, not excused. Build the
  **permanent** anchor here too — loci where every folded read witnessed the whole footprint must
  agree with production forever; that is what replaces the stage-1 differential this plan retires.
  *Depends:* C4. *Source:* spec §3.
- ☐ **D2 — the dump tool and the fixtures that must fail a wrong implementation.**
  `examples/ng_generic_loci_dump.rs`, following `ng_ssr_loci_dump.rs`. Six new fixtures, because
  **no inherited test exercises the defect**: a read adaptor-masked over part of a multi-base record
  (must be `Observed`, not `Complete`); an interior `N` and a ref-skip (no observation, counted); a
  record widened past an already-expired read (bases must not have grown); an indel whose anchor base
  was masked (the residual, pinned at one base); two read groups on one allele (two rows summing to
  the single-group total); and a deletion at a region boundary (support must match a single-region
  walk). *Depends:* D1. *Source:* spec §12, §13.
- ☐ **D3 — the measurements.** Two numbers, both deliverables rather than by-products. **The size of
  production's defect:** how many loci, reads and reference bases production credits to reads that
  never sequenced them — the number that turns the indel-deficit hypothesis into a result or kills
  it. **Throughput:** the dump over one human chromosome against production's `pileup` on the same
  CRAM. There is no candidate-only fallback, so a bad number is a performance problem to solve, not a
  design to reconsider — and the allele-list growth (spec §7) is the first place to look.
  *Depends:* D2. *Source:* spec §7, §13.

> **Checkpoint D: the generic path mints loci, the defect is measured, and the cost is known.**
> Pause for review — and for the decision whether the read-group grain stays per-`@RG`.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | the stage-1 differential still green **on reads whose events tile the footprint**; failing elsewhere by design, with each failure enumerated |
| B | one-read-group row counts unchanged; output byte-identical across repeated runs (the sort) |
| C | loci from two adjacent regions concatenate coordinate-sorted with no gap at the join; `records_outside_region` zero-sum across neighbours |
| D1 | every divergence in one of five named classes; the complete-reads differential green as a permanent anchor |
| D2 | six fixtures, each written so a span-derived or fill-preserving implementation fails it |
| D3 | the defect size reported; throughput measured against production's `pileup` |

## Out of scope (next work)

- **Step 7's use of partial observations** — freebayes' `1/k` prefix/suffix scheme, decided but not
  built here (spec §10).
- **The cohort merge**, which is what consumes these loci; **the windowed statistics**, which slide
  over their per-position depth; **parallelism**, deferred whole.
- **A second generic generator** — the active-region or haplotype-window definition of a locus, which
  the `LocusGenerator` trait exists to let sit beside this one.
