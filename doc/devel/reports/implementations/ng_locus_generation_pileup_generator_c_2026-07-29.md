# ng generic locus generator — the generator, Milestone C: the region walk

**Date:** 2026-07-29 · **Plan:**
[locus_generation_pileup_generator.md](../../ng/impl_plan/locus_generation_pileup_generator.md)
steps C1–C4 · **Spec:** [locus_generation_pileup.md](../../ng/spec/locus_generation_pileup.md)
§2, §7, §8 · **Arch:** [locus_generation_pileup.md](../../ng/arch/locus_generation_pileup.md)
§1.1, §2.1, §2.2

Implementation report for Milestone C of plan 3 of 3. Four commits, one per step.

## 1. Plan

Wrap the walk in a generator: its knobs and counts (C1), the region walk with its
query halo and its stop rule (C2), one chain-id allocator across regions with its
counters folded as deltas (C3), and the dispatcher wiring (C4).

## 2. Assumptions

Three, all recorded in the commits rather than discovered later:

1. **The read-query reference is a caller-supplied factory, not the walk's
   accessor.** `SampleReads::reads_in_region` gives each of a sample's k files its
   own raw accessor on purpose — they are stateful readers and sharing one makes k
   streams share a file cursor. Arch §1.1's field list predates that signature, so
   the generator takes a factory *beside* its own accessor rather than deriving
   one from it.
2. **A read the preparer declines (`Ok(None)`) is skipped and not tallied.** No v1
   preparer returns it.
3. **`begin_segment` is the right place to end an abandoned walk.** The trait says
   it cannot fail, and ending a walk cannot.

## 3. Changes made

### C1 — the generator's state, and the one knob that is not production's (`a5c0203`)

`PileupGeneratorConfig` (five knobs, four of them production's `pub const`s **by
name**), `PileupGeneratorCounts` (seven `RunSummary` fields plus ng's two), and
the struct they hang off.

**The ceiling is the substance.** A `ReadCoverage` run is two `u16`s minted
through `LocusLen::from_positions`, which **saturates**; production's
`--max-record-span` is an unbounded `u32` with a default of 5,000. Inheriting the
knob by name would inherit a silent truncation — a partially-witnessed record
wider than 65,535 reporting a *smaller* `positions_covered` and no error, at
exactly the long-deletion loci this generator exists to get right. `check()`
rejects it at construction; A4's `debug_assert` in `coverage_of` stays as the
statement of the invariant.

`a_span_one_past_the_ceiling_is_where_a_coverage_run_starts_lying` asserts the
reason rather than describing it. Production's default staying inside the ceiling
is a **`const` assertion**, not a test: it compares two constants, so a runtime
assert on it is a test that cannot fail.

### C2 — the region walk, with a halo that is queried and stopped (`d23e2b0`)

`begin_segment` records the region and opens nothing; the first `next_locus`
opens the query and builds the walk; records are dropped on their **anchor** and
tallied.

**The halo.** The query is `[region.start, region.end + max_record_span]`. Without
it a record anchored inside the region is emitted — by the right region, with
every counter reading zero — carrying part of its evidence.

**The stop.** `PileupWalker::stopping_after(pos)`: stop once `walker_pos > pos`
**and** nothing anchored at or before `pos` is still open, flushing exactly as
end-of-input does. `run` is untouched, so the stage-1 differential still compares
like with like.

The read stream is ng's own `PreparedRegionReads`: it prepares each read and
**sheds** fatal errors into a shared cell, because the walker consumes an
infallible `Iterator<Item = PreparedRead>`. Production hit the same seam and
solved it the same way (`ErrorSheddingAdapter`); ng keeps its own because it sheds
two error types into one and because a locus generator has no business importing
from the CLI.

### C3 — one allocator across regions (`c89d596`)

The allocator is lent to each walk (`adopting_chain_ids`) and taken back
(`into_chain_ids`), then `reset()`. `fold_region_walk` folds the walk's summary:
per-walk fields add, `active_reads_high_water` takes the max, and the two
allocator-derived counts are **deltas** against the allocator's counters as they
stood when the walk took it — because `reset()` preserves them and `summary()`
assigns them, so summing region summaries triangular-sums both.

`end_walk` is the only place the walk is cleared, and is reached from all three
ends: the walk draining, `begin_segment` abandoning it, and a fatal walker error.

### C4 — the generic slot is filled (`52d99c5`)

`impl LocusGenerator<()> for PileupGenerator`, delegating to the inherent methods.
`GeneratorSet` needed no change — `new` already takes the slots — so the
end-to-end test fills the slot itself and drives the **public**
`SampleLocusObservationsIterator`.

It also answers the FIXME that module left for this step: the iterator's field
order now drops `generators` before `reads`, because the generic generator is the
first to hold a region stream across `next_locus` calls.

## 4. Validation

| | |
|---|---|
| `cargo fmt --check` / `clippy --all-targets --all-features -D warnings` | clean |
| `cargo test --lib` | **2709 passed** (2684 at the start of the milestone) |
| `cargo test --lib ng::locus_generation::pileup` | 188 passed, 1 ignored |
| `cargo doc --no-deps` | 12 unresolved links, all pre-existing, none in this module |
| soak, host-native, `--profile soak`, 5,000 cases | 10 passed — the fold is undisturbed |

`cargo test --all-targets --all-features` is not run: it is known-red for an
unrelated reason (`benches/psp_writer_perf.rs:386` panics).

## 5. Deviations from the plan

1. **`LocusGenerationError::Walker` landed in C2, not C4** — C2 is the first step
   that can produce a walker error.
2. **`records_emitted` is not mirrored** into `PileupGeneratorCounts`, against arch
   §1.1's "the first seven mirror production's field for field" (it is seven there
   and eight in `RunSummary`). Loci *kept* is emissions minus
   `records_outside_region`, and the kept count is already
   `LocusCounts::loci_emitted`.
3. **`reference` is `Arc<R>` and `R: RefSeq` tightened to `RawRefSeq`** — the walk
   must own a `RefSeq` and must not rebuild the accessor per segment, and one type
   parameter then serves both the walk and the query factory.
4. **`preparer` + `prep_scratch` became one `Rc<RefCell<ReadPreparation>>`** — the
   stream cannot borrow a sibling field, and the shared handle is what keeps the
   scratch off the per-region stream.
5. **A read-preparation failure maps to `LocusGenerationError::Reference`**,
   matched exhaustively so a new `ReadPrepError` variant is a compile error.
6. **C3's baseline is snapshotted when the walk takes the allocator**, not at
   `begin_segment`, which now opens nothing.

## 6. The test this milestone shipped that could not fail

C2's first stop-rule test asserted that a record anchored at 99 still spans to
139 — but **a record's footprint is fixed when it opens**, so an early flush never
shortens it, and the assertion held with the second half of the stop rule
deleted. Only the mutation said so. It was rewritten around the one thing an
early stop does destroy: a **widen** anchored past the boundary, from a read that
a plain query already returns, so the halo is not the variable.

That is the seventh instance of this pattern on the branch and the second in a
test written in the same session that caught it.

## 7. Review, and what it changed

Five category agents, each in its own git worktree
([report](../reviews/ng_locus_generation_pileup_generator_c_2026-07-29.md)):
**2 Blockers, 6 Majors**, fixes in `94758d7`. Both Blockers were found
independently by more than one agent.

**Both Blockers were about the generator's *ending*, not its walking.** A shed
stream error outlived the region that shed it — reported against the next,
healthy region after it had emitted every one of its loci, or never reported at
all. And neither `Ok(None)` nor `Err` was terminal: the next call re-opened the
query and re-walked, admitting every read twice and handing one fragment two
chain ids, which is the corruption the run-lifetime allocator exists to prevent
reached from the other direction. `GeneratorSet` happened to shield both, so
they were latent through the dispatcher and live through the generator's own
public API — the API all of the module's tests drive.

**Five properties were load-bearing and pinned by nothing:** the halo's width (a
half-width halo passed every test), the stop rule's comparison (`>=` drops the
locus anchored exactly at `region.end` — a one-base hole at every boundary), two
of the fold's three rules, `end_walk` on the error path, and spec §7's lazy read
stream — an agent collected the query into a `Vec` and **the whole library stayed
green, parity included**.

**The cost, measured:** `T(k) ≈ 81.5 ms + 0.12 ms × k` over a 200 kb contig at
depth 3, flat from k=10 to k=500. A region costs about **290 loci of walking** to
set up, so regions under ~300 bp cost more to open than to walk. Peak bytes are
flat at 257,621 whether the region is 400 bp or 200 kb — the depth-shaped
footprint, measured rather than asserted. **The stop rule is worth 8.7× total
wall at 500 regions while changing not one emitted locus.**

**Two more tests that could not fail** (§6 and the review's own replacement for
it), bringing the branch's count to nine.

## 8. Checkpoint C — what is open

None of these blocks Milestone D.

1. ~~No generator error carries the region it fired in.~~ **Done — owner
   approved, 2026-07-29 (`a17f37b`).** Every generator-raised variant now carries
   a `GenomeRegion`, `Reads` split into `OpenReadQuery` + `Reads`, and **none of
   the four has a `#[from]`** — with a blanket conversion a bare `?` compiles and
   silently produces an error with no region, which is the state the change
   exists to make unreachable. See §9.
2. **`PileupGeneratorCounts` is unreachable once the generator is boxed into a
   `GeneratorSlot`** — nine counters with no reader that is not a test, and it is
   what makes a zero-knob walk *totally* silent. Needs either a counts method on
   the `LocusGenerator` trait, a fold into `LocusCounts`, or a caller-held sink.
3. **The read-query accessor factory is called once per file per region** — spec
   §8's ~564k-opens trap shape, in the one accessor the generator cannot hold for
   the run because `reads_in_region` requires a per-file one. Documented on the
   field; no non-test caller exists yet, so it is an invitation rather than a
   measured cost. D2's dump tool is the first caller that decides it.
4. **Nothing logs.** No `tracing` anywhere in the generator: the clamp, a
   declined read and a shed error are all counted or silent. Production's own
   allocator warns, so the convention exists in the copied code.
5. **Two counters are defined and never incremented** —
   `reads_silent_over_footprint` (needs a per-active-read "ever contributed" flag
   in the walk) and reads a preparer declines (no field at all). Spec §13's
   read-accounting assertion in D2 is what forces both.
6. **The iterator's field order is load-bearing with no test** (C4). The
   reviewer found and tested a cheap enforcement — an explicit `Drop` that
   replaces `generators` — which was left unapplied because the drop-order
   guarantee is a language rule, not a runtime one.
7. **Nine of the ten names `pileup/mod.rs` re-exports `pub` have no consumer**
   anywhere in the tree; only `WalkerError` is used. The module reasons
   carefully about `pub(crate)` for its *vocabulary* re-exports and never applied
   the same test to its own surface.
8. **The arch inventory still places `PileupGenerator` in `mod.rs`** and lists
   neither `generator.rs` nor `mock_reference.rs`.
9. **Checkpoint A's and B's remaining items**, and plan 2's four, stand.

## 9. The region on the error (owner, 2026-07-29)

Checkpoint C's first open item, decided and applied.

**Shape.** Four variants carry a `GenomeRegion` —`OpenReadQuery`, `Reads`,
`Reference`, `Walker` — and `TypedRegion` deliberately does not, the region
*stream* being what failed there. A `region()` accessor answers the question
without matching, so a consumer does not rot when a variant is added. Splitting
`Reads` into open-vs-mid-stream came in the same pass, as the review asked: "the
index query for this region could not be opened" and "the record stream broke
40 kb in" are different operational problems that rendered identically.

**No `#[from]` on any of the four, and that is the enforcement.** A blanket
conversion means a bare `?` compiles and yields an error with no region — the
state this change exists to make unreachable. Removing it turned every affected
`?` into a compile error, which is how the eight attachment sites were found
rather than remembered: four in the pileup generator, four in the STR generator,
which returns through the same shared type.

**Which region gets attached:** the one the generator was working over — the
segment, not the halo-widened span it queried, because the segment is the unit a
caller can act on. The one exception says so where it attaches: the STR read
fetch is a free function handed a tract-plus-margin span and knows no other.

`GenomeRegion` gained a `Display` (`contig 3:940-1100`). It says *contig 3* and
not *chr4* because the type holds an id, not a name — rendering one as the other
is a lossy translation this branch has already deleted once (A0).

### It caught a test of mine that could not fail — the tenth on this branch

`an_abandoned_regions_shed_error_is_not_charged_to_the_next_region`, written for
the review's own Blocker in `94758d7`, called `begin_segment` twice in a row —
and **`begin_segment` opens nothing**. So the "abandoned" region never ran, never
shed anything, and the error the test caught came from the *second* region's own
read. It failed when the fix was reverted, which is why it looked sound; what it
actually pinned was "a shed error is reported at all", which another test already
covers.

The region assertion is what exposed it: the error named contig 0 where the test
said contig 1. Rebuilt so the failing read is one the walker reaches **while an
earlier read is still producing loci** — the only shape in which a region can be
abandoned with its failure latched and unreported — it now reproduces reviewer
probe A exactly: with the fix reverted, chr1 emits loci while chr2's error waits
to be charged to it.

**Validation:** fmt clean, clippy `--all-targets --all-features -D warnings`
clean, `cargo test --lib` **2718 passed**, `cargo doc --no-deps` still 12
pre-existing unresolved links, host-native soak at 5,000 cases green.
