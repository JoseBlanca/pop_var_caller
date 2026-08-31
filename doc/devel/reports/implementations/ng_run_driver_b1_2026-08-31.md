# ng direct mode, step B1 — the walker as an observation source

**Date:** 2026-08-31. **Branch:** `main`. **Plan:**
[`../../ng/impl_plan/run_driver_direct_mode.md`](../../ng/impl_plan/run_driver_direct_mode.md)
step B1. **Spec:** [`../../ng/spec/run_streaming.md`](../../ng/spec/run_streaming.md) §3.4, §8.
**Architecture:** [`../../ng/arch/run_streaming.md`](../../ng/arch/run_streaming.md) §2, §5.
**Modules:** `src/ng/run/walker.rs` (new), `src/ng/run/mod.rs`,
`src/ng/locus_generation/mod.rs`.

---

## What landed

`AlignmentFilesWalker` — one sample's alignment files behind the cohort merge's source
interface. The merge asks each of its k samples "what did you see next?", one observation at a
time, in coordinate order, forward only, for the whole run; this is direct mode's answer to that
question, and the psp reader will be the other.

**It is a wrapper, and deliberately a thin one.** `SampleLocusObservationsIterator`
([`locus_generation/mod.rs:921`](../../../../src/ng/locus_generation/mod.rs)) already drives the
typed-region generators over one sample's reads and yields loci in genome order. Everything the
walker adds is one of the two things the merge's trait asks for that a plain iterator cannot
give.

- **A failure that names the sample and where the walk had reached** —
  `RunError::SourceFailed`. The merge adds nothing to a source's error and passes it through, so
  the error has to locate itself. In a run over a thousand samples, *reading a region failed*
  names neither the individual to look at nor the place to look.
- **The offer of a spent record back for reuse.** `next_observation(spare)` takes it, and this
  walker drops it. That is what the trait's own text calls the contract — "the spare is an offer
  and not an obligation" — and refilling it is step G1, which is where the measurement that
  motivates it lives (92% of the merge's frees on 63 tomato accessions over 100 kb). What B1
  buys G1 is that the hook is answered by a real implementation instead of by the blanket one,
  because the blanket implementation is the thing G1 cannot change.

`RunSegments` is the region stream a run hands its walkers: a borrow of the segmentation's own
list, in genome order. Every sample of a run reads the same object, which is what makes "k
samples over one segmentation" true rather than "k lists that happen to agree".

## Three decisions the plan left open

### The error carries how far the walk got, and that is not always a position

Arch §5 sketches `SourceFailed { sample, at: GenomePosition, source }`. **The code ships
`reached: WalkProgress` instead**, an enum of `NothingYet` and `After(GenomePosition)`, because a
source that fails on its very first draw has no position and must not invent one. A run that said
"failed at contig 0:1" when nothing had been read would send an operator to an innocent locus.

The two render as adverbial phrases inside the sentence, so the message reads either way:

```
sample zeta: reading its observations failed after contig 0:13
sample zeta: reading its observations failed before its first observation
```

This is recorded as a small deviation that keeps the intent, not a design change: it adds the
case the sketched type cannot express and changes nothing else.

### The walker is not an `Iterator`, and it cannot be

Every iterator of one sample's observations is *already* a source, through the blanket
implementation at
[`observation_cache.rs:98`](../../../../src/ng/run/cohort_merge/observation_cache.rs) that drops
the spare and calls `next`. A walker that was also an `Iterator` would therefore be a source
twice over, and Rust refuses the overlap. The walk is still reachable as an iterator one level
down — it *is* `SampleLocusObservationsIterator` — which is exactly what B2's oracle drives.

### A region stream that cannot fail says so in its type

`SampleLocusObservationsIterator` takes `Result<TypedRegion, E>` because the catalog's own reader
is fallible. A run's segments are not: they were read once at `Segmentation::build`, and every
catalog failure was reported there. So `RunSegments`'s item is
`Result<TypedRegion, Infallible>`, and `locus_generation` gains one impl —
`From<Infallible> for LocusGenerationError`, whose body is `match never {}` — to admit it. The
alternative was to reuse `RepeatCatalogError` as an error that is never raised, which would have
been a lie in the type signature.

## Verification

Every command below was run in the dev container, from the repository root.

| check | result |
|---|---|
| `cargo test --lib ng::run` | 343 passed (329 before this step, 14 new) |
| `cargo test --lib` | 5,783 passed, 13 ignored (5,769 before) |
| `cargo fmt --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo doc --no-deps` | 26 unresolved links, 23 redundant link targets — the standing baseline, unchanged |

**The fixtures are shaped against the blind spot this step has.** A walk over one segment of one
sample cannot tell *advances to the next segment* from *hands back everything it has*, and it
cannot see the spare at all. So:

- the walk fixture is **three segments emitting 2, 0 and 3 loci**, the middle one empty and the
  third on a second contig — a walker that began only its first segment yields two, and one that
  ignored segment boundaries cannot place its loci on the right contigs;
- every locus spans **two** bases and starts one base later than the last, so *the end of the
  last observation* is a different number from its start, from the first observation's end, and
  from the segment's own end — the three things `reached` could have been by accident;
- the failure fixture walks **two samples** whose walks fail at different places, so the sample
  named in the error is a claim rather than a coincidence.

## What the reviews changed

Three reviews ran in parallel, each in its own worktree over the step's diff: correctness with
mutation testing, design fidelity against the spec and architecture, and one whose only job is the
artefact a person sees.

### The correctness review: 21 mutations, 20 caught

Every targeted defect in `next_observation` (dropping the spare, not updating `reached`, taking
`region.start` for `region.end`, pinning a contig, setting `reached` on the error path, swallowing
a failure as end-of-stream), in `WalkProgress`'s `Display` (both arms, the field order), in
`RunSegments` (skip the first, stop at a contig change, hand one out twice, hand them out
backwards), in the constructor's sample name, in the counts forwarding and in the `Infallible`
impl was killed by a named test or by a compile error.

**One survived, and it is the one this step cannot pin.** A walker that stashes every offered
record and never frees one passes all fourteen tests: what they assert is that the spare does not
come back out *as an observation*, not that it is released. There is no pool to count while the
walker drops, so the test belongs to the step that starts keeping records — recorded on the plan's
G1, which is where a pool first exists to bound.

**It also found a defect in the tests rather than the code.** The draining loops were `while let`,
so a walker that yielded for ever took the test binary down with signal 9 after about four minutes
and named nothing. That is exactly the shape the spare mutation produces. Every draining loop here
is now bounded and says which test is repeating itself.

### Three claims in the prose were wrong, and one of them three times over

- **`reached` open-coded a rule the crate keeps in one place.** It built the position from
  `region.end`; the canonical `SampleLocusObservations::reach_position` is `end.max(start)`, and
  its own documentation names this mistake — `GenomeRegion` has public fields and no constructor,
  so an inverted region read straight off `end` puts an observation's reach before its own first
  base. The merge's cache keys on the canonical call. Now so does this.
- **A doc heading claimed the opposite of the code beneath it.** Arch §2's *a failure leaves the
  source live* is what lets a cover be made again; the wrapped iterator latches `done`, so a
  failed walk is **spent**. The comment asserted compliance and then described the deviation. It
  now says which it is, why nothing reaches it today — the cache propagates the error without
  marking the source spent, and both drivers abandon the cache rather than retry — and what it
  would cost if something did: cohort loci built without that sample, wrong genotypes rather than
  an error.
- **The justification for printing contigs by index was false in three ways.** It said the names
  are unreachable (a run's `Segmentation` carries the catalog's contig table, names and all), that
  a lookup would put a table in every walker (one `&Segmentation` is shared by every walker by
  construction), and that `contig 0:13` is what every other position in ng prints (ng's other
  genome *positions* print `contig N position P`; `contig N:start-end` is a *region*). The index
  stays — a bare position with no reference beside it cannot do better — but the reasons are now
  the true ones, and the spelling is the established one.

### The message now says which of its two coordinates is which

A rendered chain carries two genome positions and the reader had no way to tell them apart:

```
sample TOM-042: reading its observations failed after contig 0:13: reference fetch over
contig 0:14-2013 failed during locus generation: fetch [2012, 4012) past ContigId(0) …
```

The first is how far this sample *succeeded* and the second is where it *broke*. Where the failure
falls inside the segment the last observation came from, that reads as one fact stated twice. Both
`WalkProgress` arms now name their own role:

```
sample TOM-042: reading its observations failed; its last complete observation ended at
contig 0 position 13: reference fetch over contig 0:14-2013 failed during locus generation
```

**One recommendation was not taken.** The review asked for a "what to do next" clause, which every
other refusal in `RunError` carries. `format_error_chain` appends each cause after a colon, so an
instruction on this line would land in the middle of the sentence, ahead of the thing it tells the
reader to act on — and every variant that ends with an instruction has no cause beneath it. The
advice belongs to whatever reports the run, and the variant's documentation now says so.

## ⚑ One finding is the owner's, and it lands on C1 rather than here

**`RunSegments<'a>` borrows the segmentation, and `AlignedFilesVariantCaller` already owns one.**
The caller holds `segmentation: Segmentation` by value. To hold one walker per sample for the whole
run — which is what the architecture says it holds — it would need `Vec<AlignmentFilesWalker<RunSegments<'?>>>`
borrowing its own field, and safe Rust cannot express that. The escapes are blocked or worse:
`Segmentation` is deliberately not `Clone`, and building a walker per `next()` breaks the
one-source-per-sample clause outright.

Nothing here is wrong today, and C1 hits it first. The shape that resolves it is an owned handle —
`RunSegments` holding an `Arc<Segmentation>` and an index — which keeps exactly the reason the
borrow exists (63 walkers share one list of 100,171 segments rather than copying it) and costs the
caller a field-type change. **Not taken here**, because it changes A1's ownership of the
segmentation, which is the owner's to settle; raised at Checkpoint B.

## What this step does not prove

Every test here drives a scripted generator that ignores the `SampleReads` it is handed: the
fixture BAM is opened and never read. So the suite proves the **adapter** — the ordering, how far
the walk reached, the failure, the segment stream — and not the walk. That split is the plan's:
step B2's oracle is the existing iterator driven directly, and it is where "the walker drives real
generators correctly across a contig change" is settled.

Two gaps are recorded rather than closed:

- **The read-filter drop tallies have no route out.** They belong to a cursor from the moment it is
  made (spec §8), cursors live inside the generators, and neither the walker nor `SampleReads`
  hands one out. `generators()` is forwarded so the per-generator counts survive; the tallies need
  something that does not exist yet. The run report step owns it.
- **The walker is neither `Send` nor `Sync`**, and not by its own choice:
  `GeneratorSlot::Generator` holds a `Box<dyn LocusGenerator<S>>` with no auto-trait bound, which
  that type documents as deliberate. So a walker cannot go under `merge_cohort_in_parallel`, and a
  walker and the merge drawing from it stay on one thread. Recorded in the type's documentation and
  in arch §2 so the pool milestone reads it rather than meeting it as a compiler error.
