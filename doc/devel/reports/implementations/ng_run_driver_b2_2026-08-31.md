# ng direct mode, step B2 — a source yields what the walk yields

**Date:** 2026-08-31. **Branch:** `main`. **Plan:**
[`../../ng/impl_plan/run_driver_direct_mode.md`](../../ng/impl_plan/run_driver_direct_mode.md)
step B2. **Spec:** [`../../ng/spec/run_streaming.md`](../../ng/spec/run_streaming.md) §12, §3.4.
**Architecture:** [`../../ng/arch/run_streaming.md`](../../ng/arch/run_streaming.md) §2, §9.
**Modules:** `src/ng/run/walker.rs` (tests only), `doc/devel/ng/arch/run_streaming.md`.

---

## What landed

The differential B1 could not run. B1's tests prove the **adapter** — the ordering, how far the
walk reached, the failure, the segment stream — and every one of them drives a scripted generator
that never touches the reads it was handed. B2 runs the **real** generic locus generator over a
real indexed BAM, twice: once through `SampleLocusObservationsIterator` driven directly, which is
the machinery that existed before this step, and once through `AlignmentFilesWalker` behind the
merge's trait. Then it compares the two, observation for observation.

**The fixture, and what it is shaped to catch.** Three thirty-base reads over a two-contig all-`A`
reference: one on `chr1` matching everywhere, one on `chr1` carrying a `C` at two positions so a
locus has both a matching and a non-matching witness, and one on `chr2`. Five segments: two
generic stretches on `chr1` separated by a satellite no generator handles, a generic stretch on
`chr2`, and one more on `chr2` no read reaches. That is **62 loci**, 3 of them carrying a
non-reference base, 4 regions handled and 1 refused permanently — all asserted, not described.

**Thirty bases is a filter, not a taste.** The first draft used ten-base reads, which
`DEFAULT_MIN_READ_LENGTH` (30) drops at admission, so the generator saw no reads, every walk came
back empty, and both comparisons passed on two empty vectors. The asserted locus count is what
says so; an `is_empty` guard would have caught it too, and the count catches more.

## §12's fourth oracle is not literally true, and the divergence is the owner's

The oracle asks that a segment walked alone emit **exactly** the observations the same span emits
inside a whole-genome walk. Built, and measured: everything is equal — bases, witness, read group,
support counts, quality sums, `placed_left` — **except the chain ids**. The `chr2` read is id 0
walked alone and id 4 walked fourth.

That is not a defect. `SequenceObservation::chain_ids` says so in its own documentation — *"an id
names a read within one walk"* — and the allocator counts up across a whole walk and survives the
per-chromosome reset, so **no implementation of a walk-scoped id can satisfy "exactly"**.

The test compares the **grouping** instead: every id replaced by the order of its first appearance.
That still catches a read split in two, two reads merged into one, or a locus that lost its
witnesses. **The spec sentence is not edited** — it is the owner's, and two things wait on them:
the wording of oracle 4, and whether the same exemption is owed to oracle 1, which asks for
byte-identical psps across worker counts and inherits the same problem once chain ids reach a file.

## Verification

| check | result |
|---|---|
| `cargo test --lib ng::run` | 349 passed (343 after B1) |
| `cargo test --lib` | 5,789 passed, 13 ignored |
| `cargo fmt --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean (exit 0) |
| `cargo doc --no-deps` | 26 unresolved links, 23 redundant link targets — the standing baseline |

## What the reviews changed

Two reviews ran in parallel over the step's diff, each in its own worktree: one mutation-testing
the differential, one on design fidelity and prose.

### The differential compared the walker with itself

**The segment-independence test's whole-walk arm came from the walker, not from the iterator**, so
both its sides carried any defect the walker had. Measured: three mutations that mangled every
yielded observation — cleared chain ids, zeroed support counts, blanked reference bases — passed
that test while failing the one above it. Its whole-walk arm now comes from the iterator, which
costs nothing and makes it a differential.

### It was also walking the one segment where nothing is carried over

The test used the `chr2` segment, on the reasoning that arriving fourth differs from starting
fresh. **It does not**: `PileupGenerator::enter_chromosome` retires the old chromosome and mints a
fresh cursor, reference window and per-chromosome walker at every contig change, so the two states
are nearly identical — and the one difference, the chain-id counter, is what the comparison
normalises away. The case with genuinely carried state is a **later segment on a contig the walk is
already inside**: there the cursor has advanced past the earlier segments and the reference window
has released the bases behind it. Both cases are tested now.

### Two mechanism claims were wrong

- **"Seventeen positions were lost to exactly this."** They were not. Spec §4.3 records them lost
  by **cutting a segment** — a cut landed 74 bases inside a 91-base deletion and the part past the
  cut was emitted by no segment. Nothing in this test cuts a segment, so a failure here means
  something else entirely, and the sentence would have sent the next reader hunting a symptom that
  cannot occur.
- **The reference comparison was open-coded.** `SequenceObservation::matches_reference` exists to
  stop exactly that, and the open-coded version is wrong in general: a partial observation's bases
  stop where its read's witness stopped, so comparing them against the whole locus's reference
  bases reports a read that matched everything it saw as non-reference. It gave the right answer
  here only because every locus is one base and every observation complete. Now goes through
  `non_reference_reads()`.

### And a number that was right for the wrong reason

`reads_admitted` is **5** on this fixture, not 3: the field counts admissions, not distinct reads,
and a read is admitted once per segment it is met in, so the two `chr1` reads count twice each.
Asserting 3 would have been asserting a count the field does not keep.

### Smaller things

The `#[expect(clippy::arc_with_non_send_sync)]` copied from the file-backed probes was unfulfilled
and failed `-D warnings`: `InMemoryRefSeq` is `Send + Sync`, so the lint never fires. Two
constants shared a name and a doc comment claimed they were one value; two comments pointed at an
"empty-walk guard" that does not exist; the fixture comment said a read disagreed with the
reference when it matches everywhere. All corrected.

Two tests were added on review advice: an analysed-but-empty generic segment (so `regions_handled`
is 4 of 5 rather than 3 of 4), and the first test of `generators()`, which nothing called.

**B1's own counts test was thinner than it looked**, and one mutation showed it: its fixture had no
unhandled region, so both refusal counters were zero in every B1 test and a tally that booked one
kind of nothing to the other would have read as correct. A satellite is in that fixture now.

## What this step does not prove

- **The spare is still pinned only negatively.** A walker that stashed every offered record and
  never freed one passes every test in the file — the one survivor of both mutation passes.
  Nothing can pin the release while there is no pool to count; recorded against G1 in the plan.
- **`RunSegments` is not in the differential's ground.** The fixture hands the walker a `Vec` of
  segments directly, so the three mutations that broke `RunSegments` were caught by B1's tests
  alone. That is the split working, not a hole.
- **No repeat tract goes through it.** The STR and bundle slots are unfilled, which is the plan's
  own scope decision.

## ⚑ A placement the architecture did not anticipate

Arch §9 says the run-level oracles of spec §12 belong in `tests/`. The segment-independence oracle
landed in `walker.rs` instead: it is about one sample's walk, which is what that file owns, and an
integration test would need the same three-read BAM to say the same thing. **Recorded rather than
quietly done** — arch §9 now carries the note, and the rule is the owner's to hold or relax.

---

## Addendum — the segments are shared, not borrowed (owner's ruling, 2026-08-31)

B1's design review found that a run could not hold both halves of what the architecture says it
holds. `AlignedFilesVariantCaller` owned its `Segmentation` by value and `RunSegments<'a>` borrowed
one, so `Vec<AlignmentFilesWalker<RunSegments<'?>>>` beside that field would have been a struct
whose walkers borrow it — self-referential, and safe Rust cannot express it. Neither escape was
open: `Segmentation` is deliberately not `Clone`, and minting a walker per draw breaks the
one-source-per-sample clause of arch §2's contract.

**The owner approved the shared handle and it is applied.** `RunSegments` holds an
`Arc<Segmentation>` and an index; `AlignedFilesVariantCaller` holds the same handle and hands one
out through `shared_segmentation`; `AlignmentFilesWalker::over` takes it. The genome-sized list is
stored once however many samples read it — 63 reference counts, not 63 copies of 100,171 segments —
and the walker type now carries **no lifetime**, which is the property that lets a run store it.

**`Arc` and not `Rc`**, though a walker is `!Send` today. What makes it `!Send` is the generator
set's `Box<dyn LocusGenerator<S>>` with no auto-trait bound, one layer down; an `Rc` here would add
a second blocker that would have to be found and removed again if that one is ever lifted.

**The test is that it compiles.** `a_run_can_hold_its_walkers_beside_the_segmentation_they_read`
builds the shape C1 needs — a struct holding an `Arc<Segmentation>` and a `Vec` of walkers over
it — so a change back to a borrow fails at the compiler here rather than three steps later at the
wiring. It also asserts the sharing is real: three walkers over one segmentation leave a strong
count of four, where three copies would leave four allocations.

`ng::run` at 350, full lib at 5,790; fmt and clippy clean by exit code, `cargo doc` unmoved.
