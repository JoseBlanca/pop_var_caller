# ng direct mode, Milestone C — the merge, fed by walkers

**Date:** 2026-08-31. **Branch:** `main`. **Plan:**
[`../../ng/impl_plan/run_driver_direct_mode.md`](../../ng/impl_plan/run_driver_direct_mode.md)
steps C1 and C2. **Spec:** [`../../ng/spec/run_streaming.md`](../../ng/spec/run_streaming.md)
§3.2, §5.1, §8; [`../../ng/spec/cohort_merge.md`](../../ng/spec/cohort_merge.md).
**Architecture:** [`../../ng/arch/run_streaming.md`](../../ng/arch/run_streaming.md) §2, §3.4, §5.
**Modules:** `src/ng/run/callers.rs`, `src/ng/run/walker.rs`, `src/ng/run/mod.rs`,
`src/ng/run/cohort_merge/mod.rs`, `src/ng/read/input/test_fixtures.rs`,
`examples/ng_open_cohort_descriptors.rs`.

**C1 and C2 landed in one commit**, named. C2 adds no production code — it is the differential
that says C1 changed nothing but the sources — and splitting them would have put an unproved
claim in the tree for one commit.

---

## What landed

`AlignedFilesVariantCaller::merge_cohort` drives the single-threaded merge over one walker per
sample. Until now the merge had only ever been fed vectors: its own fixtures hand it observations
in memory, and so do both of its oracles.

Two supporting pieces, both in `walker.rs`:

- **`WalkReference`** — the reference the walk fetches bases from, with its `.fai` parsed once for
  the run. Every accessor built from it is its own, because `WindowedRefSeq` holds an open
  per-contig reader and is `Send` and deliberately not `Sync`; what is shared is the index and the
  contig table.
- **`generic_path_generators`** — the generic slot filled, both repeat-tract slots refused as
  **unbuilt**. A run over ground with tracts in it is not wrong, it is short, and the tally says
  by how much. That distinction now has a test.

And a fourth construction refusal: a reference read from a `.fai` alone holds no bases, so the
walk has nothing to fetch a reference allele from. Refused before a single alignment file opens,
under the standing rule that a run wrong at the door should be told so at the door. **It refuses
nothing a real run does** — every arm of `read_reference_verifying_or_creating_fai` keeps the
FASTA's path beside the geometry it read from the index.

## The fixture the run's tests were built on modelled a shape no run holds

`fixture_reference(false)` reads a bare `ReferenceSource::Fai`, which carries **no `fasta_path`**.
Nothing in a run ever holds one: the batteries-included read keeps the FASTA's path beside the
index's geometry, verifying on a background thread. So A1's and A2's fixtures had a reference that
could be checked against a cohort's headers and never read from.

`fixture_reference_from_its_index` is the real shape — geometry from the index, **no digests until
the background read is joined**, bases reachable throughout — and the run's tests use it.
`fixture_reference(false)` stays right for tests about a reference that genuinely has nothing
behind it, which is now what the no-bases refusal's own test wants.

## ⚑ The descriptor refusal was wrong in the unsafe direction, and it was this step that made it so

Checkpoint A measured 2 descriptors an alignment file — one for the file's reader, one for the
reference accessor its cursor holds — and the refusal has budgeted that since A2.

**A locus generator holds two more, per sample.** One accessor for the walk's own REF fetches and
one for the read preparer; each opens a reader on the FASTA at its first fetch and holds it for
the run. Re-measured with `examples/ng_open_cohort_descriptors.rs`, extended for this, on the same
63 tomato accessions:

| point | descriptors |
|---|---|
| reference open, no sample open | 3 |
| every sample's files open, nothing decoded | 4 |
| one cursor per sample, walked over one region | 130 |
| and the two accessors a generator holds per sample | **256** |

So a walking run holds **253 descriptors for 63 files over 63 samples**, where the refusal budgeted
158. **A run could pass the check and then die at `EMFILE`** — the exact failure the check exists
to prevent. The arithmetic now has two terms, 2 a file and 2 a sample, and the message shows both:

```
this run needs 36 open files and this process may open 4: 1 alignment files at 2 each,
1 samples at 2 more each for the reference bases their walks read, and 32 for the reference,
the repeat catalog and the output. Raise the limit with: ulimit -n 36, …
```

The test that pinned the old shape was named `the_descriptor_count_grows_with_files_not_with_samples`.
It now grows with both, and the test says so.

## Verification

| check | result |
|---|---|
| `cargo test --lib ng::run` | 361 passed (350 before this milestone) |
| `cargo test --lib` | 5,801 passed, 13 ignored |
| `cargo fmt --check` | clean (exit 0) |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean (exit 0) |
| `cargo doc --no-deps` | 26 unresolved links, 23 redundant link targets — the standing baseline |

**C2's differential** compares the walker-fed merge against the same observations fed from memory,
and against `merge_cohort_serially` — the merge's undivided reference implementation, the one its
parallel driver is checked against. Both go through the merge's own `render`, which destructures
`RegionOutcome` so that a field it gains has to be answered for rather than dropping silently out
of the comparison. The oracle opens its own `SampleReads` and drives the iterator directly, so
neither `AlignmentFilesWalker` nor `RunSegments` is on both sides.

**What the oracle does share, and cannot not share, is `WalkReference` and
`generic_path_generators`** — so a defect in either moves both answers together. What pins those is
the positive assertion on two named positions, not the differential.

## What the reviews changed

Two reviews ran in parallel over the milestone's diff, each in its own worktree: one mutation
testing, one on design fidelity, the artefact and the prose.

### Six mutations survived, and five now do not

The correctness review wrote 13 mutations and 6 lived. Tests were added for five of them:

- **The run's merge parameters were not used by anything.** All three — the min-alt floor, the span
  bound, the building-region width — could be replaced by unrelated constants with every test
  green. A user's threshold could have been silently ignored. Now: one sample carrying its variant
  on **one** of three reads is below the shipped floor of two and above a floor of one, so the same
  reads give one locus at one setting and two at the other. Re-injected and caught.
- **The no-bases refusal was not pinned to its position.** Moving it after every file opens left
  everything green, while the test named for the ordering still passed. Now checked the way its
  three siblings are — a bases-less reference *and* an unindexed BAM, asserting the cheap refusal
  wins. Re-injected and caught.
- **The two kinds of refusal were interchangeable.** Swapping the repeat-tract slots from *not
  implemented* to *out of scope* was invisible, and it is the difference between a run report
  saying *this caller has not built that yet* and *that ground will never be called*. Re-injected
  and caught.
- **Sample order and the silent sample** — see the next section.

**One survivor is not a defect and the test now says so.** Handing the merge the segments instead
of the analysed regions changes no output, and it should not: the merge is built so that how the
ground is divided cannot change what comes out, and its two drivers are checked against each other
for exactly that. What the swap costs is work — the same 20,000 observations take 5.4 ms over one
region and 184 ms over a thousand, 34 times for the same answer — and no assertion on the output
can see it. Recorded rather than papered over with a test that would not have worked.

### A test asserted an invariant the merge does not hold

`every_sample_is_represented_at_every_cohort_locus` asserted one `per_sample` row per sample of the
run. **Measured by the reviewer: a sample with no observations over a locus gets no row at all** —
`SampleMembers` says so in its own documentation, and identity is carried by `SampleSupport::sample`
rather than by position. The test passed only because all three fixture samples happen to cover
both positions, and it would have broken on the first fixture where one did not, while claiming the
merge had changed. It now asserts what is true and what it was written to check: the silent sample
appears under its own run index, with evidence, all of it the reference allele.

### Three claims in the prose were wrong

- **"Three accessors are built per sample"** — two are *held* per sample, and the third is a
  factory called once per file per chromosome.
- **The microsecond split was a misquote.** The measurement says 189 µs is the whole per-accessor
  cost, sharing the index brings it to 52 (34 of which is cloning the contig table), and sharing
  both leaves ~18. The text had added 188 and 34 into a total that appears in no measurement — and
  omitted the mechanism that makes sharing worth anything, which is that the parse is paid at every
  contig open rather than once per accessor.
- **Two test justifications were inverted**, each claiming the test ruled out something a sibling
  test already ruled out, or something it could not see.

### And one recommendation could not be taken

The design review asked for `merge_cohort` to be narrowed to `pub(crate)`, since arch §6 says
everything but the three iterators is crate-private. **It fails the build**: with no consumer
outside tests, dead-code analysis rejects it under `-D warnings` — which is exactly what
`callers.rs`'s own note has said since A1. What did narrow is everything reachable only through it,
`WalkReference` and `generic_path_generators`. The note now records that the block is mechanical
rather than a preference, and names the consumer that will lift it.

## ⚑ Two things wait on the owner

**`merge_cohort` destroys what a run report will need.** The walkers go into the observation cache,
the cache owns them and has no `into_sources`, so when the merge returns the run has dropped every
walker's locus tally, its generators' per-slot counts — the nine that explain a covered region
emitting nothing — and its cursors' read-filter tallies, along with the assembly-check outcome the
run computed at construction. F3 needs all four. The fix is a signature the next step is already
editing, so it is worth shaping once rather than twice.

**A run cannot set its locus generator's settings.** `PileupGeneratorConfig::default()` is
hard-coded, fixing all five knobs at production's constants, while every other knob — the calling
loop, candidate selection, the merge, the read filters — is an explicit argument to `open`. Two of
the five are the depth axis, and `max_active_reads`'s own documentation records the old value
silently refusing 19,725 reads on one ~130× tomato chromosome. The design review's recommendation
is `AlignmentInputs`, beside `read_filters`, since both answer how a run turns bytes into evidence.
`RunError::LocusGeneratorSettings` exists and cannot fire until this is settled.

## What this milestone does not do

- **It does not yield one record at a time.** `merge_cohort` consumes the run and returns every
  cohort locus at once, where spec §5.1 bounds direct mode at `callers in flight × one cohort
  locus` plus the frontier. That is not a cost this step added — `merge_cohort_through_cache`
  already accumulates one outcome per building region — and calling lands *inside* that driver, so
  the shape survives the next step. The pool milestone is where loci start being released singly.
- **No repeat tract goes through it.** Both tract slots are unfilled, which is the plan's own scope
  decision, and the run's tally reports them as this caller's gap rather than as ground refused.
