# ng — the alignment cursor, Milestone D: implementation report

*Plan: [impl_plan/alignment_cursor.md](../../ng/impl_plan/alignment_cursor.md).
Design: [spec](../../ng/spec/alignment_cursor.md), [arch](../../ng/arch/alignment_cursor.md).
Evidence: [perf review](../reviews/perf_ng-generic-pileup_2026-07-31.md).
Branch `ng-generic-perf`.*

**Milestone D points the callers at the cursor.** Everything before it built a cursor nothing
a user runs went through: at Checkpoint C the generic generator still called
`SampleReads::reads_in_region` once per region, so the probe measured the old path whichever
way the cursor behaved. D is where that stops, and it is therefore where the end-to-end
number is owed.

---

## The step boundary: D1 and D2 landed together

**Agreed with the owner before any code was written.** The plan makes them two steps — D1 the
walker holding the cursor, D2 the generator wiring it — and they cannot be separated, for a
reason that is in the dependency graph rather than in the work.

`PileupGenerator` builds a `RegionWalk` per region: `open_walk` opens the query, `genome_walk::run`
builds the walker, `end_walk` drops it. Nothing else in `src/` builds a `PileupWalker` for
more than one region. So D1's deliverable — `PileupWalker::move_to_region` — has **no caller
until D2 replaces `open_walk`**. Committing D1 alone ships a method exercised only by its own
tests and carrying a `#[allow(dead_code)]`, which is the shape B1 shipped and a reviewer found
unreachable by putting `panic!` in the body while all thirteen of its tests passed.

The counters are the second reason, and the stronger one. The largest part of D1 is not the
plumbing — six sites touch `self.reads` — but deciding, field by field, what a per-region reset
keeps. Two of those decisions are only checkable against `PileupGeneratorCounts::fold_region_walk`,
which is D2's call site. See *The trap* below.

The same coupling has now hit B1/B2 and C1/C2+C3. D3 stays its own step: the STR generator is a
separate call site with its own dump baseline.

---

## What changed

| file | change |
|---|---|
| [genome_walk.rs](../../../../src/ng/locus_generation/pileup/genome_walk.rs) | `RegionReadSource` (the trait a long-lived walker needs of its source), `LookAhead::move_to_region`, `PileupWalker::move_to_region`, `PileupWalker::chain_id_counters`, `WalkerState::begin_region`. `stopping_after` deleted — `move_to_region` sets the bound and nothing else called it. Three tests. |
| [active_read_set.rs](../../../../src/ng/locus_generation/pileup/active_read_set.rs) | `ActiveReads::begin_region` — `reset` with `silent_exits` and `next_read_id` zeroed and the emptiness assertion dropped. |
| [generator.rs](../../../../src/ng/locus_generation/pileup/generator.rs) | `PreparedRegionReads` → `PreparedSampleReads`, holding a `SampleCursor` for a chromosome and implementing `RegionReadSource`; the halo widening moved into it. `RegionWalk` split into `ChromosomeWalk` (the walker and its contig) and `RegionWalk` (the region and the counter baseline). `open_walk` → `open_walk` + `enter_chromosome`. `MakeReference` dropped from the type parameters. Three tests. |
| [read/input/mod.rs](../../../../src/ng/read/input/mod.rs) | `IngestError::Cursor`, so a cursor failure can reach `LocusGenerationError::OpenReadQuery` and `::Reads` by the routes those already have. |
| [benches/ng_generic_pileup_perf.rs](../../../../benches/ng_generic_pileup_perf.rs) | the `Driver`'s generator type loses a parameter; its factory stops being boxed at the call site. |

`examples/ng_generic_loci_dump.rs` and `examples/ng_generic_walk_probe.rs` needed **no edit** —
they pass a factory by value and the constructor boxes it.

---

## The trap, and how it is pinned

ng shares one chain-id allocator across every region of a chromosome, so two fragments of two
regions never carry the same id. `ChainIdAllocator::reset` exists for that: it drops
`pending_mates` and `active_count` and **preserves `next_id` and the three counters**.
`fold_region_walk` then folds two of those counters as *deltas* against the value they held
when the region opened, and its own doc records this corruption having happened before —
triangular-summed into a plausible `u64`, with `active_reads_high_water` surviving because it
is a max, "which is what would make the corruption look selective enough to rationalise".

A per-region reset written as `WalkerState::new(config)` is one keystroke away and compiles.

**And the mirror image sits one field away.** `ActiveReads::reset` *preserves* `silent_exits`
as a run total — which was true while every region got a fresh set — and `fold_region_walk`
sums that one **per region**. A walker that lives for a chromosome shares one active set, so
the reset has to zero it. That is why `begin_region` exists beside `reset` rather than reusing
it.

`WalkerState::begin_region` destructures exhaustively, so a field added to the struct is a
compile error until someone decides which side of the line it falls on. The decisions:

| kept | cleared |
|---|---|
| `chain_ids` — `reset()`ed, never replaced | `summary`, per region because the fold sums it |
| `config` | `active_reads` via `begin_region`, tally and all |
| the four hoisted scratch buffers, whose capacity is the point of hoisting them | `open_records` — drained without `finalise`, as dropping the walker did |
| | `chrom_id`, `walker_pos`, `last_admitted_chrom_id`, `last_admitted_locus` |

**`last_admitted_locus` is the one that fires on real data rather than in a test.** Consecutive
regions overlap: each is asked for a halo past its end while the next is asked from its own
start, so a region is regularly served a read an earlier one already admitted. Carried across,
the coordinate-order check in `admit_read` rejects the next region's first read as going
backwards.

### Mutation results

Six mutations, each killing at least one test:

| mutation | tests that fail |
|---|---|
| drop `*last_admitted_locus = None` | `a_region_that_admitted_reads_past_the_next_regions_start_still_walks_it`, `a_reused_walker_answers_a_region_exactly_as_a_fresh_one_does` |
| replace the allocator instead of `reset`ing it | `the_chain_id_allocators_counters_survive_a_region_boundary` |
| drop `*summary = RunSummary::default()` | `a_reused_walker_answers_a_region_exactly_as_a_fresh_one_does` |
| `active_reads.reset()` instead of `begin_region()` | that oracle, plus `each_regions_silent_reads_are_folded_once_and_not_again_at_the_next_region` and three abandoned-walk tests |
| drop `forget_lookahead` from `LookAhead::move_to_region` | the oracle and the out-of-order test |
| never re-mint the cursor at a chromosome change | both chromosome-boundary tests |

The oracle is worth naming: **a walker pointed at a second region must be indistinguishable
from a fresh one pointed at the same region.** It is the only test that covers every field
`begin_region` decides about at once, and four of the six mutations land on it.

---

## Deviations from the plan and the architecture

**1. The reposition is a trait, not a bound on the walker.** Arch §2.2 says `move_to_region` is
"forwarded" through the walker; bounding `PileupWalker<I, F>` on a repositionable `I` would
break `genome_walk::run` for the stage-1 differential and every unit test that hands it a
`Vec`. So `RegionReadSource: Iterator<Item = PreparedRead>` bounds one extra `impl` block
instead, with an associated `Error` so `genome_walk.rs` never names a locus-generation error.

**2. The halo widening moved into `PreparedSampleReads`.** It was in `open_walk`; the walker's
`move_to_region` takes the **segment**, because the segment is what a failure is attributed to
and the walker has no business knowing why the span is wider. `PreparedSampleReads::query`
widens for the cursor. Deciding which span to ask for is still the caller's and the caller is
still this generator — spec §2's rule is unchanged, the code for it moved one type over.

**3. `make_reference` is boxed, not deleted.** The plan says "`make_reference` deleted — which
drops a type parameter". The consequence holds — `PileupGenerator<R, MakeReference, P>` is now
`PileupGenerator<R, P>` — but the factory itself has to stay, because `SampleReads::cursor`
takes one and a sample's k cursors interleave over the same coordinates: one accessor between
them would be one file position and one sliding window, its eviction driven by whichever asked
last. The two alternatives were worse. A shared `Arc<R>` (which does implement `RawRefSeq`)
gives up the per-file accessor for k > 1 *and* has the walk's REF fetches and the filter's
mismatch fetches share one window, where an ahead-fetch's eviction can throw away what the walk
still needs. Keeping the type parameter leaves the plan's stated consequence unmet. So:
`Box<dyn FnMut() -> R>`, one virtual call per file per chromosome, on nothing that is a hot
path any more. **The cost is a `'static` bound** on the constructor's factory argument, which
every caller already satisfies by capturing owned paths and shared indexes.

Perf-review L2 is closed either way: the factory was called once per file per *region* —
9,820 FASTA opens on chromosome 21, each parsing a 2,580-record `.fai` — and is now called once
per file per chromosome.

**4. `stopping_after` deleted.** It was a consuming builder with one caller, and
`move_to_region` sets the same field. Leaving it would be a second way to bound a walk.

---

## What D1+D2 does not do

- **The STR generator is untouched** — `ssr.rs:375` still calls `reads_in_region`. That is D3.
- **Nothing is deleted from the old path.** `SampleReads::reads_in_region`, `RegionReads`, the
  pool and `readers_opened` all still exist and are still used by the tools Milestone F moves.
  Two ways to read a file coexist by design until F.
- **Single-threaded.** One cursor per generator, minted at each chromosome boundary. The
  fan-out is a later plan.

---

## Verification

**The real-data anchor, chromosome 21 of HG002 30×** (`examples/ng_generic_walk_probe`):

```
loci=236081 observations=251786 reads_admitted=54709
```

— exactly the committed baseline, digit for digit.

**Both dumps byte-identical**, `ng_generic_loci_dump` and `ng_ssr_loci_dump` on the same
chromosome, against binaries built from `ee0c94b`: 251,792 and 4,406 lines, `cmp` clean. The
STR dump is a control here — D touches only the generic generator, so a difference there would
have meant something shared had moved.

**`cargo test --lib ng::` 1,555 passed** (1,549 at `ee0c94b`, plus six). `cargo clippy
--all-targets --all-features -- -D warnings` clean, `cargo fmt --check` clean, `cargo test
--examples` green.

### The number D owes

The end-to-end walk, three runs each, same host, release:

| | seconds | peak RSS |
|---|---:|---:|
| `ee0c94b`, a query per region | 4.521 / 4.501 / 4.512 | 20.45 MB |
| D1+D2, one cursor per chromosome | 1.873 / 1.878 / 1.878 | 20.80 MB |

**2.41×, with byte-identical output.**

**Read it with its region shape, as Checkpoint C's audit requires.** This is the typed-region
walk over chromosome 21 — 102,938 generic regions averaging 392 bases, each widened by 5,000 —
walked forward. It is *not* the read path's ~23×, which is what `ng_cursor_vs_query` measures
with no walk attached: the walk itself is now most of the remaining time, so the read path's
saving is diluted by everything it does not touch. Both numbers are real and they measure
different things.

It does reconcile the two figures Checkpoint C left unreconciled. The spec's own retention
prototype put the *walk* at 5.18 s → 2.69 s, which is 1.9×; this is 2.41× on the same shape of
work, so the ~23× read-path figure and the ~2× walk figure were never in conflict — one is a
component of the other.

**Memory on this fixture: +1.7 %.** 20.45 → 20.80 MB, about 360 KB, which is the order the
spec's arithmetic predicts for one cursor's kept reads at 30× depth (§10 estimates ~0.5 MB).

### Memory again, on the shape this fixture cannot produce

Review raised a Blocker against that number: the fixture is tandem-repeat-*targeted*, covering
0.64 % of positions (spec §1), and perf-review L2 warns by name that a **shared** reference
accessor walking a contig grows monotonically — `RawChromReader::fetch` extends its window
whenever the gap is under 64 KiB and only `evict_before` shrinks it, which nothing in
`src/ng/read/` calls. This change moves the mismatch filter's accessor from *per region* to
*per chromosome*, so on dense data it would extend rather than reposition, and L2 puts that at
~250 MB on chromosome 1. A sparse fixture has gaps wide enough to reposition, so it cannot show
the term at all. That is a correct reading of L2, and it was worth measuring rather than
arguing.

Measured, on a synthetic **contiguously covered 20 Mb contig at 30×** — the shape L2 is about —
with the fixture build in its own process so its gigabytes do not mask the walk:

| walked span | `ee0c94b` | D1+D2 | delta |
|---|---:|---:|---:|
| 2.5 Mb | 8.4 MB | 13.8 MB | +5.3 MB |
| 5 Mb | 10.9 MB | 18.8 MB | +7.9 MB |
| 10 Mb | 15.8 MB | 23.4 MB | +7.6 MB |
| 20 Mb | 25.6 MB | 29.0 MB | +3.4 MB |

**The delta does not scale with contig length** — it is non-monotone and *smallest* at the
largest span, which is allocator high-water noise, not a term proportional to the walk. Fitted,
both sides grow at essentially the same ~1 byte per base, and the cursor adds a roughly fixed
~6 MB. So the change does not add a contig-scaled buffer, and the Blocker is refuted on
magnitude: at 20 Mb dense 30× it costs **+3.4 MB on a 25.6 MB baseline**, and the walk is
**1.32× faster** there too (122.5 s → 93.1 s).

**⚠ What the measurement did surface is worth more than the finding it refuted, and it is not
D's.** *Both* sides grow at ~1 byte per base. That contig-scaled term already exists at
`ee0c94b`: `PileupGenerator`'s own `reference: Arc<R>` (the walk's REF fetches) and the
preparer's accessor are both **run**-lifetime and neither is ever evicted — `Arc` cannot give
the `&mut self` `evict_before` needs. On a dense chromosome 1 that is the ~250 MB L2 predicts,
today, before this change. ng's whole-genome memory claim — 30.1 MB for 18.5 M loci — was taken
on the same sparse fixture and inherits the same blind spot. **This is a pre-existing finding
for the owner, not a D fix**, and it belongs with the deferred per-chromosome reference registry
(spec §12).

The handoff's threshold — stop and ask if peak RSS rises materially above the ~24 MB baseline —
is not crossed by this change on either fixture.
