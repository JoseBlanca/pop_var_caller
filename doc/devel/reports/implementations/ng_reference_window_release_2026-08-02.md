# ng — releasing the reference bases a walk has passed

*Branch `ng-generic-perf`, on top of `451c78e`. Found at Checkpoint D of
[the alignment-cursor plan](../../ng/impl_plan/alignment_cursor.md), but **not part of it**:
the defect predates the cursor and is fixed here on its own.*

---

## What was wrong

A reference reader in ng is a **sliding window over the FASTA**, not a whole chromosome. It
extends forward whenever it is asked for something near what it already holds, and it shrinks
only when someone tells it to.

Nobody told it. Both locus generators read the reference forward across a chromosome and never
released any of it, so a run held **one byte in memory for every base it had walked**.

Measured directly — a 150-base window fetched every 5 bases, which is a read pileup's shape,
across 2 Mb:

```
without releasing: 1,999,995 bytes resident
with releasing:            150 bytes resident
```

On human chromosome 1 that is around 250 MB, against a walk that otherwise peaks near 25 MB.

**Three readers were doing it, because they cannot be one.** The walk asks *what base is at this
position*; the read preparer asks for the window under each read it left-aligns; and each input
file's read filter asks for the window under each read it mismatch-checks. Each is a stateful
reader with its own position in the file, so sharing one would give several consumers one
position between them.

### Why it was never seen

Every measurement behind these generators used the tandem-repeat-targeted HG002 file. Its
coverage is sparse — 0.64 % of positions — so consecutive reads are far enough apart that the
reader **jumps** instead of extending. That is the one access pattern that cannot grow. The
whole-genome figure this caller is sold on, 30 MB for 18.5 M loci, was taken on that fixture.

This is not the alignment cursor's doing. It is equally present at `ee0c94b`, before any of that
work.

---

## The fix

**Eviction already existed and was already used** — the region-typing walk calls
`evict_before` after every window it fetches. The generators never did, and could not: the
method took `&mut self`, and a generator holds its reader as `Arc<R>`, because the walker must
*own* a reference and a shared handle is the only owned thing that does not rebuild the reader.
An `Arc` hands out `&self`.

Nothing had to be stretched to allow it. The window already lives behind a `RefCell`, because
*fetching* mutates it through `&self` too. So `EvictableRefSeq::evict_before` now takes `&self`,
and `Arc<T>` implements the trait by forwarding.

| file | change |
|---|---|
| [ref_seq.rs](../../../../src/ng/ref_seq.rs) | `evict_before(&self)`; `EvictableRefSeq for Arc<T>`; `resident_bases()` on the trait, so the bound is a test rather than an argument. |
| [raw_chrom_reader.rs](../../../../src/ng/raw_chrom_reader.rs) | `resident_bases()` promoted from a test-only helper. |
| [filtering.rs](../../../../src/ng/read/filtering.rs) | `ReadFilter::reference()` — shared, because releasing is all a caller needs it for. |
| [cursor.rs](../../../../src/ng/read/input/cursor.rs), [sample_cursor.rs](../../../../src/ng/read/input/sample_cursor.rs) | `evict_reference_before` / `resident_reference_bases`, reaching the filter's reader. |
| [read/mod.rs](../../../../src/ng/read/mod.rs), [left_align.rs](../../../../src/ng/read/left_align.rs) | the same pair on `ReadPreparer`, defaulting to nothing for a preparer that reads no reference. |
| [pileup/generator.rs](../../../../src/ng/locus_generation/pileup/generator.rs), [ssr.rs](../../../../src/ng/locus_generation/ssr.rs) | one release per region, to all three readers. |

### Where the release happens, and how far back it keeps

Once per region, before the region starts, less a margin.

Everything any reader will ask for from that point lies at or after the region's start minus one
of two things: a read overlapping the region may begin before it, and a record's footprint may
reach back to the read that opened it. The generic generator's margin is `max_record_span`,
which bounds both; the STR generator's is one flank, which is the furthest either of its readers
looks back from a repeat's own start.

**A wrong margin costs a re-read and never an answer.** Releasing is a hint: a base asked for
after it was released is simply read again. That is what lets the margin be generous rather than
exact, and it is why this could be got right without a proof.

---

## Verification

**Peak memory on a densely covered synthetic 20 Mb chromosome at 30×**, walked in 400-base
regions — the shape the real fixture cannot produce:

| walked span | `ee0c94b` | with the cursor (`451c78e`) | releasing |
|---|---:|---:|---:|
| 5 Mb | 10.9 MB | 18.8 MB | **3.7 MB** |
| 10 Mb | 15.8 MB | 23.4 MB | **3.4 MB** |
| 20 Mb | 25.6 MB | 29.0 MB | **3.5 MB** |

**Flat.** The term that grew with the length of the walk is gone, which is the result that
matters more than the eightfold drop at 20 Mb: extrapolated to human chromosome 1 the old
numbers pass 250 MB and these do not move.

**On the real fixture**, chromosome 21 of HG002 at 30×:

| | before | after |
|---|---:|---:|
| `ng_ssr_loci_dump` peak | 56.3 MB | **22.2 MB** |
| `ng_generic_walk_probe` peak | 21.8 MB | 21.2 MB |

The STR dump more than halves even here, because its walk covers whole contigs rather than only
the covered stretches.

**Nothing moved in the output.** Both dumps are byte-identical to binaries built from `ee0c94b`,
and the probe prints `loci=236081 observations=251786 reads_admitted=54709`. Walk time is
unchanged at 1.86–1.87 s, so the re-reads the margin allows for cost nothing measurable.

`cargo test --lib ng::` **1,570 passed**; clippy with warnings as errors, `cargo fmt --check`
and `cargo test --examples` clean.

### Mutations, each killed

| mutation | test that fails |
|---|---|
| `WindowedRefSeq::evict_before` does nothing | `a_forward_walk_holds_the_whole_span_unless_it_releases`, `a_shared_handle_can_release` |
| the generic generator stops releasing | `every_region_releases_the_reference_behind_it` |
| the STR generator stops releasing | `every_repeat_releases_the_reference_behind_it` |

**The two generator tests use a recording reference rather than a real window**, and that was a
correction during the work: the first version asserted on `resident_reference_bases` with an
in-memory reference, which holds no window and reports zero however badly the release is wired —
a test that could not fail, which is the failure mode this branch keeps hitting. Asking *what
was released, and when* is falsifiable; asking *how much is held* was not.

---

## What this does not fix

**The margin is per region, not per position.** The window still holds a region plus its halo
plus the margin — a few kilobases here — rather than a single read's width. Releasing per read
would be tighter and is not obviously worth it: the remaining term does not grow with the walk,
which was the problem.

**`OpenReference`, the bases CRAM decodes against, is untouched** and was never part of this. It
holds one chromosome at a time and drops it at a transition, and a BAM run never opens it at all.
