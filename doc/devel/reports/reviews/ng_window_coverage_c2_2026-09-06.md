# window coverage — C2 review: reliability, errors, and the numbers the step's prose claims

**Date:** 2026-09-06 (the design half), 2026-09-07 (this half)
**Reviewed:** `f2e7f31d feat(ng): C2 — the accumulator in the merge's per-sample window`
**Implementation report:** [ng_window_coverage_c2_2026-09-06.md](../implementations/ng_window_coverage_c2_2026-09-06.md)
**Fixes applied in:** the commit that follows C3 on branch `ng-window-coverage`

*Two reviews ran on C2. The design half — naming, idiom, module structure — ran in the same
session as the step and its fixes are inside `f2e7f31d` itself; the step's report lists them. **The
correctness half was still running when that session ended and wrote nothing**, so it was re-run
against `f2e7f31d` from a detached checkout. This file is the second half. The first half's
per-category output was written under `tmp/`, which is not in the repository — a mistake this file
does not repeat.*

## Verdict

**C2's central claim holds: every covered position of every sample reaches its accumulator exactly
once, and the two covers and the two evictors agree.** The per-sample cursor is a record *start*
rather than an index, so eviction cannot invalidate it; drawing one record past the reach cannot
double-count; and a ground fetch that fails part-way through a cover that crosses contigs leaves a
cache a retry completes without counting anything twice.

**What is wrong is which contig the buffer is left on.** A cover reads each contig it added records
on, ascending — which the accumulator requires — and the step's prose then says the region's own
contig comes last. Those are the same order only while no contig in the list sorts after the
region's, and one regularly does: a sample is drawn one record *past* the reach, that record can be
the first of the next contig, and the next contig sorts after. The cover then ends holding contig
*n + 1*'s bases while the builder about to run owns ground on contig *n*. It is invisible today
only because the accessor that would read it has no live caller until C4.

Three further defects sit under it, all silent, and all found by running rather than by reading.

## Findings

### 1. Major — the cover leaves the reference buffer on a later contig, over ground the builder does not own

`observation_cache.rs`, `read_the_ground_and_measure_coverage_over_it` and
`contigs_this_cover_must_read`.

**The scenario.** One sample with records at `contig 0:40` and `contig 1:5`; cover `contig 0:40-50`.
The second record is one past the reach, so it is drawn and held, unobserved — which puts contig 1
into the list of contigs to read. Sorted ascending, contig 1 is read last, and `fetch_the_ground`
refills the single buffer with contig 1's bases.

**Measured.** A probe that covers that fixture and then asks for the base at `contig 0:45` — the
middle of the region the cover was asked for — gets `None` where the reference has `A`.

**Why it matters before C4.** The serial driver's round is evict, cover, `with_observations`, so
the buffer a builder reads is the one this cover left. There is a smaller cost today as well:
`WindowedRefSeq` rebuilds its reader on any contig change, so the crossing cover opens the FASTA
for contig *n + 1* and the next cover opens it again for contig *n*.

**The fix.** Observe the contigs ascending, as the accumulator demands, then fetch the region's own
ground again where the last contig read was not the region's. One extra fetch, in the covers that
crossed and no others.

**The alternative, and why it was not taken.** Deferring the records on a contig after the
region's — leaving them to the cover whose own region reaches that contig — needs no second fetch,
but it rests on the driver never evicting such a record before it has been observed. The refetch
assumes nothing.

### 2. Major — the cursor cannot tell two records that share a start apart, and drops the second for ever

`observation_cache.rs`, `WindowCoverageInProgress::first_not_yet_observed` against `draw_next`'s
ordering guard.

The cursor is a record start and the not-yet-observed suffix partitions on `start <= seen`; the
source guard as committed was `start >= previous`. So two records at one start were one record to
the cursor. Inside a single cover both were still observed, because the observation walks by index;
**split across two covers the second never reached the accumulator at all** — no panic, no counter,
a sample simply short of one record's positions.

**Measured.** With the cursor at the first of two records sharing a start,
`summaries_not_yet_observed` over three summaries returns one, not two.

**The fix.** The guard becomes strict — `start > previous` — so what the cursor cannot represent
never enters the window. No existing fixture has two records at one start.

### 3. Major — the ground on a contig the cover is *leaving* is the last summary's reach, which is the largest only if records are disjoint

`observation_cache.rs`, `ground_this_cover_holds_on`.

C2's design review replaced a scan with two binary searches, arguing that reach is monotone across
a window of disjoint ascending records. On the region's **own** contig that holds without needing
disjointness: every held record starting at or before the chain widened the chain to its own reach,
and at most one record begins past it. **Off the region's contig there is no chain reach to lean
on**, and `held.last().reach()` is the largest only where reaches are monotone — which nothing in
the cache checks, since `draw_next` compares starts.

**Measured.** A sample holding a tract record over `contig 0:11-60` followed by a record at
`contig 0:12` panics when a later cover reads contig 0 as a contig it is leaving:

```
the cover read no reference base at GenomePosition { contig: ContigId(0), position: Position(13) },
where this sample holds a record — the ground a cover reads is computed from the records it holds,
so this is that computation and not the data
```

The panic's own message says the fault is in that computation, and it is right.

**Reachability: CONFIRMED as a panic, PLAUSIBLE as something a real store reaches.** The generic
mint cannot produce overlapping records, so this needs a source that does — which `draw_next`'s own
comment already anticipates for a run over stored files ("a source out of coordinate order is then
a fact about the file rather than a bug in this crate").

**The fix.** Take the largest reach where there is no chain reach to lean on. The scan is over one
contig's summaries in the covers that leave a contig, not over the window on every cover, so it
does not restore the cost the design review removed.

### 4. Major — the accumulator's order rule was debug-only, so the same source fills the wrong windows silently in release

`window_coverage/accumulator.rs`, `observe`.

Its doc said the check was `debug_assert!` and that whether it should hold in release "is a
question for the step that supplies the caller (plan step C2)". C2 supplied the caller and did not
answer it. **This repository ships with debug assertions off**, so with finding 3's fix in place —
which stops the missing-base panic — the same overlapping source reaches `observe` at position 12
after position 60 and nothing fires. A position behind the frontier finalises no centre and joins
the buffer out of place: the windows around it average over the wrong positions, `finalised` stops
being ascending, and `window_coverage_at`'s binary search is then searching an unsorted slice.

**The fix.** A release assert, one comparison per covered position, and the doc answers the
question rather than restating it.

**Left to the owner, and stated in the code.** The complete fix for 3 and 4 is to make the cache's
guard the disjointness the window is actually read by — `start > previous.reach()`. It was not
applied because `parallel.rs`'s
`a_builder_that_panics_leaves_the_cache_in_the_callers_hands_and_advanced` **deliberately** passes
an overlapping pair through this guard so that `build_region`'s own disjointness assertion is the
one that fires. Tightening it changes what that test is about, which is a decision rather than a
fix. Measured: with the disjointness guard, exactly two fixtures overlap and three `should_panic`
expectations need their substring changed.

### 5. Minor — the evidence C2's report calls "the strongest here" does not discriminate

The commit message and the report both say "the seven existing tests that failed when the
accumulator was first wired in are the strongest evidence here". **Reverting the multi-contig read
on the committed tree fails no pre-existing test** — only C2's own
`a_cover_that_crosses_a_contig_observes_the_records_left_on_the_one_it_leaves`. The seven were real
of an intermediate tree that had the accumulator without the contig `break` in the observation
pass; with that break a record on another contig is silently skipped rather than panicking on a
missing base. One test guards the multi-contig read today, and it is C2's own.

### 6. Minor — a source that mixes built and kept records misaligns evidence with summaries

`held_observations` and `held_bodies` are each pushed only on their own `Drawn` variant, so they are
parallel to `held_summaries` only for a source that answers uniformly. A source that mixed them
would hand the summary at index *i* another record's evidence. `depth.rs`'s release assert catches
the region mismatch and panics, so it is loud — but the message blames the store for what is a
same-index desynchronisation. **Not fixed:** it is a contract to write on `ObservationSource`, not
a line to change.

## Numbers checked

| claim | measured | how |
|---|---|---|
| `cargo test --lib --all-features` 6,358 / 0 / 15 | **6,358 / 0 / 15** ✓ | on `f2e7f31d` |
| "(6,353 before this step)" | **6,353** ✓ | on `9a6af0c6` |
| `cohort_merge` 278 passed | **278 / 0** ✓ | `-- ng::run::cohort_merge` |
| `cohort_merge::observation_cache` 47 passed | **47 / 0** ✓ | |
| "Five tests added … the whole increase" | **5** ✓ | 6,358 − 6,353 |
| **`PROJECT_STATUS.md`'s "Ten tests." for C2** | **five** ✗ | ten is C1's number, copied into C2's bullet |
| "no clippy diagnostic in this step's code" | **3 warnings**, all `needless_lifetimes`, none in changed code ✓ | |
| **the report's "`observation_cache.rs` 5" after the step** | **4 after, 5 before** ✗ | the commit message has this right; the report states the *before* number as the after |
| "about 1.2 million summary reads a cover at 3,000 samples … the ends are 6,000" | arithmetic right, **two different bases** ✗ | 3,000 × 200 = 600,000 for one contig against 6,000; 1.2 M is two contigs, whose ends are 12,000 |
| **`observation_cache.rs`'s "about one record in a thousand"** | **1 in 871 on the slice, 1 in 1,221 on the wider store** ✗ | the same commit gives 1 in 900 in `depth.rs`; this copy was left at a superseded figure |
| histogram "80.2 kB a sample" | ✓ | 50 × 401 × 4 B = 80,200 |
| **the report's "Changes made" names `read_the_ground_and_observe_what_this_cover_added`, `observe_its_new_records_on`, `window_at`, `for_each_reported_depth_of`** | **none of the four exists** ✗ | the design review renamed all four and the report's later section records two of the renames; the earlier section was not updated |
| the report's **Review:** link | **the file did not exist** ✗ | it is this one |
| the commit message's "its per-category file lands in `tmp/review_2026-09-06_window-coverage-c2/`" | **not in the repository** ✗ | `tmp/` is gitignored, so what the commit calls "the first thing to read before C3" cannot be read from a fresh checkout |

**One figure the documents do not give, and the project's range rule asks for:** at 3,000 samples
the accumulators cost about **241 MB** of histogram permanently (80.2 kB × 3,000) and up to
**786 MB** more transiently, while each sample's first 10,000 windows are held back before its depth
axis is fitted (a `Vec` reaching a capacity of 16,384 × 16 B = 262 kB, freed as soon as the width is
set). Plan step D3 prices the milestone; this is the shape of what it will find.

## Looked at and found sound

- **Exactly once, both directions.** The cursor is a start and not an index, so eviction cannot
  invalidate it; the not-yet-observed suffix is re-partitioned each cover. Two covers over the same
  ground give byte-identical windows to one.
- **The order the accumulator sees.** Within a cover, contigs ascend and each sample's held
  summaries ascend; across covers, a source that has yielded a contig-*n+1* record can never yield
  another on contig *n*, and a held contig-*n+1* record can never widen a chain reach on contig *n*.
  No regression could be built except through the non-disjoint source of finding 3.
- **A failure part-way through a cover.** `covered_to` advances only after the fetch succeeds,
  `reference_bases_from` is cleared before every fetch, and the cursor is set only after the depth
  rule returns `Ok`. A retry after a failed second-contig fetch counts nothing twice — measured.
- **Serial and parallel do not drift on the windows themselves.** Neither the VCF oracle nor the
  drivers' agreement test can see a window, because nothing read one at C2 — so "no VCF byte moved"
  is not evidence about this measurement. Compared directly over three covers and two samples of 600
  records, the two covers finalise the same windows and the two evictors drop the same ones. **Both
  probes are kept**, as the only check that the module's standing oracle holds for what C2 added.
- **Eviction keeps the finalised windows bounded.** One body for both evictors, with an exhaustive
  destructure, so a field added to the window fails to compile rather than leaking. At the defaults
  a sample retains about 750 centres between evictions — 18 kB, 54 MB across 3,000 samples.
- **`finish` has no caller in the merge**, so the last contig's final half-window of every sample is
  never finalised. That is C5's scheduled work and not a C2 defect; a contig change does finalise the
  previous contig's tail.

---

## What was applied, and where

Findings 2, 3 and 4 are fixed in the commit this file names at the top, with three regression
tests: a source with two records at one start is refused; a source whose records overlap is
refused **naming the order and not the ground**, which is what says both the ground fix and the
release assert landed; and the two probes the review recommended keeping, which compare the two
covers' windows and the two evictors' directly — the only check that the module's standing oracle
holds for what C2 added, since no output carries a window yet.

**Finding 1 is fixed in C3's commit instead**, because C3's own review found it independently and
C3's look-ahead is what turns it from a corner into the ordinary case at every contig boundary.

Every prose correction above is marked in C2's implementation report **in place**, beside the
sentence it corrects, rather than edited away.

Validation after the fixes, in the container: `cargo test --lib --all-features` 6,370 passed, 0
failed, 15 ignored; three `needless_lifetimes` clippy warnings, all predating this branch;
`observation_cache.rs` still at its four pre-existing rustfmt hunks and `accumulator.rs` clean; the
standing oracle unmoved at 2,311 records, sha256 `84ad19c2…`, on both routes.
