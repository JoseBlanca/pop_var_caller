# The per-column mate-overlap sort: skip it when no pair is there

**A mate pair is present at far fewer than half the columns, so the sort is skippable, and
skipping it is worth −4.4 % of the walk at 130× and −1.4 % at the 30× target.** It is a win
at three of four depths measured and never a loss: −4.38 % on the 130× tomato CRAM, −1.39 %
on 30× HG002 chr1, −0.93 % on 300× HG002 chr21, and nothing at all at 10× (+0.01 %, ranges
overlap). Every emitted byte is unchanged and `mate_overlap_positions` is identical on all
five fixture/contig combinations. **Worth its complexity**: the live logic is eleven lines
and one `BinaryHeap<Reverse<u32>>` field, and a debug-build assertion runs the old all-pairs
scan on every column the skip claims, so a wrong skip fails 16 tests rather than silently
changing a call.

The one thing to know before adopting it: the gain is *not* monotone in depth. It comes from
the columns where **no** pair is present, and that fraction falls as depth rises — at 300×
a pair is present at most columns and the skip stops firing.

## How often is a mate pair actually present?

`RunSummary::mate_overlap_positions` counts pairs reconciled, one per (column, pair), so
dividing by loci gives an **upper bound** on the fraction of columns holding a pair.

| fixture | depth | loci | `mate_overlap_positions` | columns with a pair |
|---|---|---:|---:|---|
| HG002 chr21 | 10× | 224,030 | 12,726 | at most **568 in 10,000** |
| tomato `SL4.0ch01`, 1 M loci | ~130× | 1,000,000 | 166,371 | at most **1,664 in 10,000** |
| HG002 chr21 | 30× | 236,081 | 39,312 | at most **1,665 in 10,000** |
| HG002 chr1 | 30× | 1,541,788 | 294,651 | at most **1,911 in 10,000** |
| HG002 chr21 | 300× | 241,038 | 410,497 | **1.70 pairs per column** — no bound below 10,000 in 10,000 |

So at both target depths **more than eight columns in ten hold no pair at all**, and today
every one of them pays a depth-sized tuple build plus a sort to find that out. The idea is
alive. At 300× it is dead by its own measurement, and that is what the 300× row of the
result table shows.

Note the tomato and HG002 columns are not comparable per-read: at 130× the tomato library
puts *fewer* pairs on a column (0.166) than the 30× human one does (0.191), because its
inserts are long relative to its reads. The tomato win is therefore the depth term (a
130-element sort) rather than a pair-scarcity term.

## Which half of the removed work is the sort?

`sample` over a `[profile.profiling]` build (`lto = false`, `codegen-units = 16` — these
shares say *where* the work is and do not transfer to release), 20 s of the 130× tomato
walk, 16,927 main-thread samples, with the skip forced off and `#[inline(never)]` wrappers
around the tuple build and the sort. Both wrappers were **reverted** before the diff below;
the raw profile is copied next to this report as `sample_attrib_mate_overlap.txt`.

| site | samples | share of main thread |
|---|---:|---:|
| `sort_chain_index` (the `sort_unstable`, all call paths) | 1,150 | **6.8 %** |
| `build_chain_index` (the `(chain_id, index)` tuple build) | 48 | **0.28 %** |

**The sort is 24× the tuple build.** The `Vec spec_from_iter` at 2.1 % in the standing
profile is therefore almost entirely someone else's allocation, not this one's.

The sort work in the profile is also *two* sorts, and only one of them is this function's:

```
Sort by top of stack, same collapsed (when >= 5):
        core::slice::sort::shared::smallsort::small_sort_general::hc0cde70314e34706  (in ng_generic_walk_probe)        684
        core::slice::sort::unstable::quicksort::quicksort::hf25bf99a5ba328e2  (in ng_generic_walk_probe)        295
        core::slice::sort::unstable::quicksort::quicksort::h78bbb40affe1b2ae  (in ng_generic_walk_probe)        258
```

`hf25bf99a5ba328e2` and `small_sort_general` sit under
`OpenPileupRecord::finalise_recycling` (620 of its 950 samples) — that is the record's own
key sort at close, which this change does not touch. `h78bbb40affe1b2ae` is the chain-index
sort. So of the ~9 % of the profile that is sorting, roughly **6.8 % is the one removed
here and the rest stays**.

## What was built

The O(1) counter, in the form the hypothesis suggested, and it was not awkward.

**`ActiveReads` keeps a min-heap of the reference positions at which each pair it holds
stops overlapping** (`pair_overlap_ends: BinaryHeap<Reverse<u32>>`). One entry is pushed in
`admit`, at the only moment a pair can be formed — the second mate's admission, where both
alignments' start and end are in hand. `may_have_mate_overlap_at(walker_pos)` pops the ends
the walker has passed and answers "non-empty". `process_position` calls it once and skips
`resolve_mate_overlap_at_pos` whole when the answer is no.

Two details make it exact rather than approximate:

- **The interval's start needs no bookkeeping.** Reads arrive in coordinate order, so the
  second mate starts at or after the first, and admission only happens once the walker has
  reached the new read's `alignment_start`. The pair is therefore already inside its overlap
  when the entry is pushed, and a read contributes to no column before it is admitted.
- **The heap needs no help from expiry.** The stored end is `min(both ends)`, which is at or
  before either read's own `alignment_end`, so the entry is always pruned by the time either
  mate expires.

The heap is cleared in `ActiveReads::reset` (chromosome boundary) and `begin_region`, both
of which restart `walker_pos` at 1; a stale end there would not lose a reconciliation but
would stop the skip firing over a whole opening megabase.

## The invariant, and where it could silently fail

**The claim: two contributors share a chain id only if they are the two mates of one pair,
and that pair was cross-linked at admission.** It rests on three properties of code this
change does not touch:

1. Chain ids are minted monotonically and **never recycled** (`chain_id_allocator.rs` module
   doc, and its test `released_ids_are_not_recycled`).
2. A second mate takes its first mate's id via `pending_mates.remove(&read.qname)` —
   `remove`, so a *third* read with the same qname mints a fresh id. At most two reads can
   ever hold one id.
3. Read ids are unique within a region, so a partner id that names an expired read cannot
   name a live one instead.

The failure mode to fear is therefore not a bug in this code but a future change to that
list: anything that lets a third read join a chain, or that puts a read into the active set
by a path other than `admit`, would make the heap incomplete and the skip would quietly drop
a reconciliation — smaller `mate_overlap_positions`, different bytes, no error.

**What pins it:**

- A `debug_assert!` in `process_position` runs `column_shares_a_chain_id` — the all-pairs
  scan the sort replaced in round 3 — on **every column the skip claims**, in every debug
  build. That is `cargo test --lib`, including the whole
  `ng_agrees_with_production_where_production_fabricated_nothing` differential.
- Five new tests: three on the predicate itself (visible at every position of an overlap,
  never visible for mates that do not meet, cleared at a region boundary) and two on the
  walk (`every_column_a_mate_pair_spans_is_reconciled_not_skipped` asserts the exact count 4
  — one staggered column plus three stacked — so losing one column still fails;
  `mates_that_never_overlap_reconcile_nothing` pins the other direction).
- **Mutation-tested.** Forcing `may_have_mate_overlap_at` to always return `false` fails
  **16 tests**, including three of the five new ones, five pre-existing mate-overlap tests,
  and four parity tests:

```
test ng::locus_generation::pileup::active_read_set::tests::a_pair_is_visible_at_every_position_its_two_alignments_share ... FAILED
test ng::locus_generation::pileup::tests::every_column_a_mate_pair_spans_is_reconciled_not_skipped ... FAILED
test ng::locus_generation::pileup::parity::ng_agrees_with_production_where_production_fabricated_nothing ... FAILED
```

The honest residual cost: one `BinaryHeap` field on `ActiveReads` (one `u32` per live
overlapping pair — at 130×, tens of entries), one push per paired read admitted, one peek
per column, and two more places that must be remembered when the active set is reset.

## The numbers

`instructions retired` from `/usr/bin/time -l`, floor-subtracted
(`PVC_PROBE_MAX_LOCI=1`, measured per fixture per binary), min of 3 runs a side,
`PVC_TRUST_REFERENCE_INDEX=1` on every run (`reference_check=trusted_unverified`).
Wall-clock is not reported: three other agents were measuring on this host throughout — one
of their `cargo build`s is visible in a `ps` taken mid-run.

| fixture | depth | baseline walk | with the skip | change |
|---|---|---:|---:|---:|
| HG002 chr21 | 10× | 10.687 G | 10.688 G | **+0.01 % (null, ranges overlap)** |
| HG002 chr1 | 30× | 110.123 G | 108.593 G | **−1.39 %** |
| tomato `SL4.0ch01`, 1 M loci | ~130× | 215.998 G | 206.540 G | **−4.38 %** |
| HG002 chr21 | 300× | 89.512 G | 88.681 G | **−0.93 %** |

Raw totals, in the order run (G instructions):

| fixture | baseline runs | skip runs | floors (base / skip) |
|---|---|---|---|
| tomato 130× | 217.316 / 217.350 / 217.602 | 208.007 / 207.857 / 208.095 | 1.318 / 1.317 |
| HG002 chr1 30× | 110.581 / 110.513 / 110.471 | 108.943 / 109.127 / 108.964 | 0.348 / 0.350 |
| HG002 chr21 300× | 91.432 / 91.409 / 91.467 | 90.620 / 90.564 / 90.603 | 1.897 / 1.883 |
| HG002 chr21 10× | 12.588 / 12.668 / 12.694 | 12.589 / 12.759 / 12.686 | 1.901 / 1.901 |

Ranges are **disjoint** on the three non-null rows (e.g. tomato: baseline min 217.316 above
skip max 208.095) and fully overlapping on the 10× row.

**Two measurement notes, both mine to own.** First, the HG002 chr1 floor I measure is
**0.348 G, not the 1.900 G the brief quotes**; 1.900 G is what a `PVC_PROBE_MAX_LOCI=1` run
costs on chr21 with the 10× and 300× BAMs, so the quoted figure looks like a chr21 floor.
Using 1.900 G instead would read −1.41 % where I report −1.39 %. The tomato floor agrees
with the brief (1.318 G measured against 1.306 G). Second, the probe binary I measured was
rebuilt after two comment edits, which shifts debuginfo and changes the binary; one
confirming run of each fixture with the final binary lands inside the measured range
(tomato 207.873 G; HG002 chr1 108.831 G, slightly better than the min I report).

Peak RSS on the 130× fixture: 391.2 MB baseline against 386.7 MB with the skip — inside the
386–391 MB run-to-run range already recorded for this walk.

### The shape against depth

The prize is *not* depth-shaped the way a per-column sort of D elements suggests, and this
is the finding to carry forward. Two terms move in opposite directions as depth rises: the
sort gets more expensive (good for the change) but the chance that some pair covers the
column rises too (bad for it). The second wins beyond ~130×:

- 10×: sort over ~10 entries — nothing to save, and 94 columns in 100 skip it. Null.
- 30×: 1,911 columns in 10,000 hold a pair. −1.39 %.
- 130×: 1,664 in 10,000 hold a pair, over a ~130-entry sort. −4.38 %, the peak.
- 300×: 1.70 pairs per column; the skip fires on the minority of columns that have none.
  −0.93 %.

## Gates

All four dumps re-run with binaries built from the final source and compared with `cmp`
against the stored baselines in
`/Users/jose/devel/pop_var_caller-ng-generic-perf/tmp/perf_review_2026-08-04_ng-generic-walk/`:

```
FINAL_GENERIC_CHR21_IDENTICAL     (251,792 lines)
FINAL_SSR_CHR21_IDENTICAL         (4,406 lines)
FINAL_GENERIC_TOM_IDENTICAL       (1,718,914 lines)
FINAL_SSR_TOM_IDENTICAL           (11,945 lines)
```

Probe counters on chr21, final binary:

```
loci=236081
observations=251786
reads_admitted=54709
mate_overlap_positions=39312
```

`mate_overlap_positions` **before and after, on both fixtures** — the direct check that the
skip never skipped a real pair:

| fixture | baseline | with the skip |
|---|---:|---:|
| tomato `SL4.0ch01`, 1 M loci (130×) | 166,371 | 166,371 |
| HG002 chr1 (30×) | 294,651 | 294,651 |
| HG002 chr21 (30×) | 39,312 | 39,312 |
| HG002 chr21 (10×) | 12,726 | 12,726 |
| HG002 chr21 (300×) | 410,497 | 410,497 |

Validation, in debug:

- `cargo test --lib` — **2,887 passed; 0 failed; 5 ignored** (2,882 on a clean tree, plus
  the five new tests).
- `cargo test --examples` — 33 targets, all `ok`.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo doc --no-deps` — 12 unresolved links, the recorded baseline; the two new intra-doc
  links resolve.

`copy_fidelity.rs`'s two pinned files (`decompose.rs`, `chain_id_allocator.rs`) are
untouched.

## The change

Diff: copied next to this report as `mate_overlap_skip.diff`, and live at
`tmp/mate_overlap/mate_overlap_skip.diff` in the worktree
`/Users/jose/devel/pop_var_caller/.claude/worktrees/agent-a1c072195b5cb91ef` (388 lines
across four files: `active_read_set.rs`, `genome_walk.rs`, `tests.rs`, and the probe, which
gained a printed `mate_overlap_positions` line so the counter is checkable from a run's
output). The live logic:

```rust
// active_read_set.rs — the field
    /// **ng's** — the last reference position at which each *pair* the set holds still has
    /// both of its alignments on the reference, smallest first.
    pair_overlap_ends: BinaryHeap<Reverse<u32>>,

// active_read_set.rs — in `admit`, inside the existing "second mate, partner is here" arm
            let overlap_end = alignment_end.min(self.reads[partner_idx].read.alignment_end);
            if overlap_end >= alignment_start {
                self.pair_overlap_ends.push(Reverse(overlap_end));
            }

// active_read_set.rs — the predicate
    pub fn may_have_mate_overlap_at(&mut self, walker_pos: u32) -> bool {
        while let Some(&Reverse(overlap_end)) = self.pair_overlap_ends.peek() {
            if overlap_end < walker_pos {
                self.pair_overlap_ends.pop();
            } else {
                return true;
            }
        }
        false
    }

// genome_walk.rs — in `process_position`
        let may_have_mate_overlap = self.active_reads.may_have_mate_overlap_at(walker_pos);
        ...
        if may_have_mate_overlap {
            resolve_mate_overlap_at_pos(contributors, &mut self.summary, &mut self.mate_overlap_buf);
        } else {
            debug_assert!(
                !column_shares_a_chain_id(contributors),
                "the mate-overlap skip fired at {}:{} on a column where two contributors \
                 share a chain id — the reconciliation was silently lost",
                self.chrom_id,
                walker_pos,
            );
        }
```

plus `pair_overlap_ends.clear()` in `reset` and `begin_region`, and
`column_shares_a_chain_id` — the all-pairs scan, live only in debug builds.

## Which numbers are mine

Everything in the tables above I measured in this worktree at `6fbbd09` + this diff. Cited
from the brief and round 3, not re-measured: the standing release profile shares, the
tomato floor of 1.306 G, and the −9.1 % that the round-3 sort itself was worth.
