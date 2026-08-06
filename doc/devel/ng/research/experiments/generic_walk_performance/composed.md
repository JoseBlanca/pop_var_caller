# The two changes in one tree: do the gains add?

**Merged as they were written, they do not. Wire them together with eight more lines, and they
do.**

Two changes were built and measured separately against `6fbbd09`. **A**, the *fast lane*: a
column with one covered base, no record open over it, no read carrying an indel and no mate
pair is answered in scalars instead of by the general fold. **B**, the *mate skip*: the
active-read set keeps a min-heap of the positions where each pair it holds stops overlapping,
and when nothing is in it the general path's depth-sized `(chain_id, index)` build and sort is
skipped whole.

| | HG002 chr1, 30× | tomato `SL4.0ch01`, ~130× |
|---|---:|---:|
| the fast lane alone | −23.81 % | −34.38 % |
| the mate skip alone | −1.70 % | −4.65 % |
| **the two merged as written** | **−24.05 %** | **−34.94 %** |
| if the two saved on disjoint work | −25.11 % | −37.43 % |
| **merged, with the heap also answering the fast lane's gate** | **−25.22 %** | **−37.23 %** |

Read the fourth row as the reference the question asks for. Two savings that touch different
work leave `(1−a)(1−b)` of the walk behind, so 23.81 % and 1.70 % compose to 25.11 %, not to
25.51 %.

**Merged as written the composition falls 1.1 points short at 30× and 2.5 points short at
130×** — one change eats most of the other. The mate skip keeps **about a fifth** of its
standalone value: on top of the fast lane it is worth −0.32 % at 30× and −0.85 % at 130×,
against −1.70 % and −4.65 % on its own. The reason is exactly the overlap the brief
anticipated: the fast lane takes 74 % of columns and those are, by its own gate, the columns
with no mate pair on them — the ones the skip was saving the sort on. The skip is left with the
quarter of columns that fall through.

**The follow-through recovers all of it.** The fast lane's gate asks "do any two contributors
share a chain id?" and settles it with a sort over the contributors' chain ids. The heap
already answers that question, exactly, in O(1). Passing its answer into `try_ordinary_column`
and running the sort only when the heap says a pair could be present brings the composition
back to the reference: **−25.22 % at 30× (0.11 points past it) and −37.23 % at 130× (0.20
points short of it)**. Now the skip keeps 93 % of its standalone value at 130× (−4.33 % on top
of the fast lane against −4.65 % alone) and slightly more than all of it at 30× (−1.86 %
against −1.70 %), because it is being applied to two copies of the same sort instead of one.

The set of columns the fast lane accepts is **unchanged** by this: `fast_columns` is 262,498 on
chr21 and 1,069,716 on tomato in every state that has a fast lane. That is what makes the
shortcut safe — the heap's `false` is exact, so it can only skip a sort whose answer was
already known.

**At the target workload — human WGS at 30× — the composed change is worth −25.2 %**, of which
the fast lane contributes 23.8 points and the mate skip 1.4.

---

## 1. Where the work was done

Worktree `/Users/jose/devel/pop_var_caller/.claude/worktrees/agent-a47148967becb854a`, detached
at `6fbbd093764662ed2496acde39424c8ee234ea1c`, left as measured. Both markers confirmed before
anything else:

```
HEAD is now at 6fbbd09 perf(ng): close a record without re-asking the fold table for what it holds
---
8
1
```

Beside this report: `composed.diff` (the tracked-file diff of the final state),
`composed_fast_column.rs` (the new module, untracked in the worktree), and the raw measurement
output — `composed_ab5.txt`, `composed_ab4.txt`, `composed_depth.txt`, `composed_confirm.txt`.

`git diff --stat` of the final state:

```
 examples/ng_generic_walk_probe.rs                 |  28 +++
 src/ng/locus_generation/pileup/active_read_set.rs | 134 +++++++++++
 src/ng/locus_generation/pileup/cigar_cursor.rs    |  84 +++++++
 src/ng/locus_generation/pileup/genome_walk.rs     | 266 +++++++++++++++++++++-
 src/ng/locus_generation/pileup/mod.rs             |  55 +++++
 src/ng/locus_generation/pileup/open_record.rs     |   2 +-
 src/ng/locus_generation/pileup/tests.rs           |  44 ++++
 7 files changed, 603 insertions(+), 10 deletions(-)
```

`copy_fidelity.rs`'s two pinned files, `decompose.rs` and `chain_id_allocator.rs`, appear in
neither diff and are untouched.

## 2. Instrument

`instructions retired` from `/usr/bin/time -l`, floor-subtracted, min of 3 runs a side, five
separately-built release binaries alternated within one script so a drift on this shared host
hits every side. `PVC_TRUST_REFERENCE_INDEX=1` throughout. Wall-clock is not reported.

**Floors are measured per binary per fixture, not cited.** They are `PVC_PROBE_MAX_LOCI=1`
runs, three rounds, interleaved with everything else. Measured:

| fixture | floor |
|---|---|
| tomato `SL4.0ch01` (49 GB CRAM) | 1.305–1.313 G |
| HG002 chr1 (30× BAM) | 0.3484–0.3500 G |
| HG002 chr21 (10×, 30×, 300× BAMs) | 1.893–1.903 G |

**The 1.900 G figure that has been circulating for the human fixture is a chr21 floor, not a
chr1 one.** On chr1 the start-up cost is 0.349 G, a fifth of it. The earlier report that flagged
this is right, and the chr21 rows of §5 below are where 1.900 G actually belongs.

## 3. The five states, on both fixtures

One script, five binaries, alternating, three rounds. Raw output in `composed_ab5.txt`; this is
its summary, verbatim:

```
== tom ==
base  runs 217.300-217.687 G  floor 1.3090-1.3127 G  walk 215.991 G
fast  runs 143.038-143.120 G  floor 1.3083-1.3103 G  walk 141.730 G  -34.38 %
skip  runs 207.262-207.523 G  floor 1.3049-1.3064 G  walk 205.957 G  -4.65 %
both  runs 141.836-141.869 G  floor 1.3058-1.3085 G  walk 140.531 G  -34.94 %
heap  runs 136.897-137.147 G  floor 1.3083-1.3106 G  walk 135.588 G  -37.23 %
== chr1 ==
base  runs 110.455-110.500 G  floor 0.3497-0.3499 G  walk 110.105 G
fast  runs 84.242-84.282 G  floor 0.3493-0.3496 G  walk 83.893 G  -23.81 %
skip  runs 108.579-108.869 G  floor 0.3484-0.3493 G  walk 108.231 G  -1.70 %
both  runs 83.976-84.009 G  floor 0.3495-0.3499 G  walk 83.626 G  -24.05 %
heap  runs 82.682-82.829 G  floor 0.3493-0.3499 G  walk 82.332 G  -25.22 %
```

`base` is a pristine build of `6fbbd09`; `fast` is A alone; `skip` is B alone; `both` is the
plain merge; `heap` is the merge with the heap wired into the fast lane's gate.

**Every pair of adjacent states has disjoint ranges on both fixtures.** The two closest are the
plain merge against the fast lane on chr1 (83.976–84.009 against 84.242–84.282) and on tomato
(141.836–141.869 against 143.038–143.120); both are separated by several times their own
widths.

### The same four states measured a second time, before the fifth existed

`composed_ab4.txt`, an earlier and completely separate run of states 1–4:

| fixture | fast | skip | both |
|---|---:|---:|---:|
| tomato ~130× | −34.40 % | −4.75 % | −35.20 % |
| HG002 chr1 30× | −23.86 % | −1.71 % | −24.03 % |

Three of the four states agree with the five-state run to within a tenth of a point. **The
plain merge is the one that does not**: its tomato minimum there was 141.415 G, below the whole
141.836–141.869 G range of the later run — a 0.3 % drift on one side, and the only place in
this review where two runs of the same binary disagree by more than their own ranges. Reading
the earlier number instead moves the shortfall from 2.5 points to 2.2 points and changes
nothing else.

### Re-measuring A and B here was worth doing

Both were re-measured in this worktree rather than cited, and both came out close to their
original reports but not identical: the fast lane −34.38 % against a reported −34.34 % at 130×
and −23.81 % against −23.87 % at 30×; the mate skip −4.65 % against −4.38 % at 130× and −1.70 %
against −1.39 % at 30×. The mate skip reads about 0.3 points better here than in its own
worktree on both fixtures.

## 4. What the fifth state is

Eight lines. `may_have_mate_overlap_at` is hoisted to the first thing `process_position` does,
above the fast-lane attempt, and its answer is passed in:

```rust
let may_have_mate_overlap = self.active_reads.may_have_mate_overlap_at(walker_pos);
```

and inside `try_ordinary_column` the chain-id sort becomes conditional:

```rust
if may_have_mate_overlap {
    scratch.chains.clear();
    scratch.chains.extend(scratch.reads.iter().map(|r| r.chain_id));
    scratch.chains.sort_unstable();
    if scratch.chains.windows(2).any(|w| w[0] == w[1]) {
        return Ok(FastColumn::Fallback);
    }
} else {
    debug_assert!( /* the same sort, in full, in debug builds */ );
}
```

The heap's `true` is an over-approximation — the pair may be silent at this base, or one mate
`N`-masked — so it is still settled by the sort. Its `false` is exact, so the columns the fast
lane accepts are the same ones, which is why `fast_columns` does not move.

The larger idea the other agent sized — handling mate overlap *inside* the fast lane, worth 9
to 13 more points of read coverage — was **not attempted**, per the brief.

## 5. The shape against depth

HG002 chr21 through the down- and up-sampled BAMs, all five binaries alternating, three rounds.
Raw output in `composed_depth.txt`. Baseline walk cost and each state's change against it:

| depth | loci | baseline walk | fast lane | mate skip | plain merge | **with the heap** | `fast_columns` |
|---|---:|---:|---:|---:|---:|---:|---:|
| 10× | 224,030 | 10.795 G | −18.78 % | −0.54 % | −18.85 % | **−19.72 %** | 276,263 |
| 30× | 236,081 | 16.014 G | −24.55 % | −1.80 % | −24.82 % | **−26.06 %** | 262,498 |
| 300× | 241,038 | 89.579 G | −6.03 % | −1.10 % | −6.19 % | **−6.83 %** | 124,293 |

Together with the two main fixtures, the composed change against depth is:

**10× −19.7 %, 30× −25.2 % (chr1) and −26.1 % (chr21), ~130× −37.2 %, 300× −6.8 %.**

**The 300× collapse is the fast lane's, not the mate skip's.** Its coverage halves there:
124,293 ordinary columns against 241,038 loci, where at 30× it is 262,498 against 236,081. A
deeper column is more likely to contain *some* read with an `I` or `D` op anywhere in its
CIGAR, and one such read disqualifies the whole column. So the change that pays 25 points at
the target depth pays 7 at 300×.

The mate skip's own depth shape reproduces the one its report describes — a peak near 130× and
a fall away on both sides — with one difference worth recording: **at 10× I measure −0.54 %
where that report measured a null (+0.01 %, ranges overlapping)**. Mine are disjoint (baseline
12.695–12.732 G against 12.635–12.639 G), so on this binary it is a small real gain rather than
nothing. It does not change the conclusion that 10× is where the skip stops mattering.

`mate_overlap_positions` is unchanged at every depth — 12,726 at 10×, 39,312 at 30×, 410,497 at
300× — which is the direct check that no reconciliation was skipped.

## 6. Gates

All four acceptance dumps, produced by binaries built from the **final composed source**, are
`cmp`-identical to the stored copies in `tmp/perf_review_2026-08-04_ng-generic-walk/`. Verbatim:

```
IDENTICAL generic chr21
IDENTICAL ssr chr21
IDENTICAL generic tomato
IDENTICAL ssr tomato
  251792 …/gate/generic_chr21.txt
    4406 …/gate/ssr_chr21.txt
 1718914 …/gate/generic_tom.txt
   11945 …/gate/ssr_tom.txt
```

Probe counters on chr21, final binary — all five exact, so both changes are still firing:

```
loci=236081
observations=251786
reads_admitted=54709
mate_overlap_positions=39312
fast_columns=262498
```

Validation, in debug, on the final source:

- `cargo test --lib` — **2,887 passed; 0 failed; 5 ignored**. 2,882 on a clean tree plus the
  mate skip's five; the fast lane adds none.
- `cargo test --examples` — 33 targets, all `ok`.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean, after the fix in §7.4.
- `cargo doc --no-deps` — 12 unresolved links, the recorded baseline; no new ones.

**The composed tree still fails loudly if the shared invariant breaks.** Forcing
`may_have_mate_overlap_at` to return `false` fails **16 tests** — the same count the mate skip's
own report recorded, so nothing was lost by merging. What changed is *which* assertion catches
it first:

```
thread '…::tests::mate_overlap_zeroes_lower_bq_contribution' panicked at
src/ng/locus_generation/pileup/fast_column.rs:243:9:
the ordinary-column path skipped its chain-id sort at 1 on a column where two contributors
share a chain id — the mate reconciliation was lost
```

That is the new assertion inside the fast lane, not the mate skip's `column_shares_a_chain_id`.
See §8.

## 7. Every conflict, and how it was resolved

`git apply` placed all hunks of both diffs by context with line offsets only — **no textual
conflict, on any file**. The merge work was semantic, in four places.

**7.1 — Two preludes at the top of `process_position`, and their order decides the answer.**
Both diffs insert at the same point. The fast lane's block returns early on 74 % of columns;
the mate skip's `may_have_mate_overlap_at` call lands below it. Applied as written that is
*correct* — pruning the heap is lazy and monotone in `walker_pos`, and both `reset` and
`begin_region` clear it — but it means the predicate is consulted only on the columns that fall
through, which is precisely what makes the skip nearly worthless in the merged build. **This is
the composition question, sitting in the diff as a line-ordering detail.** Resolved by hoisting
the call above the fast-lane attempt and passing its answer into `try_ordinary_column`; the
duplicate call in the general path was removed.

**7.2 — `let walker_pos = self.walker_pos;`.** The mate skip's hunk depends on that binding,
which in the fast lane's version of the file sits *below* the fast-lane block. Moved to the
first line of `process_position` so both paths read it.

**7.3 — The probe's stdout key set.** The mate skip adds `mate_overlap_positions` to
`ProbeReport` and to the key list its own test pins; the fast lane prints `fast_columns` to
**stderr** specifically so that key set stays as it was. They do not collide, and both survive
— the gate output above shows all five counters.

**7.4 — `try_ordinary_column` now takes eight arguments**, which trips
`clippy::too_many_arguments` at `-D warnings`:

```
error: this function has too many arguments (8/7)
   --> src/ng/locus_generation/pileup/fast_column.rs:141:1
```

Resolved with `#[allow(clippy::too_many_arguments)]` and a doc note. The alternatives are worse:
a context struct would rename the caller's own fields without hiding anything, and the change
that would genuinely shorten the list — moving the predicate out to `genome_walk.rs` — undoes
the one thing `fast_column.rs` exists to do. It is a lint suppression the merge introduced, and
it is the only one.

An attribute cannot change codegen, but it does change the binary, so the final clippy-clean
build was re-measured against the baseline (`composed_confirm.txt`): tomato 137.191–137.398 G
and chr1 82.921–82.974 G, against the measured 136.897–137.147 G and 82.682–82.829 G. Both land
0.1–0.2 % above the range recorded in §3 — host drift between two runs an hour apart, not the
attribute — and neither moves any conclusion at the 2.5-point scale this report works at.

## 8. What the merge made newly fragile

**The mate skip's debug assertion now covers a quarter of the columns it used to.**
`column_shares_a_chain_id` runs only on columns that reach the general path, and the fast lane
takes 74 % of them. That coverage is not lost — it moved into the new assertion in
`fast_column.rs`, which the mutation test in §6 shows is what now catches a broken predicate
first — but the two assertions are now each pinning half the walk, and neither one alone pins
all of it.

**One invariant, three readers.** `may_have_mate_overlap_at`'s claim that a `false` is exact was
load-bearing for one thing before: whether the general path runs its reconciliation. It is now
load-bearing for two more: whether the fast lane runs its sort, and — through that — which
columns the fast lane accepts. A false negative used to drop a reconciliation; it can now also
make the fast lane emit a locus for a column that should have fallen back. The failure is still
caught in debug, and still silent in release. The list of properties it rests on is the one the
mate skip's report gives (chain ids never recycled, a second mate takes its first mate's id via
`remove`, read ids unique within a region) — that list did not grow, but what depends on it did.

**The heap is consulted at every column again, which the plain merge had quietly stopped
doing.** Under the plain merge it was peeked only on fall-through columns, so entries
accumulated across long ordinary stretches and were popped in bursts. Correct either way, but
the fifth state restores the original amortisation, and anything that later moves the call back
below the fast-lane attempt would silently undo both that and the gain.

**Everything §5 of the fast lane's own report says still applies**, unchanged by the merge: the
four closed-form re-derivations of general-path rules, the one-step emission delay shared
between `sealed`, `close_aged_records_into`, `flush_chromosome_into`, `begin_region` and
`reached_stop`, and the two measurement knobs (`PVC_FAST_COLUMN=0`, `PVC_COLUMN_CENSUS=1`) still
in the diff and still inside the measured numbers on both sides.

Peak RSS was not re-measured here. Both agents reported it neutral for their own change and
neither change alters allocation behaviour, but that is cited, not checked.

## 9. Which numbers are mine

Measured in this worktree at `6fbbd09` plus the diffs: every instruction count, every
percentage, every floor, every counter, the depth table, the gate results, the validation
counts, and the mutation-test count of 16.

Cited and not re-measured: the census frequencies (7,898 and 7,789 ordinary columns in 10,000),
the read-coverage estimate for handling mate overlap inside the fast lane (9–13 points), the
profile shares attributing 6.8 % of the main thread to `sort_chain_index`, and the two original
reports' own A/B figures where §3 compares against them.
