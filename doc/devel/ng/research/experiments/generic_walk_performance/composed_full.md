# Four changes in one tree: the full stack, measured across the depth sweep

Follow-on to `composed.md`, which composed the fast lane (A) and the mate-overlap skip (B)
and found the plain merge worth less than the sum until the heap was wired into the fast
lane's gate. This adds the ordered-active-set work (`ordered_active_set.md`) on top: **stage
1** (the contributor carries its active-set index), and **stages 2+3** (the active set becomes
an ordered queue, and the fold table's per-record collect-and-sort is deleted).

## The answer

**The full stack is worth −28.5 % at the 30× target and −42.5 % at ~130×, and it lifts the
worst point of the sweep from −6.7 % to −18.0 %.** Walk instructions retired, floor-subtracted,
against a pristine build of `6fbbd09`.

| depth / fixture | A+B (`composed.md`) | + stage 1 | **+ stages 2+3 — the full stack** |
|---|---:|---:|---:|
| 10× HG002 chr21 | −19.71 % | −19.77 % | **−20.95 %** |
| 30× HG002 chr21 | −25.88 % | −26.44 % | **−29.28 %** |
| **30× HG002 chr1 — the target** | −25.10 % | −25.72 % | **−28.49 %** |
| ~130× tomato `SL4.0ch01` | −37.09 % | −37.79 % | **−42.53 %** |
| 300× HG002 chr21 | −6.69 % | −9.27 % | **−18.02 %** |

**The sweep is flatter, but it is not flat.** A+B alone spans 30.4 points across the sweep
(−6.7 % at 300× to −37.1 % at ~130×); the full stack spans 24.5 (−18.0 % to −42.5 %). What
changed is almost entirely the deep end: **300× nearly tripled, from −6.7 % to −18.0 %**, while
10× barely moved (−19.7 % to −21.0 %). The peak is still ~130× and the trough is still 300×.

**The two lines of attack are complementary in exactly the way the shapes predicted, and the
overlap is smaller than the column-coverage argument suggests.** The fast lane answers 74 % of
columns at 30× and never reaches the general path, so a first guess is that at most a quarter of
stages 2+3 survives. Measured, **58 % survives at 30× and 94 % at 300×**:

| | on its own (`ordered_active_set.md`) | on top of A+B here | survives |
|---|---:|---:|---:|
| stage 1 @ 30× chr1 | −1.98 % | −0.83 % | 42 % |
| stage 1 @ ~130× | −2.86 % | −1.12 % | 39 % |
| stage 1 @ 300× | −3.33 % | −2.76 % | 83 % |
| stages 2+3 @ 30× chr1 | −6.46 % | −3.72 % | 58 % |
| stages 2+3 @ ~130× | −8.75 % | −7.62 % | 87 % |
| stages 2+3 @ 300× | −10.31 % | −9.65 % | 94 % |

Two reasons the survival is higher than 26 %, and both are measured rather than argued. The
fast lane covers 74 % of *columns* but only 64.9 % of *reads folded* at 30× — the ordinary
column is shallower than average — and stage 3's saving is a sort per record close, which
scales with the depth of the columns that reach the general path. And the fast lane's own
coverage falls with depth (74 % of columns at 30×, 33 % at 300×), so what it leaves behind at
300× is nearly everything, and nearly all of it deep.

**Against the disjoint-work reference the stack is still short, by less than before.** Two
savings on different work leave `(1−a)(1−b)`: A+B at −25.10 % and all three stages at −8.31 %
compose to −31.33 % at 30× chr1, measured −28.49 %, **2.8 points short**. At ~130× the reference
is −44.24 % against a measured −42.53 %, **1.7 points short**.

---

## 1. The numbers

One script, four binaries, alternating, three rounds per fixture, floors measured per binary per
fixture. Raw output in `composed_full_sweep.txt`; this is its summary, verbatim:

```
== c21_10x ==
base  runs 12.706-12.715 G  floor 1.8992 G  walk 10.807 G  rss 17.6-18.0 MB
heap  runs 10.576-10.581 G  floor 1.8990 G  walk 8.677 G  -19.71 %  fast_columns=276263
s1    runs 10.563-10.574 G  floor 1.8922 G  walk 8.671 G  -19.77 %  fast_columns=276263
full  runs 10.439-10.449 G  floor 1.8958 G  walk 8.543 G  -20.95 %  fast_columns=276263
== c21_30x ==
base  runs 17.917-17.928 G  floor 1.9001 G  walk 16.017 G  rss 17.9-20.0 MB
heap  runs 13.771-13.780 G  floor 1.8997 G  walk 11.871 G  -25.88 %  fast_columns=262498
s1    runs 13.682-13.688 G  floor 1.9003 G  walk 11.782 G  -26.44 %  fast_columns=262498
full  runs 13.228-13.249 G  floor 1.9005 G  walk 11.327 G  -29.28 %  fast_columns=262498
== chr1_30x ==
base  runs 110.498-110.501 G  floor 0.3492 G  walk 110.149 G  rss 17.9-20.5 MB
heap  runs 82.856-82.930 G  floor 0.3493 G  walk 82.507 G  -25.10 %  fast_columns=1746198
s1    runs 82.169-82.316 G  floor 0.3493 G  walk 81.820 G  -25.72 %  fast_columns=1746198
full  runs 79.122-79.213 G  floor 0.3496 G  walk 78.773 G  -28.49 %  fast_columns=1746198
== tom_130x ==
base  runs 217.443-217.602 G  floor 1.3089 G  walk 216.134 G  rss 368.9-371.8 MB
heap  runs 137.278-137.470 G  floor 1.3066 G  walk 135.972 G  -37.09 %  fast_columns=1069716
s1    runs 135.758-135.776 G  floor 1.3093 G  walk 134.449 G  -37.79 %  fast_columns=1069716
full  runs 125.509-125.529 G  floor 1.3072 G  walk 124.202 G  -42.53 %  fast_columns=1069716
== c21_300x ==
base  runs 91.447-91.490 G  floor 1.8958 G  walk 89.552 G  rss 19.6-21.1 MB
heap  runs 85.450-85.452 G  floor 1.8924 G  walk 83.557 G  -6.69 %  fast_columns=124293
s1    runs 83.148-83.196 G  floor 1.8949 G  walk 81.253 G  -9.27 %  fast_columns=124293
full  runs 75.307-75.314 G  floor 1.8956 G  walk 73.411 G  -18.02 %  fast_columns=124293
```

`base` is a pristine build of `6fbbd09`; `heap` is the state `composed.md` recommends (A + B +
the heap answering the fast lane's gate); `s1` adds stage 1; `full` adds stages 2 and 3.

**Every adjacent pair has disjoint ranges on every fixture**, with one exception: at 10× the
`heap` and `s1` ranges overlap (10.576–10.581 against 10.563–10.574), so **stage 1 is a null
result at 10×** and is reported as one.

**Peak RSS is neutral in these measurements.** At 30× chr1 the baseline spans 17.9–20.5 MB and
the full stack 17.8–20.2 MB — overlapping, so no claim either way. The −4 % to −11 % that
`ordered_active_set.md` reports at 30× was measured against a tighter baseline (21.00–21.14 MB)
than the one this host gave me today, and I do not reproduce it as a distinct effect.

**300× is measured on chr21 here**, where `ordered_active_set.md` used chr1. The two 300× rows
are therefore not directly comparable, and the disjoint-work arithmetic above is given only for
30× and ~130×, where the fixtures match.

## 2. Gates — and two of them moved, exactly as forecast

**All four dumps keep every line count** (251,792 / 4,406 / 1,718,914 / 11,945). The two SSR
dumps are `cmp`-identical. The two generic dumps differ by **exactly one line each**, and it is
the header:

```
=== chr21 generic ===
5c5
< # record_widen_events=423 column_depth_truncations=0 regions_in=205875 regions_handled=102938 loci_emitted=236081
---
> # record_widen_events=425 column_depth_truncations=0 regions_in=205875 regions_handled=102938 loci_emitted=236081
=== tomato generic ===
5c5
< # record_widen_events=622 column_depth_truncations=0 regions_in=527599 regions_handled=263800 loci_emitted=1711775
---
> # record_widen_events=621 column_depth_truncations=0 regions_in=527599 regions_handled=263800 loci_emitted=1711775
```

`diff` reports two changed lines per dump and nothing else. **423→425 and 622→621 are the exact
values `ordered_active_set.md` reports, and every locus line is identical.** The merge did not
make the damage worse.

Probe counters on chr21, final binary — all five exact:

```
loci=236081
observations=251786
reads_admitted=54709
mate_overlap_positions=39312
fast_columns=262498
```

Validation, in debug:

- `cargo test --lib` — **2,885 passed; 2 failed; 5 ignored**. The two failures are
  `parity::ng_agrees_with_production_where_production_fabricated_nothing` and
  `parity::every_divergence_from_production_is_one_of_the_six_named_classes` — the same two
  stage-2 tests, and **no others**. A+B alone was 2,885 + those two green, i.e. 2,887/0.
- `cargo test --examples` — 33 targets, all `ok`.
- `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo doc --no-deps` — 12 unresolved links, the recorded baseline.

### The one place the emitted loci genuinely move

At 300× — and only there — the walk's output differs, in the way `ordered_active_set.md`'s §6
already flagged for an owner decision:

| | baseline | full stack |
|---|---:|---:|
| `loci` (chr21 300×) | 241,038 | **241,030** |
| `mate_overlap_positions` (chr21 300×) | 410,497 | **410,496** |

`column_depth_truncations` is 0 at 10×, 30× and ~130× and 909 at 300×, so this is the depth-cap
sampling change and nothing else: reordering contributors changes which reads
`contributors.truncate(cap)` keeps. **`mate_overlap_positions` is unchanged at every depth where
the cap does not fire** — 12,726 at 10×, 39,312 at 30× chr21, 294,651 at 30× chr1, 166,371 at
~130× — so the one lost reconciliation at 300× is a consequence of the cap keeping different
reads, not of the mate-overlap heap. **My merge reproduces that divergence and adds nothing to
it.**

## 3. The conflicts, and how each was decided

Stage 1 applied over the composed tree by context, offsets only (`+55`, `+97`, `+248`). Stage 3
applied the same way to `genome_walk.rs` and `open_record.rs`. **`active_read_set.rs` was a
genuine textual conflict** — `git apply` refused it, because stage 3 deletes the `ahash::AHashMap`
import that the mate skip's `admit` hunk sits on.

**How it was resolved: stage 3's file was taken whole, and the mate skip's five pieces were put
back on top of it by hand.** That direction rather than the reverse, because stage 3 rewrites the
container and the mate skip only hangs one field and one method off it.

**3.1 — Imports.** `use std::cmp::Reverse;` and
`use std::collections::{BinaryHeap, VecDeque};` — the mate skip's two and stage 3's one, merged
into one `use`.

**3.2 — `admit`, the only substantive rewrite.** The mate skip pushed its heap entry inside the
`if let Some(partner) = … && let Some(partner_idx) = self.by_read_id.get(&partner)` arm. Stage 3
deletes `by_read_id` and finds the partner with
`self.reads.binary_search_by_key(&partner, |entry| entry.read_id)`. The push now hangs off the
binary search's `Ok(partner_index)`, and `alignment_start`/`alignment_end` are read off
`active.read` before `push_back` moves it — stage 3 moved the read into the queue at a different
point in the function than the code the mate skip was written against.

**3.3 — The heap and `min_alignment_end`, which the brief asked about.** They coexist without
interacting, and the field doc now says why: **`min_alignment_end` is a minimum over *reads*'
`alignment_end` and gates the expiry scan; `pair_overlap_ends` is a minimum over *pairs*'
overlap ends and gates the reconciliation.** Neither is derivable from the other — a set can
have many live reads and no live pair, or one pair whose overlap ends long before either mate
does. Both are cleared in `reset` and `begin_region`. `flush_all` resets `min_alignment_end` and
**not** the heap, which is how the mate skip shipped: `flush_chromosome_into` calls `reset`
immediately afterwards, and a stale heap can only stop the skip firing, never lose a
reconciliation. Left as measured rather than "tidied", so the composed state is the state the
numbers describe.

**3.4 — `swap_remove` → `remove` in `resolve_mate_overlap_at_pos`, on a path my heap now gates.**
This is safe, and the reason is worth stating rather than assuming: **the gate skips that
function only on columns where no two contributors share a chain id, and on such a column the
function removes nothing.** So the contributor ordering stage 2 needs is preserved on exactly
the columns where it could have mattered, and the `remove` is reached on every column the gate
lets through. The `debug_assert!` the mate skip carries — the all-pairs scan, run on every column
the skip claims — is what pins that, and it is still armed.

**3.5 — The three mate-skip tests** were re-added to stage 3's `active_read_set.rs` test module
unchanged; they use the pristine `paired_read` helper, which stage 3 does not touch.

## 4. What this round made newly fragile

**Three overlapping facts about read ends now live on `ActiveReads`**: the queue's own ascending
`read_id` order, `min_alignment_end`, and `pair_overlap_ends`. Each has one writer and a clear
set of reset sites, but a future change to admission or expiry must now keep three invariants
where before this merge each change kept one. The queue's order is load-bearing for stage 3's
speed, `min_alignment_end` for the expiry guard's correctness, and the heap for the
reconciliation's correctness.

**The census got slower, and only measurement runs pay it.** The `PVC_COLUMN_CENSUS=1` block
calls `get_by_read_id` once per contributor per column, which stage 2 turned from a hash lookup
into a binary search — O(D log D) per column where it was O(D). It is off by default and outside
every number in this report except §6's coverage census.

**Everything `composed.md` §8 lists still applies**, unchanged: the mate skip's assertion covers
only the columns that reach the general path, the fast lane's covers the rest, and
`may_have_mate_overlap_at`'s "a false is exact" claim is load-bearing in three places.

## 5. What I recommend

1. **Take the whole stack at 30× and below.** −28.5 % on the target workload, every locus line
   identical, one header counter moved on each of two dumps.
2. **The two decisions `ordered_active_set.md` raises are still the blockers**, and neither is a
   performance question: the mate-overlap tie-break's divergence from production (the two red
   `parity::` tests), and the depth-cap sampling change at coverages that reach
   `max_snp_column_depth`. My merge changes neither and reproduces both exactly.
3. **Do not build the `matches_only` refinement** — §6.
4. **Do not take the fast lane's sort removal** — §7.

## 6. The optional refinement, sized — and its premise does not hold

The suggestion was that at 300× the fast lane's coverage halves because `matches_only`
disqualifies a whole column when **any** read carries an `I` or `D` op **anywhere**, and that an
exact per-read test (no `D` op at all, plus no indel anchored here) would keep the coverage.

**Measured, on chr21 at 300×**, with the fast lane switched off so the census sees every column
(`PVC_FAST_COLUMN=0 PVC_COLUMN_CENSUS=1`), verbatim:

```
census_columns=376464
census_columns_ordinary=134559
census_contributors=61783163
census_contributors_ordinary=11340198
census_reject_record_already_open=2605
census_reject_indel_event=2151
census_reject_read_has_deletion=149970
census_reject_mate_overlap=186847
census_reject_depth_cap=909
census_reject_multi_read_group=0
census_reject_read_has_indel=170392
census_columns_simple=124542
census_contributors_simple=10310913
census_columns_simple_with_mate=205527
census_contributors_simple_with_mate=27703662
```

`columns_simple` is what the built predicate admits (124,542, against `fast_columns=124,293`
actually taken); `columns_ordinary` is what the exact predicate would admit (134,559). Read as
the share of **reads folded**, which is what costs:

| chr21 300× | reads in a fast column | of every 10,000 |
|---|---:|---:|
| the built predicate (`matches_only`) | 10,310,913 | 1,669 |
| the exact predicate proposed | 11,340,198 | 1,835 |
| if mate overlap were handled too | 27,703,662 | **4,483** |

**The refinement is worth 1.7 points of read coverage at 300×.** At 30× chr21 the same
comparison is 6,485 in 10,000 against 7,020 — 5.3 points — so it is worth *more* at 30× than at
depth, the opposite of the premise.

**What actually halves the coverage at 300× is mate overlap, not `matches_only`.** It
disqualifies 186,847 of 376,464 columns (49.6 %), against 170,392 (45.3 %) for "some read
carries an indel op somewhere". Handling it would take read coverage from 1,669 in 10,000 to
4,483 — **a gain of 28 points, sixteen times the refinement's**. So the refinement is not the
lever at 300×; the lever is the one both earlier reports declined, and it still carries the cost
they named (a second copy of the reconciliation tie-breaks). **Not built, and not recommended.**

## 7. One cut the merge enabled, measured and declined

With stage 2 the active set iterates in ascending `read_id`, so the fast lane's own
`sort_unstable_by_key(read_id)` — which exists to fix `q_sum`'s summation order — is sorting
already-sorted input, once per ordinary column. Deleting it (replaced by a `debug_assert!` that
the input is ordered) is **byte-identical** — the same two header lines and nothing else, 2,885
passed / 2 failed — and worth:

```
== c21_30x ==   full walk 11.346 G   without the sort 11.304 G   -0.37 %
== chr1_30x ==  full walk 78.753 G   without the sort 78.608 G   -0.18 %
== tom_130x ==  full walk 124.274 G  without the sort 123.727 G  -0.44 %
== c21_300x ==  full walk 73.406 G   without the sort 73.348 G   -0.08 %
```

Ranges disjoint on all four, so the gain is real and it is between a fifth and half a percent.

**Reverted, and the reason is not the size.** It would make this module's `q_sum` summation
order — and so the emitted bytes — depend on a container choice in `active_read_set.rs`, where
today it depends only on a sort this module owns. `ordered_active_set.md` is careful to keep
that property for its own `FoldedReads` (§4: *"correctness does not depend on ascending
arrival … only the speed does"*), and this cut would give it up for half a percent. The patch is
kept as `composed_full_fast_column_sortremoved.rs` so the decision can be revisited, and the
declined sort now carries a comment saying what it costs.

## 8. Files and worktree

Worktree `/Users/jose/devel/pop_var_caller/.claude/worktrees/agent-a47148967becb854a`, detached
at `6fbbd093764662ed2496acde39424c8ee234ea1c`, **left in the full-stack state** (A + B + heap
gate + stages 1, 2, 3; the §7 cut reverted).

```
 examples/ng_generic_walk_probe.rs                 |  28 ++
 src/ng/locus_generation/pileup/active_read_set.rs | 346 +++++++++++++++++-----
 src/ng/locus_generation/pileup/cigar_cursor.rs    |  84 ++++++
 src/ng/locus_generation/pileup/genome_walk.rs     | 283 +++++++++++++++++-
 src/ng/locus_generation/pileup/mod.rs             |  55 ++++
 src/ng/locus_generation/pileup/open_record.rs     | 281 +++++++++++++-----
 src/ng/locus_generation/pileup/tests.rs           |  44 +++
 7 files changed, 958 insertions(+), 163 deletions(-)
```

| file | what it is |
|---|---|
| `composed_full.diff` | the full stack, tracked files, over `6fbbd09` |
| `composed_full_fast_column.rs` | the new module (untracked in the worktree) |
| `composed_full_fast_column_sortremoved.rs` | the §7 variant, declined |
| `composed_full_sweep.txt` | raw output of the four-state depth sweep |
| `composed_full_sortremoval.txt` | raw output of the §7 A/B |
| `composed.diff`, `composed*.txt` | the A+B stack and its measurements (`composed.md`) |

`copy_fidelity.rs`'s two pinned files, `decompose.rs` and `chain_id_allocator.rs`, are untouched
by all four changes.

## 9. Which numbers are mine

Measured here: every row of §1, the gate diffs and counters of §2, the validation counts, the
census of §6, and the A/B of §7.

Cited from `ordered_active_set.md` and not re-measured: the three stages' standalone figures
(−1.98 / −2.86 / −3.33 % for stage 1; −6.46 / −8.75 / −10.31 % for stages 2+3 on top of it;
−8.31 / −11.36 / −13.30 % for all three), the +23 %-at-300× reproduction of the reverted attempt,
and its peak-RSS reduction at 30×. Their 300× measurements are on chr1; mine are on chr21.
