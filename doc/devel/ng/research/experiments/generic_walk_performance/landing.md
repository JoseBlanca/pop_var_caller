# Five changes in one tree, measured and ready to commit

Worktree `/Users/jose/devel/pop_var_caller/.claude/worktrees/agent-a46a0a3c7f7aef843`,
detached at `6fbbd093764662ed2496acde39424c8ee234ea1c`, left in the assembled five-change
state. Nothing is committed and nothing is pushed. Both required markers were confirmed
before anything was applied: `spans_only_its_anchors` appears 8 times in `cigar_cursor.rs`,
`fn finalise_recycling` once in `open_record.rs`.

Diff beside this file as `landing.diff` (11 tracked files, 2,156 insertions, 325 deletions),
plus the two new modules the diff cannot carry — `landing_fast_column.rs` and
`landing_read_sampling.rs`, which belong at
`src/ng/locus_generation/pileup/{fast_column,read_sampling}.rs`.

---

## The answer

**The four performance changes are worth −28.6 % at the 30× target and −42.5 % at ~130×.
Adding the behaviour change gives back half a percent: −28.2 % and −42.2 %.** Walk
instructions retired, start-up floor subtracted, minimum of three runs a side, the three
binaries alternated inside one script.

| depth / fixture | performance changes only | **all five — what the tree holds** |
|---|---:|---:|
| 10× HG002 chr21 | −21.00 % | **−20.67 %** |
| 30× HG002 chr21 | −29.21 % | **−28.95 %** |
| **30× HG002 chr1 — the target** | −28.55 % | **−28.23 %** |
| ~130× tomato `SL4.0ch01` | −42.48 % | **−42.15 %** |
| 300× HG002 chr21 | −18.01 % | **−17.32 %** |

The behaviour change — capping evidence per position instead of refusing whole reads at the
door — costs **+0.41 %, +0.37 %, +0.44 %, +0.57 % and +0.84 %** across those five rows. It
is the price of not discarding reads, and it is under one percent everywhere.

**The performance side reproduces `composed_full.md` to within a tenth of a percent on all
five fixtures** (−20.95 / −29.28 / −28.49 / −42.53 / −18.02 there, against −21.00 / −29.21 /
−28.55 / −42.48 / −18.01 here), measured today on a fresh build of the same commit. The
behaviour side reproduces `depth_cap.md`'s +0.37 / +0.14 / +0.48 within half a percent.

**The gates held.** All four acceptance dumps keep their exact line counts; the two SSR
dumps are byte-identical; the two generic dumps differ by exactly one line each, and it is
the header counter `record_widen_events` moving 423 → 425 on chr21 and 622 → 621 on tomato —
the two values `ordered_active_set.md` forecast. **No locus line moves on any of the four
fixtures.** All five probe counters are exact.

### Three things surprised me

**1. The two changes really did interact, in a place nobody looked, and one test caught it.**
The scalar path for the ordinary column emits its locus and returns *before* the block that
the depth-cap change added to count positions the hold ceiling left short of the cap. So on
a fast column the owner's test of success — "is any position left short of the cap while
reads covering it exist?" — could not fire, and the heap of ceiling losses was never pruned.
`no_position_is_short_of_the_cap_while_reads_covering_it_exist` failed on the first build of
the merged tree and is what found it. **Fixed** by falling back to the general path whenever
that heap holds anything; §3.3 has the reasoning and the cost. With the fix, forcing the
ceiling down to 4,096 at the deep tomato spot reports `positions_short_of_cap=346`, which is
exactly what `depth_cap.md`'s own sweep reports for that ceiling.

**2. One of the two red tests came back green, and it should have.** The brief expected both
`parity::` tests to stay red. Only one does. The other,
`ng_agrees_with_production_where_production_fabricated_nothing`, failed on the four-change
tree at a **capped column** — case 7, where two walkers kept different reads (mapping
qualities 32 against 23) — and the depth-cap change's exclusion of capped cases covers
exactly that. §6 gives the exact counts and the verification against the four-change tree.

**3. The "88,351 rows changed" figure is wrong by three orders of magnitude, and the
behaviour change removes the effect entirely.** Re-measured the way `depth_cap.md` §5
re-measured its own: **88 loci** of 240,912 differ in evidence, not 88,351 rows. And with
the depth cap in the tree first, the ordered active set changes **zero** loci at 300×. §7.

---

## 1. The numbers

One script, three binaries, alternating, three rounds per fixture, floors measured per binary
per fixture with `PVC_PROBE_MAX_LOCI=1`, `PVC_TRUST_REFERENCE_INDEX=1` throughout. Wall clock
is not used anywhere: this host has 6 high-performance and 12 low-energy cores. Raw output in
`landing_sweep.txt` (10×, 30× chr21, 30× chr1, 300×) and `landing_sweep_tom130x.txt` (~130×).

```
== c21_10x ==
base  runs 12.691-12.704 G  floor 1.8913 G  walk 10.799 G  rss 17.4-17.5 MB
four  runs 10.429-10.449 G  floor 1.8980 G  walk 8.531 G  rss 17.7-19.4 MB  -21.00 %
five  runs 10.470-10.486 G  floor 1.9031 G  walk 8.567 G  rss 17.5-17.7 MB  -20.67 %
     five vs four: +0.41 %
== c21_30x ==
base  runs 17.903-17.922 G  floor 1.8975 G  walk 16.005 G  rss 17.5-17.6 MB
four  runs 13.230-13.246 G  floor 1.9006 G  walk 11.330 G  rss 17.7-17.8 MB  -29.21 %
five  runs 13.274-13.277 G  floor 1.9027 G  walk 11.372 G  rss 17.8-20.1 MB  -28.95 %
     five vs four: +0.37 %
== chr1_30x ==
base  runs 110.499-110.525 G  floor 0.3494 G  walk 110.149 G  rss 17.8-20.1 MB
four  runs 79.052-79.161 G  floor 0.3497 G  walk 78.703 G  rss 19.9-20.1 MB  -28.55 %
five  runs 79.403-79.414 G  floor 0.3501 G  walk 79.053 G  rss 19.9-20.3 MB  -28.23 %
     five vs four: +0.44 %
== tom_130x ==
base  runs 217.188-217.361 G  floor 1.3059 G  walk 215.883 G  rss 368.2-375.0 MB
four  runs 125.489-125.551 G  floor 1.3076 G  walk 124.182 G  rss 368.7-371.7 MB  -42.48 %
five  runs 126.192-126.223 G  floor 1.3091 G  walk 124.883 G  rss 369.1-374.4 MB  -42.15 %
     five vs four: +0.57 %
== c21_300x ==
base  runs 91.446-91.454 G  floor 1.8953 G  walk 89.550 G  rss 19.3-20.2 MB
four  runs 75.320-75.331 G  floor 1.8957 G  walk 73.424 G  rss 20.0-21.1 MB  -18.01 %
five  runs 75.938-75.980 G  floor 1.8958 G  walk 74.042 G  rss 19.9-21.1 MB  -17.32 %
     five vs four: +0.84 %
```

`base` is a pristine build of `6fbbd09`; `four` is the four performance changes; `five` adds
the depth-cap behaviour change. **Every adjacent pair has disjoint run ranges on every
fixture**, so every percentage above is a real difference rather than run-to-run scatter.

**Peak RSS is neutral at every depth in the sweep.** At ~130×, where the walk holds the most,
the three ranges are 368.2–375.0, 368.7–371.7 and 369.1–374.4 MB — mutually overlapping, so
no claim either way. At 30× and below the ranges sit inside 17–21 MB for all three.

**Peak RSS at the deep spot, where the raised ceiling actually holds more reads,** is in §5
with the rest of that measurement: the walk holds 10,747 reads there instead of the 4,096 the
old ceiling allowed, and the peak does not move measurably.

**A fixture warning, because it cost me a full sweep.** The ~130× tomato fixture is
`benchmarks/tomato_big_cram/DRR000741.p1.cram`, **not** the tomato CRAM the gates use
(`benchmarks/ssr_tomato1/crams/SRR5079860.p1.bench.cram`). The gate CRAM is a sparse
population sample: over the same contig and the same 1 M loci it admits 61,753 reads against
1,501,890, and its walk is 28.7 G instructions against 215.9 G. Both walk `SL4.0ch01` against
the same reference, so the mistake is silent in everything except the counters. The gate
CRAM's rows are still in `landing_sweep.txt` under `tom_130x` and are **superseded** by
`landing_sweep_tom130x.txt`.

## 2. Gates

**All four dumps keep their exact line counts** — 251,792 / 4,406 / 1,718,914 / 11,945. The
two SSR dumps are `cmp`-identical to the stored copies in
`tmp/perf_review_2026-08-04_ng-generic-walk/round3/`. The two generic dumps differ by exactly
one line each, verbatim:

```
DIFFERS         g_chr21    251792 lines
5c5
< # record_widen_events=423 column_depth_truncations=0 regions_in=205875 regions_handled=102938 loci_emitted=236081
---
> # record_widen_events=425 column_depth_truncations=0 regions_in=205875 regions_handled=102938 loci_emitted=236081
BYTE-IDENTICAL  s_chr21      4406 lines
DIFFERS         g_tom   1718914 lines
5c5
< # record_widen_events=622 column_depth_truncations=0 regions_in=527599 regions_handled=263800 loci_emitted=1711775
---
> # record_widen_events=621 column_depth_truncations=0 regions_in=527599 regions_handled=263800 loci_emitted=1711775
BYTE-IDENTICAL  s_tom     11945 lines
```

`diff` reports one changed line per dump and nothing else, so **every locus line is
identical** on all four fixtures. 423 → 425 and 622 → 621 are the values
`ordered_active_set.md` and `composed_full.md` both predict.

Probe counters on chr21 at 30×, final binary — all five exact:

```
loci=236081
observations=251786
reads_admitted=54709
mate_overlap_positions=39312
fast_columns=262498
```

The release binary rebuilds byte-for-byte identical after every intermediate tree state I
visited while measuring, which is how I know the tree left behind is the one these numbers
describe.

## 3. The conflicts, and how each was decided

I merged by taking the four performance changes' tree as one side, a pristine
`6fbbd09 + depth_cap.diff` tree as the other, and their common ancestor as the base, then
running a three-way merge per file. **Two of the five shared files conflicted; three merged
with no conflict at all** (`mod.rs`, `tests.rs`, `examples/ng_generic_walk_probe.rs`). The
four files only the depth-cap change touches — `chain_id_allocator.rs`, `copy_fidelity.rs`,
`generator.rs`, `parity.rs` — applied cleanly, as did `read_sampling.rs` (new).

### 3.1 `active_read_set.rs` — one conflict, and it is the same one both prior merges hit

The depth-cap change adds three methods at the end of the `impl` and keeps `swap_remove`; the
ordered-active-set change deletes `swap_remove` and the `read_id → index` hash map it
maintained, because the queue is now ordered by `read_id` and a binary search replaces the
map.

**Decided: keep the ordered queue's shape and rewrite the depth-cap methods against it.**

- `worst_sampling_key` — taken unchanged. It scans the set and hashes as it goes; a
  `VecDeque` iterates exactly as the `Vec` it replaced did.
- `evict_by_read_id` — rewritten. It finds the read by `binary_search_by_key(read_id)`
  instead of the deleted hash map, and removes it with `VecDeque::remove` instead of
  `swap_remove`. **The removal has to preserve order**: three separate lookups — this
  method's own search, `admit`'s partner lookup, and `get_by_read_id` — are binary searches
  that a swap would silently break, and "silently" is the word that matters, because they
  would start missing reads rather than failing.
- `note_exit` — dropped. It existed to tally a silent exit from inside `swap_remove`'s
  caller; the ordered `expire_passed` already does that tally inline.

**Two stale-value questions, both answered by "stale is conservative, so leave it".**
`min_alignment_end` is a lower bound on the earliest end still held, and evicting a read can
only push the true minimum later — a stale value makes the expiry scan run when it need not,
which is slower and never wrong. Recomputing it would put an O(held) pass on the eviction
path for nothing. The same holds for the mate-overlap heap: an evicted read's pair entry
stays in it, and a stale entry can only stop the skip firing, never lose a reconciliation.
Both are stated in the doc comment rather than left for the next reader to re-derive.

### 3.2 `genome_walk.rs` — two conflicts, both adjacency

- **The region-reset destructure.** Both changes add a field to an exhaustive destructure
  that exists to make a new field a compile error until someone decides whether it is
  region-scoped. Merged by keeping both fields and both reset lines.
- **Two independent blocks inserted at the same point after the depth cap.** The depth-cap
  change adds the block that prunes ceiling losses and raises `positions_short_of_cap`; the
  performance work adds the measurement-only column census. They share no state. Merged by
  keeping both, with the ceiling-loss block first because it reads the same `depth` and `cap`
  the lines above compute.

### 3.3 The conflict that was not textual, and the one line of new code in this merge

The scalar path for the ordinary column returns as soon as it has built its locus. That
return sits **above** the ceiling-loss block. Two consequences, and the first is the one that
matters:

- `positions_short_of_cap` — the counter the whole behaviour change exists to drive to
  zero — could not fire on a fast column;
- the heap of ceiling losses was never pruned on a fast column, so it could grow across a run
  of them.

**Decided: while that heap holds anything, take the general path.** The gate becomes
`self.fast_column_enabled && self.ceiling_losses_by_end.is_empty()`.

It costs one `is_empty()` per column. The heap receives a push only when the set is at the
ceiling, and at the shipping ceiling of 32,768 that happens on **no fixture measured** —
`reads_shed_at_admission` and `reads_evicted_at_ceiling` are both 0 at 10×, 30×, ~130×, 300×
and at the deep tomato spot. So on a normal run the branch is false for the whole walk and
the fast lane's coverage is unchanged: `fast_columns` is 262,498 on chr21 at 30× and
1,069,716 at ~130×, both identical to the four-change stack's.

The test that found this is `no_position_is_short_of_the_cap_while_reads_covering_it_exist`,
which asserts the counter *fires* when the ceiling is put below the cap. Confirmed on real
data as well as in the test: at the deep tomato spot with the ceiling forced back to 4,096,
the assembled tree reports `positions_short_of_cap=346`, matching `depth_cap.md`'s sweep row
for that ceiling exactly, and `positions_short_of_cap=0` at 16,384 and above.

## 4. Re-taking the census figure the raised ceiling was supposed to move

**The 1,203 figure did not move, and it never could have — it is not measuring what
`depth_cap.md` §8 says it measures.**

The scalar path's first test compares the number of reads the walk is holding against the
8,000-read cap on how many reads one position may fold. It is deliberately conservative from
above: the reads that contribute at a position are a subset of the reads held, so a column
that passes the test cannot possibly be truncated by that cap. **That is still correct under
the raised ceiling** — the test can now reject, and rejecting means falling back to the
general path, which is always right.

`common_column.md` reports 1,203 columns at ~130× as the gap between what the census's
predicate admits and what the scalar path actually took, and `depth_cap.md` §8 attributes
that gap to this length test and predicts it will move. Re-measured on the assembled tree
over the same window (`landing_census.txt`):

```
census_columns_simple=1070919
fast_columns=1069716
active_reads_high_water=193
```

The gap is **1,203, unchanged**. And the third line is why: over that window the walk never
holds more than **193** reads, so a test against 8,000 cannot reject a single column there.
The 1,203 is the other conservatism the same report names — the indel test asking the whole
active set rather than only the reads that contribute here.

**Where the length test does fire, it costs 139 columns in 95,914.** At the deep spot on
`SL4.0ch01`, where the set reaches 10,747 held reads and 10,744 contributors at one position:

```
census_columns_simple=95914
fast_columns=95775
active_reads_high_water=10747
column_depth_high_water=10744
```

139 columns of 95,914 is **1.4 in 1,000**, and that gap still contains the indel
conservatism as well, so 139 is an upper bound on what the length test costs. It costs
nothing measurable in instructions — §5 shows the whole stack is *further* ahead at the deep
spot than at ordinary depth, not behind.

## 5. The mate-overlap heap and the ordered fold table, measured where the walk holds ten thousand reads

**Both get better as the walk holds more reads, not worse.** `depth_cap.md` §8 reasoned that
neither would be a problem and asked for a measurement; here it is.

The comparison is the depth-cap change alone against the depth-cap change plus the four
performance changes, over the same window — 100,000 loci starting at position 33,000,000 of
`SL4.0ch01`, reached with `PVC_PROBE_FROM_BP`. Floors were taken with the same
`PVC_PROBE_FROM_BP`, so the CRAM decode in front of the window is common to run and floor and
subtracts out. Minimum of three runs a side. Raw output in `landing_deepspot.txt` and
`landing_ceilingsweep.txt`.

| hold ceiling | reads actually held | walk, depth cap alone | walk, plus the four | difference |
|---:|---:|---:|---:|---:|
| 4,096 | 4,096 | 105.546 G | 98.316 G | **−6.85 %** |
| 8,192 | 8,192 | 96.796 G | 89.560 G | **−7.48 %** |
| 16,384 | 10,747 | 85.842 G | 78.369 G | **−8.71 %** |
| 32,768 (shipping) | 10,747 | 85.768 G | 78.390 G | **−8.60 %** |

The performance stack's advantage **grows** from 6.9 % to 8.6 % as the ceiling goes from
4,096 to what the data actually needs. If the min-heap of pair-overlap ends or the ordered
fold table degraded with the number of held reads, this column would shrink. It does the
opposite, and the two rows where the ceiling no longer binds (16,384 and 32,768, both holding
10,747) agree with each other to 0.03 G.

**Peak RSS at the deep spot does not move measurably either.** At the shipping ceiling the
depth-cap-alone binary spans 379.9–392.8 MB and the assembled tree 377.8–378.6 MB; at a
ceiling of 4,096 the same two span 374.2–375.2 and 374.1–379.7 MB. The ranges overlap in
every pairing, so the raised ceiling costs a few megabytes at most on this fixture — which is
what `depth_cap.md` measured too (it reports 10.7 MB, about 1.6 KB per extra held read). The
worst case it names — a completely full 32,768-read set costing about 53 MB, roughly 46 MB
more than the old ceiling's worst case — is **cited, not re-measured**: nothing on any real
fixture comes near a third of the ceiling.

## 6. Tests, lint and docs

Validation in debug, on the assembled tree:

- `cargo test --lib` — **2,893 passed; 1 failed; 5 ignored**.
- `cargo test --examples` — **33 targets, all `ok`, 0 failed**.
- `cargo clippy --all-targets --all-features -- -D warnings` — **clean**.
- `cargo doc --no-deps` — **12 unresolved links, exactly the recorded baseline**, none of
  them in any file this merge touches. (`ClassicStutterModel`, `SsrSegment`,
  `SsrSegmentCriteria::bundle_threshold`, `em::genotype_prior`, `is_close`, `scan_windowed`,
  and six document paths.)
- `cargo test --release` was **not** run: it is red on a clean tree.

### The one failure, and the one that came back

**Failing:** `parity::every_divergence_from_production_is_one_of_the_six_named_classes`, on
seed `0x5eed0001` case 18, and only on the run-summary counter comparison:

```
  left: Some(SummaryCounters { reads_admitted: 15, records_emitted: 76, record_widen_events: 4, mate_overlap_positions: 5, ... })
 right: Some(SummaryCounters { reads_admitted: 15, records_emitted: 76, record_widen_events: 3, mate_overlap_positions: 5, ... })
```

That is `record_widen_events` differing by one on a case where **no cap fires**
(`column_depth_truncations: 0` on both sides) — the same divergence the gate dumps show as
423 → 425, and the one the owner accepted and will name in the spec.

**Passing, where the brief expected it to fail:**
`parity::ng_agrees_with_production_where_production_fabricated_nothing`. Verified rather than
assumed: I rebuilt the four-change tree and ran it, and it fails there —

```
test result: FAILED. 2885 passed; 2 failed; 5 ignored; ...
    ng::locus_generation::pileup::parity::every_divergence_from_production_is_one_of_the_six_named_classes
    ng::locus_generation::pileup::parity::ng_agrees_with_production_where_production_fabricated_nothing
```

— on seed `0x5eed0001` **case 7**, where the two walkers emit the same locus with different
reads surviving (mapping quality 32 against 23, chain ids 9 against 10). That is a capped
column, and the depth-cap change's `Case::caps_can_fire` skips exactly those cases. **So the
two test edits compose: the exclusion covers one of the two failures and leaves the other,
which is a widen counter on an uncapped column and has nothing to do with caps.**

The counts add up with nothing lost: the four-change tree runs 2,887 tests, the assembled
tree 2,894, and the difference of 7 is exactly the depth-cap change's seven new tests.

### One thing that is red and was red before

`cargo fmt --check` reports **14 formatting differences** on the assembled tree against
**7 on pristine `6fbbd09`**. All seven new ones are in code the four performance changes
introduced — two in `fast_column.rs`, three in `genome_walk.rs`, two shifted blocks in
`open_record.rs`. **Left alone deliberately**: reformatting would change the binary every
number in this report was measured on, for no behaviour change. It is listed as its own
commit in §8.

## 7. Re-checking "88,351 of 341,094 rows"

**The figure overstates the change by three orders of magnitude, and with the behaviour
change in the tree the change disappears entirely.**

`ordered_active_set.md` §6 reports that the ordered active set moved 88,351 of 341,094 rows
on a 300× chr21 dump. `depth_cap.md` §5 found that the equivalent raw count on its own
change, 98,922 rows, collapses to 538 loci once chain-id renumbering is discounted — chain
ids are minted in admission order and numbered per file, so a handful of changed admissions
renumbers every id after them. Re-measured the same way, on the same fixture, both with and
without the behaviour change in the tree (`landing_recheck.txt`):

**Without the behaviour change** — pristine `6fbbd09` against the four performance changes:

```
removed lines: 88356
added lines:   88373
diff hunks:    61917
loci in both: 240912
loci only in A: 126
loci only in B: 118
loci whose evidence differs (chain-id column dropped): 88
loci differing only in chain-id numbering: 67009
```

The raw line counts reproduce `ordered_active_set.md`'s exactly. **Once the chain-id column
is dropped, 88 loci differ** — of 240,912 present in both dumps, that is **4 loci in 10,000**.
A further 67,009 loci differ *only* in the numbering of their chain ids, which is what
inflates the row count: 10 fewer reads admitted over the whole contig (532,761 → 532,751)
renumbers every chain id minted after them. Another 244 loci are present in one dump and not
the other. **So the figure a perf review states as 88,351 changed rows is 88 changed loci
plus 244 that appear or disappear — 332 in 241,000, or about one in 700.**

**With the behaviour change first** — the depth-cap change alone against the depth-cap change
plus the four performance changes:

```
removed lines: 1
added lines:   1
diff hunks:    1
loci in both: 241028
loci only in A: 0
loci only in B: 0
loci whose evidence differs (chain-id column dropped): 0
loci differing only in chain-id numbering: 0
```

The one changed line is the header: `record_widen_events` 1,634 → 1,627. **Not one locus
line moves, and not one chain id is renumbered, at 300×.** That is `depth_cap.md` §8's claim
confirmed on real data: once a capped position keeps the reads with the smallest sampling
keys, the kept set is a function of the reads alone, so reordering the container the walk
holds them in cannot change it. The test that pins the property,
`the_kept_set_does_not_depend_on_the_order_the_set_holds_reads_in`, is passing.

**This is the strongest single argument for the commit order in §8**: land the behaviour
change before the ordered active set and the ordered active set costs no emitted evidence at
any depth, only a header counter.

## 8. Suggested commit sequence

Five commits, each of which builds and passes its own gate. The order differs from the order
the changes were built in, and the reason is §7: putting the behaviour change fourth rather
than fifth means the last commit changes no emitted evidence anywhere.

| # | one-line message | gate this commit must pass |
|---|---|---|
| 1 | `perf(ng): skip the mate-overlap reconciliation when no pair is held` | all four dumps byte-identical; `cargo test --lib` all green |
| 2 | `perf(ng): answer the ordinary column in scalars` | all four dumps byte-identical; all green; `fast_columns=262498` on chr21 at 30× |
| 3 | `perf(ng): the contributor carries its active-set index` | all four dumps byte-identical; all green |
| 4 | `fix(ng): cap evidence per position instead of refusing reads at the door` | all four dumps byte-identical; all green, including 7 new tests; `positions_short_of_cap` 351 → 0 on the whole ~130× tomato contig |
| 5 | `perf(ng): hold the active set and the fold table in read order` | every locus line identical on all four dumps **and at 300×**; `record_widen_events` 423 → 425 and 622 → 621; one named `parity::` test red |

Notes the owner will need while doing it:

- **Commits 1 and 2 are one pair in effect but two commits in fact.** Each is byte-identical
  on its own. The min-heap that commit 1 adds is also what answers commit 2's own gate, so
  commit 2 written against a tree without commit 1 would be worth measurably less — that is
  `composed.md`'s finding, not a new one.
- **Commit 4 must carry the fast-lane fallback** described in §3.3 — one condition on the
  scalar path's gate. Without it commit 4's own new test fails.
- **Commit 5 must not be split.** Separated, the ordered fold table without the ordered
  active set is a 23 % regression at 300× (`ordered_active_set.md`). They land together or
  not at all.
- **A sixth, optional commit: `style(ng): cargo fmt the pileup module`.** The assembled tree
  has 14 `cargo fmt --check` differences against 7 on `6fbbd09`; the seven new ones all come
  from commits 1–3. Running `cargo fmt` on those files changes no behaviour and would need
  only a rebuild, not a re-measurement. I left the tree unformatted so it is exactly what the
  numbers above describe.

## 9. Which numbers are mine

**Measured here, on this host, today:** every row of §1 and §2, including a fresh
measurement of the four-change stack so the two columns are comparable; the census figures of
§4; the deep-spot and ceiling-sweep tables of §5; every test, lint and doc count in §6,
including the four-change baseline; and both re-checks in §7.

**Cited and not re-measured:** `depth_cap.md`'s whole-contig result that
`positions_short_of_cap` goes from 351 to 0 with 834,832 read-positions of evidence recovered
(a 13-minute run per state); its worst-case memory arithmetic for a full 32,768-read set
(~53 MB, ~46 MB above the old ceiling); `composed.md`'s finding that the scalar column path is
worth less without the mate-overlap heap feeding its gate; and `ordered_active_set.md`'s
finding that the ordered fold table without the ordered active set costs +23 % at 300×.

## 10. Files

| file | what it is |
|---|---|
| `landing.diff` | the assembled five-change tree, tracked files, over `6fbbd09` |
| `landing_fast_column.rs` | new module → `src/ng/locus_generation/pileup/fast_column.rs` |
| `landing_read_sampling.rs` | new module → `src/ng/locus_generation/pileup/read_sampling.rs` |
| `landing_sweep.txt` | raw depth sweep, 10× / 30× chr21 / 30× chr1 / 300× (its `tom_130x` rows use the gate CRAM and are superseded) |
| `landing_sweep_tom130x.txt` | raw ~130× sweep on `DRR000741.p1.cram` |
| `landing_deepspot.txt` | raw deep-spot A/B at the shipping ceiling |
| `landing_ceilingsweep.txt` | raw ceiling sweep at the deep spot, 4,096 → 32,768 |
| `landing_census.txt` | raw ordinary-column census, ~130× window and deep-spot window |
| `landing_recheck.txt` | raw 300× dumps compared, both with and without the behaviour change |

Scripts in the worktree under `tmp/landing/`: `sweep.sh`, `gates.sh`, `deepspot.sh`,
`ceilingsweep.sh`, `census.sh`, `recheck.sh`, and the two awk summarisers.
