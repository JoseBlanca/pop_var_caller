# Capping depth per position instead of throwing reads away at the door

Worktree `/Users/jose/devel/pop_var_caller/.claude/worktrees/agent-a0e8c78f4a1afebee`,
detached at `6fbbd09`. Diff beside this file as `depth_cap.diff` (10 files, 1,372
insertions, 167 deletions). Every number below is measured on this host unless it says
*cited*; commands and raw output are in the worktree's `tmp/depth_cap/`.

---

## The owner's test of success

**Is any position left short of the cap while reads covering it exist?**

On the whole ~130× tomato chromosome `SL4.0ch01` — 86,139,952 loci, 113.6 M reads, the
deepest real fixture available:

| | positions short of the cap | reads missing from them |
|---|---|---|
| **before** (ceiling 4,096, refuse at the door) | **351** | **834,832** |
| **after** (ceiling 32,768, sample per position) | **0** | **0** |

Those 351 positions are not a rounding error in what they lost: 834,832 read-positions of
evidence, an average of 2,378 reads missing from each. They sit in one place, around
33,037,543 on that contig, where the pile-up is about a hundred times the local typical
depth — and the walk was answering it by refusing 19,725 reads outright, each of which
then contributed at **no** position at all.

On 30× human chr21, 300× human chr21 and 30× human chr1 the number was already zero and
still is: the ceiling never binds there.

**The property has a boundary and it is worth stating plainly.** "No position is short"
holds *up to the ceiling*. Push coverage past what the walk will hold and reads start being
dropped again. The ceiling is 32,768 held reads; the deepest position measured on real data
needs 10,747. Where it starts to fail is measured, not guessed — the table in §2 shows
4,096 losing 346 positions, 8,192 losing 66, and 16,384 and above losing none.

---

## 1. What was wrong, in one paragraph

Three caps existed and only the wrong one was doing any work. `max_snp_column_depth`
(8,000 reads at a position with no indel) could never fire, because the walk held at most
4,096 reads. `max_indel_column_depth` (250) fired. And `max_active_reads` (4,096) refused
reads **at admission** — before the read was ever decomposed, so it reached no position and
nothing downstream could say which loci it would have covered. The effective rule was
therefore *throw whole reads away at the door, and never cap a position*, which is the
behaviour the owner called wrong: positions ended up with less coverage than the BAM had
for them.

A fourth ceiling sat underneath: the chain-id allocator's map of first mates waiting for a
partner had a hard, fatal cap of 10,000 entries. That constant was out of reach while the
walk held 4,096 reads and became reachable the moment it held more — **measured at 11,384
on the tomato deep spot** — so raising the hold ceiling without raising it too would have
converted a bad sample into an aborted run.

---

## 2. Part A — the ceiling stays, but high enough that it shapes nothing

The ceiling was **not** removed (owner: *"still with a high enough cap, otherwise we could
run out of memory"*). `DEFAULT_MAX_ACTIVE_READS` goes from **4,096 to 32,768**, and
`MAX_PENDING_MATES` from **10,000 to 1,000,000**.

### Where 32,768 comes from — depth against cost, both measured

Swept at the deep spot itself (`SL4.0ch01`, 100,000 loci from 33,000,000, ~130× CRAM,
`PVC_PROBE_MAX_ACTIVE_READS` overriding the ceiling):

| ceiling | reads refused | reads evicted | **positions short** | reads missing | held peak | deepest column | pending mates peak | peak RSS |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4,096 | 6,795 | 13,030 | **346** | 709,752 | 4,096 | 4,096 | 8,937 | 392.7 MB |
| 8,192 | 649 | 4,430 | **66** | 51,506 | 8,192 | 8,192 | 11,094 | 395.0 MB |
| 16,384 | 0 | 0 | **0** | 0 | 10,747 | 10,744 | 11,384 | 403.5 MB |
| 32,768 | 0 | 0 | **0** | 0 | 10,747 | 10,744 | 11,384 | 403.4 MB |
| 65,536 | 0 | 0 | **0** | 0 | 10,747 | 10,744 | 11,384 | 403.0 MB |

Read the "held peak" column first: it stops rising at **10,747** once the ceiling is out of
the way, so 10,747 is what this data actually needs. **The deepest position holds 10,744
contributors**, not the 12,792 the earlier note cited — that figure is not reproduced here,
and 10,744 is what the walk measures after read preparation on this contig.

16,384 would be enough for this sample. 32,768 is one doubling above it, chosen so a sample
somewhat deeper than the deepest one on hand does not silently start losing positions — the
failure mode is silent by nature, which is why `positions_short_of_cap` now exists.

### What the ceiling costs in memory

At the deep spot, holding 10,747 reads instead of 4,096 cost **10.7 MB** of peak RSS
(392.7 → 403.4 MB) — **about 1.6 KB per held read**, which is the read's sequence, its base
qualities, its CIGAR and its cursor.

**The worst case is the number to weigh, because that is what an out-of-memory failure
is.** A completely full 32,768-read ceiling costs about **53 MB** over an empty set,
against about **6.6 MB** for the old 4,096. So the change raises the walk's worst-case
resident footprint by roughly **46 MB per concurrent walk**. Nothing on any real fixture
comes near it: the largest set ever held here was a third of the ceiling, for a few
thousand positions of one chromosome.

On the **typical** case it costs nothing measurable. Whole-contig peak RSS on the ~130×
tomato chromosome, verbatim:

```
STATE 1 (ceiling 4,096)   386383872  maximum resident set size
STATE 2 (ceiling 32,768)  386400256  maximum resident set size
STATE 3 (shipping)        386580480  maximum resident set size
```

386.4 → 386.6 MB, a difference of 0.05 %. The walk's peak on this fixture is dominated by
the CRAM decode holding a reference contig, not by the active set.

---

## 3. Part B — which reads a capped position keeps

The old rule was `contributors.truncate(cap)`: keep the first `cap` contributors in whatever
order the active set happened to hold them. Two faults, and the second is the one the
owner's request is about.

- **It was not a fact about the reads.** The order is a permutation produced by
  `swap_remove`, so a change to how the set stores reads moved 88,351 of 341,094 emitted
  rows on a 300× walk without any read changing (*cited*, `ordered_active_set.md` §6).
- **It preferred reads that reach a position from the left.** Reads arrive sorted by
  alignment start, so wherever arrival order survives, "first `cap`" means
  "leftmost-starting `cap`".

The new rule: every read gets a 64-bit **sampling key** — FNV-1a over its query name plus
its mate role, run through a splitmix64 finaliser — and a cap keeps the `cap` smallest keys.
New module `src/ng/locus_generation/pileup/read_sampling.rs`.

Of the three properties asked for, **two hold exactly and one approximately**:

- **Deterministic across runs and processes** — exactly. The hash is written out rather than
  taken from `ahash` or `DefaultHasher`, both of which are free to seed per process.
  `parity::ng_emits_the_same_bytes_in_a_second_process` still passes, and the literal hash
  values are pinned in a test so changing the rule has to be deliberate.
- **Unbiased with respect to alignment start** — exactly, not approximately. The key is a
  function of the query name, which carries no positional information. Measured in §5.
- **Stable between adjacent positions** — approximately. The kept set is the `cap` smallest
  keys among the reads covering the position, so two neighbouring positions covered by the
  same reads keep the same reads, and the set changes only as reads enter and leave. What it
  does not survive is a change in the cap itself, or the ceiling evicting a read a later
  position would have wanted — and the second only happens in a region deep enough to fill
  the walk, which no fixture measured here reaches at the shipping ceiling.

**The key is not the chain id**, though the chain id is already on every contributor and
would need no lookup. The chain id is minted in admission order, so it is a fact about the
walk's history rather than about the read; hashing it would hide the ordering dependence
rather than remove it. The cost of using the query name is one hash-map lookup per
contributor **at capped positions only** — 909 of 251,792 positions on the 300× fixture, and
none at 30×.

---

## 4. The ceiling now evicts fairly rather than refusing arrivals

The owner asked for this if it were free. **It is free, and it is built.**

When the set is full, the arriving read and the reads already held are compared on the
*same* sampling key a capped position uses. If the arrival beats the worst held read, that
read is given back and the arrival takes its place; otherwise the arrival is refused.

**Nothing is paid for it on the common path.** No key is stored on a read and no index is
kept in key order: the keys are computed by an O(held) scan inside the branch that only runs
when the set is already full — the same condition that used to make the walk refuse the read
outright. Measured (§6): the sampling rule and the eviction branch together are within noise
of zero at both 30× and ~130×, in both directions.

**What eviction has to unwind, and why it turned out to be nothing.** A read given back may
already have folded into open records. That is not a new state: it is exactly the state
every read reaches when the walker passes its end while a record it folded into is still
open, and `refold_live_reads` has always handled it by skipping a folded read that is no
longer live. Its chain id stays minted (ids are unique `u64`s that are never reused, so a
spent one is harmless), its pending-mate entry is governed by `mate_lookup_window` exactly
as an expiry's is, and its mate cross-link dangles exactly as an expiry's does. **Eviction
adds no case to the fold; it makes an ordinary case arrive sooner.** The only consequence is
that the evicted read's witness stops being updated if the record later widens — again,
precisely what an ordinary expiry means. So the fallback the brief offered — evict only
reads that have not yet contributed — was not needed, and would have been nearly useless
anyway: at the moment the set fills, the only reads that have not yet contributed are the
ones that arrived at the very position being processed.

At the shipping ceiling, on every fixture measured, **eviction never fires**
(`reads_evicted_at_ceiling = 0` everywhere). It is graceful degradation for data deeper than
anything on hand, not a mechanism the ordinary path uses.

**Does it help?** Where the ceiling binds, it moves the count a little and the fairness a
lot. At the deep spot with the ceiling forced back to 4,096:

| | reads refused | reads evicted | positions short | reads missing |
|---|---:|---:|---:|---:|
| refuse at the door | 19,725 | 0 | 351 | 834,832 |
| evict fairly | 6,795 | 13,030 | 346 | 709,752 |

Fair eviction alone does **not** fix the problem — a ceiling below the depth leaves
positions short whichever read it drops. Raising the ceiling is what fixes it. Eviction is
what makes the degradation above the ceiling graceful instead of ordered by arrival time.

---

## 5. What changed in the output, and whether the sample is fair

### Where no cap fires: nothing changed

All four acceptance dumps are byte-identical to the stored copies, and the probe counters
are exact:

```
BYTE-IDENTICAL  gen_chr21       251792 lines
BYTE-IDENTICAL  ssr_chr21         4406 lines
BYTE-IDENTICAL  gen_tom        1718914 lines
BYTE-IDENTICAL  ssr_tom          11945 lines
loci=236081  observations=251786  reads_admitted=54709
```

### Where a cap fires: 2 loci in 1,000 changed their evidence

300× HG002 chr21, `ng_generic_loci_dump`, pristine `6fbbd09` against this branch:

| | |
|---|---|
| emitted loci, either dump | 241,201 |
| **loci whose evidence differs** | **538** |
| loci differing only in chain-id numbering | 75,297 |
| loci present in only one dump | 173 / 163 |
| positions the cap acted on | 909, both sides |

**The raw line diff badly overstates the change and should not be quoted.** It reads 98,922
rows removed and 98,916 added; drop the chain-id column and it collapses to **538 loci**.
Chain ids are minted in admission order and are file-scoped, so 14 fewer admissions over the
whole contig renumber every chain id after them. (The 14 come from capped columns dropping a
different deletion-carrying read, which changes a record footprint, which changes where a
region's walk stops — `record_widen_events` moved 1,640 → 1,634.) **The circulating
"88,351 of 341,094 rows" figure for the ordered-active-set change is likely inflated the
same way, and is worth re-checking before it is used to judge that work.**

Depth barely moves, which is what should happen — a capped position folds exactly `cap`
reads either way:

```
loci in both dumps: 240865
  depth unchanged: 240827
  depth gained:    15
  depth lost:      23
range of change: -5 .. +31
```

### The fairness comparison: the sample is now fair, and the old one was not

At every capped position the walk counts the reads it *saw* and the reads it *kept*, each
split by whether the read began left of the position or exactly at it. A read that begins
exactly at the position is the one an arrival-ordered prefix throws away first.

300× HG002 chr21, over the 909 capped positions:

| | reads seen | of those, beginning **at** the position | reads kept | of those, beginning **at** the position |
|---|---:|---:|---:|---:|
| population | 276,419 | 1,873 | — | — |
| **old rule** (prefix) | 276,419 | 1,873 | 227,250 | **85** |
| **new rule** (sampled) | 276,419 | 1,873 | 227,250 | **1,555** |

A fair draw of 227,250 from 276,419 would keep 1,540 of the 1,873. The new rule kept
**1,555 — within 1 % of fair**. The old rule kept **85, about a twentieth of fair**: of every
100 reads that began exactly at a capped position, the old cap kept about 5 and the new cap
keeps about 83.

The same census over the whole ~130× tomato chromosome, about 3,510 capped positions: the
old rule kept **438** of the 28,054 at-position reads it saw; the new rule keeps **13,179**
of 27,955, against about 14,700 for a fair draw. (The pooled expectation is only
approximate here, because the kept fraction varies a lot between columns of very different
depth; the chr21 figure above, where it does not, is the cleaner one.)

---

## 6. Cost

### Instructions retired, floor-subtracted, min of 3, binaries alternated

`PVC_TRUST_REFERENCE_INDEX=1` throughout; floors measured **per binary per fixture**
(`PVC_PROBE_MAX_LOCI=1`), interleaved with the runs. Wall clock is inadmissible on this
host.

| fixture | floor (base / new) | walk, `6fbbd09` | walk, this branch | change |
|---|---|---:|---:|---:|
| HG002 chr1, 30× | 0.3495 / 0.3498 G | 110.140 G | 110.543 G | **+0.37 %** |
| tomato `SL4.0ch01`, ~130×, 1 M loci | 1.3082 / 1.3087 G | 216.158 G | 216.460 G | **+0.14 %** |
| HG002 chr21, 300× | 1.8944 / 1.8975 G | 89.575 G | 90.002 G | **+0.48 %** |

All three are below this host's ~1 % run-to-run spread but consistent in sign across rounds,
so treat them as a real half-percent rather than as noise.

**None of it is the sampling machinery.** Isolated by comparing a binary that carries every
new counter and the raised ceiling but keeps production's two rules, against the shipping
binary — the difference between them is only the sampling rule and the eviction branch:

| fixture | rules-only change |
|---|---:|
| HG002 chr1, 30× | +0.10 % |
| tomato ~130×, 1 M loci | −0.23 % |

Opposite signs, both inside noise: **the code that decides which reads survive a cap costs
nothing to the reads it never touches.** What the half-percent buys is the per-position
counters and a slightly larger walker state.

### The whole ~130× chromosome, where the ceiling actually binds

| | instructions | peak RSS | loci |
|---|---:|---:|---:|
| STATE 1 — refuse at 4,096, prefix cap (today) | 18.088 T | 386.38 MB | 86,139,952 |
| STATE 2 — refuse at 32,768, prefix cap | 18.524 T | 386.40 MB | 86,139,952 |
| STATE 3 — evict fairly at 32,768, sampled cap (shipping) | 18.493 T | 386.58 MB | 86,139,875 |

**+2.24 % from state 1 to state 3, and essentially all of it is state 1 → state 2.** That
step is the walk doing the work it used to skip: it admits 19,719 more reads, holds up to
10,747 instead of 4,096, and folds 2.64 M contributors at capped positions instead of
1.74 M. It is the price of not discarding evidence, concentrated in one deep region of one
chromosome. State 2 → state 3, the sampling rule, is **−0.17 %**. (Single runs each at this
scale — 13 minutes apiece — so read the sign, not the last digit.)

---

## 7. Tests that moved, and what that says about ng versus production

Validation used, in debug: `cargo test --lib` **2,889 pass** (2,882 on a clean tree, plus
the seven new ones), `cargo test --examples` **all pass**, `cargo clippy --all-targets
--all-features -- -D warnings` **clean**, `cargo doc --no-deps` **12 unresolved links —
exactly the baseline, none of them new**. `cargo test --release` was not run: it is red on
a clean tree.

`cargo fmt` was applied to the changed files only. `6fbbd09` is not itself fmt-clean —
`cigar_cursor.rs` and `open_record.rs` have pre-existing violations, and reformatting them
was reverted so they stay out of this diff.

### Eight tests moved. Four are the change being asserted; four are the boundary.

| test | what happened |
|---|---|
| `walker_vocabulary_tests::the_copied_active_reads_cap_is_still_productions` | **now asserts the two constants differ.** Its own doc said what to do when they were deliberately allowed to differ: *"that test is what says so."* Both values are pinned, so a retune on either side is a decision rather than drift. |
| `generator::tests::the_default_knobs_are_productions_five_constants` | four of the five still read production's constants by name; the fifth is pinned at ng's 32,768 beside production's 4,096. |
| `tests::column_depth_cap_keeps_first_n_of_admission_order` | **renamed** `column_depth_cap_keeps_the_smallest_sampling_keys`. It no longer expects a fixed pair of reads: it works the expected answer out from the reads themselves, which *is* the property. |
| `genome_walk::tests::reads_past_the_active_read_cap_are_shed_and_the_walk_survives` | asserts the arithmetic closes (admitted + refused = seen; evictions = admissions past the ceiling) rather than a literal split, because the split now depends on the reads. |
| `genome_walk::tests::a_slot_freed_by_an_expired_read_admits_the_next_one` | asserts the late read produced its own locus, which is what "the ceiling counts residents, not admissions" actually means and the sampling rule cannot move. |
| `parity::ng_agrees_with_production_where_production_fabricated_nothing` | **skips capped cases and counts them**, asserting the count is at least a tenth of all cases. |
| `parity::every_divergence_from_production_is_one_of_the_six_named_classes` | same skip, plus class 3 is now measured off ng's walk alone on the skipped cases. |
| `parity::the_generator_exercises_what_the_port_can_break` | unchanged and still green — the generator still draws tiny caps one case in three. |

Seven new tests carry what the parity exclusion gives up, plus the two properties the whole
change exists for:

- `read_sampling::tests::*` (four) — the key ignores alignment start, the two mates draw
  separately, the hash values are pinned, and 10,000 sequencer-style names spread uniformly.
- `genome_walk::tests::no_position_is_short_of_the_cap_while_reads_covering_it_exist` — and
  it also asserts the counter *fires* when the ceiling is put below the cap, because a
  counter that reads zero whatever the configuration is measuring nothing.
- `genome_walk::tests::the_kept_set_does_not_depend_on_the_order_the_set_holds_reads_in` —
  the same reads offered in two arrival orders produce the same loci.
- `genome_walk::tests::the_ceiling_keeps_the_smallest_sampling_keys` — with a fixture whose
  fair answer is deliberately not the first-come answer.

### Is ng's divergence from production still describable by its six named classes?

**Yes on an uncapped column, and on a capped one the honest answer is that it is no longer a
divergence at all — it is two walkers seeing different evidence.** The six classes describe
two walkers that saw the same reads and rendered them differently. At a capped column ng
keeps the `cap` smallest sampling keys and production keeps a prefix, so the bases, the
counts, the footprint a dropped deletion would have widened, and even how many records exist
all differ. Filing that as a seventh class would give every class an escape hatch: a real
transcription slip could hide inside a class that has to excuse an arbitrary difference.

So the differential states the boundary instead of blurring it, and states it in code —
`Case::caps_can_fire`. The excluded cases are counted and the count is asserted non-trivial,
so a generator that quietly stopped drawing tiny caps would fail rather than turn the
exclusion into an exclusion of nothing. **One real piece of coverage had to be rescued:**
class 3 (a per-locus counter production has no counterpart to) fires on
`reads_discarded_by_cap`, which only a capped column fills — so on the skipped cases ng's
walk is still run, alone, and the class is counted off it. What is no longer checked there is
a comparison, because there is none to make.

---

## 8. What this makes fragile, and the four changes landing separately

Built on `6fbbd09` alone; none of the four in-flight performance changes are merged. Three
of them are touched by a larger hold ceiling, and here is what this implies for each —
**stated, not built.**

- **The scalar fast column.** Its gate is `if active_reads.len() > max_snp_column_depth ||
  active_reads.is_empty() { Fallback }` — the *active set's* length against the 8,000-read
  SNP cap, deliberately conservative from above. While the walk held 4,096 reads that test
  could never reject anything. It now can: on the tomato deep spot the set holds 10,747, so
  the gate falls back for the first time, over a few thousand positions of one chromosome.
  That is the right behaviour and costs nothing measurable, but the comment explaining the
  bound as free should be re-read against a walk that can exceed it, and the census figure
  it quotes (1,203 columns lost to this test at ~130×) was taken under the old ceiling and
  will move.
- **The fold table as an ordered `Vec`** and **the mate-overlap min-heap on the active set.**
  Both scale with the number of held reads. The heap holds one entry per admitted mate pair
  and is popped by position, so a set eight times larger makes it eight times larger in the
  worst case; the fold table's `locate()` is a binary search whose depth grows with the log
  of the held count, so 10,747 reads costs about 1.4 comparisons more than 4,096. Neither
  looks like a problem; neither has been measured against a walk that holds this many.
- **The ordered active set** is the change this work exists to decouple from. Under the old
  prefix rule, restoring admission order changed which reads a capped position keeps —
  88,351 rows on a 300× dump (*cited*, and see §5 on why that figure is probably inflated).
  Under the sampling rule the kept set is a function of the reads alone, so the coupling is
  **removed rather than re-baselined**:
  `the_kept_set_does_not_depend_on_the_order_the_set_holds_reads_in` is the test that says so.

Two things this branch makes fragile in its own right:

- **The sampling hash is now a format.** It decides which reads reach a `.psp` at every
  capped position. Changing it — including swapping it for a library hash — changes the
  output. The pinned literal values in `read_sampling.rs` are the tripwire.
- **`MAX_PENDING_MATES` at 1,000,000 is untested at its own ceiling.** The measured peak is
  11,384, so there is a factor of 88 of headroom, but the constant is now large enough that
  reaching it means the map has genuinely grown to a million entries — tens of megabytes —
  before anything says so. `pending_mates_high_water` is reported on every run, so the
  headroom is a number somebody has read rather than one somebody has reasoned about.

Also changed, because it had become untrue: the allocator's one-shot high-water warning said
*"the run will fail with ActiveReadsExhausted"*. ng has not failed that way since the walk
started giving reads back instead of aborting; the message would have told an operator to
expect a crash and hidden the thing that does happen. It now names
`positions_short_of_cap`, which is the number to act on.

`copy_fidelity.rs` releases `chain_id_allocator.rs`, leaving `decompose.rs` as the last file
still byte-for-byte production's. The release table records what moved and why.

---

## 9. What could not be done

- **The 12,792-read position was not reproduced.** The brief cites 12,792 reads over one
  position near 33,037,565 on `SL4.0ch01`. With the ceiling out of the way the deepest column
  measured there is **10,744** contributors, from an active set peaking at 10,747. The
  discrepancy is not explained here; 10,744 is what this walk measures after read preparation,
  and it is the number the ceiling was chosen against.
- **Reaching the deep spot needed a new knob.** No `PVC_PROBE_MAX_LOCI` prefix gets near
  33 Mb, and a whole-contig walk is 13 minutes. `PVC_PROBE_FROM_BP=n` drops every region
  ending before position `n`. It does **not** skip the reads in front of it, so a run with it
  set is not comparable with one without on time or instructions — only on the counters at
  the place it reaches, which is what it was used for.
