# Sketch 4 — a block-filling generator on one side of the cohort merge

**Status:** throwaway experiment, complete. Code is a diff, not a branch to merge.
**Date:** 2026-08-06. **Tree:** `1e5ffa8`, with sketch 1's producer diff reapplied.
**Plan:** `doc/devel/ng/impl_plan/locus_stream_shape_experiments.md`, §10's open question.

---

## The answer, in one sentence

**Yes — a block-filling generator makes the in-memory fold pay: it cuts the merge from
1.642 to 0.239 × 10⁹ instructions, a factor of 6.9, worth 5.3 % of the whole path; but nine
tenths of that saving comes from the block carrying a four-byte per-locus keep column and
only one tenth from folding it, and against today's shipped records-everywhere the merge
contributes one twentieth of a 16.1 % end-to-end gap that is otherwise all producer.**

---

## Vocabulary, defined once

- A **locus** is one covered reference base of one sample — what the generator emits.
- The **merge** is the k-way, position-by-position join across samples. A position is
  **variable** if the maximum, over samples, of that sample's summed non-reference
  observations reaches 2. Every arm applies the identical rule.
- A **record** is today's owned `SampleLocusObservations`: a boxed slice of reference bases,
  a `Vec` of observations, and per observation an allele string and a `Vec` of chain ids.
- A **block** is many loci with each field held as one array — sketch 1's `LocusBlock`.
- The **keep column** is the block's new `locus_nonref_obs`: one `u32` per locus, the summed
  `num_obs` of the observations whose bases differ from that locus's reference bases. It is
  what this sketch had to add, and §3 is about what it cost.

---

## 1. The five states, and why there are five

The plan asks for four. I built five, because the fourth confounds two things.

| state | producer | merge | mode name |
|---|---|---|---|
| **A** | records — today's shipped generator | one owned record of lookahead per sample, O(N) head scan | `rec-records` |
| **B** | records | per-sample owned queues, a light column derived per record, folded over a window | `rec-fold` |
| **C** | blocks | every locus refilled into a per-sample scratch record, then the same O(N) head scan | `blk-records` |
| **D** | blocks, keep column on | the block's keep column folded across samples; only variable positions read the heavy arrays | `blk-fold` |
| **E** | blocks, keep column on | O(N) head scan again, but the **keep rule reads the keep column**, so a locus becomes a record only when its position is kept | `blk-records-sum` |

**E is the arm that separates two claims D confounds:** *folding a cheap column across
samples* and *the block carrying a cheap column at all*. Without E, the whole of C → D would
have been reported as the fold's, and 90 % of it is not.

Each producer also has a **floor** — the same walk with every locus dropped and no merge at
all — so the merge's own contribution is a difference, not a share of a large number:

| floor | what it is | mode name |
|---|---|---|
| record floor | every generator driven to exhaustion, every record dropped | `rec-drop` |
| block floor | every generator driven to exhaustion, every block cleared unread, **no keep column** | `blk-drop` |
| block floor + column | the same, **keep column maintained** | `blk-drop-sum` |

`blk-drop-sum − blk-drop` is the keep column's price and nothing else.

**A and B are re-measurements, not citations.** Different worktree, different binaries. They
land 1.8 % above sketch 2's ng-fed figures for a reason §6 gives and quantifies.

---

## 2. The result

10 tomato bench CRAMs, 300 regions of `ssr_regions.bed`, 2,830,932 loci, 309,018 covered
positions, **3,060 of them variable — 1 position in 101**. Single-threaded throughout.
Instructions are `instructions retired` from `/usr/bin/time -l`, **minimum of five runs**,
all eight modes alternated inside one script; run-to-run spread was 0.07–0.20 %. All numbers
**measured** unless marked cited.

| state | instructions, whole path | vs A | **the merge alone** (state − its floor) | merge's share of the path | peak RSS | loci ever materialised | bytes copied at the handoff | lines of source |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **A** records + records | 27.034 × 10⁹ | — | **0.448 × 10⁹** | 1.66 % | 260.1 MB | 2,830,932 records | 282.0 MB | **87** |
| **B** records + fold | 27.007 × 10⁹ | −0.10 % | **0.420 × 10⁹** | 1.56 % | 268.5 MB | 2,830,932 records | 282.0 MB | 158 |
| **C** blocks + records | 23.932 × 10⁹ | −11.5 % | **1.642 × 10⁹** | 6.86 % | 262.3 MB | 2,830,932 scratch refills | 350.2 MB | 190 |
| **E** blocks + records, keep column | 22.820 × 10⁹ | −15.6 % | **0.384 × 10⁹** | 1.68 % | 262.5 MB | **28,718** scratch refills | 240.0 MB | 202 |
| **D** blocks + fold | **22.674 × 10⁹** | **−16.1 %** | **0.239 × 10⁹** | **1.05 %** | 262.4 MB | **0** | 238.2 MB | 165 |

Every state materialises the same 28,718 merged-evidence objects (3,060 positions × the
samples covering each), from pooled storage, so that column is omitted.

### Read three things off this

**With a record producer the fold is still worth nothing, exactly as sketch 2 found.**
0.420 × 10⁹ against 0.448 × 10⁹ is a 6 % difference against a run-to-run spread of
±0.035 × 10⁹ on the totals — the two are indistinguishable, and which is nominally ahead
flips between the `full` and `digest` consumers. Sketch 2's finding reproduces.

**With a block producer the fold is worth a factor of 6.9 at the merge** — 0.239 × 10⁹
against C's 1.642 × 10⁹ — and 5.3 % of the whole path (22.674 against 23.932 × 10⁹). This is
the answer to the plan's question, and it is a real number.

**But the mechanism is not the fold.** Splitting C → D at E:

| step | what changes | instructions at the merge | share of C → D |
|---|---|---:|---:|
| C → E | the block carries a keep column, so 100 loci in 101 never become records | 1.642 → 0.384 × 10⁹ | **1.258 × 10⁹ = 90 %** |
| E → D | the per-position O(N) head scan becomes a dense window fold | 0.384 → 0.239 × 10⁹ | 0.146 × 10⁹ = 10 % |

**Nine tenths of the prize is a four-byte column, not a columnar merge.** E keeps the merge
record-shaped — the same O(N) head scan as A, reading records — and still banks 90 % of it.

### The producer, not the merge, is where the path changed

State A against state D is 4.360 × 10⁹ instructions. Attributed by the floors:

| | |
|---|---:|
| producer: record floor → block floor with keep column (26.586 → 22.435 × 10⁹) | 4.151 × 10⁹ — **95.2 %** |
| merge: 0.448 → 0.239 × 10⁹ | 0.209 × 10⁹ — **4.8 %** |

So the merge's contribution to the combined win is **one twentieth**. The merge does get
47 % cheaper; it is just small.

---

## 3. What the block had to grow, and what it cost — the design consequence

**The block as sketch 1 built it cannot be folded.** Its per-locus arrays are contig, start,
end, reads-without-observation, reads-discarded-by-cap, and two offset arrays. Nothing in
that set answers *"is this locus worth materialising"*. Answering it means, for every
observation of the locus, slicing `obs_bases` through its offsets and comparing against the
locus's reference bases — which touches the heaviest array in the block and is exactly the
work the fold exists to avoid.

Sketch 2's `.psp` fold never met this, because a `.psp` file stores its alleles
reference-first and hands the non-reference observation count over as a column it already
has. **ng's observations carry raw bases and no reference marker.** So:

> **The merge's needs shape the generator's output layout.** A block that is only a
> producer-side implementation detail is not enough; the block must carry a summary the
> merge can read, and that summary has to be *decided by the merge and computed by the
> producer*.

**What it cost to compute — measured.** The keep column is accumulated inside
`LocusBlock::push_observation`: compare the observation's bases against the open locus's
reference bases, and on a difference add `num_obs` into the locus's slot.

| | |
|---|---:|
| block floor without the column | 22.291 × 10⁹ |
| block floor with the column | 22.435 × 10⁹ |
| **the column costs** | **0.145 × 10⁹ = 0.65 % of the block walk** |
| per locus (2,830,932 loci) | **51 instructions** |
| per observation (2,858,788 observations) | **51 instructions** |
| allocations it adds (dhat, 60 regions) | **100**, out of 3.28 million |

**51 instructions for a one-byte comparison is an upper bound, and a loose one.** This
implementation re-derives something the walk already knows: it rebuilds the reference slice
from an offset it stored two lines earlier, calls `last_mut()` on a `Vec` with a bounds
check, and does a slice comparison. A producer-side implementation would pass a
`matches_reference` flag the emit site already has. I did not build that, because it would
have meant editing the walk's hot files and this sketch's whole point is what the *block*
must carry. **Read the 51 as "the naive version costs a fifth of what the fold saves at the
merge, and it can be made cheaper."** Net of it, D still beats C by 1.258 × 10⁹.

**What it cost to carry — measured.** Four bytes per locus: block payload rose from
225,482,113 to 236,805,841 bytes over 2,830,932 loci, which is 4.000 bytes per locus to the
third decimal. Peak resident memory rose 0.2 MB (262.1 → 262.3 MB), 0.08 %.

**What it cost in source.** 29 code lines, all inside `block.rs`. Nothing in
`fast_column.rs`, `open_record.rs` or `genome_walk.rs` — the walk's three hottest files —
changed by one line for it. That is the one part of this that came out better than expected.

---

## 4. Memory, allocations and bytes moved

**Peak resident memory — there is a small loss, and it is in the producer, not the merge.**
Minimum of five, `/usr/bin/time -l`:

| state | peak RSS | vs A |
|---|---:|---:|
| A records + records | 260.1 MB | — |
| B records + fold | 268.5 MB | **+3.2 %** |
| C blocks + records | 262.3 MB | +0.8 % |
| E blocks + records, keep column | 262.5 MB | +0.9 % |
| D blocks + fold | 262.4 MB | **+0.9 %** |

Two readings. **The block producer costs 2.2 MB more resident than the record producer**
(259.9 against 262.1 MB at the floors) — ten samples each holding one block whose payload
high-water is 164 kB, plus array slack. That matches sketch 1's finding that the columnar
arms are 0.2–0.4 % *higher*, in the same direction and a little larger here.

**But the block fold does not pay the record fold's buffer cost.** State B needs per-sample
owned queues of records and costs 8.4 MB over state A; state D needs nothing, because the
block *is* the buffer, and costs 0.1 MB over its own floor. So the fold's memory penalty —
which sketch 2 measured at 4.8 % of peak heap on `.psp` and 10 MB of RSS on ng records —
**disappears when the producer emits blocks.** D is 6.1 MB *below* B.

**Allocations — dhat, `--features dhat-heap --target-dir target-dhat`, 60 regions:**

| state | allocations | total churn | peak heap |
|---|---:|---:|---:|
| record floor | 5,126,114 | 1,028.1 MB | 219.601 MB |
| A | 5,126,247 (+133) | 1,028.2 MB | 219.601 MB |
| B | 5,126,565 (+451) | 1,033.6 MB | 219.602 MB |
| block floor | 3,279,019 | 911.6 MB | 219.602 MB |
| block floor + keep column | 3,279,119 (+100) | 911.7 MB | 219.602 MB |
| C | 3,279,293 (+274) | 911.6 MB | 219.603 MB |
| E | 3,279,393 (+274) | 911.8 MB | 219.603 MB |
| D | 3,279,386 (+267) | 911.8 MB | 219.602 MB |

**The record floor, A and B reproduce sketch 2's ng allocation counts to the unit** —
5,126,114, +133 and +451. That is the cross-experiment check that this harness is sketch 2's
harness.

**The block producer removes 36 % of the process's allocations** (5,126,114 → 3,279,019) and
11 % of its churn. **Every merge shape contributes between 133 and 451 allocations out of
3.3 million.** As in sketch 2, there is nothing at the merge for a layout change to remove;
what removes allocations is the producer.

**Peak heap is 219.60 MB in all eight modes**, identical to five decimal places. It is the
reference window and the read decoder; neither the block nor any merge shape is visible in
it.

**Bytes copied at the handoff** — measured, 300 regions. The two producer figures use
different fixed-part conventions (the record figure counts 96 bytes per observation of owned
struct; the block figure counts the block's own 28-per-locus / 49-per-observation payload
accounting), so compare them as magnitudes, not to the byte:

| state | producer → consumer | consumer's own copy | merged evidence | total |
|---|---:|---:|---:|---:|
| A, B | 280.6 MB into records | 0 | 1.42 MB | **282.0 MB** |
| C | 225.5 MB into blocks | 123.3 MB into scratch records | 1.42 MB | **350.2 MB** |
| E | 236.8 MB into blocks | 1.77 MB into scratch records | 1.42 MB | **240.0 MB** |
| D | 236.8 MB into blocks | 0 | 1.42 MB | **238.2 MB** |

**State C is the worst state on every axis**, and that is the honest shape of "blocks in the
producer, records at the merge, no keep column": it writes the block and then copies all of
it back out, one locus at a time, to answer a question the block could have answered in four
bytes. **E's 1.77 MB against C's 123.3 MB is the same 1-in-101 skew, in bytes.**

---

## 5. Correctness

**All five states, both floors, and the accounting mode agree bit for bit.** The digest is
sketch 2's: FNV-1a over the raw bit patterns of the merged evidence — position, contributing
samples, and per allele the sequence bytes, observation count, `q_sum` bits, forward count,
placed-left, MAPQ sum and sum-of-squares, and the chain-id list — then a
genotype-likelihood + allele-frequency EM whose log-likelihood bits are folded in too.

```
rec-records  rec-fold  blk-records  blk-records-sum  blk-fold
  digest = 0xf46f8f924ee468aa   in all five
  loglik_acc = -56602.839011    in all five
  loci_produced = 2,830,932   positions_seen = 309,018   positions_kept = 3,060
  merge_objects = 28,718      merge_bytes = 1,418,502
```

**Tolerance allowed: none, and none was needed.** The merge performs no floating-point
arithmetic and the EM sees identical inputs in identical order.

**And the digest is sketch 2's own value.** Sketch 2's ng-fed arms reported
`0xf46f8f924ee468aa` and `-56602.839011` on the same fixture in a different worktree with
different code. Three of the four states here reach it by a route sketch 2 did not have.

**The digest is block-boundary invariant.** Across block budgets from 4 KiB to 1 MiB — 58,702
blocks down to 3,000 — every merging mode produces the same digest (§7).

**The shipped walk is where sketch 1 left it.** `cargo test --lib`: **2,893 passed, 1 failed**,
the failure being `parity::every_divergence_from_production_is_one_of_the_six_named_classes`
— the accepted clean-tree divergence, and the exact baseline stated in the brief. (`cargo
test --lib --release` reports 2,884 / 10 because ten of those tests assert on debug
assertions; that is not a regression and not a comparison anyone should make.)

---

## 6. The control — what arm A is carrying that sketch 2's arm A was not

Sketch 1's producer installs a runtime `LocusSink` enum in the walk. Arm A goes through it
even though it always chooses records. **That is not free, and both A and B are measured
with it.** To size it I reverted `src/` to clean `1e5ffa8`, rebuilt sketch 2's own binary in
this worktree, and measured the same fixture:

| | clean `1e5ffa8` library | this sketch's library | tax |
|---|---:|---:|---:|
| record floor (`walk` / `rec-drop`) | 26.114 × 10⁹ | 26.586 × 10⁹ | **+1.81 %** |
| records (A) | 26.552 × 10⁹ | 27.034 × 10⁹ | **+1.81 %** |
| fold (B) | 26.514 × 10⁹ | 27.007 × 10⁹ | +1.86 % |

Two things follow.

**Sketch 2's ng-fed measurements reproduce here to within 0.1 %** — its `walk` 26.128, its
`records` 26.565, its `fold` 26.542 × 10⁹ against 26.114, 26.552, 26.514 measured now. The
merge shares reproduce too: 0.438 and 0.400 × 10⁹ here against sketch 2's 0.437 and 0.414.
Nothing about this host or fixture has drifted.

**The four-state table overstates the block producer's win by 1.8 points.** Against a clean
record path (26.552 × 10⁹), state D's 22.674 × 10⁹ is **−14.6 %**, not −16.1 %. Both sides
carry the branch — a blocks-only build would be slightly cheaper than 22.674 too — so 14.6 %
is a floor and 16.1 % a ceiling on the same quantity. **The real production version has no
branch at all**, because a shipped pipeline picks one shape at compile time; the enum is a
sketch artefact and sketch 1 said so.

I report the taxed table as primary because the brief asks for the two producers compared
inside one experiment, which requires one binary.

---

## 7. Two sweeps that test whether the answer is an artefact of the fixture

### Sample count — the fold's prize does **not** grow with it

The brief expected it to, since the fold's saving is supposed to come from many samples
agreeing a position is invariant. **It does not, on this fixture.** Minimum of three, 300
regions:

| samples | merge alone, A | merge alone, D | A ÷ D | merge alone, C | C ÷ D | end to end, D vs A | positions kept |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 0.121 × 10⁹ | 0.062 × 10⁹ | 1.96× | 0.480 × 10⁹ | 7.8× | −15.5 % | 1,615 |
| 10 | 0.421 × 10⁹ | 0.233 × 10⁹ | 1.80× | 1.640 × 10⁹ | 7.0× | −16.1 % | 3,060 |
| 20 | 0.910 × 10⁹ | 0.538 × 10⁹ | 1.69× | 3.407 × 10⁹ | 6.3× | −15.9 % | 3,844 |
| 40 | 1.898 × 10⁹ | 1.099 × 10⁹ | 1.73× | 6.729 × 10⁹ | 6.1× | −16.2 % | 4,359 |

**The end-to-end gap is flat at 16 % from 3 to 40 samples**, and the fold's advantage over the
record merge stays between 1.7× and 2.0× with no trend. What *does* move is C ÷ D, and it
moves the wrong way: **7.8× at three samples down to 6.1× at forty.** The reason is in the
last column — going
from 3 samples to 40 raises kept positions only from 1,615 to 4,359 while raising the loci
walked 13-fold, so the *fraction* kept falls, but the **evidence objects each kept position
demands rise linearly with the cohort** (4,635 → 160,926). The materialisation the keep
column lets you skip grows slower than the materialisation you cannot skip.

**So the fold's prize is set by how few positions are variable, not by how many samples
agree.** At 1 position in 101 it is what it is; a cohort ten times larger does not improve it.

### Block size — the boundary constraint costs 0.3 %, and only when blocks are tiny

At the default 256 KiB budget this fixture's ~1 kb regions each fit in one block, so every
sample holds the whole region and sketch 2's watermark rule — *the merge advances only to
the minimum, over samples, of the last position in each currently-held block* — is never
exercised. Shrinking the budget forces a region across many blocks. Minimum of three:

| block budget | blocks | D `blk-fold` | D's merge alone | E `blk-records-sum` | C `blk-records` |
|---:|---:|---:|---:|---:|---:|
| 4 KiB | 58,702 | 22.736 × 10⁹ | 0.285 × 10⁹ | 22.842 × 10⁹ | 23.947 × 10⁹ |
| 16 KiB | 16,348 | 22.678 × 10⁹ | 0.243 × 10⁹ | 22.839 × 10⁹ | 23.961 × 10⁹ |
| 64 KiB | 5,682 | 22.666 × 10⁹ | 0.236 × 10⁹ | 22.834 × 10⁹ | 23.945 × 10⁹ |
| 256 KiB | 3,000 | 22.661 × 10⁹ | 0.221 × 10⁹ | 22.855 × 10⁹ | 23.941 × 10⁹ |
| 1 MiB | 3,000 | 22.680 × 10⁹ | 0.235 × 10⁹ | 22.849 × 10⁹ | 23.936 × 10⁹ |

At 4 KiB each region spans about twenty blocks per sample and the watermark advances twenty
times instead of once. **The fold's merge cost rises 29 % (0.221 → 0.285 × 10⁹) and the whole
path rises 0.33 %.** Everything from 16 KiB up is flat, matching sketch 1's finding that
instructions are flat across block size — with the addition that **the cohort fold has a
weak floor the single-sample pre-pass did not**, and it is below 16 KiB.

**Every block size produces the same digest**, which is the property that matters: the merge
is not sensitive to where the producer chose to cut.

---

## 8. What the fixture is, and where the numbers do not transfer

| | |
|---|---:|
| loci | 2,830,932 |
| observations | 2,858,788 (1.01 per locus) |
| reads at a position, `Σ num_obs` | 12,884,986 |
| **depth** | **4.55 reads at a position per covered base** |
| covered positions across the cohort | 309,018 |
| variable positions | 3,060 — **1 in 101** |

**This is a shallow, tandem-repeat-targeted fixture, and that inflates the producer's share
of the win.** Sketch 1 measured the block layout worth **−12.1 % at 4.9×, −7.8 % at 20.7×,
−2.0 % at 98×** (cited), because the saving is a fixed ~1,344 instructions per covered base
and what changes with depth is the denominator. At 4.55× my 15.6 % producer saving sits
right where sketch 1's depth curve says it should, and it must not be read as the number at
the 30× whole-genome target.

**The merge's own numbers transfer better than the producer's**, because the merge's cost
scales with kept positions and cohort size, not depth. But its *share* of the path would
shrink at 30×, since the denominator grows and the merge does not. If sketch 1's 7.8 % at 30×
is the producer's real saving there, then a merge contributing 0.209 × 10⁹ of a 4.360 × 10⁹
combined gap here would contribute a larger fraction of a smaller gap — that is a
plausible-sounding extrapolation and **I did not measure it**, so it is not a number.

The 1-in-101 variable rate is not directly comparable with the "4 positions in 100" the plan
cites for production's grouped keep rule, nor with sketch 2's 11.8 % on its 50-sample `.psp`
cohort: this is a simplified per-position rule on a different cohort.

---

## 9. What building it was actually like

**Half a day, and the surprise was which arm was hard.**

**Getting to a starting line took about forty minutes and nothing fought back.** Sketch 1's
producer diff applied clean to `1e5ffa8`; sketch 2's ng harness applied clean beside it; both
were built and reproducing their own reported figures within the hour. That is a good deal
better than these things usually go, and it is a credit to how the two earlier sketches
packaged their code.

**The design point hit in the first ten minutes, exactly where the brief predicted.** I
opened `sketch1_src_block.rs` to find the column to fold and there is none — the block's
per-locus arrays are positions and offsets, and the only way to ask *"is this locus
variable"* is to walk `obs_bases` through two offset arrays and compare. I had the answer to
the sketch's most interesting question before I had written a line: **the block as sketch 1
designed it cannot be folded, and the merge's requirement propagates into the producer's
output layout.**

**Adding the column was the easiest 29 lines in the sketch**, and I expected it to be the
hardest. It lives entirely in `block.rs`: remember where the open locus's reference bases
start, compare each observation's bases against them, accumulate. No hot walk file changed.
I had braced for the "one more field, kept in step by hand" tax sketch 1 reported for
`finalise_into_columns`, and it did not arrive — because the column is derived from data the
block already holds, not pushed through the walk.

**Arm D needed one piece of reasoning I would not want a maintainer to have to re-derive.**
The round structure is sketch 2's — load mutably, fold and materialise from shared views,
drop the views, advance mutably — and it works here for the same reason: N distinct
generators can all be borrowed shared at once because each borrow is of a different
generator. But refilling is only legal for a sample whose block is fully consumed, and a
sample holding a partial block cannot be refilled, so the loop can stall. It does not stall,
and the reason is not obvious: **the sample that sets the watermark is always the one that
gets fully consumed**, because the watermark is the minimum over samples of each block's last
position. I convinced myself of that on paper before writing the loop. Had I not, it would
have deadlocked on the first multi-block fixture — which, as §7 says, the default block size
never produces, so it would have gone unnoticed until the block sweep.

**It needed one new library line too**: `PileupGenerator::block_sink()`, thirteen lines,
because sketch 1 only exposes `block_sink_mut()` and the fold needs N shared borrows at one
moment. Small, but it is the shape of the thing — **the merge's borrowing pattern is a
requirement on the generator's interface, not just on its layout.**

**Arm E was the genuine surprise and took thirty minutes.** I built it only to attribute
D's win, expecting it to lose badly. It banks 90 % of the fold's saving, needs **no round
structure at all** — it borrows one sample's block at a time, so there is no moment where
every sample's block must be held — and its merge loop is the same O(N) head scan as arm A.
The one thing it wants from the block is three per-locus array reads. **If the block carries
the keep column, the merge does not have to become columnar to collect nearly all of the
prize**, which is sketch 1's *"no consumer ever needs to see columns"* arriving a second time
from a different direction.

**Arm C was easy to write and is simply the wrong shape.** It compiled first time and it is
the state the architecture would land in by accident: blocks in the producer because the
producer wanted them, records at the merge because records read better, and no column between
them — so every locus is copied out of the block to answer a question worth four bytes. It
is the most expensive of the five on instructions at the merge, on bytes copied, and on
objects materialised.

**The digests caught nothing, and I would not have trusted a number without them.** Every arm
agreed on its first run. That is a statement about how carefully the arms were written, not
about how hard they are to get wrong — arm D indexes three per-locus arrays and two
per-observation index spaces, all `usize`, and nothing in the type system distinguishes them.
Sketch 1 said the same about arm B, and having now written the cohort version of it, I would
put it more strongly: **the fold's correctness rests entirely on a digest that a production
pipeline would not have.**

**Lines of source.** Code lines only, blanks and comments stripped, counting each arm's merge
loop plus the helpers only it uses:

| | consumer lines | producer lines |
|---|---:|---:|
| A records + records | **87** | 0 |
| B records + fold | 158 | 0 |
| C blocks + records | 190 | 830 (sketch 1) + 13 (`block_sink`) |
| E blocks + records, keep column | 202 | 843 + **29** (keep column) |
| D blocks + fold | 165 | 843 + 29 |

A is 87 lines and every one of them says what it means. D is 165 — less than twice A, and
much less than the 3.9× sketch 2's `.psp` fold cost, because the block is denser and better
behaved than a `.psp` block: one locus per covered base means a dense window instead of a
sorted k-way merge, and a cursor into a block instead of a queue that has to be drained.
**The fold got cheaper to write when the producer got columnar.** That is a real and slightly
counter-intuitive result, and it is worth more than the 10 % of instructions E → D buys.

---

## 10. The decision this leaves the owner

The open question in §10 of the plan is whether a block-filling generator is worth **830
lines in the walk's three hottest files** for **7.8 % of a covered base at 30×**. This sketch
was asked whether the merge changes that trade.

**It changes it by 4.8 %.** Of the 4.360 × 10⁹ instructions between records-everywhere and
blocks-plus-fold on this fixture, 4.151 × 10⁹ is the producer and **0.209 × 10⁹ is the
merge**. The merge saving is real, reproducible and 6.9-fold in its own terms — and it is one
twentieth of the combined figure. **If 830 lines was not worth 7.8 % before, the merge does
not make it worth it.**

**Three things do change, and none of them is the instruction count.**

1. **The block must carry a keep column.** This is not optional and it is not a
   producer-side detail: the merge decides what the column is, the producer computes it, and
   without it the block is *worse* than records at the merge (state C, 11.5 % better than A
   end to end but with the merge itself 3.7× more expensive than A's). Cost: 29 lines in
   `block.rs`, 4 bytes per locus, 51 instructions per locus in a naive form that could be
   made much cheaper. **If blocks are adopted, this must be in the design from the start**,
   because retrofitting it is how state C happens.

2. **The merge does not have to become columnar.** Arm E — record-shaped, O(N) head scan,
   reading only the keep column — collects 90 % of the fold's saving in 202 lines with no
   round structure and no shared-borrow phase. Sketch 1 concluded that no consumer needs to
   see columns; the cohort merge is a second, independent instance of it. **The remaining
   10 % is 0.6 % of the whole path and costs the round structure, the watermark, the dense
   window, and the block-size floor below 16 KiB.** I would not spend it.

   **E is 202 lines against D's 165, and I am still recommending E.** 66 of E's lines are
   `LocusScratch` — the buffer that exists purely so the merge stays record-shaped. That is
   the plan's §2 preference bought explicitly, in one place, in code that copies thirteen
   fields and does nothing clever. D's 165 lines are 37 fewer and nearly all of them are
   cursor and window bookkeeping over untyped `usize` indices. **Fewer lines is not the same
   as less to get wrong**, and this is a case where the two point opposite ways.

3. **The fold's memory penalty goes away.** Sketch 2 measured the fold arm 4.8 % above the
   record arms on `.psp` peak heap and 10 MB above on ng RSS, and called it a loss. With a
   block producer it is 0.1 MB, because the block is already the buffer the fold needed. What
   remains is the block producer's own 2.2 MB (0.9 %), which is sketch 1's number and points
   the same way: **peak memory is not an argument for or against this, in either direction.**

**My opinion, labelled as one.** The measured case for blocks did not improve enough to move
me. What did move is the *shape* of the thing: the two states worth having are A (records
everywhere, 87 lines, nothing to get wrong) and E (blocks with a keep column, a record-shaped
merge, no columnar consumer anywhere), and **D is not worth the round structure over E**. If
the pre-pass has to fit inside the walk's current budget — which §10 says is the real
question — then the answer is E, not D, and the extra thing to specify is the keep column,
not a columnar merge.

**What this sketch does not settle.** Whether 7.8 % at 30× is worth 830 lines. That was
sketch 1's question, it is still open, and the merge has now been measured well enough to say
it does not decide it.

---

## Files

Beside this report:

| file | what it is |
|---|---|
| `sketch4_code.diff` | the whole sketch against `1e5ffa8` — sketch 1's producer, the keep column, `block_sink()`, and the eight-mode harness |
| `sketch4_src_block.rs` | `LocusBlock` / `LocusSink` with the keep column |
| `sketch4_example_producer_merge.rs` | the five states, two floors and the accounting mode |
| `sketch4_instructions_raw.txt` | main sweep: 5 reps × 2 consumers × 8 modes |
| `sketch4_control_clean_tree_raw.txt` | §6's control: sketch 2's binary on a reverted library |
| `sketch4_nsweep_raw.txt` | §7's sample-count sweep, N = 3, 10, 20, 40 |
| `sketch4_blocksweep_raw.txt` | §7's block-budget sweep, 4 KiB to 1 MiB, with per-run digests |
| `sketch4_kept_by_n_raw.txt` | variable positions and evidence objects by sample count |
| `sketch4_counts_raw.txt` | the accounting run — bytes, objects, digests. Not a timing run |
| `sketch4_dhat_raw.txt`, `sketch4_dhat_*.json` | heap profiles, 60 regions, all eight modes |
| `sketch4_reduced.txt` | every figure quoted above, reduced from the raw files |
| `sketch4_scripts/` | every script that produced the above |

To reproduce any state:

```
PVC_TRUST_REFERENCE_INDEX=1 /usr/bin/time -l \
  target/release/examples/sketch4_producer_merge blk-fold full \
    $HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa \
    benchmarks/ssr_tomato1/ssr_regions.bed 300 \
    benchmarks/ssr_tomato1/crams/*.bench.cram
```

Modes: `rec-drop`, `rec-records`, `rec-fold`, `blk-drop`, `blk-drop-sum`, `blk-records`,
`blk-records-sum`, `blk-fold`, `stats`. `PVC_SKETCH4_BLOCK_KB` sets the block budget.

---

## Instrument and its limits

- **`instructions retired` from `/usr/bin/time -l`**, minimum of five runs, all modes
  alternated inside one script, single-threaded throughout — no rayon, no threads, in any
  mode. Run-to-run spread 0.07–0.20 %.
- **The merge's own cost is an ablation**, never a share: each state minus its own floor,
  where the floor is the same walk with every locus dropped.
- **No wall clock is quoted and none informs any conclusion.** The field is recorded in the
  raw files and ignored.
- **No sampling profile was taken and none is quoted.** Every attribution here is an
  instruction-count ablation or a static allocation tally from dhat, both per-process.
- **Peak resident memory** from the same `/usr/bin/time -l`.
- **`PVC_TRUST_REFERENCE_INDEX=1` throughout.** No run here is compared against one that
  verified the FASTA.
- `--features dhat-heap` always with `--target-dir target-dhat`.
- Cited, not measured here: sketch 1's depth curve (−12.1 % / −7.8 % / −2.0 % at 4.9× /
  20.7× / 98×), its 1,344-instruction per-base term and 830-line producer count, sketch 2's
  ng-fed figures (which §6 re-measures anyway), and the plan's "4 positions in 100".
  Everything else is measured in this worktree.
