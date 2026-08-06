# Sketch 2 — the k-way cohort merge: does it want columns?

**Status:** throwaway experiment, complete. Code is a diff, not a branch to merge.
**Date:** 2026-08-06. **Tree:** `1e5ffa8`.
**Plan:** `doc/devel/ng/impl_plan/locus_stream_shape_experiments.md` §4, Sketch 2.

---

## The answer, in one sentence

**The merge wants columns only when a columnar file is already on the other side
of it: reading `.psp`, folding a cheap column first halves the instructions and
cuts allocations 414-fold (19,571,302 down to 47,315); reading ng's generator
directly, the
merge is 1.6 % of the path and the two shapes are indistinguishable — so the
merge is not a reason to make ng's locus stream columnar.**

---

## What was built

Two binaries, both throwaway, both in the diff beside this file.

`examples/sketch2_merge_shape.rs` — the **`.psp`-fed** configuration. Three arms
over 50 tomato `.psp` files, single-threaded, all covering byte-identically the
same 955 kb of covered genome:

| arm | what it is |
|---|---|
| `merger` | production's own `PerPositionMerger` (`src/var_calling/per_position_merger.rs:145`) driven by the shipped `OwnedRecordsIter`. The literal default architecture. |
| `records` | the same shape written fresh: one owned `PileupRecord` of lookahead per sample in a reused array, an O(N) scan over the heads at each output position, no per-position allocation. |
| `fold` | one `BlockColumnReader` per sample; derive the light columns from each loaded block, fold them across samples over the watermark window, materialise only the variable positions. |

I ran **both** `merger` and `records` because `PerPositionMerger` turned out to be
public and reusable, and because keeping them apart separates *the record shape*
from *production's implementation of it*. That distinction is worth 388 MB of
churn — see "What production's merger costs beyond its shape" below.

`examples/sketch2_ng_merge.rs` — the **ng-fed** configuration, which is the one
the plan's question is actually about. Ten tomato benchmark CRAMs, ten
`PileupGenerator`s driven in lockstep over 300 regions, no `.psp` anywhere:

| mode | what it is |
|---|---|
| `walk` | every generator driven to exhaustion, every locus dropped. No merge. This is the **producer floor** — the share of the path no merge shape can change. |
| `records` | one owned `SampleLocusObservations` of lookahead per sample, O(N) head scan. |
| `fold` | per-sample owned buffers of loci, a light column derived as they arrive, folded across samples, only variable positions materialised. |

### Correctness

Every arm feeds the same sink, which folds each merged position into an FNV-1a
digest over raw bit patterns — position, contributing samples, and per allele the
sequence bytes, observation count, `q_sum` bits, forward count, placed-left,
placed-start, MAPQ sum and MAPQ sum-of-squares, and the chain-id list — then runs
a genotype-likelihood + allele-frequency EM over it and folds the resulting
log-likelihood's bits in too.

**Tolerance allowed: none.** The merge performs no floating-point arithmetic (it
only moves `f64`s), and the EM sees identical inputs in identical order in every
arm, so exact bit equality is achievable and was demanded. All three `.psp` arms
agree on `0x55827295ee08398e` and `loglik_acc = -4621584.027974`; all three ng
modes that merge agree on `0xf46f8f924ee468aa` and `-56602.839011`. The per-arm
counters (positions kept, evidence objects, bytes copied) match to the unit.

### Fixtures, and which `.psp` files

`benchmarks/tomato1/results/ours/cohort/psp/`, first 50 files, **written
2026-07-07**. These are the regenerated ones: their headers carry `kind = "snp"`
and the REF alleles carry empty chain-id lists, so the 96.6 % of chain ids that
the merger discards are already absent from the file. `tmp/aligned_psp/` was
**not** used — those were written 2026-06-02, predate the `kind` header field,
and are the stale set the earlier session measured as 17–20 % heavier.

ng configuration: `benchmarks/ssr_tomato1/crams/*.bench.cram` with
`ssr_regions.bed`. **These CRAMs are tandem-repeat-targeted, not whole-genome** —
coverage is roughly 1 kb islands scattered along the chromosome. Any per-base
number below is from that mode and should not be compared with the whole-genome
per-base baseline in the plan's §3.

### Instrument

`instructions retired` from `/usr/bin/time -l`, three runs a side, arms
alternated within one script, **single-threaded throughout** (no rayon, no
threads of any kind, in any arm). Run-to-run spread was under 0.1 % for the
`.psp` arms and under 0.2 % for the ng arms; medians are quoted.

Floor for the `.psp` arms: a run that opens all 50 files and parses all 50 block
indexes, then stops before merging anything — measured per arm per rep, roughly
0.27 × 10⁹ instructions. For the ng arms the `walk` mode *is* the floor, so no
subtraction is applied and the merge cost is quoted as the delta.

Heap: `--features dhat-heap --target-dir target-dhat`, always.

Wall-clock is recorded in the raw files but no conclusion rests on it.

---

## Result 1 — the `.psp`-fed merge

50 samples, 200 blocks, 954,851 covered positions, 112,419 of them variable
(11.8 %) at production's `--min-alt-obs-per-sample` default of 2, 46,701,870
sample-loci. All numbers **measured**.

| arm | instructions per merged position | per sample-locus | peak RSS | peak heap (dhat) | allocations (dhat) | bytes copied at the handoff | positions materialised as objects | lines of source | what it cannot do |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `merger` | 967,449 | 2,328 | 151.1 MB | 80.73 MB | 19,765,312 | 4,638 MB | 46,701,870 records + 5,519,822 evidence | 65 | build nothing per covered base — it makes a `Vec<Option<PileupRecord>>` for every one of them, 88 % discarded |
| `records` | 954,946 | 2,299 | 160.1 MB | 80.72 MB | 19,571,302 | 4,638 MB | 46,701,870 records + 5,519,822 evidence | 107 | see a position's cohort evidence before every sample's record for it has been built |
| `fold` | **478,254** | **1,151** | 149.6 MB | 84.59 MB | **47,315** | **252 MB** | **5,519,822 evidence** | 251 | hand anything downstream that it did not decide to materialise *in advance*, from the light columns alone; and hold a borrowed view across a block advance |

Instructions are floor-subtracted medians with the full consumer (digest + EM).
Allocations, peak heap and total churn come from a dhat run at 40 blocks (9.4 M
sample-loci), which is why they are not on the same scale as the instruction
column; the ratios are what matter and they are stable. Peak RSS is
`/usr/bin/time -l`'s maximum resident set at the 200-block scale and is noisy to
about ±10 % — dhat's peak heap is the trustworthy memory number.

**The instruction gap is 2.00×** (107.35 × 10⁹ against 53.76 × 10⁹). Without the
EM — digest consumer only, which isolates the merge — **it is 2.74×** (84.32 × 10⁹
against 30.72 × 10⁹). The EM itself costs 23.04 × 10⁹ instructions in *both* arms,
to the second decimal place, which is the check that both arms really did hand
the consumer the same work.

**The allocation gap is 414×**: 19,571,302 allocations against 47,315, for the
same output. Bytes churned fall from 1.69 GB to 0.81 GB.

**The byte gap at the handoff is 18.4×.** Both record arms copy 4,386 MB into
`PileupRecord`s — every covered base of every sample, sequence bytes plus chain
ids plus the fixed-width support statistics — and then copy 252 MB again into the
merged evidence for the 11.8 % of positions that survive. The fold arm copies
only that second 252 MB. The 4,386 MB is not the merge's fault in any narrow
sense; it is the price of the merge *wanting records*, and it is the number the
plan asked to have measured.

**Peak heap goes the other way, and by a small amount: the fold arm is 4.8 %
higher** (84.59 MB against 80.72 MB). The 3.87 MB difference is almost exactly
the fold's own scratch: three `u32` arrays — absolute position, summed non-REF
observations, cumulative allele offset — over 4,592 records per block for 50
samples is 2.75 MB, plus the dense fold window and vector slack. Both shapes
already hold one decoded block per sample, because the record iterator decodes a
whole block before yielding its first row; that is why they otherwise tie.

### What production's merger costs beyond its shape

`PerPositionMerger::next()` builds a fresh `vec![None; n_samples]` for every
output position (`per_position_merger.rs:306`). At 40 blocks that is 193,774
positions × 50 samples × 40 bytes = 387.5 MB — and the measured churn difference
between `merger` and `records` is 387.6 MB (2.079 GB against 1.692 GB), with
194,010 extra allocations. In instructions it is worth 1.4 %: 108.74 × 10⁹ against
107.35 × 10⁹.

So of the 2× gap, essentially none is production's implementation. It is the
record shape.

---

## Result 2 — the ng-fed merge, which is the question that decides this

10 CRAMs, 300 regions, 2,830,932 loci produced, 309,018 covered positions, 3,060
of them variable (1.0 %). All numbers **measured**.

| mode | instructions (median of 3) | merge's share of the path | peak RSS | peak heap (dhat) | allocations (dhat) |
|---|---:|---:|---:|---:|---:|
| `walk` (producer floor, no merge) | 26.128 × 10⁹ | — | 272.4 MB | 219.5878 MB | 5,126,114 |
| `records` | 26.565 × 10⁹ | 0.437 × 10⁹ = **1.7 %** | 272.7 MB | 219.5878 MB | 5,126,247 |
| `fold` | 26.542 × 10⁹ | 0.414 × 10⁹ = **1.6 %** | 282.7 MB | 219.5887 MB | 5,126,565 |

The dhat rows are from a 60-region run (546,178 loci); instructions and RSS are
from the 300-region runs.

Three things to read off this.

**The merge is 1.6–1.7 % of the path, and the two shapes are within noise of each
other** (0.414 × 10⁹ against 0.437 × 10⁹, with a run-to-run spread of about
0.03 × 10⁹). Whichever is nominally ahead flips between the two consumers, which
is the honest way of saying they are the same.

**The merge shape is invisible in the allocation profile.** The producer makes
5,126,114 allocations walking those loci — 9.4 per locus, which is the
thirteen-allocations-per-covered-base shape the plan cites. On top of that, the
record merge adds **133** and the fold merge adds **451**. Out of 5.13 million.
There is nothing here for a layout change to remove.

**The fold arm costs 10 MB more resident** (282.7 MB against 272.7 MB) for its
per-sample owned buffers, and buys nothing back.

### Why the fold cannot win here, stated mechanically

In the `.psp` configuration the fold reads a decoded block it can borrow from, so
"materialise late" means the 88 % of positions that get discarded are never built
as objects at all. In the ng configuration there is no block to borrow from and,
more decisively, **`PileupGenerator::next_locus` is the only emitter the
generator has, and it returns a fully owned `SampleLocusObservations`** — there is
no API that yields any per-locus summary, non-reference observation count
included, without building the whole record first. (The buffer-recycling path
`OpenPileupRecord::finalise_recycling` is private and recycles the fold table and
bucket list only; the reference bytes still leave with every locus.)

So the fold arm's light column has to be *derived from records that were already
built*, and its "materialise late" saves only a copy it was never going to make —
the record arm checks the keep rule before copying too. What is left is pure
overhead: buffering, a light-column derivation, and a window fold.

That overhead is small, which is why the two arms tie rather than the fold
losing. But the direction of the result is what matters: **the `.psp`
configuration's 2× is a property of the file being columnar, not of the merge
being a merge.**

---

## What the per-sample cursor cost when it could not borrow

This is the part that outlives the experiment, so it is reported in the terms the
plan asked for.

`BlockColumnReader::columns()` returns a `BlockColumns<'_>` borrowing the
reader's own reusable decode buffers (`src/psp/reader.rs:1380`). A view cannot
survive the next `load_current()` or `advance()`. That is production's
self-referential cursor — and `BlockColumnReader` itself exists because of it:
its own doc says it owns its `PspReader` "so a cohort producer can hold one per
sample in a `Vec` without a self-referential borrow" (`reader.rs:1234-1237`).

**I solved it with neither a reference count nor a per-sample owned copy: a round
structure.** N distinct readers can all be borrowed *shared* at the same moment,
because each borrow is of a different reader. So the loop is three phases:

1. **mutable** — every live sample loads one block, and its light columns are
   derived into scratch that lives outside the reader vector;
2. **shared** — every sample's `BlockColumns` is taken at once into a
   `Vec<Option<BlockColumns<'_>>>`, the fold runs across all of them, and the
   variable positions are materialised straight out of the columns;
3. **mutable** — the views are dropped, and every sample whose block is drained
   is advanced.

Cost of the constraint, in three parts:

- **Zero bytes copied.** Re-taking the views each round is free — `columns()`
  only assembles a struct of slices. This is the part that is *better* than what
  production does, which copies row ranges out of columns into new columns
  (`append_range`), the 49 %-of-peak item.
- **3.87 MB of per-sample scratch**, measured, and the reason the fold arm's peak
  heap is 4.8 % above the record arms'. The light columns cannot live inside the
  reader's buffers, so they live beside them.
- **An architectural bound that does not go away.** Because no view can be held
  while any reader advances, the merge can only progress to the cohort watermark
  — the minimum, over samples, of the last position of the block each sample
  currently has loaded. The window the fold works on is set by whichever sample
  has the shortest loaded block, not by anything the algorithm wants. Production
  hit the same wall and answered it with `Arc`-shared segments plus a
  decode-once straddler cache; the round structure answers it by simply not
  advancing until every sample is ready, which is simpler and, at these block
  sizes, adequate. It would stop being adequate if block boundaries were badly
  misaligned across samples. On this fixture they are not: across two files,
  178 of 200 distinct block boundaries are shared by both, a union-to-per-file
  ratio of 1.1×.

In the ng-fed configuration the same constraint has a different and blunter
answer: **there is nothing to borrow**, so the fold arm's buffer is a queue of
records it owns outright. That is "a per-sample owned buffer" in its purest form,
and it is what the shape degenerates to whenever the producer emits records.

---

## What I changed as I built it, and why

**Added a third arm.** The plan describes two. On finding `PerPositionMerger`
public and reusable, I ran production's merger *and* a fresh implementation of
the same shape, because otherwise a per-position `Vec` allocation inside
production's merger would have been charged to "records" as a property of the
shape. It is not: it is 1.4 % of instructions and 388 MB of churn, and separating
it makes the remaining 2× unambiguously about the shape.

**Replaced the fold's sorted k-way merge with a dense array over the watermark
window.** Production's `CohortSpanFold::fold_sample_light` merges sorted position
arrays pairwise, N times, because `.psp` rows are sparse relative to the genome.
ng's stream is one locus per covered reference base, so a `Vec<u32>` indexed by
`position − window_start` is enough, and it removes the merge entirely. This is a
small result in its own right: **when the locus stream is dense, the fold does
not need a sorted merge.**

**Then had to apply the same change to the ng arm, and it mattered.** My first
ng fold used a sorted union with binary-search insertion, which made the ng fold
arm look 69 % more expensive than the record arm (0.699 × 10⁹ against
0.415 × 10⁹). With the dense window it is 0.414 × 10⁹ — a tie. I would have
reported a false result. The tied numbers are the ones quoted; the superseded run
is not.

**Simplified the keep rule.** Production groups positions by overlapping reach
and keeps or drops a whole group, so an indel's neighbours ride along. Here a
position is variable iff the maximum, over samples, of its summed non-reference
observations reaches 2. Every arm applies the identical rule, so the comparison
is sound, but the 11.8 % variable-position rate is not directly comparable with
the "4 positions in 100" the plan cites for production's grouped rule on a
different cohort.

---

## What building each arm was actually like

The default position is records everywhere with columns only in the file, so this
section is the one that has to be honest.

**`merger` took about twenty minutes and I got it right first time.** It is
fifteen lines of wiring: open the readers, check the chromosome agreement, turn
each into an owning record iterator, hand the vector to `PerPositionMerger`, and
loop. There is nothing to get wrong because the merger owns the invariant.

**`records` took about an hour, and the only difficulty was borrow-checker
plumbing** — needing `cur[si].head.as_ref()` and `cur[si].pull(...)` in the same
loop body forced an index loop instead of an iterator. The algorithm is three
O(N) scans and I could hold all of it in my head. It was right the first time it
compiled.

**`fold` took most of a day, and I got it wrong twice in ways only careful
reading caught.** Naming them precisely, because this is what the "lines of
source" column is standing in for:

- The first version re-scanned each sample's block from its round-start cursor
  for *every* kept position, which is quadratic in the window. Correct output,
  and it would have been reported as a fold arm slower than it really is. The fix
  needed a second per-sample cursor threaded through the shared phase — a cursor
  that cannot live in the same structure as the round cursor, because one is
  borrowed immutably while the other is advanced.
- The window has to be capped, because the position key packs the chromosome id
  into the high 32 bits and a chromosome change makes the key jump by 2³²; a
  dense array over that window would ask for terabytes. The cap then interacts
  with the watermark and the advance rule, and getting all three to make forward
  progress in every case — capped window, uncovered gap, exhausted sample,
  partially consumed block — is the part I would not want to maintain.

The offset arithmetic is the tax the plan predicted. Deriving absolute positions
from `delta_pos`, cumulative allele offsets from `n_alleles`, and then slicing
two ragged CSR columns with `offsets[j]..offsets[j+1]` is re-derived in this
sketch at three places, in 251 lines against the merger's 65. That is the same
shape of cost the plan records as "re-derived at five call sites" in production,
and I reproduced it faithfully without meaning to.

**So: the fold arm is twice as fast on `.psp` and I found it materially harder to
write correctly.** It is 3.9× the source of production's merger and 2.3× the
source of the hand-written record arm, and both of its bugs were in the cursor
and window bookkeeping that only exists because a view cannot outlive a block.

---

## What this does and does not settle

**Settles:** the merge is not a reason to make ng's generic locus stream
columnar. Fed by ng's generator the merge is 1.6 % of the path, both shapes cost
the same, and the merge contributes 451 allocations out of 5.13 million.
Whichever shape is easier to read wins, and by the plan's §2 that is records.

**Also settles, in the other direction:** if the pipeline *does* persist to
`.psp` and read it back — which the plan says is usually done — then the merge
side of that read should fold before it materialises. It halves the instructions,
cuts allocations 414-fold, and cuts bytes copied at the handoff by 18×,
for 4.8 % more peak heap and roughly 190 extra lines. That is a decision about
the `.psp` read path, not about the locus stream's shape.

**Does not settle:** whether the generator itself should emit columns. Everything
above says the *merge* would not benefit, because the merge is not where the cost
is; whether the *producer* benefits is Sketch 1's question, and this sketch's
`walk` mode (26.13 × 10⁹ instructions for 2.83 M loci, no merge at all) is
evidence only that the producer is where the cost is.

**Not measured, and worth saying:** the two-phase column-selective decode —
inflating only the heavy columns the fold's keep mask actually wants — is a
capability only the fold arm can have, and it is `pub(crate)`
(`BlockColumnReader::decode_current_two_phase`), so this sketch could not reach
it. Production has measured it; I did not re-derive it, and the fold arm's 2× is
therefore a floor rather than a ceiling for what the columnar `.psp` read path
can do.

---

## Files

| file | what it is |
|---|---|
| `sketch2_code.diff` | both sketch binaries, as a patch |
| `sketch2_instructions_raw.txt` | `.psp` arms: 3 reps × 2 consumers × 3 arms, with floors |
| `sketch2_counts_raw.txt` | `.psp` arms: the accounting run (counters on — not a timing run) |
| `sketch2_dhat_raw.txt`, `dhat_{merger,records,fold}.json` | `.psp` arms: heap profiles |
| `sketch2_ng_instructions_raw.txt` | ng modes: 3 reps × 2 consumers × 3 modes |
| `sketch2_ng_counts_raw.txt` | ng modes: the accounting run |
| `sketch2_ng_dhat_raw.txt`, `dhat_ng_{walk,records,fold}.json` | ng modes: heap profiles |
| `measure_psp.sh`, `measure_ng.sh`, `dhat.sh`, `dhat_ng.sh`, `counts.sh`, `ng_counts.sh`, `agree.sh` | every script that produced the above |
