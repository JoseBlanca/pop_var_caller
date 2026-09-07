# ng — window coverage: each sample's depth and GC around a locus, and the histogram that turns them into a copy number

*Status: draft, 2026-09-06 — the design settled with the owner in the conversation preparing the
hidden-duplication filter's port. No code yet. This is the first of two documents: it produces the
numbers, and [`hidden_paralog_filter.md`](hidden_paralog_filter.md) consumes them. Implementation
plan: [`../impl_plan/window_coverage.md`](../impl_plan/window_coverage.md). No architecture
document: the type blocks in §3 are the code shape the plan cites.*

---

## 1. What this is

**One computation, run once per sample inside the calling pass, with two consumers.** Given one
sample's records in coordinate order, it emits for every covered position the **mean read depth**
and the **GC fraction** of the 500-base window centred on that position. The two consumers are:

- **a per-locus observation** — at every locus the caller writes, each sample's pair at the
  locus's first base, handed to the filter beside the record;
- **a per-sample histogram** — every emitted pair, binned by GC and by depth, complete when the
  pass ends, from which the filter fits what one copy's depth looks like in that sample.

The pair is the *measurement* and the histogram is the *yardstick*; the filter divides one by the
other to get a copy number, so they have to be the same quantity measured the same way. Producing
both from one accumulator in one place is what this document exists to arrange.

**Where it runs is the whole of the design.** The merge's observation cache draws every record of
every sample in coordinate order, in both modes — direct mode from the alignment files, psp mode
from stored files — and that draw is the one place every position of every sample passes
([`draw_next`, `observation_cache.rs:902`](../../../../src/ng/run/cohort_merge/observation_cache.rs)).
The accumulator lives there. Nothing is stored in the psp and nothing is recomputed in a second
pass over the reads.

### 1.1 Goals

1. **Every written record carries, per sample, the window's mean depth and GC at the locus** — or
   an explicit *absent* where the sample had no usable window there.
2. **Every sample's histogram is complete at the end of the calling pass**, in both modes.
3. **One definition of the window, in one type**, so the measurement and the yardstick cannot
   drift apart.
4. **No psp format change.** The psp's head already carries the number this reads at most
   positions; the rest is decoded from bodies the run would otherwise skip, at a cost bounded by
   how much of the genome is repeat tract or indel.
5. **Mode-independent**: direct mode and psp mode run the same code on the same numbers and
   produce identical pairs and histograms.
6. **Memory per sample stated and bounded**, charged against the calling stage's per-sample
   budget ([`run_streaming.md`](run_streaming.md) §7.2).

### 1.2 Non-goals, and what this does not do

- **It scores nothing and fits nothing.** The coverage model, the likelihood ratio, the verdict —
  all [`hidden_paralog_filter.md`](hidden_paralog_filter.md).
- **It does not put a histogram in the psp trailer.** Considered and set aside (§8): a trailer
  gives psp mode a source direct mode cannot have, so the filter would run through two paths.
- **It changes no call and no VCF byte.** With the filter off, a run's output is byte-identical to
  today's (§10).
- **It does not fill the repeat-bundle hole.** `SsrBundle` and `Satellite` regions emit no loci
  today ([`locus_generation.md`](locus_generation.md) §11), so their positions are absent from
  every window. That is consistent — absent from the measurement and from the yardstick alike —
  and this document does not change it.

### 1.3 Vocabulary

- **Covered position** — a reference position at which the sample has a record, whose reference
  base is not `N`. Only covered positions enter a window; everything else is simply not there.
- **Window** — the covered positions within 250 bases either side of a centre, on one contig.
  Width 500 is production's `--gc-window-bp` default and is inherited, not re-measured.
- **Window mean depth** — the sum of per-position depth over the window's covered positions,
  divided by their count. **Window GC** — the fraction of those positions whose reference base
  is `G` or `C`.
- **Observation depth** — at a position, the number of reads whose observation covers it.
  Reads that produced no observation (`reads_without_observation`) and reads a depth cap
  discarded are not in it (§3.1). It is below read depth, consistently on both sides.
- **The head count** — `reads_compared_with_reference`, the fixed field in every psp record's
  head ([`record.rs:138-144`](../../../../src/ng/psp/record.rs)): reads whose whole sequence over
  the locus was compared against the reference.

---

## 2. Where it sits

**The cache is mode-blind, and that is what makes this one design rather than two.** Each sample
has one forward reader for the whole run; the cache draws them ahead to cover a building region,
holds what it drew for the builders, and drops what nothing can reach
([`cohort_merge.md`](cohort_merge.md) §6.4). What a draw yields is a `Drawn`
([`observation_cache.rs:167`](../../../../src/ng/run/cohort_merge/observation_cache.rs)): in direct
mode `Built`, the whole record the walker just minted; in psp mode `Kept`, the head's summary with
the body left in the source, buildable on demand through `source.build`
([`observation_cache.rs:710`](../../../../src/ng/run/cohort_merge/observation_cache.rs)). Both
shapes carry the record's region and the head count; only one carries the body.

**What consumes the result.** A cohort locus assembled by a builder
([`CohortObservation`, `build.rs:1055`](../../../../src/ng/run/cohort_merge/build.rs)) is called and
turned into a record on the merge thread, in genome order
([`call_cohort_from_sources_handing_each_record_over`, `callers.rs:799`](../../../../src/ng/run/callers.rs)).
The per-sample pair joins the evidence gathered at that moment
([`SampleEvidenceForOutput`, `assemble.rs:48`](../../../../src/ng/vcf/assemble.rs)) and travels
beside the record to the run's sink. The histograms come out of the cache when the merge returns
its sources.

---

## 3. The design

### 3.1 The depth at a position — decision: observation depth, from the head where the head has it

**Per position, the depth is the number of reads whose observation covers that position.**
Production computes the same thing — the sum of `num_obs` over every allele at a record
([`pileup_to_psp.rs:101`](../../../../src/pileup/per_sample/pileup_to_psp.rs)) — and the filter
was validated on it.

Where the number comes from depends on the record's span, which every `Drawn` carries:

- **A record spanning one base: the head count.** At a single-base locus no read's witness can
  stop *inside* the locus, so every observation is whole and the head count equals the sum of
  `num_obs`. No body is decoded. **This equality is asserted, not assumed**: the plan's first
  measurement walks a real store and checks it at every single-base record (§10).
- **A record spanning more than one base: the body.** A deletion-widened generic record and a
  repeat tract both have partial witnesses — reads whose evidence ends inside the locus — and the
  head count excludes them by design (its own doc: at such a locus it sits "far below the number
  of reads that covered the ground",
  [`locus_generation/mod.rs:244`](../../../../src/ng/locus_generation/mod.rs)). So the body is built and
  [`num_obs_along_locus()`](../../../../src/ng/locus_generation/mod.rs) gives the depth at each
  position, spreading a partial witness over exactly the positions it witnessed. In direct mode
  the body is already in hand; in psp mode it is built through the source, which
  [`cohort_merge_psp_path.md`](cohort_merge_psp_path.md) §3.1 says any thread may do from the
  bytes alone. The predicate is on the head's region, so choosing which bodies to build is free.

**Which positions a record speaks for — decided by kind, because of how the walk emits.** The
generic walk emits one record per covered position, and a record widened by a deletion is
anchored at its first base while the interior positions have records of their own (production's
walker, copied unchanged: `locus_generation/pileup/`). So a **generic** record contributes depth at
its **anchor only** — its first position — or the interior of every deletion would be counted
twice. A **repeat-tract** record is the only record over its ground (region typing partitions the
reference, and a tract is one segment), so it contributes at **every position of its span**. The
kind is read off the body where one is built, and a single-base record is generic by construction.

**What is not in the depth, on both sides alike:** reads that covered a locus and produced no
observation, and reads a depth cap discarded. Neither counter is ever set on the generic path (no
assignment outside the STR generator; the per-position cap is 8,000 reads,
[`walker/mod.rs:83`](../../../../src/pileup/walker/mod.rs)); both are set at tracts. So a tract's
depth is observation depth, below its read depth by the reads that anchored no border, and the
shortfall grows with tract length. It is consistent — the same number trains the yardstick and
makes the measurement — but a tract window is not directly comparable to a generic one, and the
filter has no way to tell which it is looking at. The filter document carries that caveat forward.

### 3.2 The GC — per sample, over the same positions

**A window's GC is averaged over the positions that contributed its depth, and no others.**
Production stores GC per sample rather than one reference GC per position for exactly this reason
([`types.rs:210`](../../../../src/var_calling/types.rs)): each sample's window spans *its* covered
positions. A GC over all reference positions in the window would describe ground that contributed
no depth, and the two would diverge wherever coverage is patchy — which is where the filter is
being asked its hardest question.

So the accumulator takes the reference base *with* the depth at each position, and the reference
supplies bases and nothing more. The cache is handed **one reference accessor for the whole run**,
minted beside the padding accessor the run already holds ([`callers.rs:799`](../../../../src/ng/run/callers.rs)
takes `padding_reference`), and fetches the bases for a cover's ground **once per cover** into a
scratch buffer every sample reads by offset ([`RefSeq::fetch_into`, `ref_seq.rs:143`](../../../../src/ng/ref_seq.rs);
[`WindowedRefSeq`, `ref_seq.rs:594`](../../../../src/ng/ref_seq.rs) is the forward-sliding,
evictable implementation the run already uses). A position whose base is `N` is not covered.

### 3.3 The window — production's sliding accumulator, copied

**The centred sliding window is production's, transcribed into ng unchanged with its tests**:
[`SlidingWindowCoverageAccumulator`, `coverage.rs:331`](../../../../src/sample_summary/coverage.rs).
Fed `(contig, position, reference base, depth)` in coordinate order, it finalises a centre once the
stream has advanced half a window past it — or the contig changes, or `finish` is called — and
emits that centre's `(GC fraction, mean depth)`. Memory is one window's worth of positions. The
rules its tests pin, and which ng inherits by transcribing them
([`coverage.rs:814-1014`](../../../../src/sample_summary/coverage.rs)): an `N` position advances
the frontier without joining any window; a window never spans contigs; a gap wider than half a
window leaves each side's centres with the positions they have; a single covered position emits a
window over itself alone.

**That last rule is the one ng changes, and it is the one addition to the accumulator.** A
window built from a handful of positions is as confident-looking as one built from five hundred,
and in ng the holes are not scattered: a tract region that emits no loci, an analysed-region edge, a
stretch no read reached. So a window finalised over fewer than `min_window_positions` covered
positions **emits an absent value** — `NaN` in both fields — and the filter skips the sample at
that locus, which its scorer already does for an absent sample. **The default is soft and is set
by measurement** (plan step D1): the distribution of covered positions per window on the tomato
slice and on HG002. Until measured, 50.

**The look-ahead is a contract the cache has to keep.** A centre at `p` is complete only when the
stream has passed `p + 250`. A builder handed a region ending at `r` therefore needs every sample
drawn to `r + 250`, where today `cover` draws to the region's end plus the reach of any observation
chaining past it ([`observation_cache.rs:761`](../../../../src/ng/run/cohort_merge/observation_cache.rs)).
**`cover` draws half a window further.** What that costs is half a window more of held summaries
per sample — about 250 at three reads a position, some tens of kilobytes at a thousand samples —
and in psp mode the kept bytes behind them, which the source's arena holds regardless
([`cohort_merge_psp_path.md`](cohort_merge_psp_path.md) §3.4 records that arena as unreleased; this
design neither fixes nor worsens that). **Getting this wrong is silent**: a region's last centres
would be finalised by the *next* cover, after their builder has run, so the builder finds no value
and the sample reads as absent at exactly the loci nearest every region boundary. The plan lands it
as its own commit with the whole-store oracle green either side (§10).

**Emitted values are held per sample until eviction.** The accumulator's ready values go into a
per-sample deque keyed by position, alongside the held summaries; `evict_before` drains it with
them; a builder reads it through the cohort window it is already handed
([`WindowedCohort`, `observation_cache.rs`](../../../../src/ng/run/cohort_merge/observation_cache.rs)),
as a third parallel view beside summaries and records.

**Truncation at edges needs no rule of its own.** Positions with no record are absent from the
window — a tract hole, an uncovered stretch, the gap between two analysed regions on one contig —
and the count in the denominator is what shrinks. The owner's ruling: treat a tract edge as
production treats a contig end. The floor above is what stops a mostly-absent window from
speaking.

### 3.4 The histogram — production's shape, with the depth axis scaled to the sample

**Every finalised window that clears §3.3's floor is folded into a per-sample count matrix
indexed by GC bin and depth bin**, production's [`CoverageByGcHistogram`, `sample_summary/mod.rs:118`](../../../../src/sample_summary/mod.rs)
— the input the filter's model fit reads and the one it was validated on. The accumulator already
does this fold ([`coverage.rs:445`](../../../../src/sample_summary/coverage.rs), `finish` returns
the histogram beside the tail).

**A window under the floor is finalised and is deliberately *not* folded.** Both casts of `NaN`
come back `0`, so folding one would pile every too-sparse window into a single cell — the lowest
GC bin's lowest depth bin — in exactly the samples and regions where coverage is thinnest, and
the yardstick would carry a spike of windows that had no depth to report. The histogram counts
the silenced windows separately (§3.5) so that a yardstick fitted from a small share of a
sample's positions is visible as such.

**The bin scheme is where ng departs, because production's does not fit the budget or the
range.** Production bins depth at 0.5× in 2,000 bins to 1,000×, 50 GC bins, four bytes a cell:
**400 kB a sample**. The calling stage's per-open-sample budget is 500 kB, of which a human
reference already spends 480 ([`run_streaming.md`](run_streaming.md) §7.2). And a fixed width
serves one end of the depth axis at a time: at three reads a position a 0.5× bin is a sixth of the
single-copy peak's own scatter, and at 300 reads the range has to reach 1,200× for a four-copy
carrier.

**Decision: 50 GC bins and 400 depth bins, with the depth bin's width set per sample from the
sample's own depth.** The accumulator buffers its first `depth_scale_windows` finalised windows —
10,000 — takes the median of their mean depths, and sets the width so that the 400 bins span ten
times that median (`width = median / 40`); it then folds the buffered windows and every later one.
A sample whose whole run finalises fewer than 10,000 windows sets the width from what it has at
`finish`; a sample with none has no histogram and is carried absent. **80.2 kB a sample** —
50 × 401 × 4 bytes — plus a **262 kB** transient while the first windows are held back — 10,000 pairs of `f64` is
160 kB of live data, and the list grows by doubling, so its capacity reaches 16,384 (measured,
2026-09-06; this paragraph previously said 120 kB, which left the budget below looking about
145 kB roomier than it is); 80 MB at a
thousand samples, 240 MB at three thousand. On tomato that fits today (108 kB + 80); on a human
reference it fits once the run's readers share one contig list, the psp path plan's step D1
(480 kB → 123). The plan measures the overflow fraction — the share of windows above the range,
which is the fit's own rejection guard — on both benchmarks (step D2).

**Soft, and marked so.** The 10,000-window scale sample, the factor of ten, and 400 bins are
starting values chosen to fit the budget and cover both ends of the depth axis; none has been
measured against the fit's accuracy. The alternative that lost: a depth axis on a log scale, which
would give the same resolution at every depth in fewer bins, and lost because the fit's three
functions assume uniform bins and would have to be rewritten rather than copied.

**Determinism holds at any thread count.** Each sample's records arrive in one fixed coordinate
order whatever the cover's schedule (the parallel cover reaches the same fixpoint by any
schedule, [`observation_cache.rs:792`](../../../../src/ng/run/cohort_merge/observation_cache.rs)),
so the first 10,000 windows, the width, and every count are the same run to run.

### 3.5 What a locus carries, and how it reaches the filter

**A cohort locus carries one pair per covering sample, read at the locus's first base.** The pair
is the window centred on that position for that sample; a sample covering the locus with no value
at that exact position — its coverage begins a base later, or the window was under the floor — is
absent. Production keys the same way, on the record at the anchor position or `NaN`
([`variant_caller.rs:338`](../../../../src/var_calling/variant_caller.rs)).

```rust
/// One sample's window at one locus. Both fields `NaN` where the sample has no usable
/// window there — compared and stored by bit pattern, never by `==`.
pub struct WindowCoverage {
    pub gc_fraction: f32,
    pub mean_depth: f32,
}

// On the cohort locus, parallel to `per_sample` (the covering samples):
//   pub window_coverage: Vec<WindowCoverage>,
// On the evidence gathered for output, dense over the run's samples:
//   pub window_coverage: WindowCoverage,   // NaN pair where the sample is absent
```

**The record itself does not change.** The pair travels *beside* the record to the run's sink —
the sink's signature gains it — so that a run with the filter off writes byte for byte what it
writes today, and a run with it on has what the filter's spill needs without re-deriving anything
from the record.

**The histograms leave the cache when the merge does.** `into_sources` already hands the readers
back so a mode can report per-sample facts
([`observation_cache.rs`](../../../../src/ng/run/cohort_merge/observation_cache.rs), `into_sources`);
the histograms come out the same way, one per sample in the run's sample order, each either a
fitted histogram or the reason there is none.

```rust
/// A sample's finished histogram, or the reason it has none. The three silences are
/// three different reports and this is the only place the difference exists: the pass
/// reached nothing of this sample; every window it did finalise was under the floor,
/// which is a reading on `min_window_positions` rather than a fault; or its windows
/// reported no positive median depth, which is. A caller that does not act on the
/// difference says `.fitted()` and gets an `Option` back.
pub enum SampleHistogram {
    Fitted(CoverageByGcHistogram),
    NoWindowFinalised,
    EveryWindowUnderTheFloor,
    MedianDepthNotPositive,
}

/// Production's `CoverageByGcHistogram` shape, carrying what the model fit reads:
/// the scheme (window width, GC bins, depth bin width and count), the counts, the
/// number of windows folded, and the number the floor silenced — a yardstick fitted
/// from a tenth of a sample's positions is one to trust less, and nothing else
/// downstream can recover that, because an absent window leaves no trace in any cell.
pub struct CoverageByGcHistogram { /* copied; §7 */ }
```

### 3.6 The accumulator's interface

```rust
/// Production's centred sliding window, transcribed, plus the floor and the
/// per-sample depth scale.
pub struct WindowCoverageAccumulator { /* private */ }

impl WindowCoverageAccumulator {
    /// `window_bp` 500, `gc_bins` 50, `depth_bins` 400, `min_window_positions` (soft, §3.3),
    /// `depth_scale_windows` 10,000 (soft, §3.4).
    pub fn new(config: WindowCoverageConfig) -> Self;
    /// One covered position, in non-decreasing `(contig, position)` order.
    pub fn observe(&mut self, contig: ContigId, position: Position, reference_base: u8, depth: u32);
    /// The next finalised centre, oldest first; `None` when none is complete yet.
    pub fn pop_ready(&mut self) -> Option<(GenomePosition, WindowCoverage)>;
    /// Finalise the tail and hand back the histogram, or the reason there is none.
    pub fn finish(self) -> (Vec<(GenomePosition, WindowCoverage)>, SampleHistogram);
}
```

---

## 4. One sample and three thousand, three reads and three hundred

- **One sample.** Nothing branches: the accumulator is per sample and the reference is fetched
  per cover whatever the cohort size.
- **Three thousand samples.** 240 MB of histograms, plus half a window of extra held summaries
  per sample. Both are charged to the per-sample budget §7.2 of `run_streaming.md` names, and
  §3.4 says which reference fits it today.
- **Three reads a position.** A window of 500 positions at three reads each is a mean over
  ~1,500 reads, which is the reason the filter windows at all (production's spec §4: per-base
  counts cannot tell one copy from two at this depth). The floor matters most here — a window
  over 50 positions at 3× is a mean over 150 reads.
- **Three hundred reads a position.** The depth axis scales with the sample (§3.4), so the
  single-copy peak and a four-copy carrier both land inside the range where production's fixed
  scheme would have overflowed and rejected the fit.

## 5. Cross-cutting concerns

- **Memory** — §3.4's histogram, the ready deque (one **24-byte** entry per held position per
  sample: §3.6's `pop_ready` hands back a `GenomePosition` beside the pair, and ng's `Position`
  is a `u64`), and the half-window of extra held summaries. Nothing grows with the genome.
- **Errors** — a reference fetch that fails is a `RunError` naming the region, as the padding
  fetch's failure already is (`RunError::PaddingBaseUnreadable`). Nothing else here can fail: an
  absent window is a value, not an error.
- **Concurrency** — the accumulator is per-sample state behind the same `&mut SampleWindow` the
  cover already holds per sample, so the parallel cover needs no new synchronisation; the
  reference slice for a cover is shared read-only.
- **Determinism** — §3.4; the standing oracle is byte-identical output at any thread count.

## 6. Traps — what will bite the coder

1. **The same number on both sides.** Feed the histogram anything but the depth §3.1 defines —
   the walker's raw depth in direct mode, say, because it is in hand — and the yardstick is on a
   different scale from the measurement. Nothing fails; every copy number is biased.
2. **Anchor only for a generic record.** Spreading a deletion record's depth over its span
   double-counts the interior, which has records of its own. A tract has none.
3. **The look-ahead** (§3.3). Own commit; oracle green before and after.
4. **`NaN` by bits.** A codec, a comparison or a test that normalises `NaN` turns *absent* into a
   number. Production's `LocusWindowCoverage` compares by bit pattern for this reason
   ([`types.rs:233`](../../../../src/var_calling/types.rs)).
5. **The first 10,000 windows are buffered, not folded** — a fold before the width is known has
   no bins to fold into. `finish` must handle the under-10,000 case and the zero case.
6. **The arena is not released** ([`cohort_merge_psp_path.md`](cohort_merge_psp_path.md) §3.4).
   Building bodies for wide records adds nothing to that, but a coder measuring memory here will
   see it and should know whose it is.
7. **`reads_compared_with_reference == Σ num_obs` at a single-base record is a claim.** It is
   asserted by the plan's first measurement; a coder who finds a counter-example has found a
   partial witness at a one-base locus, and the rule in §3.1 changes to "build every body".

## 7. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| the sliding window and its histogram fold | [`SlidingWindowCoverageAccumulator`, `coverage.rs:331-540`](../../../../src/sample_summary/coverage.rs) and its tests `:780-1050` | transcribed into `src/ng/window_coverage/`, production untouched (the freeze rule); the floor and the per-sample width are the only additions |
| the histogram type | [`CoverageByGcHistogram`, `sample_summary/mod.rs:118`](../../../../src/sample_summary/mod.rs) | copied; fields the model fit does not read dropped — `callable_positions` (production's own reader is `var_calling::diversity`) and `n_skipped_tiles`, which the sliding model can only ever write as zero; `windows_under_the_floor` added (§3.5) |
| depth per position from a record | [`num_obs_along_locus`, `locus_generation/mod.rs:67`](../../../../src/ng/locus_generation/mod.rs) | called as-is on bodies spanning more than one base |
| the head count | [`LocusSummary::reads_compared_with_reference`, `observation_cache.rs`](../../../../src/ng/run/cohort_merge/observation_cache.rs) | read as-is at single-base records |
| the draw, hold, evict rhythm | [`draw_next`, `evict_before`, `cover`](../../../../src/ng/run/cohort_merge/observation_cache.rs) | the accumulator is observed at the draw, its values held and evicted with the summaries; `cover` gains the look-ahead |
| the reference | [`RefSeq`, `WindowedRefSeq`, `ref_seq.rs`](../../../../src/ng/ref_seq.rs) | one accessor for the run, one fetch per cover |
| the window definition's rationale | [`hidden_paralog_pileup_window_coverage.md`](../../architecture/hidden_paralog_pileup_window_coverage.md) §2–3 | the "one computation, two consumers" argument, carried whole |

**The parity oracle** is a whole-store recomputation (§10), not production's stored columns:
production's per-position depth is raw depth, ng's is the head count, so the two would differ at
exactly the records §3.1 builds bodies for.

## 8. Deferred, with a recommended home

- **A histogram in the psp trailer**, so that a run restricted to a few regions can use a
  whole-genome yardstick, and so psp mode can skip the fit's third pass. Home:
  [`psp_file_format.md`](psp_file_format.md) §3.4, which reserved the trailer for it. **Only if a
  measurement forces it**, because it adds a second path the filter has to run through — direct
  mode has no trailer — and the owner's ruling is that one path serving both modes outranks it.
  `replace_trailer` means it can be added to existing psps without rewriting them.
- **Depth for `SsrBundle` and `Satellite` ground**, so those positions stop being holes in
  every window that overlaps them. Home: [`locus_generation.md`](locus_generation.md) §11, which
  already names it.
- **A GC window at the fragment scale**, which production's spec §9 marks as a refinement owed.
  Home: the filter's own later work; the accumulator takes its width as a parameter.

## 9. Resolved decisions & open questions

- **Where the computation runs — resolved: the cache's draw, both modes.** Against a stage in
  front of the merge (the owner's first proposal): the keep rule is a cohort decision, so a
  per-sample reader cannot filter to variable loci, and a reader over all files that could is the
  merge; and the merge evaluates that rule inside its parallel builders, which a serial stage in
  front would pull back onto one thread.
- **Where the depth comes from — resolved: head at single-base records, body elsewhere** (§3.1).
  Against the head everywhere: the tract undercount is large and systematic. Against the body
  everywhere: it is the skip the psp path exists to avoid.
- **GC per sample — resolved** (§3.2). Against one shared reference GC stream: the pairing with
  depth breaks where coverage is patchy.
- **Histogram at calling time, not psp-write time — resolved** (§1.2, §8). What it costs is the
  filter's third pass and the histograms' memory during pass one; what it removes is a format
  change, a format-version story, and a second driver whose divergence from the first would have
  been silent.
- **OPEN — the floor's default.** Leaning 50 of 500; **settled by plan step D1's distribution**,
  confirm before the value is written into a default.
- **OPEN — the bin scheme's three constants.** Leaning as §3.4; **settled by plan step D2's
  overflow fraction on both benchmarks**.

## 10. How we know it works

- **The transcribed tests pass unchanged** — production's own suite for the sliding window,
  including the contig-boundary and large-gap cases.
- **The whole-store recomputation oracle.** A probe walks one psp end to end, building every
  body, computes each position's depth by §3.1's rule and each window by §3.3's, and records
  every position's pair and the finished histogram. A calling run over the same psp must produce,
  at every written record, per-sample pairs **bit-identical** to the probe's at that position, and
  histograms **byte-identical** to the probe's. This is the oracle for the look-ahead commit and
  for the depth rule.
- **The single-base equality**, asserted by the same probe: at every record spanning one base,
  the head count equals the body's `num_obs` sum, over a real store.
- **Mode equivalence with the pair attached**: direct mode and psp mode over the same cohort
  produce identical pairs and identical histograms
  ([`scripts/ng_mode_equivalence_oracle.sh`](../../../../scripts/ng_mode_equivalence_oracle.sh)
  extended to compare them).
- **Nothing moved with the filter off**: the VCF is byte-identical to the run before this work,
  on the run fixtures and the tomato slice.
- **Thread-count invariance**: the pairs and histograms are identical at one thread and at
  sixteen.
