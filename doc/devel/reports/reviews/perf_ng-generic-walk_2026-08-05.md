# Performance Review: ng-generic-walk (structural)
**Date:** 2026-08-05
**Reviewer:** rust-performance-review skill (orchestrator)
**Scope:** ng's generic (SNP/indel) locus generator — the organisation of the walk, not another pass over its loops
**Verdict:** Apply the listed wins — the first two need an owner decision first
**Hot-path evidence:** a per-visit cost census on two fixtures, a resolved profile attribution, four independently built and measured changes, and a four-state composition measured across a 10× / 30× / 130× / 300× depth sweep

---

## 1. Scope and constraints

**What was reviewed.** How the walk is organised: `src/ng/locus_generation/pileup/` —
[genome_walk.rs](../../../../src/ng/locus_generation/pileup/genome_walk.rs),
[open_record.rs](../../../../src/ng/locus_generation/pileup/open_record.rs),
[active_read_set.rs](../../../../src/ng/locus_generation/pileup/active_read_set.rs),
[cigar_cursor.rs](../../../../src/ng/locus_generation/pileup/cigar_cursor.rs) — and the read
input path they drive. This was not a fresh six-category sweep; the 2026-08-04 review did that,
and this round was asked for **one structural idea or an honest statement that there isn't
one**.

**Reviewed against.** `6fbbd093764662ed2496acde39424c8ee234ea1c` on branch `ng-generic-perf`.
Every agent detached its own worktree at that SHA and confirmed two branch-only markers before
measuring.

**Targets and hardware.** Human whole-genome sequencing, one sample, ~30×. The generator is
single-threaded; the parallel fan-out is out of scope. One locus per covered reference base is
settled design — emitting only candidate sites is not available and was not proposed. Host: 6
high-performance cores and 12 low-energy cores, 128 kB L1 data cache per performance core.

**Hot-path evidence available.** The standing profile from 2026-08-04 (`sample` over a
`[profile.profiling]` build, 3 M loci of the ~130× tomato CRAM), plus, new this round: a
counted census of what one read-at-one-position costs; a resolution of the profile's largest
unattributed line; four changes built, gated and measured with `instructions retired`; and a
composition of all four measured across four depths. Audit trail in gitignored
`tmp/perf_review_2026-08-05_ng-generic-structure/` — `census.md`, `common_column.md`,
`mate_overlap_sort.md`, `ordered_active_set.md`, `composed.md`, `composed_full.md`, with every
diff, every raw measurement file and the rejected variants beside them.

**Deliberately out of scope.** Designing the parallel fan-out; the SSR/STR generator except
where it shares the cursor; `src/pileup/` and `src/var_calling/` (frozen production).
`copy_fidelity.rs` still pins `decompose.rs` and `chain_id_allocator.rs` byte-identical to
production, so findings against those two are reported, not proposed.

**Categories dispatched.** Five agents in five worktrees, each measuring one hypothesis rather
than running one category checklist: a cost census; specialising the common column; restoring
the walk's lost read order; removing the per-column mate-overlap sort; and composing the
results.

---

## 2. Verdict

**Apply the listed wins.** There is a structural idea, it is large, and it has been built four
times over and measured together:

| depth | fixture | full stack against `6fbbd09` |
|---|---|---:|
| 10× | HG002 chr21 | **−20.95 %** |
| 30× | HG002 chr21 | **−29.28 %** |
| **30×** | **HG002 chr1 — the target workload** | **−28.49 %** |
| ~130× | tomato `SL4.0ch01` | **−42.53 %** |
| 300× | HG002 chr21 | **−18.02 %** |

Every locus line of all four acceptance dumps is identical at every depth below the one where
the per-column depth cap fires. `cargo test --lib` 2,885 pass and 2 fail; clippy `-D warnings`
clean; `cargo doc` at its 12-link baseline.

**The structural finding, in one sentence: the walk runs a fully general machine at every one
of 86 million covered bases, and the general case occurs on fewer than 71 columns in 10,000.**
No column in four million touched more than one open record. A record widens on 14 columns in
10,000; some read carries an insertion or deletion on 24; the record is wider than one
reference base on 33. Mate overlap is the only common exception at one column in ten. Around
that rarity the walk maintains a hash map keyed by read identity, sorts a depth-sized list
twice per covered base, and re-derives per read what a scalar already knew.

**Why the headroom was there to take.** One read at one position costs **1,300 retired
instructions on a BAM and 1,800 on a CRAM**, against 300–600 for the work it performs — pull at
most two events from a cursor, compare one base against one allele, add six scalars, store one
map entry. The gap is not memory: the entire column, at either depth, fits inside a 128 kB L1
data cache (10.2 kB at 30×, 50.4 kB at ~130×, 114 kB at the deepest column seen). **The
instructions are being executed, not stalled on.** That also explains why three struct-narrowing
attempts in earlier rounds measured null, and predicts a fourth would too.

**Two decisions are yours and neither is a performance question.** Both come from restoring the
read order, and both are stated in full in §5 (H2) with their measurements:

1. **ng's mate-overlap tie-break diverges from production's** when quality, mate role and
   alignment start all tie, because its last resort is which read came first in the list. Two
   comparison tests fail. It cannot fire on real data — a genuine mate pair has opposite roles —
   but it is a divergence outside the six named classes, which is a spec call.
2. **At coverages deep enough to trip the per-column depth cap, the emitted loci change.**
   Measured on a 300× chr21 dump: 88,351 of 341,094 rows differ, from 909 capped columns. The
   cap keeps the first N reads, which under the restored order means the N leftmost-starting
   ones rather than a scrambled subset. **At 30× and ~130× the cap never fires and every locus
   line is identical.**

**My recommendation: take all four changes.** The two decisions above are worth making rather
than working around. The first is a tie-break nobody chose — it is an artefact of a container's
removal strategy, and a spec that names it will be more honest than one that inherits it. The
second replaces an arbitrary subsample with a stated one; the question worth asking is not
whether the bytes moved but whether the calls do, and that needs a variant-calling comparison
rather than a walk measurement.

**One thing the measurement discipline caught that the acceptance gates did not.** The
specialised path finishes its result at the base it walks, while the general path holds the
record open one step longer. That one step is invisible to a completed walk and visible to one
that aborts or is abandoned — all four dumps were `cmp`-identical while two tests failed. The
dumps are necessary and are not sufficient.

---

## 3. Measurement plan

Ordered by what unblocks what.

1. **Decide the two questions in §2 before landing the second half of the stack.** Nothing
   below depends on how they are answered, but the ordered-active-set changes cannot land until
   they are. **Threshold:** none; these are owner decisions.
2. **Land the two byte-identical changes now** — specialising the ordinary column and skipping
   the mate-overlap sort, with the eight lines that make them compose (H1, H4, H5). Together
   they are −25.2 % at the 30× target with every dump `cmp`-identical and every counter exact.
   **Threshold:** already met.
3. **Then re-profile before queueing more work.** The walk is about 28 % shorter at the target
   and 43 % shorter at ~130×; the profile has moved again, as it did after each of the last two
   rounds.
4. **If a 300×-class workload is ever a target, the largest remaining lever is handling mate
   overlap inside the specialised path** — measured at 28 points of read coverage there, against
   1.7 for the refinement that looked more obvious. It was declined twice on complexity (a
   second copy of the reconciliation tie-breaks) and should stay declined until someone needs
   it. **Threshold:** it must be worth more than a second implementation of a rule that already
   exists once.
5. **Get a real 30× whole-genome BAM.** Carried unchanged from 2026-08-04. The human fixture is
   tandem-repeat-targeted, and two instructions in five of a `chr1` run go on traversing 613,682
   typed regions to reach 1,541,788 covered bases — 2.5 covered bases per region. Every finding
   here was checked against a second fixture with a different format for that reason.

**The reproducible commands.**

```
TREF=$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa
BIG=/Users/jose/devel/pop_var_caller/benchmarks/tomato_big_cram/DRR000741.p1.cram
HREF=$HOME/genomes/h_sapiens/gca_grch38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
BAM=/Users/jose/devel/pop_var_caller/benchmarks/ssr_hg002/bam/30x/HG002_TR_v1.0.1_Tier_30x.bam

PVC_TRUST_REFERENCE_INDEX=1 PVC_PROBE_MAX_LOCI=1000000 /usr/bin/time -l \
  ./target/release/examples/ng_generic_walk_probe "$TREF" "$BIG" SL4.0ch01
PVC_TRUST_REFERENCE_INDEX=1 /usr/bin/time -l \
  ./target/release/examples/ng_generic_walk_probe "$HREF" "$BAM" chr1
```

**The depth sweep is not optional and now has four points**, since three changes last round won
at one depth and lost at another: `benchmarks/ssr_hg002/bam/{10x,30x,300x}/` plus the ~130×
tomato CRAM.

**⚠ A start-up floor that has been quoted wrong.** `PVC_PROBE_MAX_LOCI=1` on human **chr1**
costs **0.349 G instructions**, not the 1.900 G that has been carried in three briefs — 1.900 G
is a **chr21** floor. Three agents measured 0.349 G independently. Using the wrong one moves a
30× result by about 1.4 % of itself.

**The correctness gate for every experiment.** All four dumps compared with `cmp`, not line
counts — 251,792 / 4,406 / 1,718,914 / 11,945 lines, md5s in
`tmp/perf_review_2026-08-04_ng-generic-walk/` — plus probe counters exact
(`loci=236081 observations=251786 reads_admitted=54709`), and now two more that say each change
is still firing rather than silently disabled: `fast_columns=262498` and
`mate_overlap_positions=39312`.

---

## 4. Build / toolchain configuration

**No change recommended, and none was needed this round.** The `[profile.*]` audit passed in
2026-08-04 and nothing about it moved. The one caveat worth repeating: `[profile.profiling]`
sets `lto = false, codegen-units = 16`, so self-time shares taken under it say *where* work is
and do not transfer to release — which is why every number in §5 is `instructions retired` from
a release build.

**Two build findings from 2026-08-04 are still open and were not touched here**: `cargo test
--release` is red on a clean tree with nine failures in three root causes, and
`benches/psp_writer_perf.rs` panics rather than producing a number. Both remain as filed.

---

## 5. Code-level findings

### Hot-path

**H1: [src/ng/locus_generation/pileup/genome_walk.rs:850](../../../../src/ng/locus_generation/pileup/genome_walk.rs#L850)
and [open_record.rs:2316](../../../../src/ng/locus_generation/pileup/open_record.rs#L2316) — the
general fold runs at every covered base, and 7,898 columns in 10,000 need none of it.**
*Confidence: High.* **Built, gated and measured: −23.81 % at 30×, −34.38 % at ~130×,
byte-identical.**

- **Hot-path evidence.** A census over every column of both fixtures, not an extrapolation:
  1,445,179 tomato columns and 2,430,138 human ones. Columns needing nothing general —
  one covered base, no record already open over it, no read carrying an indel, no mate pair, one
  read group — are **7,898 in 10,000 at ~130× and 7,789 in 10,000 at 30×**. Weighted by reads
  folded, which is what costs: 7,721 and 6,997 in 10,000. **No column in either run touched more
  than one open record** — 0.997 and 0.972 affected records per column — so the `Vec<u32>`
  allocated per column at [open_record.rs:2325](../../../../src/ng/locus_generation/pileup/open_record.rs#L2325),
  its membership scan and its sort are machinery for a case that did not occur once in four
  million.
- **Mechanism.** For an ordinary base the answer is a handful of scalars — one allele bucket,
  one observation, a quality term, a strand bit, two mapping-quality moments, a complete
  witness. The walk instead builds a `SmallVec` of 40-byte events per read, pushes a
  104-byte contribution, opens a record, fetches a reference byte, builds a haplotype string,
  canonicalises a witness set, hashes the read into a per-record map, and at close collects that
  map, sorts it, and re-derives each read's witness and reference agreement.
- **What was built.** `fast_column.rs`, 322 lines; the general path is untouched and takes every
  column that does not qualify. The gate is decided before any per-read work: no open record
  covers the base, the active set is within the column cap, and no read's CIGAR carries an
  insertion or deletion. Mate overlap is settled after the pass. It fires on **74.0 % of columns
  at ~130× and 71.9 % at 30×**.
- **Byte-identity, and where it was at risk.** Quality sums are `f64`, so summation order is
  load-bearing: the specialised path sorts its own compact buffer by read identity and
  accumulates in that order, which is the same sequence of additions the general path gets from
  its sort at close. Reference agreement, the observation key and the observation sort each
  reduce to a closed form for a one-base record, and each was checked against the code rather
  than assumed.
- **The subtle piece, and the gates did not catch it.** The specialised path finishes its result
  at the base it walks; the general path holds a one-base record open one step longer. That step
  is observable to a walk that **aborts** (a reference fetch past the contig end loses whatever
  is still open) or is **abandoned** (a consumer that stops pulling stops the walk a position
  later). All four dumps were `cmp`-identical while two tests failed. Fixed with a one-slot hold
  that `reached_stop` also has to see — and getting *that* wrong lost 91,572 of chr21's 102,938
  regions' out-of-region records in the first draft.
- **Complexity cost, stated honestly.** Two code paths through the hottest code in the walk that
  must agree, and four general-path rules re-derived as closed forms with nothing in the type
  system tying them together: if what makes an observation's identity changes, the specialised
  path silently keeps the old rule. A debug-build assertion compares its scalar read against the
  general cursor query for every read at every ordinary column, which pins the read but not the
  closed forms. The one-step hold is a shared invariant across five call sites and two of the
  five were wrong in the first draft.
- Audit: `common_column.md`, `common_column_fastlane.diff`, `common_column_fast_column.rs`.

**H2: [src/ng/locus_generation/pileup/active_read_set.rs:249](../../../../src/ng/locus_generation/pileup/active_read_set.rs#L249)
and [open_record.rs:403](../../../../src/ng/locus_generation/pileup/open_record.rs#L403) — the
walk destroys the read order it is given and buys it back by sorting, twice per covered base.**
*Confidence: High.* **Built and measured: −6.46 % at 30×, −8.75 % at ~130×, −10.31 % at 300×.
Not byte-identical — see the two decisions.**

- **Hot-path evidence.** The two sorts are **9.6 % of the walk at ~130× and 4.7 % at 30×**,
  priced in release instructions by running each twice and differencing. They sort essentially
  the same list — mean 85.3 and 85.6 elements at ~130×, 17.2 and 17.7 at 30× — once each per
  covered base. `small_sort_general` (4.4 % of the standing profile) is the sort at record
  close; `Vec spec_from_iter` (2.1 %) is the collect that feeds it, not a decoder buffer as its
  neighbours in the ranking suggest.
- **Mechanism.** Reads are admitted in ascending identity. `ActiveReads::expire_passed` removes
  with `swap_remove`, so from the first expiry the active set is a permutation of admission
  order. Everything downstream pays: the per-record fold table must be a hash map keyed by read
  identity; closing a record must sort that map to fix the `f64` summation order; the active set
  carries a second hash map from read identity to slot, maintained on every admission and every
  expiry and queried once per (contributor × affected record); and `swap_remove` moves a large
  `ActiveRead` by value on every expiry.
- **The measurement that settles it.** The same fold-table rewrite built **without** the
  ordering reproduces the 2026-05-12 revert: **−3.18 % at 30× and +23.15 % at 300×** (history
  recorded −3.2 % and +16.4 %). With the ordering, the same code is **−6.46 % / −8.75 % /
  −10.31 %**. A 33-point swing at high depth from arrival order alone. **The old regression was
  caused by the arrival order and by nothing else.**
- **What was built, and my proposed shape was wrong.** A slab with a side list of live slots
  measured **worse** (+1.05 % at 30×, +1.59 % at ~130×): it pays a dependent load on ~260 M read
  accesses per million bases to save ~3 M hash operations. What worked is a `VecDeque` — the
  push at admission *is* the ordering guarantee, and reads leave in very nearly the order they
  arrived, so the removal point is at or near the front — plus a `min_alignment_end` field that
  skips the expiry scan entirely at the three positions in four where nothing can expire. The
  fold table becomes an ordered `Vec` whose append case is one comparison against the last key,
  and the sort at close is deleted.
- **What is load-bearing is the opposite of what one assumes.** Correctness does **not** depend
  on ascending arrival: the container inserts in sorted position whatever order reads arrive in,
  so the summation order is fixed by construction. Only the speed depends on it — ascending
  arrival is what makes the insert a push. That is why the +23 % row is a performance result and
  not a wrong answer.
- **Decision 1 — the mate-overlap tie-break.** After quality, mate role and alignment start all
  tie, `pick_agree_keeper` / `pick_overlap_loser` in
  [genome_walk.rs](../../../../src/ng/locus_generation/pileup/genome_walk.rs) fall back to which
  contributor has the smaller index in the list — which is exactly what reordering changes. Two
  comparison-against-production tests fail. On real data it cannot fire: a genuine mate pair has
  opposite mate roles, so the comparison resolves a step earlier. The synthetic fixtures build
  reads that tie all the way down. **This is a divergence outside the six named classes, and
  naming it — or making the tie-break order-independent — is a spec decision.**
- **Decision 2 — the depth cap keeps different reads.** Measured, not inferred, on a 300× chr21
  dump: **88,351 of 341,094 rows differ, from 909 capped columns**; the *set* of loci barely
  moves (one in a thousand). `contributors.truncate(cap)` keeps the first N, which under
  ascending order means the N leftmost-starting reads rather than an arbitrary subset — and a
  capped read carrying a deletion opens and widens nothing, so record footprints move with it.
  The code's own comment already flagged the corner. **At 30× and ~130× the cap never fires
  (`column_depth_truncations=0` on every fixture), so the target workload is unaffected and
  every locus line is identical.**
- **Third consequence, smaller.** The `record_widen_events` run counter moves by ±2 on the gate
  fixtures (423→425 on chr21, 622→621 on tomato) — whether a record opens at full width or opens
  narrow and then widens depends on which contributor is processed first. **Every locus line is
  identical**; nothing else reads that counter.
- **Complexity cost.** Three overlapping facts about read ends now live on the active set: the
  queue's own order, the expiry guard's minimum, and (with H4) the pair-overlap heap. Each has
  one writer and clear reset sites, but a future change to admission or expiry must keep three
  invariants where it used to keep one.
- Audit: `ordered_active_set.md`, `stage2.diff`, `stage3.diff`, and `stage1_plus_3.diff` — the
  deliberate reproduction of the old regression.

**H3: [src/ng/locus_generation/pileup/open_record.rs:2429](../../../../src/ng/locus_generation/pileup/open_record.rs#L2429)
— the fold looks up by hash the read the contributor was built from.** *Confidence: High.*
**Built and measured: −1.98 % at 30×, −2.86 % at ~130×, −3.33 % at 300×, byte-identical, 2,882
tests green.**

`ReadContribution` is built directly from an `&ActiveRead` in
[genome_walk.rs:860](../../../../src/ng/locus_generation/pileup/genome_walk.rs#L860), and the
fold then finds that same read again through the secondary hash index, once per affected record.
Carrying the slot index removes the lookup. **It stands alone** — independent of H1, H2 and H4,
byte-identical, and worth 2–3 % on its own at every depth measured except 10×, where it is null.
*Complexity:* one `u32` on the contribution, and one invariant promoted from true-but-unused to
load-bearing — the index is valid only while the active set is untouched since the list was
built, which it already was, and the accessor panics rather than returning `None` so a violation
stops the walk instead of dropping a read.

**H4: [src/ng/locus_generation/pileup/genome_walk.rs:1146](../../../../src/ng/locus_generation/pileup/genome_walk.rs#L1146)
— a depth-sized sort runs at every covered base to answer a question admission already
settled.** *Confidence: High.* **Built and measured: −1.70 % at 30×, −4.65 % at ~130×,
byte-identical.**

- **Hot-path evidence.** `resolve_mate_overlap_at_pos` builds one tuple per contributor, sorts
  it, and returns immediately if no two neighbours share a chain identity. **More than eight
  columns in ten hold no mate pair at either target depth.** Attributed with an
  `#[inline(never)]` wrapper: the sort is 6.8 % of the profiling build's main thread, the tuple
  build 0.28 % — **the sort is 24 times the build**.
- **Mechanism.** Two contributors share a chain identity only if they are the two mates of one
  pair, and whether a pair's alignments overlap is decidable when the second mate is admitted —
  both alignments are in hand and the cross-link is already made there. The active set keeps a
  min-heap of the positions where each pair it holds stops overlapping; when nothing is in it,
  neither the build nor the sort happens. The heap's "no" is exact; its "yes" is an
  over-approximation settled by the sort as before.
- **What pins it.** A debug-build assertion runs the old all-pairs scan on **every column the
  skip claims**. Forcing the predicate to always answer "no" fails 16 tests including the
  production comparison. Five new tests, including one asserting an exact reconciliation count
  so losing a single column still fails.
- **The shape against depth, and it is not monotone.** The gain comes from columns with no pair,
  and that fraction falls as coverage rises: nothing at 10×, −1.70 % at 30×, **−4.65 % at ~130×
  (the peak)**, −1.10 % at 300×.
- *Complexity:* one heap field, one push per paired read admitted, one peek per column, and two
  reset sites to remember.
- Audit: `mate_overlap_sort.md`, `mate_overlap_skip.diff`.

**H5: composing H1 and H4 needs eight more lines, or one eats four fifths of the other.**
*Confidence: High.* **Measured: the plain merge is 1.1 points short at 30× and 2.5 points short
at ~130×; wired together it lands on the arithmetic reference for two independent savings.**

Merged as written, the mate-overlap skip keeps about a fifth of its value — the specialised path
takes 74 % of columns and those are, by its own gate, precisely the no-pair columns the skip was
saving the sort on. The fix: hoist the heap's answer above the specialised path's attempt and
pass it in, so the chain-identity sort inside that path runs only when a pair could be present.
The skip then keeps 93–109 % of its standalone value, because it is applied to two copies of the
same sort instead of one. **The set of columns the specialised path accepts does not move**
(`fast_columns` identical), which is what makes the shortcut safe. *Complexity:* one invariant —
"a 'no' from the heap is exact" — now has three readers instead of one, and through the third it
decides which columns the specialised path accepts. Audit: `composed.md`.

### Note

**N1 — a per-read cursor hint is a regression, and it is now measured rather than assumed.** The
CIGAR cursor is stateless and rescans from the first operation at every position, while the walk
queries each read at monotonically increasing positions. A remembered operation index costs
**+1.44 % at ~130× and +0.85 % at 30×**, ranges disjoint. An Illumina CIGAR under this cursor is
one to three operations, so the scan it skips is one or two integer comparisons while the hint
costs a cell load, a bounds check, a running variable through the loop body and a cell store.
Reverted; recorded so it is not re-run.

**N2 — the profile's largest unattributed line is resolved, and 40 % of it is in a file that
cannot be edited.** `<deduplicated_symbol>` (8.7 %) is three linker-merged functions: 49 % the
mate-overlap sort (H4), 11 % a hash inside the CRAM codec, and **40 % `ChainIdAllocator`'s sweep
of its whole pending-mate map on every read admission** — 954 of 1,060 samples in
`evict_stale_pending`, **3.6 % of the walk, and a term that grows as the square of depth**.
`chain_id_allocator.rs` is pinned byte-identical to production by `copy_fidelity.rs`, so this is
reported rather than proposed. It is the largest single thing this review found and did not
take.

**N3 — the working set fits in L1, so no struct-size change in this scope can pay.** One column
touches 10.2 kB at 30×, 50.4 kB at ~130× and 114 kB at the deepest column seen, against 128 kB
of L1 data cache per performance core. This closes the question from the other side after three
null narrowing experiments in earlier rounds, and predicts a fourth.

**N4 — the sort inside the specialised path can be deleted once the read order is restored, and
should not be.** It is byte-identical and worth −0.18 % to −0.44 %, ranges disjoint. Declined
because it would make that module's quality-sum order — and so the emitted bytes — depend on a
container choice in another file, where today it depends on a sort the module owns. The patch is
kept so the decision can be revisited.

**N5 — the obvious refinement to the specialised path's gate is worth 1.7 points and its premise
is wrong.** Its coverage halves at 300× because any read carrying an indel anywhere in its CIGAR
disqualifies the whole column; making that test exact per read would raise read coverage from
1,669 to 1,835 in 10,000 there — and it is worth *more* at 30× (5.3 points) than at depth, the
opposite of the reasoning that motivated it. What actually halves coverage at 300× is **mate
overlap** (49.6 % of columns against 45.3 %), and handling it inside the specialised path would
be worth 28 points — sixteen times the refinement. Both declined; the second carries a second
copy of the reconciliation tie-breaks.

**N6 — the 2026-08-04 non-wins stand and were not re-proposed.** Emitting only candidate sites
(forbidden by design); the fold table as a sorted `Vec` built the old way (now explained rather
than merely recorded — see H2); struct narrowing three times; three fold-table size estimators;
mimalloc; the lazy record type; caller-side region coalescing; the tandem-repeat scan at 0.17 %
of the walk on real breadth.

**N7 — allocation count is still not the currency**, and the exchange rate held everywhere it
was used this round: ~340 instructions per allocation, so a million allocations is ~1.35 % of
the walk.

---

## 6. Out-of-scope observations

- **⚠ The chain-id allocator's pending-mate sweep is 3.6 % of the walk and grows as depth
  squared** (N2). It runs `retain` over the whole map on every read admission. The file is
  pinned byte-identical to production, so taking it needs the same owner call that released
  `cigar_cursor.rs` and `raw_chrom_reader.rs` — and unlike those, the same waste is in
  production's walk.
- **The human fixture spends two instructions in five getting to the reads, not walking them** —
  613,682 typed regions for 1,541,788 covered bases on `chr1`, 2.5 covered bases per region,
  about 68,600 instructions per region traversed. That is a property of a tandem-repeat-targeted
  BAM rather than of the walk, and it is why every finding here was checked on a second fixture.
  It will disappear with a real whole-genome alignment; until then, region setup is 40 % of that
  fixture's instructions and any ratio taken from it alone is wrong.
- **The acceptance dumps cannot see a walk that stops early.** Two of the four changes were
  caught by tests after the dumps were green (H1's one-step hold, H2's tie-break). Worth a
  standing note in the gate description rather than rediscovering it each round.

## 7. What's already good

- **The cost census is repeatable because the measuring tool reports the dispatcher's counters,
  not its own.** Five agents this round proved their variant did the same work rather than less
  of it, from `loci`, `observations` and `reads_admitted` alone — and the two new counters
  (`fast_columns`, `mate_overlap_positions`) extend the same idea to "is this change even
  firing", which is the failure a timing cannot show.
- **`copy_fidelity.rs` made the divergence decisions explicit instead of accidental.** Two files
  are still pinned, and the one finding this review could not act on (N2) is in one of them —
  which is the mechanism working, not failing.
- **Two of the three per-record containers were already recycled and one was already sized from
  the last record to close**, so the fold table's remaining cost was its hashing rather than its
  allocation — which is why H2's rewrite is a container swap and not a memory-management
  project.

---

## Author responses

Owner decisions taken 2026-08-05, in conversation, the same day the review was written.

| finding | response |
|---|---|
| **Decision 1** — the mate-overlap tie-break's divergence from production (H2) | **Accepted as recommended: name the divergence in the spec.** The tie-break's last resort is an artefact of a container's removal strategy, not a chosen rule; preserving it would mean keeping `swap_remove` so that an arbitrary tie stays arbitrary in the same way. |
| **Decision 2** — the depth cap keeping different reads (H2) | **Rejected, and the underlying rule called wrong** — see below. Not "re-baseline the output", but "the cap is the wrong shape and should be fixed". |
| **The `copy_fidelity.rs` pin** (N2, §6) | **Lifted entirely.** Owner: *"You might copy any file from production and then change it. ng is destined to replace production in the near future."* So N2 — the chain-id allocator's pending-mate sweep, 3.6 % of the walk and growing as depth² — is unblocked, and production is not to be fixed alongside it. |
| **H1, H3, H4, H5** — the byte-identical half | Unblocked; no decision required. |

### Decision 2, restated: the cap is wrong, and the review had it half right

Owner: *"Ideally we should cap all positions at the same capped depth. … what it is wrong is to
leave positions with less coverage because we have discarded reads that cover it."*

**The review reported the symptom against the wrong cap.** Three facts, established after the
decision and recorded here because they change what the fix is:

- **The per-position cap does not lose coverage.** `contributors.truncate(cap)` leaves exactly
  `cap` contributors at a capped position, and a read it drops there still folds at every other
  position it covers.
- **`max_snp_column_depth` = 8,000 can never fire**, because the walk holds at most
  `max_active_reads` = 4,096 reads, so a position cannot gather 8,000 contributors.
  `PileupGeneratorConfig`'s own doc comment at
  [generator.rs:113-115](../../../../src/ng/locus_generation/pileup/generator.rs#L113) says so.
- **So the effective SNP rule is the admission shed, which discards whole reads.** A read refused
  at the door is never decomposed and contributes at **no** position. On the ~130× tomato
  `SL4.0ch01` that refused 19,725 reads of 113,629,764 and pinned `active_reads_high_water` at
  exactly 4,096; one position on that contig has 12,792 reads over it against a local typical of
  86–133. `max_indel_column_depth` = 250 is the only per-position cap that fires at all — 909
  positions on a 300× chr21 walk.

**The fix, in flight, in two parts.** (1) Raise the hold ceiling past the deepest real column, and
with it the chain-id allocator's hard 10,000-entry ceiling on first mates awaiting a partner,
which is what actually binds — that file is now editable. **The ceiling stays** (owner, same
conversation: *"yes but still with a high enough cap, otherwise we could run out of memory"*), so
the target is a number high enough that it never shapes the evidence and low enough to bound
memory, chosen against measured depth and measured peak RSS. (2) Choose which reads a capped
position keeps by a deterministic function of the read rather than a prefix of the container's
order. **Part 2 dissolves the finding that
prompted this decision**: once the kept subset is independent of how the active set stores reads,
H2's reordering changes no output at any depth, and the last blocker on the second half of the
stack disappears. Measured in `tmp/perf_review_2026-08-05_ng-generic-structure/depth_cap.md`.

### The assembled state, and a correction to H2's own evidence

All five changes are in one tree, measured against a pristine build of `6fbbd09`, uncommitted —
`tmp/perf_review_2026-08-05_ng-generic-structure/landing.md`, with `landing.diff` and a suggested
five-commit sequence.

| depth | four performance changes | **all five, with the cap fix** |
|---|---:|---:|
| 10× chr21 | −21.00 % | **−20.67 %** |
| 30× chr21 | −29.21 % | **−28.95 %** |
| **30× chr1 — the target** | −28.55 % | **−28.23 %** |
| ~130× tomato | −42.48 % | **−42.15 %** |
| 300× chr21 | −18.01 % | **−17.32 %** |

The behaviour change costs 0.37–0.84 %. Peak RSS neutral at every depth, including at the deep
spot where the walk now holds 10,747 reads instead of 4,096. All four dumps keep exact line
counts, both SSR dumps are byte-identical, both generic dumps differ by the one
`record_widen_events` header line, and **no locus line moves**. `cargo test --lib` 2,893 pass and
1 fails — `every_divergence_from_production_is_one_of_the_six_named_classes`, on
`record_widen_events` differing by one on an **uncapped** case, which is Decision 1's divergence
and nothing else.

**⚠ H2's row-count evidence was inflated by three orders of magnitude, and the effect is now
zero.** This report and `ordered_active_set.md` state that restoring the read order changed
88,351 of 341,094 rows on a 300× dump. The raw line counts reproduce exactly — and collapse, once
the chain-id column is dropped, to **88 loci of 240,912 whose evidence differs, plus 244 present
in only one dump**: about one locus in 700. Chain ids are minted in admission order and numbered
per file, so ten fewer admissions over a contig renumber every id after them and inflate a line
diff by 67,009 loci. **And with the cap fix in the tree first, the ordered active set changes
zero loci at 300×** — one header line and nothing else — because a capped position's kept reads
are now a function of the reads alone. That is why the suggested sequence lands the behaviour
change fourth and the ordering fifth.

**Two things the assembly caught that no single change could.** The scalar column path returned
before the block that counts positions the ceiling left short, so on an ordinary column the
success counter could not fire and the ceiling-loss bookkeeping was never pruned — caught by the
cap fix's own new test on the first build, fixed with one condition on the scalar path's gate.
And `depth_cap.md` §8's attribution of 1,203 lost columns to the scalar path's length test is
wrong: that window never holds more than 193 reads, so a test against 8,000 cannot reject there.

## Author response convention

Address each finding by its identifier (H1, N4, …) with one of: `applied in <commit>` /
`experiment shows no gain — closing` / `disputed because …` / `deferred to <issue>` /
`won't fix because …`. Four of the five hot-path findings arrive already built, gated and
measured.
