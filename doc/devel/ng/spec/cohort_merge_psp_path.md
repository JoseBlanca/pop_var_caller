# ng — the cohort merge's psp path: summaries first, evidence on demand

*Status: design spec draft, 2026-09-04. No code yet — this settles the design of the work
[`../impl_plan/run_driver_psp_mode.md`](../impl_plan/run_driver_psp_mode.md) reserves a slot for
under "psp-mode performance". It builds on that plan's finished state: `call-from-psps` exists,
and its VCF equals direct mode's byte for byte (that plan's Milestone F). **This spec moves no
byte of output** — it changes what a psp-mode run decodes, holds and allocates to produce the
same bytes. Two of its citations land with the `ng-psp-mode` branch:
[`psp_head_compared_reads.md`](psp_head_compared_reads.md) (the head field this design cannot
run without) and `src/ng/run/psp_source.rs` (the source it upgrades); line numbers to that
branch are as of `d2c7113e`. Parents:
[`cohort_merge.md`](cohort_merge.md) (the merge, path-agnostic),
[`run_streaming.md`](run_streaming.md) §3.3 and §10 (the three-step read and the deferred
"cheap question" seam), [`psp_file_format.md`](psp_file_format.md) §4.3 (the record head).*

---

## 1. What this is

**Today a psp-mode run builds every stored observation into a full record before the merge looks
at it, and about 99 in 100 of those records are then discarded unread.** The first psp source
([`psp_source.rs:31-38`](../../../../src/ng/run/psp_source.rs)) is an adapter over the
build-everything walk, deliberately: the merge's source interface has one method and it returns
a whole built observation
([`ObservationSource::next_observation`,
`observation_cache.rs:70-90`](../../../../src/ng/run/cohort_merge/observation_cache.rs)), so a
source behind it has no way to say anything cheaper.
[`run_streaming.md`](run_streaming.md) §10 records this exact gap — *"the built merge cannot
express that"* — and defers the fix to this document.

**The design: every stored observation reaches the merge in two parts, and the second is built
only on demand.** The *summary* — where the record sits, how far it reaches, its non-reference
read count and its compared-read count — is read off the record's head without decoding the
body. The *evidence* — sequences, qualities, read groups, chain ids — is built only for the
records that overlap a locus the cohort decided to keep. Everything the merge decides before
assembly (grouping, the span verdict, the keep verdict) it decides from summaries; only
assembly (§4.4 of [`cohort_merge.md`](cohort_merge.md)) touches evidence.

**What that is worth, measured on ng's own files (§5):** a walk that skips the bodies a run
does not need is **2.57× faster** than one building everything at the rate a 63-sample cohort
actually needs — one record in eight — and **2.66×** on a high-depth sample at its own rate of
one in fifty-three. Block decompression happens either way and is what remains. The memory side is unmeasured and
expected to be the larger win at scale: a built record is a heap object with several
allocations; its raw encoded bytes are about 5 a record at that depth on disk, and what this
design holds per sample over the in-flight ground is summaries and raw bytes, with built
records existing only for the roughly one locus in a hundred that survives.

### 1.1 Goals

1. **The same VCF, byte for byte.** Direct mode and psp mode already agree
   ([`run_driver_psp_mode.md`](../impl_plan/run_driver_psp_mode.md) F2); every change here keeps
   that oracle green, and the decode-everything source survives as the test oracle the two-phase
   source is compared against, record for record — production's eager-decode pattern
   ([`sample_reader.rs:20-26`](../../../../src/var_calling/sample_reader.rs)).
2. **One merge, not two.** The merge's grouping, verdicts, ownership and organiser stay one code
   path for both modes. What varies is the cached item's state, not the algorithm (§3.1) — the
   same shape that kept the calling loop single when the psp source arrived
   (`run_driver_psp_mode.md` D2).
3. **Bounded across the committed range** — one sample to several thousand, three reads a
   position to several hundred. Every retained quantity in §3.3 is a formula in sample count and
   in-flight ground, and §6 states what each end of each axis does to it.

### 1.2 Non-goals, and what this does not do

- **No new subcommand, no new flag.** `call-from-psps` exists and its interface does not change.
- **No change to any verdict.** The keep rule, the span bound, grouping, projection,
  unification — all fixed by [`cohort_merge.md`](cohort_merge.md). This document changes when
  evidence is materialised, never what is decided.
- **No change to the psp format.** Its one format prerequisite — the head carrying the keep
  rule's denominator — is [`psp_head_compared_reads.md`](psp_head_compared_reads.md)'s, and
  landed on main on 2026-09-04 (`47d4c7e1`), costing 3.9 % of the compressed file at 10 reads
  a position and 8.5 % at 280 (`ebf8d2f5`, re-measured in `fe9df2a3`).
- **It does not parallelise the merge.** The organiser's serial cover phase and the 1.4×-on-8-
  threads wall are measured facts of the direct path
  ([`cohort_merge.md`](cohort_merge.md) §6.2); whether the psp path's decode belongs on several
  threads is settled by measurement *after* this design is built and timed (§5), not designed
  here on suspicion.

## 2. Where it sits, and why the head can answer the merge

The run shape is fixed upstream: reading one sample's evidence is three steps — a stored
summary, the cheap per-position numbers, the full evidence — and *"the merge puts every
sample's answers together and decides which positions are worth calling"*
([`run_streaming.md`](run_streaming.md) §3.3; its step 1, a per-stretch skip without
decompression, was removed from the format on 2026-08-25 and is dead text — the cheap read is
step 2, after decompression).

**The head's numbers are the verdict's inputs, not an approximation of them.** The merge's keep
verdict calls, per member record, `non_reference_and_compared_reads()` and feeds the pair to
`MinAltReads::reached_by`
([`close.rs:665,713`](../../../../src/ng/run/cohort_merge/close.rs);
[`mod.rs:526`](../../../../src/ng/run/cohort_merge/mod.rs)). The head's two count fields are
defined as that function's two return values, derived by the writer at encode time and checked
against the rebuilt body at decode time
([`psp_head_compared_reads.md`](psp_head_compared_reads.md) §3). So a fold over summaries and a
fold over built records read identical numbers by construction — the byte-identity of the two
modes' VCFs needs no tolerance argument, only the record-equality oracle of goal 1.

Grouping needs only each record's region, which the head also carries.

**The locus's kind is not part of a summary, and an earlier draft of this section was wrong to
imply it.** That draft said the width verdict and the never-mix assertion "need the record's
kind before anything is assembled", which reads as though every observation must answer for its
kind. Neither does. `max_cohort_locus_span` governs generic loci and not repeat tracts, whose
span the reference fixes — but that verdict is passed **once on the closed locus**, from the
kind the opening observation carried
([`close.rs:114-124`](../../../../src/ng/run/cohort_merge/close.rs)), not once per member. The
only per-member read is a release assertion that a locus never mixes the two, guarding
something segments already guarantee. So a summary is three facts and a region, and carrying a
kind in it would have cost every observation a field to serve a check that cannot fire.

**What a run whose evidence stays undecoded does for that one per-locus kind is C's to settle,
and the owner's ruling of 2026-09-04 (`fe9df2a3`) already supplies the answer it will use**: a
locus's kind is the kind of the typed region it falls in, typed regions come from the reference
and the catalog before any read is looked at, and every psp records the segmentation inputs its
typing used against a run that refuses a cohort disagreeing with its own. So a coordinate is
enough to look one up — once per locus, not once per record. The tag stays at the end of the
body, where the body decoder needs it to know whether a tract's motif and flanks follow.

The one thing summaries cannot answer is assembly, which is the point.

## 3. The design

### 3.1 The cached observation has two states, and the merge stops caring which

The observation cache's item today is a built `SampleLocusObservations`. It becomes a
two-state item:

- **summary available, evidence built** — direct mode's only state: the walker minted the whole
  record, and its summary is read off it (the same derivation call the verdict makes today,
  [`close.rs:665`](../../../../src/ng/run/cohort_merge/close.rs));
- **summary available, evidence stored** — the psp state: the summary came from the head, and
  the evidence is raw bytes the source still holds. Asking for the evidence builds it, once,
  and the item moves to the first state.

Everything before assembly — the fold, grouping, both verdicts, ownership — reads summaries.
Assembly asks for evidence, and in direct mode that ask is free. This is
[`cohort_merge.md`](cohort_merge.md) §1's *"one step in the direct path and two in the psp
path"*, made a type instead of a sentence.

**Building is per record and pure**, so nothing about which records get built, or when, can
move the output: the set of kept loci is a function of the summaries
([`cohort_merge.md`](cohort_merge.md) §9), and a built record equals the record the
decode-everything walk would have produced — the oracle of goal 1 checks exactly that.

### 3.2 The psp source: one reader, two cursors, and builds in coordinate order

The upgraded source replaces the decode-everything adapter. Per sample, over one forward
reader:

- **The summary cursor** walks record heads, applying each record's chain-id changes to its
  live set as every reader must
  ([`record.rs:1881-1902`](../../../../src/ng/psp/record.rs)), and hands each record's summary
  plus its raw bytes — head and unbuilt body — into the retained window. It never builds a
  body. This is the measured 0.141 s walk.
- **The retained window** is the raw bytes of every record between the two cursors, per
  sample — a FIFO the cover advances and eviction drains, the same rhythm the cache already
  has ([`observation_cache.rs:113-118`](../../../../src/ng/run/cohort_merge/observation_cache.rs)).
- **Building a body needs nothing but its own bytes**, so any thread may do it, at any time,
  in any order. That is a property of the format rather than a convenience of this design: a
  psp body writes every count absolutely and every observation's read list from zero rather
  than as a difference from the record before it, **which is the reason a skipped body costs
  nothing in the first place** ([`record.rs:828-846`](../../../../src/ng/psp/record.rs)). A
  reader that never saw the records in between has missed nothing it needs.

**So there is no build cursor, and an earlier draft of this section invented one.** It had a
second cursor trailing the first, replaying each skipped record's chain-id changes into its own
live set before building — and argued at length why replaying beat storing a snapshot per
record. Both were answers to a problem the format does not pose. **The state a reader carries
across records lives in the head, not the body**, and every reader parses every head whether or
not it wants the record ([`record.rs:1881-1902`](../../../../src/ng/psp/record.rs)); the one
thing the body decoder is handed from outside is the live set the head walk already produced.

**Which resolves the concurrency question this design would otherwise have had.** The merge
hands several builders one shared window at a time and they run together, so a build that
needed mutable state — a cursor to advance, a live set to replay into — would need either a
lock per sample or a serial build phase between deciding and assembling, and the second would
put decoding, which is 43% of a psp-mode run at 63 samples (§5), on one thread. Neither is
needed: the retained bytes are immutable, a body is built from them alone, and builders share
nothing but read-only memory. **The one piece of state a run must still carry is the live set,
and the head walk carries it** — so what the window retains beside each record's bytes is the
set as the head walk left it, at the one place it is asked for.

⚠ **This rests on the body staying self-contained, and today that is under-tested because the
chain ids are not yet written.** `encode_record_body` drops them — a record read back has empty
lists where it had ids ([`record.rs:830-835`](../../../../src/ng/psp/record.rs)) — so the live
set is inert in every file this caller currently produces, and the residual read list a body
derives from it is trivial. When the encoding's Milestone E writes them, the body's own
guarantee is what must continue to hold; if it does not, this section's argument fails and the
lock-or-serial-phase question returns. **The test that would catch it is the record-equality
oracle of §1.1 goal 1 run on a file whose records carry chain ids**, and it does not exist
because such a file cannot yet be written.

**What this replaces:** the adapter's refusal of head-only records
(`ObservationBodyNotBuilt`, [`psp_source.rs:95-112`](../../../../src/ng/run/psp_source.rs))
exists because the current interface cannot carry a deferred body; in this design a head-only
record is the normal case and that variant goes, its job taken by the build-order refusal
above. The reader-side machinery both cursors need exists: the head-only walk and the
build-some walk are `RecordIter` and `building_only_where`
([`walk.rs:215`](../../../../src/ng/psp/walk.rs)).

### 3.3 What is retained, and when it goes

Per sample, between eviction and the cover frontier:

> summaries (a fixed few dozen bytes each) + the records' raw encoded bytes + one live set per
> cursor

and built records only where assembly asked. The raw bytes replace the built records the cache
holds today over the same ground, and a built record costs a multiple of its raw bytes that
nobody has measured — [`cohort_merge.md`](cohort_merge.md) §8 marks the built size unmeasured,
and the raw size is about 5 bytes a record compressed at three reads a position
([`psp_record_encoding.md`](psp_record_encoding.md) §11; the decompressed encoded size is also
unmeasured — the plan measures both). Eviction is unchanged in shape: when the organiser
releases ground, the source drops the retained bytes behind the window's left edge.

**The alternative that lost: two readers per sample** — a head-only reader ahead and a
building reader behind, no retained window at all. It loses on arithmetic: each reader
decompresses every block, and decompression is 0.104 s of the 0.141 s walk, so the second pass
gives back most of what the skip won — about 0.27 s a sample against the single-reader design's
roughly 0.17 s, on the corpus above. It would win only if the retained window's memory turned
out to matter more than the time, which §6's formulas say it does not at either end.

### 3.4 Where the retained bytes live — settled, and not the implementer's after all

**One growing byte arena per sample, records appended and addressed by range.** An earlier
draft left the choice open between that and a box per record, leaning towards the boxes and
saying the arena could wait for an allocator profile. **That leaning was wrong, and the reason
is the measurement this design is bought with.**

The 2.57× a skipping walk gives (§5) is the speed of a walk that reads each head, advances past
the body, and **keeps nothing** — no allocation anywhere in it. A box per retained record puts
an allocation and a copy back on every record in the window, which is the per-record cost the
skip exists to remove: what would be saved is the record's several allocations, and what would
be spent is one, on every record rather than on the one in eight that gets built. The arena
spends no allocation per record at all — the append is a memcpy into a buffer that is already
long enough after the first window — so it keeps the shape the measurement was taken on.

**The constraint that was already right stands: eviction must return memory.** An arena that
only grows is the cache leak this module exists to avoid. Since the window advances in
coordinate order and drops a prefix, what an arena needs is the same prefix drain the record
window already does — the bytes behind the window's left edge go, and what survives moves down.

*Unmeasured, and the first thing to look at if the memory is wrong: whether the drain's move
costs more than it saves at large windows. It is the same move `held_observations` already
makes, over bytes rather than records.*

## 4. The run-level companions

Two items travel with this work because they live in the same objects, not because the
two-phase read needs them:

- **One contig list for the run.** On a human reference 357 kB of an open sample's 480 kB is
  a copy of the reference's contigs, identical in every sample — 1.07 GB of the same list at
  three thousand ([`run_streaming.md`](run_streaming.md) §10, with the owner's 2026-08-30
  ruling: the *file* keeps carrying its own list; what the *readers* work from is the run's).
  `PspVariantCaller::open` already reads every header, so it checks the lists agree and hands
  every reader one shared list.
- **The spare offer, taken.** The merge hands every source a spent record for reuse and the
  psp adapter currently drops the offer
  ([`psp_source.rs:23-29`](../../../../src/ng/run/psp_source.rs);
  [`observation_cache.rs:58-69`](../../../../src/ng/run/cohort_merge/observation_cache.rs)
  names a decoder as the hook's best customer). With builds now rare, the offer matters most
  exactly where records are large: at depth, a built record's buffers are refilled instead of
  reallocated. Adopted only if measured — the merge's 39% Linux allocator share is the number
  that says it might pay ([`../research/cohort_merge_parallel_cost_plan.md`](../research/cohort_merge_parallel_cost_plan.md) §2).

## 5. Performance: what is known, what is measured first, what is gated

**Measured 2026-09-04, on ng's own store — the plan's Milestone A, and it replaces the
prototype figures this section carried.** A psp-mode run had never been timed; both corners of
the range now have been, against direct mode over the same ground, the same cohort and the same
parallel path on the same day:

| | 63 tomato accessions, 200 kb, ~3 reads a position | HG002, 76 kb, ~280 reads a position, one sample |
|---|---|---|
| calling phase, direct mode | 13.91 s | 3.34 s |
| calling phase, psp mode | 2.05 s | 0.05 s |
| reading the records back | 881 ms (42.9%) | 47.6 ms (91.2%) |
| assembling the loci | 392 ms (19.1%) | 2.3 ms (4.4%) |
| genotyping them | 713 ms (34.7%) | 1.4 ms (2.7%) |

**Calling from stored files is 6.8× faster than calling from alignments at the cohort corner
and 67× at the single high-depth sample**, and the run's centre of gravity moves: at 63
accessions, reading records back is 43% of the work where decoding reads is 81% in direct mode,
and the merge plus the genotyper rise from 18% of a direct run to 54% of a psp one.

**What the skip is worth, and the prediction it refutes.** This section used to say the skip was
2.06× and that *its value shrinks with depth*, because the chain-id changes ride in the head and
grow from 0.432 bytes a position at 11.4 reads to 6.42 at 293
([`psp_record_encoding.md`](psp_record_encoding.md) §6). **Measured on ng's own files, it does
not shrink** ([`examples/ng_psp_skip_value.rs`](../../../../examples/ng_psp_skip_value.rs)):

| store | reads a record | bodies kept | skipping against building everything |
|---|---|---|---|
| tomato accession | 8.5 | 1 in 8 — **the rate a 63-sample cohort actually needs** | **2.57×** |
| tomato accession | 8.5 | 1 in 16 | 2.90× |
| tomato accession | 8.5 | 1 in 100 | 2.91× |
| HG002 | 281.9 | 1 in 53 — the rate one sample actually needs | **2.66×** |

**The keep rate is per record, not per locus, and that correction matters more than the depth
axis.** A kept cohort locus needs *every covering sample's* record built, so at 63 accessions
24,538 kept loci ask for roughly 1.5 million of the 12.1 million records drawn — **one record in
eight, not one in a hundred**. The design was sized against the locus rate in an earlier draft;
at the true rate the skip is still 2.57×, which is why it survives the correction.

**So the saving this design can take is about a quarter of the calling phase at the cohort
corner** — 43% of the work running 2.57× faster — **and over half of it at the single
high-depth sample**, where reading records back is 91% of a very short call. Both are estimates
composed from a bare walk's speed-up applied to a run's measured share, not end-to-end
measurements of the built design; Milestone E is what replaces them.

**Known of the merge:** its parallel driver gives 1.4× on 8 threads because the organiser's
cover runs serially between rounds — tolerable in direct mode behind a generator 14–23× its
cost, and much more exposed in psp mode, where the generator is gone
([`cohort_merge.md`](cohort_merge.md) §6.2).

**What is still gated on measurement:** per-sample cover parallelism, overlapping the cover with
building, and the rest of
[`../research/cohort_merge_parallel_cost_plan.md`](../research/cohort_merge_parallel_cost_plan.md)'s
psp half fold into this work **only if the share of a run they could recover says so**, and
[`run_streaming.md`](run_streaming.md) §11 question 7 is where the conclusion is owed either
way.

## 6. Degradation at the edges

- **One sample.** The two-phase read stands on its own: the fold is one sample's summaries and
  the keep rule is unchanged at k = 1 ([`cohort_merge.md`](cohort_merge.md) §7.2). What is
  saved is the same 2× walk; what is retained is one sample's window.
- **Three thousand samples.** Retained bytes scale as `samples × window`, the same product the
  built cache already pays with a larger per-record constant; open readers are
  `samples × 123 kB` once the contig list is shared (against 480 kB before —
  [`psp_file_format.md`](psp_file_format.md) §5.2). The plan reports peak resident beside wall
  at the sweep's top end.
- **Three reads a position.** The skip's best corner: bodies dominate records, one in a
  hundred is built.
- **Three hundred reads.** Predicted to be the skip's worst corner, twice over — the head is
  most of the record, and error alone clears the keep rule's floor at about 4 positions in 100
  so more loci reach assembly. **Measured, it is not** (§5): on HG002 at 281.9 reads a record,
  keeping the 1 record in 53 that sample's 1,439 kept loci actually need, the skipping walk is
  **2.66×** a full one, against 2.57× on tomato at its own true keep rate. What does grow at
  depth is the body: 40.97 bytes a record against tomato's 23.17, so each body skipped is worth
  more, which offsets the head's growth rather than being swamped by it. The design's floor
  argument stands unchanged either way — it degrades to what the decode-everything source
  already is, plus one head walk.

## 7. Cross-cutting concerns

**Errors.** Two new refusals, both naming the sample: a backwards build ask (§3.2), and a build
ask for ground already evicted — both are merge bugs surfacing, not file damage, and say so.
File damage keeps its existing shapes (the bounded body decode, the head-body checks).

**Concurrency.** Unchanged: one reader per sample, drawn by the organiser; builders read, never
draw. Nothing here adds a lock, and the source stays single-threaded per sample so that the
gated cover-parallelism work (§5), if it comes, parallelises *across* samples.

**Memory.** §3.3's formula, measured in the plan with the dhat probes that already exist
([`examples/dhat_psp_reader.rs`](../../../../examples/dhat_psp_reader.rs)).

## 8. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| the head walk and the skip | `RecordIter`, `building_only_where` ([`walk.rs:168,215`](../../../../src/ng/psp/walk.rs)) | the summary cursor is the head walk, which also keeps the live set; a body is built from its own retained bytes and needs no second walk |
| the bounded body build | `decode_the_body_of` ([`record.rs:1984`](../../../../src/ng/psp/record.rs)) | called once per built record, unchanged |
| the verdict's numbers | `non_reference_and_compared_reads` → `reached_by` ([`close.rs:665,713`](../../../../src/ng/run/cohort_merge/close.rs)) | the summary carries the same pair; the verdict code does not change |
| the source seam and its errors | `PspObservationSource`, `PspSourceError` ([`psp_source.rs`](../../../../src/ng/run/psp_source.rs)) | upgraded in place; the decode-everything form survives as the record-equality oracle |
| cover / evict rhythm | `ObservationCache` ([`observation_cache.rs:113`](../../../../src/ng/run/cohort_merge/observation_cache.rs)) | the retained window advances and drains on the same calls |
| the mode-equivalence oracle | `run_driver_psp_mode.md` F2 | rerun green after every milestone — the definition of "moves no byte" |

## 9. Deferred, with a recommended home

- **Parallelising the psp path's cover or decode** — gated on §5's timing; its conclusion is
  owed to [`run_streaming.md`](run_streaming.md) §11 question 7, and if built it is its own
  plan.
- **Dropping the record's stored reference bases and re-fetching from the run's reference**
  ([`run_streaming.md`](run_streaming.md) §11 question 4) — untimed leaning, untouched here;
  builds get rarer under this design, which moves that trade and is worth saying when it is
  finally timed.
- **The callers-in-flight default** — [`run_streaming.md`](run_streaming.md) §11 question 2's
  open half, unchanged by this work.

## 10. Open questions

1. **The retained window's storage shape** — §3.4, the implementer's, leaning stated.
2. **Is the spare offer worth taking?** — §4; settled by the allocator share in the plan's
   before/after profile, not by argument.
3. **Does anything structural get built for scaling?** — §5; settled by the end-to-end shares,
   answered into [`run_streaming.md`](run_streaming.md) §11 question 7.
