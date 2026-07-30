# ng — the generic locus generator: a non-STR stretch → many loci

*Status: design spec, 2026-07-26. **No code yet.** The second
[`LocusGenerator`](locus_generation.md) and the last loci-minting mechanism — it consumes the
`Generic` region kind and produces one locus per covered reference position. Inherits the locus
type, the contract, the dispatch and the error model from
[locus_generation.md](locus_generation.md); read that first. It is the only consumer of read
preparation ([read_preparation.md](read_preparation.md)), so it is also what validates that step.
Reuse target: production's `src/pileup/walker/`. Naming: **STR** in prose, `ssr` in code.
Code-facing companion: [`../arch/locus_generation_pileup.md`](../arch/locus_generation_pileup.md)
(the types, the signatures, and the reconciliation table). Build order:
[prerequisites](../impl_plan/locus_generation_pileup_prerequisites.md) →
[the port](../impl_plan/locus_generation_pileup_port.md) →
[the generator](../impl_plan/locus_generation_pileup_generator.md).*

---

## 1. Scope — goals, non-goals, and what it does not do

**What it does.** For one stretch of genome that region typing ruled non-repetitive, walk it base
by base over the sample's reads and, at every position where a read said anything, emit a locus:
the distinct sequences the reads showed there, with their support.

This is the *data-defined* half of the locus stream ([ng_proposal.md](ng_proposal.md) §1). An STR
locus is minted from the reference before a read is opened; a generic locus does not exist until
the reads are read.

**Goals.**

1. Fill `GeneratorSet`'s `generic` slot, which today holds
   `GeneratorSlot::Unfilled(UnhandledReason::NotImplemented)`
   ([locus_generation/mod.rs:377](../../../../src/ng/locus_generation/mod.rs#L377)) — so ng can
   mint a SNP or indel locus at all, and the shared per-locus core downstream has a generic input.
2. Get production's walker into an in-memory generator ng owns outright, **provably without
   changing what it computes** (§3), so that whatever ng later changes is a measured delta and not
   an accident.
3. Settle the four things `locus_generation.md` §11 hands to this step: `chain_ids`, the
   read-position-bias fields, mate-overlap reconciliation, and the shape of the parity oracle (§6,
   §5, §3).
4. Say what a per-position record cannot say, and add it: which reads saw only part of a locus
   (§6).

**Non-goals — could be goals, deliberately are not.**

- **A different definition of "locus".** `ng_proposal.md` §1 names the definition itself as a
  swappable axis — one position, an active-region window, a haplotype window grown to a fixpoint.
  This generator is the **per-position** implementation, the one with a parity oracle. The others
  are siblings in the same folder, later (§10).
- **Local reassembly.** Settled, not deferred: the generic path has no GATK-style local
  reassembly, because production already beats GATK without one.
- **Calling.** Evidence, not genotypes — no candidate alleles, no likelihood, no QUAL.

**What it does not do.** No cohort merge; no windowed depth/GC statistics; no `SsrBundle`
generator; no parallelism; no BAQ. All five are deferred elsewhere with a home
([locus_generation.md](locus_generation.md) §11, [read_preparation.md](read_preparation.md) §10)
and none is re-opened here.

---

## 2. Where it sits — one region, one walk, one locus per position

The dispatcher hands this generator a `Generic` region and the sample's reads. Its `begin_segment`
starts a walk; each `next_locus` returns the next closed record, or `None` when the walk is done.
The `Generic` region kind carries no payload, so the segment type is `()` and the region arrives
through `begin_segment` alone.

```
GenomeRegion ─▶ SampleReads::reads_in_region ─▶ ReadPreparer ─▶ the walk ─▶ locus, locus, locus…
                (AlignedRead)                    (PreparedRead)   (§4)
```

**One read *query* per segment — and a query is a lazy stream, never a batch.** This is worth
stating flatly because the opposite reading is alarming and wrong. `reads_in_region` returns a pull
iterator ([read/input/mod.rs:589](../../../../src/ng/read/input/mod.rs#L589),
`Iterator<Item = Result<AlignedRead, IngestError>>`); it collects nothing. The walk pulls a read
only when it reaches that read's start, and drops it once it passes the read's `alignment_end`
([active_read_set.rs:149](../../../../src/pileup/walker/active_read_set.rs#L149)). **What is
resident is the reads overlapping the current position — order (local depth × read length), hard
-bounded by `max_active_reads` — never the region's reads.** A 100 kb `Generic` region is walked
with the same footprint as a 1 kb one; the region's length sets how long the walk runs, not how
much it holds. Holding a segment's reads at once would indeed be unaffordable, which is why the
walker consumes an iterator and why that property has to survive the port (§7).

*Per-query cost, and one piece of it that matters here.* The STR generator issues ~10⁶ region
queries and drove three fixes for it on `main` — the CRAM container cache (`b918fb6`), resolving a
record's footprint once per decode rather than once per query (`70003c2`), and sharing one
reference reader across queries (`4bc3ef9`, which also closes `locus_generation.md` §8 item 1's
"Arc gap" by forwarding `RefSeq`/`RawRefSeq`/`ContigTable` through `Arc<T>`). **One query per region
means this generator pays almost none of that.** *(All three are merged into this branch as of
2026-07-27, together with the read-group work — see the arch doc's branch-state note.)*

*`4bc3ef9`'s reposition threshold is not a hazard for this generator — checked, 2026-07-27.*
`RawChromReader::fetch` repositions instead of extending when a window lies more than
`REPOSITION_GAP = FILE_READ_CHUNK = 64 KiB` clear of the buffered one, in either direction. The
pileup lands on the *extend* side by construction — a contiguous forward walk, which is the case the
threshold exists to preserve — and on the *reposition* side only where it genuinely jumps: a
coverage gap the walker skips ([driver.rs:549](../../../../src/pileup/walker/driver.rs#L549)), the
STR or satellite regions lying between two `Generic` ones, or a BED-restricted run. Repositioning is
the wanted answer in all three. Two margins make it comfortable rather than marginal: record
widening cannot approach the threshold (`max_record_span` is 5000), and read preparation's fetches
lag the walker's by at most a read length, so even a shared accessor would extend rather than
alternate repositions — and they are not shared, each holding its own. The one thing that *would*
hurt is rebuilding the accessor per segment (§8).

**A locus per covered position, REF-only ones included.** Production emits a record at every
position at least one read observed, including positions where every read matched the reference
([walker/tests.rs:158](../../../../src/pileup/walker/tests.rs#L158)); uncovered positions emit
nothing ([:466](../../../../src/pileup/walker/tests.rs#L466)). ng keeps that, because a cohort
merge must be able to tell "this sample is REF here" from "this sample has no coverage here", and
because per-position depth is what the windowed statistics slide over
([locus_generation.md](locus_generation.md) §11).

**Why it has to be every position, not just the variant-carrying ones — settled (owner,
2026-07-26).** ng splits the per-sample stage from the cohort stage, so a sample's evidence is
gathered *before* anyone knows which positions will turn out to carry a candidate allele. If this
generator emitted only the positions where *this* sample showed something non-reference, then a
position where another sample carries a candidate would have no data for this one — and the SNP
call at that position could not be made. The per-position record is what makes the two-stage split
possible; a candidate-only stream would force the cohort back into the per-sample pass.

**The cost, stated plainly.** One `SampleLocusObservations` per covered base — order 10⁹ loci on a
resequenced human sample, each owning a `GenomeRegion`, a `Box<[u8]>` of reference bases, and a
`Vec<SequenceObservation>`. The streaming contract keeps one resident at a time
(`locus_generation.md` §4), so this is a throughput cost, not a memory one — and there is no
fallback design to fall back to (§7).

**Region boundaries: clamp on the anchor, and query with a halo.** The read query returns every read
*overlapping* the region, including reads starting before it, so the walk emits records just outside
the region's bounds. Drop those on the record's anchor position, exactly as production's region
driver does (`(start..=end).contains(&record.pos)`,
[pileup_to_psp.rs:271](../../../../src/pileup/per_sample/pileup_to_psp.rs#L271)). Because typed
regions tile the genome gap-free and disjointly, every record is emitted by exactly one region — no
duplicates, no holes.

**The evidence behind those records is a separate question, and the naive query gets it wrong.** A
record anchored inside a region can have a footprint reaching up to `max_record_span` (5000) past
the region's end — a long deletion does exactly that. Reads that fold into it may lie **entirely
beyond the region**, and a query for "reads overlapping the region" never returns them. The record
is then emitted once, by the right region, with **part of its support missing** — and no counter
notices, because the record itself is not lost. So the query is
`[region.start, region.end + max_record_span]`: reads overlapping *that*. The extra reads are walked
and their records clamped away unless anchored in the region, so the halo costs bounded work at each
boundary and buys correctness at exactly the loci — long indels next to a repeat — where it matters
most.

**But the halo must be *stopped*, not just queried.** The walk has no right bound: it runs until the
reads are exhausted and the active set is empty, so a naive halo walks all 5000 extra positions at
full depth, finalises every record in them, and throws them away at the clamp. With regions tiling
the genome that tax can exceed the region interiors. **Stop as soon as `walker_pos > region.end` and
no open record is anchored at or before `region.end`** — nothing anchored later can survive the
clamp, so the halo then costs only what the boundary-crossing records actually need.

*What the halo does not change, because it is not a boundary effect:* the interior positions of a
wide record get no locus of their own, in any region. An event inside an open record's footprint
folds into that record rather than opening its own (`find_overlapping`,
[open_record.rs:303](../../../../src/pileup/walker/open_record.rs#L303)), so "a locus per covered
position" has always meant "per position not already inside a wider record". That is production's
behaviour and `locus_generation.md` §3 already records the overlapping-loci caveat it creates.

---

## 3. How the port happens — copy the walker, prove it identical, then change it

**The rule this step is written under:** production is frozen — ng reads it, reuses what is already
`pub`, and copies into `src/ng/` otherwise. **As it turns out this step needs no exception at all**
(§3's decision). The owner granted one on 2026-07-26 — a visibility-only lift of
`src/pileup/walker/` — and it is no longer needed; it is recorded here only so nobody re-derives the
need for it. One structural fact decides the rest of this section.

**The walker is not entangled with `.psp`.** Grepping every file under `src/pileup/walker/` finds no
dependency — no `Write`, no `serde`, no encoder, no `psp::` import. (`column` matches, in
`max_snp_column_depth` and friends, and `errors.rs:42` names `.psp` in a message: vocabulary, not
coupling.) The two real couplings are one-line policies, not machinery: `finalise()` sets
`windowed_gc`/`windowed_coverage` to `f32::NAN` for the writer's seam to fill
([open_record.rs:180](../../../../src/pileup/walker/open_record.rs#L180)), and it drops the REF
bucket's chain ids ([:150](../../../../src/pileup/walker/open_record.rs#L150)). ng re-takes both
(§6); neither is an obstacle.

**It is, however, mostly sealed** — which is why there was never a middle route between copying it
and driving it whole. `walker/mod.rs:15-22` declares `active_read_set`, `chain_id_allocator`,
`cigar_cursor`, `decompose`, `driver`, `open_record` and `errors` **private**; the re-exports are
`run`, `PileupWalker`, `RunSummary`, `WalkerError` and `DEFAULT_MAX_ACTIVE_READS`; and two pieces
this step cares about are private `fn`s *inside* those private modules (`base_in_adaptor`,
[cigar_cursor.rs:85](../../../../src/pileup/walker/cigar_cursor.rs#L85);
`resolve_mate_overlap_at_pos`, [driver.rs:648](../../../../src/pileup/walker/driver.rs#L648)). Not
all of it is sealed — `indel_norm` is `pub(crate)` and ng's read preparation already calls it, which
is why it is absent from the copy inventory below. So `run` is the only copy-free way in: the
black-box option, rejected next.

**Rejected outright: driving production's walker as a black box.** Calling the public
`run(reads, fetcher, config)` and adapting each `PileupRecord` needs no lift at all — but it
**cannot fill the locus type**. Three of `SampleLocusObservations`' six fields have no recoverable
source in a `PileupRecord`: `read_witness` is gone by `finalise()`, which consumes the per-read
fold state, and `reads_without_observation` / `reads_discarded_by_cap` exist only as run-level
totals on `RunSummary`, never per record. Worse, the thing `read_witness` exists to prevent is
already baked in: a read covering only part of a record's footprint has the rest filled from the
reference ([open_record.rs:522-531](../../../../src/pileup/walker/open_record.rs#L522),
[:568-573](../../../../src/pileup/walker/open_record.rs#L568)), folded as if it had witnessed bases
it never saw (§6).

**Chosen: copy the whole walker into ng, and touch production not at all (owner, 2026-07-28).**
ng needs the read to carry its library (§6), so ng needs its own read type — and every module that
looked reusable names `PreparedRead` in its signatures, so that type reaches all of them anyway.

| what | where it lands |
|---|---|
| `driver.rs`, `open_record.rs`, `cigar_cursor.rs`, `decompose.rs`, `active_read_set.rs`, `chain_id_allocator.rs`, `errors.rs` (~5,500 lines) | **copied verbatim** into `src/ng/locus_generation/pileup/`, then changed. One rename on the way in — `driver.rs` → `genome_walk.rs`, the only one of the seven named for a role rather than for what it owns (arch, *Module home*) |
| `PreparedRead`, `MateRole`, `ReadLengthError` | **copied into ng** and extended there with `read_group` (§6) |
| `CigarOp`, `PileupRecord`, `AlleleObservation`, `AlleleSupportStats` | reused as-is — already `pub`, and unchanged by ng |

**What it buys.** *Production is edited zero times* — no visibility lift, no field on a frozen type,
no relocation of `ReadGroupId` to reach it. The freeze holds absolutely, which also means no
production review and no risk to the shipping caller. And because ng is a lab, the read type will
change again (BAQ if it returns, the re-align mode); each of those would otherwise be a fresh edit
to frozen code. Against that, ~5,500 duplicated lines — which cannot drift, because the source they
came from does not move.

**Transcribe first, change second — now over the whole walker.** All seven files land verbatim,
still emitting `PileupRecord`, and are proven to compute exactly what production computes before a
line is edited. That is the rule that paid three times on this branch — `delimit_read` (200,000
randomized cases, zero divergences), `left_align_indels`, the SplitMix64/reservoir port. Copying
everything also makes the differential **complete** — it tests every line ng will actually run, with
no half taken on trust.

**The parity oracle, in two stages.**

*Stage 1 — the port is identical.* A `#[cfg(test)]` differential in the `delimit_parity` /
`left_align_parity` shape already used twice in `src/ng/`: build one `Vec<PreparedRead>`, hand it
to `crate::pileup::walker::run` and to the ng copy's `run` with the same `WalkerConfig` and the
same reference fetcher, and assert the two `Result<PileupRecord, WalkerError>` streams are equal
element for element, plus `RunSummary`. "Byte-identical" is well defined here because the output
type is still production's, and `PileupRecord`'s hand-written `PartialEq` compares the two `f32`s
**by bits** ([pileup_record.rs:208](../../../../src/pileup_record.rs#L208)), so the `NaN`
placeholders compare equal and the comparison is total. Run it on the walker's own fixtures
(`MockFasta`, `snp_read`, `paired_snp_reads`,
[walker/tests.rs:32-152](../../../../src/pileup/walker/tests.rs#L25)) and at real scale on GIAB
HG002 and a tomato CRAM, under the `PVC_PARITY_CASES` convention.

*Stage 2 — the ng-shaped output is a projection plus a named list of divergences.* Once the copy
emits `SampleLocusObservations`, parity runs through a total `project(PileupRecord) ->
SampleLocusObservations` in the test module. **An earlier draft claimed exactly two divergence
classes; there are six**, and naming them all is the point — an unlisted class gets triaged as a
listed one, which would contaminate the measurement §12 exists to produce.

The count has moved twice, and how the sixth was found is worth one sentence, because the same
method is the only one that will find a seventh: **D1 built the permanent anchor, and the anchor
failed on its first run.** Not a count arrived at by reading the two walkers — a count arrived at
by asserting equality where equality was believed to hold and being shown a locus where it did
not.

| divergence | why | how it is checked |
|---|---|---|
| a read that did not witness the whole footprint becomes an observation carrying an `Observed` witness, whose bases are what it saw, where production folded a reference-filled haplotype into REF | production fabricates the unwitnessed bases (§4, §6) | the divergent loci are enumerated and each hand-verified against the read's CIGAR; **this count is the deliverable** |
| an allele supported from several read groups is several observations | the group joins the key (§6); `project` cannot split a production observation | observation counts reconcile: summing a locus's observations by `(bases, coverage)` must reproduce production's per-allele totals |
| `reads_without_observation` / `reads_discarded_by_cap` are non-zero | production keeps neither per record | asserted against a hand-counted fixture, not against production |
| production emits a REF bucket with `num_obs == 0`; ng emits no such observation | production creates `alleles[0]` at record open regardless ([open_record.rs:110-118](../../../../src/pileup/walker/open_record.rs#L110)); ng derives observations from reads that folded | the projection drops zero-support observations before comparing, and that drop is asserted to be the only one |
| observation **order** | production's is bucket-creation order; ng's is whatever the fold yields | ng **sorts** before emitting (below), so the comparison is against a sorted projection |
| **production's stale widen**: a read folded before the record widened keeps a **reference tail** production appended on its behalf, where ng's observation says what the read witnessed | see below — it is the one class that is *production being wrong*, not the two walkers describing the same evidence differently | `stale_widen_shape` ([parity.rs:2107](../../../../src/ng/locus_generation/pileup/parity.rs#L2107)): every production observation ng does not have must be some ng observation's bases plus a reference tail, with the per-`bases` evidence unmoved |

**The sixth class, and the mechanism, because the obvious version of it is wrong.** It is tempting
to say production's `widen` "re-folds nothing". It does not: `process_position` re-folds every
**contributor** whose events overlap an affected record, subtract-then-add, and production's own
comment argues that appending is equivalent to re-folding *"modulo the new bytes never being
event-modified by this read"*
([open_record.rs:405-410](../../../../src/pileup/walker/open_record.rs#L405)).

**The reads that "modulo" excludes are the class.** `widen` appends the new reference bases to
every allele's `seq` ([open_record.rs:416](../../../../src/pileup/walker/open_record.rs#L416)), and
records not affected at a step are deliberately not re-folded, because folding again would
double-count ([open_record.rs:658-661](../../../../src/pileup/walker/open_record.rs#L658)). So a
read that is **not a contributor at the widening step** keeps that appended tail unrevised. Two
kinds are: a read still live but silent there (an `N`, an adaptor mask, a ref skip), and a read
that has already expired.

ng answers each differently, and only the first is a function production lacks:
`refold_live_reads` ([open_record.rs:1441](../../../../src/ng/locus_generation/pileup/open_record.rs#L1441))
re-places the **live non-contributors** — it skips contributors precisely because the fold loop is
about to re-fold them — and an **expired** read cannot be reached by either walker, its cursor
having gone with it, so ng leaves its observation saying what it saw where production leaves the
tail.

**It is not class 1 and must not be filed there.** Class 1 is "a read did not witness the whole
footprint". Here the read may have witnessed every position of the final footprint and production
is wrong about it anyway. Filing it under class 1 would put reads production **mis-folded** into
the count of reads production **credited with bases they never sequenced** — and §13.2 asks for
those two as separate three-number measurements.

**Observation order has to be specified, not left to the fold.** ng's observations are derived from a per-read map
(`folded_reads` is an `AHashMap`, [open_record.rs:90](../../../../src/pileup/walker/open_record.rs#L90)),
whose iteration order is seeded per process — so emitting in fold order would make the output
**non-deterministic run to run**, contradicting §7 and §13's byte-identity check. Sort by
`(bases, read_witness, read_group)` before emitting. The STR generator already sorts for exactly
this reason ([ssr.rs:1055-1075](../../../../src/ng/locus_generation/ssr.rs#L1055)); this is the same
requirement on the same shared type.

*What stage 1 cannot cover, and it is the important half:* it runs on the **unchanged** copy, so it
proves the transcription and nothing about the two changes that reshape the output most — the
witnessed-extent rule and the group split. Those are stage 2's, and stage 2 has no byte-identity to
lean on. That asymmetry is why §12 requires stage 1 to be shown *discriminating* before it is retired,
and why a **permanent anchor** has to replace it — some differential that keeps holding after stage 1
is gone.

**What that anchor is, corrected: a fixture, not a filter.** This spec used to say it was the
narrower *complete-reads* differential — loci where every folded read witnessed the whole footprint
had to agree with production forever. **That is false, and class 6 above is why:** a read can
witness every position of the final footprint and production still be wrong about it, so filtering
on complete reads does not select an agreeing set. D1 found the counter-example on the anchor's
first run.

The anchor is therefore built on a fixture where the *cause* is absent rather than a filter over one
where it is not. `generate_uniform_events` gives every read on a contig **one shared event set**
(one CIGAR, one start; bases, qualities, MAPQ, strand and pairing still vary), so every widening
event is an event of every live read, every live read is a contributor at every widen, and none is
left stale. Records still widen — an earlier draft asserted they would not and was contradicted at
7 widens on its first case — so what the fixture removes is the *staleness*, not the widen.

The property to state, then, is **"on a fabrication-free fixture the two walkers agree at every
locus"**, and it is worth knowing that the complete-reads filter is retained on top of it and
currently excludes **nothing**: measured, 216,203 of 216,203 loci qualify, because uniform events
leave no way to be blind over part of a footprint. The filter is a tripwire for fixture drift, not
the thing doing the work — and a reader of the test should not have to guess which.

*One thing the harness must get right:* feed **the same** `PreparedRead` stream to both walkers.
ng's read preparation uppercases the reference where production does not
([read_preparation.md](read_preparation.md) §6), so preparing separately would make the two
streams differ on a soft-masked reference for reasons that are step 2's, not this step's. Parity
here is about the walk.

*And one thing that is safely identical:* the walker's reference contract is already canonical —
`MultiChromRefFetcher::fetch` returns "uppercase ASCII over `{A,C,G,T,N}`, canonicalised by the
fetcher implementation" ([fasta/mod.rs:117](../../../../src/fasta/mod.rs#L117)) — which is exactly
ng's `RefSeq` view. The step-2 case divergence does **not** recur here, and the adapter from
`RefSeq` to `MultiChromRefFetcher` is a shim with no semantic content.

---

## 4. The walk — what ng changes, and what it must not lose

The walk itself is production's and is not restated here: the code is the specification, and the
port is verbatim before it is anything else (§3). What follows is only the part a reader of that
code would not infer — the one behaviour ng changes, and the decisions inside the walk that a
careless port silently drops.

**Five decisions a rewrite loses without noticing**, each a wrong answer rather than a crash:
indels are dropped at the first or last CIGAR op; a ref-skip emits nothing and lets both flanks
emit independently; an indel's base quality is a `min` over a padded window, because deleted bases
carry no quality of their own; adaptor masking and `N` **silence** a base rather than flagging it;
and the column cap depends on the column's content (250 with an indel present, else 8000),
truncating in admission order. Every one is in the code being copied; they are listed so a reviewer
can check the copy kept them.

**The one invariant worth stating in prose, because it is the least obvious and was once a real
bug:** each (record, read) pair folds exactly **once over the record's lifetime**, not once per
position — a six-base footprint would otherwise count a spanning read six times. The mechanism is
subtract-then-add against `FoldedReadState`. It is the single thing most likely to be lost.

**ng's rule: nothing is ever written into an observation that its read did not witness.** That is one
rule, and it replaces both the fill and the widen-extension:

- **The haplotype builder does not fill.** `apply_events_to_ref_into` emits reference bytes for every
  offset no event covered ([:522-531](../../../../src/pileup/walker/open_record.rs#L522),
  [:568-573](../../../../src/pileup/walker/open_record.rs#L568)); ng's emits only what the events
  cover, and reports the extent it covered. **Reads whose events tile the footprint come out byte
  for byte as before** — there are no gaps to fill — so the complete class stays parity-comparable
  and only the fabricating cases move.
- **`widen` extends the REF bucket only.** `alleles[0]` is the record's own reference sequence and
  genuinely grows; the other buckets hold what reads witnessed and never grow. A live read re-folds
  against the wider window and lands wherever its bases put it; an expired read keeps a bucket whose
  bases already describe exactly what it saw. This is what makes the rule implementable: production
  cannot express it because the bases live on the shared bucket and `FoldedReadState` holds none of
  its own.
- **The witnessed extent is stored in absolute reference coordinates**, not relative to the
  then-current footprint, on `FoldedReadState` — which today holds only
  `{allele_index, contribution, chain_id}`
  ([open_record.rs:100-104](../../../../src/pileup/walker/open_record.rs#L100)). Coverage is then
  resolved **once, at `finalise()`**, by comparing that extent against the *final* footprint. A read
  `Complete` when it folded becomes `Observed` after a widen with nothing about the read having
  changed, and no re-fold is needed to notice — the read may be long gone, since `expire_passed`
  touches no open record.

**And ng's new "no observation" path sits outside that protection, which is where it will go wrong.**
The non-contiguous case (§6) is reached at *every* position the record is affected at, so a counter
incremented there multiplies by the footprint length — the same bug, on the one path with no
inherited test. It must be a per-record set of read ids, not a counter. Worse, a read can fold
contiguously and *become* non-contiguous when the window widens right across an interior gap: the
natural `continue` skips the subtract-then-add above it, leaving a live contribution in a bucket for
a read that now has no observation. That breaks `chain_ids.len() <= num_obs`, silently, and only on
multi-base records.

*Why absolute coordinates and not a bucket-time tag:* the bucket a read folds into is chosen from
its bases, which are fixed at fold time; its coverage is relative to a footprint that is not. Keeping
the two on different clocks is what lets bases stay immutable while coverage stays truthful.

**One base still comes from the reference, and it is worth naming rather than glossing.** An
insertion or deletion emits its **anchor** base from `ref_seq` when no `Match` already emitted it
([open_record.rs:546](../../../../src/pileup/walker/open_record.rs#L546),
[:556](../../../../src/pileup/walker/open_record.rs#L556)). Normally the `Match` is there and the
read's own base wins, so nothing is borrowed. The corner is a read whose base at the anchor was
dropped — `N` or adaptor-masked — while its indel at that position was still emitted: it witnessed
*the indel* but not *the anchor base's identity*, and the builder supplies the reference base for
that one position. **Recorded as a known residual, not fixed:** it is one base inside an event the
read genuinely witnessed, the alternative (discarding an observed indel over a masked anchor) loses
more than it saves, and the acceptance dump can count how often it fires before anyone pays to
change it.

**This converges the two generators, which is the outcome to protect.** The STR path already carries
witnessed-only bases — a read running off mid-tract yields `bases = b"CACA"`, not the full tract
([ssr.rs:886-894](../../../../src/ng/locus_generation/ssr.rs#L886); it tags that `PartialLeft(4)`
today, `Observed { offset_in_locus: 0, positions_covered: 4 }` after the reshape). Before this change the generic
path would have filled the same shared field with reference-padded bases: **one field, two
incompatible meanings.** After it, both mean "what the read witnessed", and that is now a
cross-generator invariant worth stating explicitly — `bases.len()` must be consistent with
`read_witness` on every observation, whichever generator minted it (§13).

**It also makes the step-7 scheme work as designed.** freebayes matches a partial's bases as a
prefix or suffix of each candidate and splits its weight `1/k` across the candidates it cannot
distinguish (§10). Fed a reference-padded partial, that machinery would match the REF candidate
exactly, find `k = 1`, and cast a **full-weight confident vote for REF** — hardening the defect
instead of softening it. Fed the witnessed bases, `AAA` from a 6-base locus is a prefix of both
`AAAAAA` and a deletion's `AAAA`, so `k = 2` and the read splits its weight, which is the whole
point of adopting the scheme.

---

## 5. Mate-overlap reconciliation — ported whole, and the half that gets lost

Read preparation is per read and pairwise-independent by design, so it pushed this downstream
([read_preparation.md](read_preparation.md) §10). It lands here, and it is ported **in full** — the
rules, the tie-breaks and the samtools constants are in the code and are not restated.

Three things that code does not say on its face:

- **Detection is by shared chain id**, not by read name at fold time — names are matched once, at
  admission. So a design that drops the chain-id allocator drops mate reconciliation with it: they
  are one mechanism.
- **Two regimes, and they differ in kind.** With an indel anchored at the position the pair
  *collapses* to a single observation; match-only, both survive and the quality is moved onto one of
  them. A port that treats them alike is wrong in one direction or the other.
- **The decision is taken in the walk and *replayed* in the fold.** The fold re-pulls each read's
  events from the cursor, which knows nothing about mates, so the outcome rides on the contributor
  as two flags. **Forget the replay and reconciliation silently applies at one position out of a
  record's whole footprint** — the failure is quiet, position-dependent, and invisible to any test
  whose records are one base wide.

One consequence for the group split (§6): a zeroed mate still counts as an observation with no
quality mass, in **its own** group's observation, and in the indel regime one group loses the read outright.
A per-group model therefore reads a slightly biased `(count, quality)` cross — the exact quantity the
split exists to serve. Recorded, not corrected: attributing a reconciled pair to one group would
invent an attribution the data does not have.

---

## 6. What a `PileupRecord` cannot say — the fields ng adds

Production's record maps onto the shared locus type almost field for field. The table is the port;
the paragraphs after it are where it is not.

| `SampleLocusObservations` | from production |
|---|---|
| `region` | `pos ..= pos + alleles[0].seq.len() - 1` — `GenomeRegion` is 1-based **inclusive** ([types.rs:79](../../../../src/ng/types.rs#L79)), so the end is the last covered position, not one past it |
| `reference_bases` | `alleles[0].seq` — the REF bucket's sequence *is* the record's reference bytes; there is no separate field ([open_record.rs:110-118](../../../../src/pileup/walker/open_record.rs#L110)) |
| `observations` | every `AlleleObservation`, REF bucket included |
| `kind` | `LocusKind::Generic` |
| `SequenceObservation.{bases, num_obs, num_fwd, q_sum, mapq_sum, mapq_sum_sq, chain_ids}` | `AlleleObservation.seq` / `AlleleSupportStats.{num_obs, fwd, q_sum, mapq_sum, mapq_sum_sq}` / `chain_ids` |

**`read_witness` — the one that changes an answer, and the one production gets wrong.** Production
folds a read into a record even when the read did not witness the whole footprint, filling what it
did not witness from the reference
([open_record.rs:522-531](../../../../src/pileup/walker/open_record.rs#L522),
[:568-573](../../../../src/pileup/walker/open_record.rs#L568)). At a six-base deletion locus a read
that saw only the first two bases is counted as a full witness of a six-base reference haplotype it
never saw, and `widen` extends the same fabrication retroactively to reads that have already left
the active set.

> **The fabrication primitive is "no event → reference base", and that is wider than "outside the
> read's span".** A read emits no event at a position that is `N`
> ([cigar_cursor.rs:272-274](../../../../src/pileup/walker/cigar_cursor.rs#L272)), adaptor-masked
> ([:278-280](../../../../src/pileup/walker/cigar_cursor.rs#L278)), inside a ref-skip, or covered by
> an indel the first/last-op rule dropped. Every one of those positions is *spanned* by the read's
> alignment and *not witnessed* by it. **A coverage tag derived from the alignment span is therefore
> blind to all four** — an adaptor-masked reverse-strand read would be tagged a complete witness of
> a haplotype it half-invented. So ng derives coverage from the **events**, which are the truth, and
> the fix lives in the haplotype builder rather than at the border (§4).

**It is not a design choice: no rationale exists anywhere, the walker spec asserts
the opposite invariant as already true, and freebayes — which this walker is derived from — has an
explicit span gate plus a down-weighted partial channel that a 2026-05-08 review scoped out of
comparison rather than read.** The mechanism, the archaeology, what GATK and bcftools do instead,
and why it should bias against calling deletions are in
[pileup_partial_coverage_ref_fill_2026-07-27.md](../../reports/research/pileup_partial_coverage_ref_fill_2026-07-27.md);
they are a finding about production, not rationale for this step, so they live there.

**What ng does.** The read's observation carries the bases it witnessed and an `Observed` coverage
tag computed from its **events** against the record footprint (§4) — a lower bound, kept as a
separate observation from the complete ones, with `complete_observations()` keeping it away from a
likelihood that would score it as a short allele. **This step records; it does not weigh.** Whether a partial is *used*, and how, is the
caller's — and that is now settled: step 7 adopts freebayes' partial-support scheme (owner,
2026-07-27; §10). The division of labour is the point of the split — a weighting needs the candidate
allele set, which does not exist until step 6.

**What step 7 will find here.** Not a list of what this step should add for it — that is step 7's
call, not this one's — but what its output already contains, so the question can be answered without
re-reading this code:

| freebayes uses | where it is in ng's output |
|---|---|
| the partial read's observed bases | `SequenceObservation.bases` |
| which end the read ran off (the prefix-vs-suffix test) | `read_witness`: `offset_in_locus == 0` is flush left (a prefix), `offset_in_locus + positions_covered == region.len()` is flush right (a suffix). An interior run is expressible too, which freebayes has no equivalent for |
| the count and quality to divide by `k` | `num_obs`, `q_sum` |
| the candidate allele set (to compute `k`) | step 6's output, not this step's |
| **read bases *outside* the locus window** | not on a single locus — see below |

**The last observation is reachable, and the REF-chain-id drop does not block it.** `assignPartialSupport`
extends a partial's sequence with the read's own bases beyond the window before testing
prefix/suffix (`Sample.cpp:377-381` → `Allele::read5p`/`read3p`, `Allele.cpp:992-1020`). Because
this generator emits **a locus per covered position** (§2), those bases are recorded at the
*neighbouring* loci, and step 7 reading a few loci either side gets them — **without needing the
read's identity in any REF observation**, because absence of a chain id *is* the information: a read
that appears in no non-reference observation at a locus was reference there, and a reference read's
bases are `reference_bases`, already on the locus. Identity is only needed where a read departs from
the reference, and that is exactly where it is kept (§6). Coverage does not need marking either: an
uncovered position emits **no locus at all** (§2), and a read that covered a locus and said nothing
is counted in `reads_without_observation`. So the outward walk is well founded — locus present plus
id absent means reference, and the read's own `read_witness` pins the border it ran off.

**The size of production's defect is this port's to measure, not to assume** — stage 2 of the parity
oracle counts it (§3, §13), and that number is what decides whether the research note's
indel-deficit hypothesis is a result or is dead.

**`ReadWitness` becomes `Complete` + one `Observed` variant — decided (owner, 2026-07-28).** Three
partial variants cannot describe what a read witnesses once the events, not the span, define it: a
read can be blind in the middle of a footprint (an interior `N`, a ref-skip) or blind at either end,
and a widened record can be wider than a read on both sides. One variant covers all of it:

```rust
pub enum ReadWitness {
    /// The read witnessed every position of the locus.
    Complete,
    /// The stretch it did witness, in **locus positions** — the axis `bases` is not on
    /// (that is allele content, in read coordinates). Derived from the read's *events*,
    /// never from its alignment span.
    Observed {
        /// Locus positions between the locus's left border and the first one witnessed.
        /// `0` = flush with the left border, i.e. a prefix constraint.
        offset_in_locus: u16,
        /// How many locus positions were witnessed, running from `offset_in_locus`.
        positions_covered: u16,
    },
}
```

**`Complete` is kept rather than folded into `Observed`.** It is the overwhelmingly common case, it
keeps `complete_observations()` a cheap equality instead of arithmetic against the footprint, and it
is exactly the STR path's "reached both borders". Prefix-versus-suffix survives as a derivation —
`offset_in_locus == 0` is flush left, `offset_in_locus + positions_covered == region.len()` is flush
right — so the STR path's "a prefix and a suffix are different constraints" is preserved, not lost.

**A non-contiguous witness yields no observation, and the read is counted.** `Observed` describes one
run, so a read blind in the middle cannot be summarised honestly and goes to
`reads_without_observation` instead. That is rare by construction — adaptor masking and the
dropped-indel rule always truncate from one side, so they stay expressible — and it has a useful
side effect: it gives that counter a real population (below).

**This is a change to a built ng module, and it is not small — a correction to an earlier estimate
in this spec.** `ReadWitness` is not `#[non_exhaustive]` and has **six** exhaustive match sites, not
one: `num_obs_along_locus` ([locus_generation/mod.rs:81-83](../../../../src/ng/locus_generation/mod.rs#L81)),
the STR generator's complete/partial tally ([ssr.rs:1015](../../../../src/ng/locus_generation/ssr.rs#L1015))
and its deterministic sort key ([ssr.rs:1072](../../../../src/ng/locus_generation/ssr.rs#L1072)),
plus four dump tools that `--all-targets` builds. The STR generator also **mints** coverage in four
places ([ssr.rs:713](../../../../src/ng/locus_generation/ssr.rs#L713), `:716`, `:889`, `:910`) — and
two of those pass `ReadWitness::PartialLeft`/`PartialRight` **as function values**, which a struct
variant cannot be, so that helper is restructured rather than retyped. `locus_generation.md` §3 and
`../arch/locus_generation.md` §1 carry the type and need the change folded in (§10).

**The read group — carried, and part of the observation's identity (owner, 2026-07-27).** The
read-group work merged on 2026-07-27 makes the group a first-class object: `AlignedRead` carries
`read_group: ReadGroupId` ([aligned_read.rs:36](../../../../src/ng/read/aligned_read.rs#L36)),
"the library it was prepared in", and the run's table maps that id to a `ReadGroup` holding the
sample, the library and the experiment ([read_groups.rs:43](../../../../src/ng/read/input/read_groups.rs#L43)).
**Today it does not survive to here, and the fix is a type ng owns (owner, 2026-07-28).**
Production's `PreparedRead` ([walker/mod.rs:236](../../../../src/pileup/walker/mod.rs#L236)) has no
such field — checked after the merge — so the group dies at the preparation boundary. A read's
library membership is a property of the read, exactly like its MAPQ or its strand; that the type
lacks it is an accident of one written before read groups existed. **So ng copies `PreparedRead`
into `src/ng/read/` and adds `read_group` there**, rather than editing production or having this
generator rebuild the group from a side channel.

**Copying the type is what makes copying the walker the right call**, not an extra cost on top of it
(§3): the four modules that looked reusable all name `PreparedRead` in their signatures, so an
ng-owned read type reaches them whatever else is decided. Two things follow, and both are
simplifications:

- **Production is not edited at all** — no field on a frozen type, no relocation of `ReadGroupId` to
  a crate-visible home, and no production→ng dependency inside the frozen caller (`src/pileup`,
  `src/psp`, `src/var_calling`, `src/vcf`, `src/pop_var_caller`), which has none today — the
  experimental `src/pop_var_caller_exp` bin already imports `crate::ng`, and is not what "frozen"
  means here. An earlier draft proposed exactly that edit; it is withdrawn.
- **The next field is free.** ng's read type will change again — BAQ if it returns, the re-align
  mode — and under the withdrawn plan each would have been another edit to frozen code.

`ReadPreparer` returns ng's `PreparedRead`, which reverses
[read_preparation.md](read_preparation.md) §3's "reuse production's as-is" — a fold-in that spec
owes (§10). *Unchanged by this:* `PreparedRead`'s home inside `pileup/walker/` is a recorded
misplacement (preparation produces it, the pileup only consumes it) deferred to the port-back. ng
copying the type neither pays that debt nor worsens it.

**Why carry it — and the strongest reason is not this path's.** `SequenceObservation` is **shared with
the STR generator**, and the STR model is *already* parameterised per sample group: the stutter
**level** is per group, the shape is per `(group, period)`, and the per-base error `ε` is a
per-group value ([ssr_cohort_mark2.md:288-289](../../specs/ssr_cohort_mark2.md)). Its glossary makes
the intent plain by listing `sample group` and `chemistry` side by side ([:529](../../specs/ssr_cohort_mark2.md)).
**And today those groups are inferred rather than known** — that spec calls them "data-driven soft
clusters", a proxy standing in for a grouping the data did not carry. The read group *is* that
grouping, declared by the sequencing run. This is also, verbatim, what motivated the read-group work
in the first place: *"emitting read-group columns from the stutter dump — the analysis that motivated
this work"* ([read_groups.md:34](../../ng/impl_plan/read_groups.md)). **Stutter depends heavily on
the library**, so the STR path is the near-term consumer, not a bystander.

The generic path's own case is real but further off: **damage is a property of the library
preparation** — ancient DNA names it, where C→T at read ends is a per-library rate — so a per-library
error model is something this caller may want. Either way the constraint is the same: it cannot be
added later from a tally that merged the libraries together.

**Why *split* the tally rather than attach counts.** Such a model needs the **allele × group cross
with its quality moments**: how many reads from library X supported this allele, at what qualities.
A per-group count beside one merged observation gives the first and loses the second. So the group
joins the dedup key: an observation is one distinct `(bases, read_witness, read_group)` combination.

**At the finest grain — per `@RG` (owner, 2026-07-28).** `ReadGroupId`
([types.rs:178](../../../../src/ng/types.rs#L178)) identifies one `@RG`, i.e. one lane; the run's
table carries the library and experiment each one belongs to. Carrying the *finest* grain is what
lets the consumer choose: library, experiment, or the read group itself, none of them foreclosed
here. The coarser ids do not exist anyway — the read-group work merged without a `LibraryId` or an
`ExperimentId` (checked at `eb2857c`) — but that is not the reason; the reason is that picking a
grain is a modelling decision and this step is not the modeller.

**And collapsing is exact, which is what makes the fine grain safe rather than merely
conservative.** Every field an aggregation touches is additive — `num_obs`, `num_fwd`, `mapq_sum`,
`mapq_sum_sq` sum; `q_sum` sums in log-error space; `chain_ids` union — and `bases` and
`read_witness` are identical across the observations being merged, since they are part of the key. So a
downstream fold from read group to library loses nothing a single-grain tally would have kept.

**What it costs.** Nothing where a sample has one read group, which is most of them — the observation count
is unchanged. Where it has several, the observations at a locus multiply by the groups covering it, and at
per-lane grain a four-lane library splits four ways. That cost is accepted as the price of not
deciding for the consumer. **If it ever bites**, the answer is to aggregate at the point of
consumption, or to mint a coarser id then — measuring it is cheap, and §13's dump is the
instrument.

**One interaction with mate overlap, small but exactly on target.** Mates are detected by shared
chain id, which is group-blind, so reconciliation itself is unaffected. But its *outcome* lands in
group observations: in the match-only regime the zeroed mate stays an observation of **its own**
group, with `num_obs = 1` and no quality mass; in the indel regime one group loses the read outright. So a
per-group model reads a slightly biased `(count, quality)` cross — which is precisely the quantity
the split exists to serve. Recorded, not corrected: the alternative is to attribute a reconciled
pair to one group, which invents an attribution the data does not have.

**Two consequences worth stating plainly.** First, `observations` is no longer a table of
distinct *sequences* — `read_witness` already made it a table of observations, and this widens the key
again. **A consumer that wants per-allele totals must aggregate over coverage *and* group**, and one
that treats each entry as an allele will now count the same allele several times. Second, the field
lands on the **shared** type, so the STR generator fills it too and its observations split the same
way. That rebaselines its fixtures — a real cost — but it is not a cost *imposed* on that path: it
is the path that wanted the split first, and today substitutes inferred clusters for it. Scheduling
belongs to whoever owns it (§10).

**`placed_left` — carried; `placed_start` — not (owner 2026-07-27, corrected 2026-07-28).**
`locus_generation.md` §11 assigns both read-position-bias counters to this generator: freebayes'
`placedLeft`/`placedStart`, how many supporting reads started strictly left of the record's anchor
and how many started exactly on it
([open_record.rs:790-791](../../../../src/pileup/walker/open_record.rs#L790)).

**`placed_left` is read**, and an earlier draft of this spec wrongly said nothing read it — the grep
behind that claim was truncated and never reached `src/vcf/`.
[vcf/qual_refine.rs:101,145,170](../../../../src/vcf/qual_refine.rs#L101) turns it into the
read-position-bias term subtracted from QUAL, live via `final_qual` into the cohort VCF writer and
the `--min-qual` gate.

**So `placed_left` is carried.** Dropping it would forfeit the ability to reproduce production's
QUAL, which outranks the tidiness argument entirely. `per_group_merger.rs` also rescales it through
the cohort merge, so it is not inert there either.

**`placed_start` is still not carried.** No model consumes it — it is merged, serialised and printed
by the `psp-to-pileup` dump, and read by nothing that computes. The YAGNI argument survives for it
alone, and it stays cheap to reverse: both counters are pure functions of the read's start against
the record anchor, and the read's span is already on `FoldedReadState` for the coverage rule (§4), so
a later consumer re-derives it at `finalise()` with no change to the fold.

**`reads_without_observation` — new, and the obvious definition is empty.** Production keeps no
per-record equivalent. "Reads considered, minus reads folded" is **identically zero**: contributors
are exactly the reads with an event at the walker position
([driver.rs:417-422](../../../../src/pileup/walker/driver.rs#L417)), every affected record's
footprint contains that position, so no contributor's window is ever empty and the fold's
`is_empty()` guard is unreachable. Considered *is* folded — worth knowing before anyone implements
the obvious thing.

Under the event-derived rule (§4) the counter means something: **a read whose witnessed
positions inside the footprint are non-contiguous** — an interior `N`, a ref-skip — cannot be
summarised by one `Observed` run, so it yields no observation and is counted here. Still counted
run-level rather than per record, because they are never contributors at all: reads whose bases over
the *whole* footprint are `N` or adaptor-masked, in
`PileupGeneratorCounts::reads_silent_over_footprint`. One consequence of that class worth knowing —
a position covered **only** by such reads opens no record, so it is indistinguishable from
uncovered.

**`reads_discarded_by_cap` — new, and "per record" needs defining because the cap is per
position.** Production counts `column_depth_truncations` on `RunSummary` — *positions* truncated,
run-wide. A read can be truncated at one position of a record's footprint and survive at another,
and if it folds at all it folds with its **whole** window, so its evidence is not subsampled. The
naive per-record count would therefore flag records whose support is complete.

The quantity the locus type actually wants — *"the support counts are a subsample, not the depth"*
(`locus_generation.md` §3) — is: **reads that had events inside this footprint and were truncated at
every position where they did, so folded nowhere.** Track the read ids truncated at any position of
an open record's footprint and count, at `finalise`, those absent from `folded_reads`. **The cap
truncates in the walk, before any record is identified** ([driver.rs:472-476](../../../../src/pileup/walker/driver.rs#L472)),
so those ids have to be plumbed into the fold and registered per affected record — a signature
change to the copied fold, not a local addition. Two cases have no clean answer and are recorded as
such: a read truncated where no record is open is unattributable to any locus, and a truncated read
carrying a deletion would have *widened* a record, so dropping it changes the footprint and hence
every other read's coverage. That is exact and bounded, and it is the same membership bookkeeping the
non-contiguous class needs (above), so it is paid once.

*Two consequences to record rather than discover.* The cap is applied before the fold, so a capped
read never influences `read_witness` — coverage describes the reads that folded, not the reads that
covered. And with the read group in the key, a truncated column yields a **group-biased** subsample,
since truncation is in admission order and says nothing about groups; a per-group model reading a
capped locus is reading a biased cross. Neither is a defect to fix here, but a model that ignores
`reads_discarded_by_cap` will be wrong in a direction it cannot see.

**`chain_ids` — carry, never for REF, and make it structural.** The memory lesson is already banked
in the code being ported: `finalise()` skips `allele_index == 0`, so REF chain ids are never
materialised ([open_record.rs:150-160](../../../../src/pileup/walker/open_record.rs#L150)); the
~96.6%-of-all-chain-ids / ~31%-of-peak-live-heap figure is what that fix removed, and stale `.psp`
files written before it still carry the cost. So the port inherits the fix — **provided it is not
lost in the flattening**. ng's `observations` is a flat list with no "index 0 is REF"
position, so the rule has to be restated in ng's own terms: **the observation whose `bases` equal
`reference_bases` carries no chain ids** — restated below, because observation splitting breaks that
wording. The alternative — dropping chain ids entirely, since ng's only consumer would be
a cohort merge that does not exist yet — was rejected because the allocator has to run regardless
(§5) and the marginal cost is one empty `Vec` per REF observation.

**The rule needs a per-read referent, not a per-observation one.** "The observation whose bases equal
`reference_bases`" named a unique observation while there was one observation per allele. It no longer does: observations
now split by coverage and by read group, so a reference-matching read can sit in a partial observation whose
bases are a *prefix* of `reference_bases` and never compare equal to it. Production's own rule is
positional (`allele_index == 0`, [open_record.rs:158](../../../../src/pileup/walker/open_record.rs#L158))
and equally unportable. **State it per read instead: a read contributes no chain id where it agreed
with the reference across everything it witnessed.** That is decidable at fold time from what the
read is, survives every split of the observations, and reduces to production's rule exactly when the observations
are one-per-allele.

**Why the REF drop is right, and not merely cheap.** Production justifies it by artifact size (~96.6%
of all chain ids, kept out of the `.psp`). The better reason is semantic: **a chain id marks which
haplotype a read came from, and the reference is the default — a default needs no tag.** Absence of
an id is therefore not missing information, it is the encoding: a read carrying no id at a locus came
from the reference haplotype there, and the reference haplotype's bases are `reference_bases`,
already on the locus. Only departures from the default have to be marked, and departures are the
rare case — which is why the saving is large. *Verified against freebayes, 2026-07-27:* it has no
chain-id analogue at all, because it never detaches an observation from its read (buckets hold
`Allele*` pointing into the read's own vector, `AlleleParser.cpp:3742`); its per-allele `readID` is
carried on reference and non-reference alleles alike but is **never used algorithmically** — all
five occurrences are debug or printing — and it does no cross-site phasing whatsoever (`GT` is
always sorted and `/`-joined, `Genotype.cpp:122`; no `PS`, no `HP`). So the drop forfeits **no
freebayes capability**.

**Nor does it conflate reference with absence at the locus level.** A position nobody covered emits
no locus at all (§2), so a locus's existence already says it was covered; reads that covered it and
said nothing are counted; and a read's own `read_witness` pins where it stopped. Chain ids
therefore carry one thing — *which non-reference haplotype* — and carry it only where there is
something to say.

**One limit on that, which bounds what step 7 may assume.** The encoding is per *locus*, not per
*read*: a read that was wholly silent over a locus — every base `N` or adaptor-masked — is a
contributor to nothing and appears in no per-locus field (it is counted only at run level, §6). So
"chain id absent at this locus" means *reference **or** silent*, and the two are not separable from
the locus alone. For step 7's outward walk that is a bounded error — it would read reference bases
for a position the read had masked — and it is bounded by how often a read is silent where a
neighbouring locus needs it, which the acceptance dump can count. Recorded rather than fixed:
separating them per read costs a per-locus membership test whose only consumer is hypothetical.

---

## 7. Config, counts, and cross-cutting concerns

**Config.** Per `locus_generation.md` §7 a generator owns its knobs and takes them at construction.
ng gets **its own constants**, starting at production's values but free to diverge — the same rule
the STR generator set for its reservoir cap. All five are production's, **inherited and never
measured by ng**; that is the map of what is safe to move.

```rust
pub struct PileupGeneratorConfig {
    /// Reads folded at a position with no indel anchored there. Production: 8000.
    pub max_snp_column_depth: u32,
    /// Reads folded at a position where any read has an indel. Production: 250.
    pub max_indel_column_depth: u32,
    /// Widest record footprint before the walk fails. Production: 5000.
    pub max_record_span: u32,
    /// How far a first mate stays available for pairing. Production: 10000.
    pub mate_lookup_window: u32,
    /// Active-read ceiling. Production: 4096.
    pub max_active_reads: u32,
}
```

**Counts.** Run-level, alongside the shared `LocusCounts`; production's `RunSummary` is the model
and most fields port straight across (`reads_admitted`, `record_widen_events`,
`mate_overlap_positions`, `chain_allocations`, `active_reads_high_water`, `mate_lookup_evictions`,
`column_depth_truncations`). Two are new: reads silent over a whole footprint (§6), and records
dropped by the region clamp (§2) — the latter should be *observably* zero-sum across neighbouring
regions, which is how the tiling argument stays checkable rather than asserted.

**Memory — everything this generator holds is bounded by depth, not by region length.** One locus
resident at a time is the shared contract; internally it holds the active read set (`max_active_reads`,
4096), the open-record table (`max_record_span` plus the walk's own locality), and **no read buffer
at all** — reads are pulled from a lazy stream and retired as the walk passes them (§2). That last
one is the property a port can quietly destroy: collecting the query into a `Vec` to make an
ownership problem go away would turn a depth-shaped footprint into a region-shaped one, and a
`Generic` region can run to hundreds of kilobases. The walker consumes an `Iterator`, and it must
keep consuming an `Iterator`.

One resident cost that is *not* the walker's and is worth knowing: on CRAM input the reader handle
now caches a decoded container (~10,000 records, measured at +34 MB peak RSS on `main`, `b918fb6`),
per worker. It is a constant, not a leak, but it is charged to this step's process.

**Throughput, and why it is worth measuring first.** A locus per covered base (§2) is the design,
and there is no candidate-only fallback to retreat to — the two-stage split requires it. So a bad
number here is a performance problem to solve, not a design to reconsider, and it is cheaper to
learn that early: time the acceptance-test dump (§13) over one human chromosome against
production's `pileup` subcommand on the same CRAM, as soon as the dump exists. Production is the
right yardstick because it does the same walk and writes an artifact per position too.

**The allele list gets longer, and that is the cost centre to watch.** Two of §4's changes multiply
buckets. Under the no-fill rule, reads witnessing different extents of a multi-base record no longer
share one REF bucket — each distinct extent is its own observation, so the baseline count tracks distinct
read starts rather than distinct alleles. And REF-only widening removes the property that makes a
live read re-fold into its *existing* bucket: production appends the same bytes to every bucket
precisely so the re-fold matches, a 25-line comment above the loop
([open_record.rs:390-415](../../../../src/pileup/walker/open_record.rs#L390)) exists to prove it, and
without it every widen leaves the old bucket behind at `num_obs == 0` and creates a new one. That
matters because `find_allele_index` is a **linear scan with a full byte compare**, run once per
(record, contributor) at every position of the footprint — and the comment above it, "records
typically carry ≤ a few alleles", stops being true exactly at the long-deletion loci the port exists
to fix. **Design for it rather than discover it:** evict `num_obs == 0` buckets at widen, and expect
to key buckets by a hash of the bytes rather than a scan. Note the counter-pressure too — a
positional `allele_index` on `FoldedReadState` is what makes eviction awkward.

**Determinism.** No sampling anywhere: the depth cap truncates in admission order, and the
mate-overlap map is `ahash` specifically so iteration is stable run to run
([driver.rs:674](../../../../src/pileup/walker/driver.rs#L674)). Output is a deterministic function
of (reference, config, reads).

**Errors.** The walker's `WalkerError` variants are all fatal and terminal for the iterator
([driver.rs:213](../../../../src/pileup/walker/driver.rs#L213)); they map onto
`LocusGenerationError` as a new variant, since none of the three existing ones
(`TypedRegion`/`Reads`/`Reference`) covers a malformed read or an exhausted chain-id space. That
is a widening of a `#[non_exhaustive]` enum in `locus_generation/mod.rs` — ng's own code, not
production's.

---

## 8. Traps — the ones that cost a day

- **Hoisting the chain-id allocator corrupts two counters, because `reset()` preserves them and
  `summary()` *assigns* them.** `reset()` deliberately keeps `counters` — they are file-scoped
  ([chain_id_allocator.rs:179-185](../../../../src/pileup/walker/chain_id_allocator.rs#L179)) — and
  `RunSummary` takes them by assignment, not addition
  ([driver.rs:249-254](../../../../src/pileup/walker/driver.rs#L249)). Production is safe because it
  builds a fresh walker, hence a fresh allocator, per region. ng shares one across regions (above),
  so every region's summary reports the **run-to-date** total, and adding them up gives a triangular
  sum — `chain_allocations` and `mate_lookup_evictions` inflated by roughly the region count, in a
  plausible-looking `u64`. Snapshot the counters at `begin_segment` and fold the delta, or read the
  allocator once at run end. (`active_reads_high_water` is a max and survives, which makes the
  corruption look selective enough to rationalise.)
- **`events_overlapping` does not clip a deletion to the window.** Matches are clipped; a `Deletion`
  comes back whole whenever its footprint intersects, so its span can run past the record end, and
  one anchored before the record can report an anchor below `record_pos`. The witnessed extent must
  be intersected with `[record_pos, record_end)` or `offset_in_locus`/`positions_covered` are wrong
  at exactly the deletion loci this port exists to fix.
- **`bases.len()` is not `positions_covered`.** An insertion's footprint is 1 reference position but
  contributes several bases; a deletion's is `deleted_len + 1` and contributes one
  ([decompose.rs:55-61](../../../../src/pileup/walker/decompose.rs#L55)). The §13 consistency check
  is *"no observation claims more positions than its events account for"*, **not** an equality — an
  implementer who makes it one will either fail it on every indel observation or "fix" it by deriving the
  extent from the byte length, which is the span-vs-events confusion all over again.
- **Build the reference accessor once, not per segment.** It is a field on the generator;
  `begin_segment` must not replace it. A fresh accessor per region throws away the sliding buffer at
  every boundary *and* re-pays `RawChromReader::for_contig` — a full `.fai` parse plus two
  `open(2)`s. That is not hypothetical: it is exactly the bug `4bc3ef9` fixed on the STR side, where
  a per-query accessor cost ~564k opens on one tomato chromosome and 14% of a cohort run. The same
  applies to the preparer's own accessor.
- **Share the allocator's *counter* across regions, and reset the rest at every region end.** ng
  walks one region at a time where production walks a chromosome. A fresh allocator per segment gives
  two different fragments the same id, and a later phasing step chains them — so the allocator lives
  on the generator. But it is **not** just a counter: it also holds `pending_mates`, `active_count`
  and the eviction window
  ([chain_id_allocator.rs:78-112](../../../../src/pileup/walker/chain_id_allocator.rs#L78)), and
  production scopes all three per chromosome. Carried blindly across regions they misbehave in two
  ways — a pending first mate from one contig can pair with a read on another (the eviction test
  compares raw positions), and `active_count` never returns to zero if a walk is abandoned before it
  drains, ending the run in `ActiveReadsExhausted`. The fix is production's own idiom: call `reset()`
  at the end of each region walk, which clears `pending_mates` and `active_count` **and preserves
  `next_id`** ([:179](../../../../src/pileup/walker/chain_id_allocator.rs#L179)). *This is a change
  in the walk, not the fold* — worth saying, because the port's change set is otherwise concentrated
  in `open_record.rs`.
- **A read overlapping two regions is admitted in both walks and gets two chain ids.** One fragment,
  two identities — the mirror of the mates-in-different-regions limitation below, and it bites
  single-end reads too. Unavoidable while regions are walked independently; recorded so a phasing
  step does not assume the id is a fragment key across the genome.
- **Mate pairs are chained only within a region.** A consequence of the one-segment contract that
  no amount of care removes: if two mates land in different `Generic` regions, each walk sees only
  one of them and they get different chain ids. *Reconciliation* (§5) is unaffected — physically
  overlapping mates are always in the same walk — but compound haplotypes break at region
  boundaries. Whether that matters is the phasing step's question (§10).
- **Left-alignment does not move a read's start**, so preparation cannot break the walker's
  coordinate-order requirement. Production calls its normalizer with end-deletion stripping off
  precisely so `reference_offset` never moves ([read_preparation.md](read_preparation.md) §5). Worth
  knowing because the walker's order check is fatal (`WalkerError::OutOfOrder`), and a preparer
  that later gains the re-align mode *would* move starts.
- **`PreparedRead::length()` is a precondition, not a validation.** The cursor indexes `seq` and
  `bq_baq` by CIGAR-derived offsets, so a CIGAR whose read-consuming length disagrees with `seq` is
  checked once at admission ([driver.rs:396](../../../../src/pileup/walker/driver.rs#L396)) and
  assumed thereafter.
- **The walker skips uncovered gaps but walks covered ones base by base.** `advance` jumps to the
  next read's start **only when the active set is empty** ([driver.rs:549](../../../../src/pileup/walker/driver.rs#L549)).
  A rewrite that jumps more eagerly loses every REF-only locus.
- **`events_at` and `events_overlapping` are deliberately not the same query.** The first returns
  events *anchored* at a position; the second returns events whose *footprint intersects* a window,
  which includes a deletion anchored before the window whose run reaches into it. Substituting one
  for the other compiles and is wrong.

---

## 9. Reuse over rewrite — the map

Every row read at the cited line, 2026-07-26.

| what | existing code | ng reuse |
|---|---|---|
| CIGAR cursor, decompose, active read set, chain-id allocator | [cigar_cursor.rs](../../../../src/pileup/walker/cigar_cursor.rs), [decompose.rs](../../../../src/pileup/walker/decompose.rs), [active_read_set.rs](../../../../src/pileup/walker/active_read_set.rs), [chain_id_allocator.rs](../../../../src/pileup/walker/chain_id_allocator.rs) | **copy** — each names `PreparedRead` in its signatures, so ng's own read type reaches all four (§3) |
| the prepared read | `PreparedRead` / `MateRole` / `ReadLengthError` [walker/mod.rs:236](../../../../src/pileup/walker/mod.rs#L236) | **copy into `src/ng/read/`** and extend with `read_group` (§6); reverses [read_preparation.md](read_preparation.md) §3's reuse-as-is |
| the walk loop and the open-record fold | [driver.rs](../../../../src/pileup/walker/driver.rs) (→ ng's `genome_walk.rs`), [open_record.rs](../../../../src/pileup/walker/open_record.rs) | **copy verbatim** into `src/ng/locus_generation/pileup/`, prove byte-identical, then change (§3) |
| `PreparedRead`, `CigarOp`, `MateRole`, `WalkerConfig`, `WalkerError`, `RunSummary` | [walker/mod.rs:236,43,188,133](../../../../src/pileup/walker/mod.rs#L236) | already `pub` — reuse as-is, no copy |
| the read input | `SampleReads::reads_in_region` ([read/input/mod.rs:508](../../../../src/ng/read/input/mod.rs#L508)) | reuse as-is; one query per region (§2) |
| `AlignedRead` → `PreparedRead` | `LeftAlignPreparer` ([read/left_align.rs:87](../../../../src/ng/read/left_align.rs#L87)) | call — this generator is step 2's only consumer |
| the read group | `AlignedRead.read_group` ([aligned_read.rs:36](../../../../src/ng/read/aligned_read.rs#L36)), `ReadGroupId` ([types.rs:178](../../../../src/ng/types.rs#L178)), the run's table ([read_groups.rs:43](../../../../src/ng/read/input/read_groups.rs#L43)) | **carry** into the observation's identity (§6). `PreparedRead` does not hold it, so the generator keeps it beside the walk rather than through it |
| the reference | `RefSeq` (canonical) → `MultiChromRefFetcher` | a shim; the two contracts already agree on canonical bytes (§3) |
| the region clamp | [pileup_to_psp.rs:271](../../../../src/pileup/per_sample/pileup_to_psp.rs#L271) | reuse the **rule**, not the code — the rest of that file is `.psp` machinery |
| the parity harness shape | [alignment/delimit_parity.rs](../../../../src/ng/alignment/delimit_parity.rs), [read/left_align_parity.rs](../../../../src/ng/read/left_align_parity.rs) | model for §3's differential |
| walker test fixtures | [walker/tests.rs:32-152](../../../../src/pileup/walker/tests.rs#L25) | `pub(crate)` under `#[cfg(test)]` — usable from ng's tests today |
| the windowed-statistics look-ahead buffer | `SampleSummaryAccumulators` ([pileup_to_psp.rs:57](../../../../src/pileup/per_sample/pileup_to_psp.rs#L57)) | **do not port** — it exists to fill two `.psp`-only fields the ng type does not have |

**Parity oracle:** production's `PileupWalker` on the same `PreparedRead` stream, stage 1 exact and
stage 2 through a projection with two named divergence classes (§3).

---

## 10. Deferred, with a recommended home

- **A second generic generator** — the active-region or haplotype-window definition of a locus
  ([ng_proposal.md](ng_proposal.md) §1 names it as a swappable axis). **Home:**
  `src/ng/locus_generation/`, beside this one; the trait exists so they sit side by side.
- **Consuming partial observations — the behaviour is decided, the work is step 7's (owner,
  2026-07-27).** ng adopts freebayes' partial-support scheme: match the bases the read actually has
  as a prefix or suffix of each candidate allele, count how many candidates it cannot tell apart
  (`k`), and give it `1/k` of an observation and `1/k` of its quality to each — **fractional depth
  included**, which is the part easiest to drop by accident (`Sample.cpp:354-420`, `:41`, `:86`,
  `DataLikelihood.cpp:76-95` — MIT, portable; worked example in the research note §4). It lands at
  step 7 because `k` cannot be computed without the candidate allele set, which is step 6's output.
  **Nothing changes in this step because of it** — §6 maps what step 7 will find in this step's
  output, including the one item (read bases outside the locus window) that lives at the neighbouring
  loci rather than on one locus. **Home: step 7**, together with the STR path's censored-observation
  question ([locus_generation_ssr.md](locus_generation_ssr.md) §8): the two are the same question on
  two paths, and `complete_observations()` is the guard on both until they are answered.
- **A fold-in to [read_preparation.md](read_preparation.md) §3** — it prints `PreparedRead`'s fields
  and records the decision to reuse production's type as-is. ng now owns its own copy with a
  `read_group` field (§6), so that spec's owner should record the reversal and that the preparer
  threads the group through from `AlignedRead`. **Home: that spec.** *(An earlier draft made this a
  prerequisite pass that edited production's `PreparedRead` in place; withdrawn — ng owns the type,
  so production is untouched.)*
- **Four changes to the shared locus type — schedule them as one pass, because they share one
  fixture rebaseline.** All four land on types both generators fill, so the STR generator's output
  moves with them; done separately they rebaseline its fixtures four times.
  1. `ReadWitness` → `Complete` + `Observed { offset_in_locus, positions_covered }` (§6). Six
     exhaustive match sites and four minting sites, two of which pass the variant as a function
     value.
  2. `SequenceObservation` gains the **read group**, which splits its observations (§6). The STR path is the
     one with a model waiting for it, so this is a hand-off rather than a fold-in — whether it then
     replaces its inferred sample groups with declared ones is the STR cohort work's call.
  3. `SequenceObservation` gains **`placed_left`** (§6). The type carries neither bias field today, so
     `locus_generation.md` §11's instruction to add both is half right: add `placed_left`, not
     `placed_start`. Left as-is, the value is computed, carried through the fold, and has nowhere to
     go at the last line of `finalise`.
  4. `locus_generation.md` §3's `reads_without_observation` wording — "reads that covered this locus
     and produced no observation at all" — is broader than the generic path fills, since the
     wholly-silent read is counted at run level (§6). Either the doc gains the caveat or the field is
     read as a lower bound.

  **Home: `locus_generation.md` §3 and `../arch/locus_generation.md` §1, plus the STR specs for (1)
  and (2).** Sequencing them before this generator's fixtures exist keeps the rebaseline to one.
- **An aggregating accessor on the shared type** — with `read_witness` and now the read group in
  the key, `observations` splits by both, and per-allele totals need a fold over both axes
  (§6). `complete_observations()` is the precedent for putting that guard on the type rather than in
  every consumer. **Home: `locus_generation.md`**, when the first consumer needs it; noted so the
  need is not rediscovered as a bug.
- **The cohort merge, the windowed depth/GC statistics, the `SsrBundle` generator, parallelism** —
  all already deferred with homes in [locus_generation.md](locus_generation.md) §11; nothing here
  changes them. This generator is the one that makes the windowed statistics *possible*, since
  per-position depth over a generic stretch is what they slide over.
- **Cross-region phasing** — chain ids do not survive a region boundary (§8). Whether compound
  haplotypes need to, and at what cost, is the phasing step's question. **Home: step 10.**
- **BAQ** — deferred sine die by [read_preparation.md](read_preparation.md) §10. Recorded here only
  because the walker is where its effect would be felt, through `bq_baq`.
- **Moving `PreparedRead` and `CigarOp` out of `pileup/walker/`** — a pre-existing misplacement,
  already on record twice. **Home: the port-back**, when production unfreezes.

---

## 11. Resolved decisions & open questions

**Resolved.**

- **Copy the whole walker into ng; edit production zero times (owner, 2026-07-28).** ng needs its own
  read type (§6), and every module that looked reusable names `PreparedRead`, so it reaches all of
  them. *Rejected:* driving production's walker as a black box — it cannot fill three fields of the
  locus type and bakes in the fabricated-REF fold. §3.
- **A locus per covered position, REF-only included (owner, 2026-07-26)** — because the per-sample
  and cohort stages are split, so a sample's evidence is gathered before anyone knows which
  positions carry candidates. Emitting only this sample's non-reference positions would leave no
  data at a position where *another* sample turns up a candidate, and that SNP could not be called.
  The candidate-only alternative is therefore ruled out, not deferred. §2.
- **This step records partial coverage; step 7 weighs it, using freebayes' `1/k` scheme (owner,
  2026-07-27).** The split is forced — `k` needs the candidate allele set, which does not exist
  until step 6 — and it leaves this step's output unchanged. §6 maps what step 7 will find here.
- **Nothing is ever written into an observation that its read did not witness (owner, 2026-07-27/28).**
  One rule replacing two mechanisms: the haplotype builder does not fill, and `widen` extends the REF
  bucket only. Three ports follow — the witnessed extent comes from the **events**, not the alignment
  span (the span is blind to `N`, adaptor-masked, ref-skipped and dropped-indel positions); it is
  stored in absolute reference coordinates on `FoldedReadState`; and coverage is resolved once at
  `finalise()` against the final footprint. §4, §6.
- **`ReadWitness` becomes `Complete` + `Observed { offset_in_locus, positions_covered }` (owner,
  2026-07-28).** One run in **locus** coordinates describes every case the events can produce —
  blind at either end, blind in the middle, or a record wider than the read on both sides. `Complete`
  is kept: it is the common case and keeps `complete_observations()` an equality test. A
  non-contiguous witness yields no observation and is counted instead. The blast radius is six exhaustive match
  sites and four minting sites, two of which pass the variant as a function value. §6, §10.
- **Chain ids: carried, never for the REF observation, as a type invariant.** Production already
  drops them; ng restates the rule in terms of its flat observation list so the flattening cannot
  lose it. Dropping them wholesale was rejected — the allocator runs regardless, since it *is* the
  mate-overlap predicate. **The rationale recorded is the semantic one, not production's
  artifact-size one (owner, 2026-07-27):** a chain id marks a read's haplotype, the reference is the
  default, and a default needs no tag — so absence *is* the encoding, and only departures are
  marked. Checked against freebayes: no chain-id analogue, `readID` never used algorithmically, no
  cross-site phasing — nothing is forfeited. §5, §6.
- **`reads_without_observation` is the cheap per-record count; the wholly-silent read is tallied at
  run level (owner, 2026-07-28).** The per-open-record counter is "reads considered for this record,
  minus reads that folded into it". It misses one class by construction: a read silent over the
  *whole* footprint — fully adaptor-masked, or all `N` — is never a contributor at any position, so
  it is never considered. Counting it per locus would need a membership test against reads that
  overlap the footprint but have already left the active set; that cost is not worth paying for a
  number nothing yet reads. The class is not lost — it goes to
  `PileupGeneratorCounts::reads_silent_over_footprint`. **The consequence is a contract statement,
  not a defect:** the per-locus value is an **honest lower bound**, and `locus_generation.md` §3
  describes the shared field more broadly ("reads that covered this locus and produced no
  observation at all"), so either that doc gains the caveat or the generic path's value is read as
  the subset it is — a fold-in (§10). Revisit only if a filter turns out to read the per-locus
  number.
- **`LocusGenerationError` gains a `Walker` variant.** None of `TypedRegion` / `Reads` / `Reference`
  describes a malformed read or an exhausted chain-id space (§7); ng's enum is `#[non_exhaustive]`,
  so the addition is source-compatible. No fork — recorded because it changes a shared ng type.
- **The read group is carried at `@RG` grain (owner, 2026-07-28).** The finest grain available, so
  the consumer picks library, experiment or read group rather than inheriting this step's guess —
  and the aggregation is exact, since every support field is additive and the merged observations share
  their `bases` and `read_witness` by construction. The cost is an observation split per lane where a sample
  has several, accepted deliberately. §6.
- **The read group is carried and joins the observation's identity (owner, 2026-07-27/28).** Two
  reasons, the stronger one the STR path's: its stutter and `ε` are already fit **per sample group**,
  and those groups are currently *inferred* ("data-driven soft clusters") because the evidence did
  not carry the real one — the read group is it. The generic path's own case, a per-library damage
  model, is real but further off. Either needs the allele × group cross *with its quality moments*,
  which only a split tally expresses. Free at one read group. *Rejected:* per-group counts beside one
  merged observation, which keeps the counts and loses exactly the moments the models need. §6.
- **`placed_left` carried; `placed_start` not.** `vcf/qual_refine.rs` turns `placed_left` into the
  read-position-bias penalty subtracted from QUAL, so dropping it would forfeit QUAL parity.
  `placed_start` is consumed by no model and is derivable later from the span §4 already keeps. §6.
- **Mate-overlap reconciliation is ported as a real algorithm**, both halves — the decision in the
  driver and the replay inside the fold. §5.
- **Region boundaries are handled by clamping on the record's anchor**, production's rule; gap-free
  disjoint tiling makes it exact. §2.
- **The chain-id allocator lives on the generator, across segments** — a per-segment allocator
  collides ids between regions. §8.

**Open (confirm before code).**


**Not to fix here.** Two standard validation commands are red independently of this work, and every
step of every ng plan currently has to except them by hand: `cargo test --all-targets
--all-features` panics in [benches/psp_writer_perf.rs:386](../../../../benches/psp_writer_perf.rs),
and `cargo doc --no-deps` fails on 11 unresolved intra-doc links. Both are tracked under
PROJECT_STATUS "Standing project-wide items".

---

## 12. What happens to production's tests

The walker.s **113** tests come with the copy — 44 end-to-end in `tests.rs` plus 69 inline across the
seven files (70 `#[test]` markers, one pair mutually exclusive by `cfg(debug_assertions)`). *(This
spec and the port plan both said "46" until 2026-07-29, when A4 counted them.)* Knowing in advance which must still pass and which must
change is how the port is kept honest — a test that changes silently is a behaviour change nobody
decided.

**At stage 1 (the verbatim copy) all 113 must pass, unmodified.** That *is* the gate. Anything that
needs touching at this stage is a transcription error, not a design change.

**At stage 2 they split three ways.**

*Must still pass, after mechanical adaptation only* (the read type moves to ng, `PileupRecord`
becomes `SampleLocusObservations`) — this is the regression floor, and it is most of the suite:
CIGAR decomposition and the cursor's parity against the `decompose` oracle; the `PreparedRead`
length checks; admission and ordering errors; adaptor masking; the column-depth cap in all four of
its cases; every mate-overlap and chain-id-pairing test; record coordinate ordering; uncovered
positions emitting nothing; insertion and deletion record shape; the widen-events counter.

*Must change, because the behaviour deliberately changed:*

| test | why |
|---|---|
| `placed_left_and_placed_start_are_per_record` | `placed_start` is no longer carried (§6); the `placed_left` half stands |
| `refold_after_widen_clears_chain_id_from_old_bucket` | `widen` now extends the REF bucket only (§4), so what a re-fold moves between buckets changes |
| `chain_ids_persist_across_chromosome_boundaries` | ng's boundary is the **region**, and the allocator is `reset()` there (§8) — the property is the same, the fixture is not |
| any test asserting allele order within a record | ng sorts observations before emitting (§3); production's order is bucket-creation order |

*Must be **added**, because nothing existing covers them* — and this is the finding worth carrying:
**no test in the suite exercises the defect the port exists to fix.** The two that come closest do
not: `deletion_record_does_not_double_count_ref_reads`'s spanning read covers the whole footprint, so
it is `Complete` either way, and `g1_walker_drops_match_observations_past_adaptor_boundary` uses
one-base records, where a read either witnesses the position or is not a contributor. **The
fabrication is untested in production, which is a large part of why it survived** — the same
"test that cannot fail" pattern this project has hit repeatedly. So the port cannot rely on the
inherited suite to catch a regression here, and owes new fixtures:

1. A multi-base record with a read **adaptor-masked over part of it** — must be `Observed`, not
   `Complete`, with its bases the length it witnessed. A span-derived implementation passes
   everything else and fails this.
2. The same with an interior `N`, and with a ref-skip — the non-contiguous case, which must yield no
   observation and be counted.
3. A record widened past an **already-expired** read — its bases must not have grown.
4. A read witnessing an indel whose own anchor base was masked — the one residual where a reference
   base is still borrowed (§4); the test pins that it is one base and no more.
5. Two read groups supporting one allele — two observations, summing to the single-group total.
6. A deletion at a region boundary whose footprint reaches into the next region — the halo case
   (§2); its support must match a single-region walk over the same span.

## 13. Acceptance test

This generator emits no variant calls, so "done" must not decay into "compiles". Three things,
in the order they should be built.

1. **The port is identical.** §3's stage-1 differential: production's walker and ng's copy over one
   `PreparedRead` stream produce equal `PileupRecord` streams and equal `RunSummary`s, on the
   walker's own fixtures and at scale on GIAB HG002 and a tomato CRAM. Zero divergences, or the
   port is not done. *The lesson from the three previous ports applies: a differential that passes
   immediately is more likely to have an inadequate generator than a correct port — the fixture
   must be shown to exercise mate overlap, adaptor masking, widening, re-folds, and the column cap,
   by mutating each and watching the differential fail.*
2. **The ng-shaped output is a projection.** §3's stage-2 assertion, with every divergent locus
   falling into one of the **six** named classes and counted, not excused. The count of
   partial-coverage divergences is **the deliverable, not a by-product**: how many loci, how many
   reads, and how many reference bases production credits to reads that never sequenced them — and,
   **separately, the same three numbers for the stale widen** (§3, class 6). That is the measurement
   that turns the indel-deficit hypothesis into a result or kills it.

   **Which reads the second triple is about — corrected, because the original wording named the
   wrong population.** It used to say "the reads `widen` extended after they had already left the
   active set". Those reads leave ng holding an `Observed` row, so the census files them under
   **class 1**, whose triple already counts them: class 6 is gated on there being *no* partial
   witness ([parity.rs:1999-2000](../../../../src/ng/locus_generation/pileup/parity.rs#L1999)).
   The population the second triple is owed for is the other kind of non-contributor — a read
   **still live** at the widening step but silent there — where ng re-places the read and
   production keeps the reference tail. **Both numbers are still owed:** the census currently
   carries class 6 as a locus count only, with no reads and no reference bases
   ([parity.rs:1699](../../../../src/ng/locus_generation/pileup/parity.rs#L1699)), so at the
   default case count it can say production mis-folds reads at 264 loci and cannot say how many
   reads or how many bases moved.
3. **A dump tool over a committed fixture** — `examples/ng_generic_loci_dump.rs`, following
   `examples/ng_ssr_loci_dump.rs`: a `#`-prefixed `key=value` counts header, then a TSV. Asserted,
   so it is a regression test and not a demo:
   - Every position covered by at least one observing read yields exactly one locus, and
     uncovered positions yield none.
   - Every fetched read is accounted for: it supports an observation, is counted in
     `reads_without_observation`, was discarded by the cap, or was silent over the whole footprint
     (§6).
   - **No observation carries a chain id for a read that agreed with the reference across
     everything it witnessed** — §6's invariant, stated per read because observation splitting means "the
     REF observation" is no longer a unique observation. Checked rather than commented.
   - **On a fixture with two read groups, an allele supported by both is two observations**, one per group,
     whose `num_obs` sum to the single-group total — the check that §6's split is computed and not
     defaulted, and the measurement the grain question needs (§11). On a one-group fixture the observation
     count must be **identical** to a run with the field ignored, which is what "free at one read
     group" has to mean in practice.
   - observations with an `Observed` witness exist and are separate from `Complete` ones with the same bases — which is
     what proves `read_witness` is computed rather than defaulted.
   - **A read blind in the middle of a footprint yields no observation and is counted**, and a read
     adaptor-masked over part of a footprint is `Observed`, not `Complete` — the two checks that
     coverage comes from the **events** and not the alignment span (§4). A span-derived
     implementation passes every other test in this list and fails these two, which is the point of
     them.
   - **A locus whose record widened past a read's end shows that read as `Observed`**, with its
     bases the length it witnessed and not the footprint's — the check that the no-fabrication rule
     survived the port, and the one production cannot pass by construction. Needs a fixture with a
     deletion long enough to widen a record past a read that has already expired.
   - **No observation claims a position the locus does not have**, and none claims zero — every
     observation with an `Observed` witness satisfies
     `offset_in_locus + positions_covered ≤ footprint` and `positions_covered > 0`, asserted
     globally, because a witness that reaches past the footprint it is measured against, or that
     witnessed nothing and still has an observation, is the shape §4's change exists to make
     impossible.

     > **Corrected 2026-07-30.** This clause used to ask for *"the consistency check between
     > `bases` and `read_witness`"*, and the dump tool's doc restated it verbatim above an
     > assertion that never reads `bases`. **The check as worded cannot be written.**
     > `positions_covered` is *derived from* the read's events, so checking it against them is
     > tautological; and §8's own trap is that no inequality relates `bases.len()` to the footprint
     > in general — an insertion adds bases without positions, a deletion positions without bases,
     > so an observation's byte count may exceed or fall short of its position count either way.
     > What is checkable is the bound on the *positions*, which is what the tool has always
     > asserted. The weakening is in the wording, not in the code: nothing was checked before that
     > is unchecked now.

     The evidence that `bases` and `read_witness` agree is carried elsewhere and by construction:
     both come from the same `apply_events_into` call, which returns the bases and the witnessed
     extent together and returns `None` rather than a run it cannot describe honestly (§6).
   - Loci from two adjacent `Generic` regions concatenate into a coordinate-sorted,
     duplicate-free stream with no gap at the join (§2).
   - Output is byte-identical across repeated runs.
