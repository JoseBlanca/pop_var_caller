# ng run driver, psp mode — implementation plan

**Status:** draft, 2026-09-03. The build order for **psp mode**: the walk stage that writes
one psp (and census) per sample, and the calling stage that reads a cohort of psps to a VCF —
the two run objects `run_streaming.md` names `SampleObservationGatherer` and
`PspVariantCaller`, and the subcommands `generate-psps` and `call-from-psps`. Design is
settled in [`run_streaming.md`](../spec/run_streaming.md) (spec) and
[`../arch/run_streaming.md`](../arch/run_streaming.md) (types & interfaces); the file format
in [`psp_file_format.md`](../spec/psp_file_format.md) and its executed plan
([`psp_file_format.md`](psp_file_format.md), milestones A–H ✅). This plan turns that design
into build order; it is **not** a place for new design — three questions it surfaced are
routed upstream below, and none blocks step 1.

This is the plan [`run_driver_direct_mode.md`](run_driver_direct_mode.md) deferred
("Everything psp … their own plan, once the encoding settles" — its Scope/Out). The encoding
settled 2026-08-30; direct mode completed 2026-09-01 and is this plan's oracle.

## ⚑ Another plan runs inside this one, and it has to land before Milestone F

[`psp_head_compared_reads.md`](psp_head_compared_reads.md) — one milestone, lettered **H** so
its steps cannot be confused with this plan's A–G — adds a sixth field to the psp record head:
`reads-compared-with-reference`, the keep rule's denominator. The rule that admits a locus asks
for `max(2, 2% of the sample's compared reads)` non-reference reads, and the head carries the
numerator alone, so a future cheap-numbers read could apply the rule at three reads a position
(where the floor decides) and not at three hundred (where the share does).
[`run_streaming.md`](../spec/run_streaming.md) §3.3 flagged the requirement and the settled head
never picked it up.

**The ordering is a hard constraint, and it comes from the format rather than from either
plan.** A head layout change costs nothing while no psp exists outside tests and costs a format
version afterwards — [`record.rs:133-137`](../../../../src/ng/psp/record.rs) says so at the
field list itself: *"It costs nothing today because no psp exists; from Milestone F it costs a
version."* **F is where this plan makes psps that people keep**, so H goes in before it.

H is independent of Milestones A–E and touches files none of them touch, so it slots wherever
it fits. **Sequence chosen 2026-09-04: A–E, then H, then F–G.** A fresh conversation picking
this plan up after Checkpoint E should read H's plan next, not F's.

---

## Scope

**In:**

- The psp header's remaining fields from spec §6.1 — analysed regions, segmentation inputs,
  the read-group table, the observation reach ceiling (`max_record_span`, owed by
  `psp_file_format.md` §11), read filters into provenance.
- `SampleObservationGatherer` (arch §3.3): one sample's alignment files in, its observation
  stream into a `PspWriter`, the census accumulator fed at the yield point, census file
  written at `finish()`.
- The `generate-psps` subcommand (name agreed 2026-08-28; recorded in
  `run_driver_direct_mode.md` Preconditions).
- A psp-backed `ObservationSource`: an adapter over `PspReader::records()`, plus lifting the
  source-agnostic calling tail out of `AlignedFilesVariantCaller` so both modes drive one
  body (spec §3.1: "the two differ only in what a source is").
- `PspVariantCaller` (arch §3.4): open every header, run every §6.2 refusal, merge the
  read-group tables into a run-wide numbering, call through the lifted tail.
- The `call-from-psps` subcommand, `--parameters`/`--defaults` exactly as direct mode.
- The oracles: mode equivalence (§12.3), worker-count invariance end to end (§12.1), file
  order (§12.6), separately-walked cohort (§12.7), analysed-but-empty round trip (§12.9),
  every refusal at construction (§12.5).

**Out (later plans — nothing dropped):**

- **The fit stage**: a run that fits parameters from census files, the census built *from a
  psp*, the §7.12 byte-for-byte census-equality oracle, and the `generate-census` migration
  command (`parameter_prepass_joint_records.md` §6.1) — a follow-on plan
  ([`parameter_prepass_runs.md`](parameter_prepass_runs.md), written 2026-09-04). This plan
  writes the census beside the psp (Milestone G) so no sample ever needs re-walking; it does
  not read one back.
- **psp-mode performance**: the cheap-numbers source method (spec §3.3 step 2, deferred at
  spec §10 / arch §8), the shared contig list (spec §7.2, §10 — 357 kB a sample), record
  leasing through the `spare` offer (Milestone G of the direct-mode plan, deferred by the
  owner 2026-09-02), and §11 q7's psp half. Correctness first; a performance plan follows
  the first measured psp-mode run.
- **Open questions the spec holds**: reference bases in the record (§11 q4 — the codec's
  current answer ships as-is), cohorts over different analysed regions (§11 q5 — §6.2's
  refusal is the behaviour), the trailer's contents, the block byte ceiling's default
  (`psp_file_format.md` §12 q2).

## Questions routed upstream — all three ruled (owner, 2026-09-03)

1. **The record count is dropped from the header** (spec §6.1 amended).
   The census identity's count travels beside the header — `PileupIdentity::of_header`
   already takes it as its own argument, supplied by `WriteStats` — so Milestone A builds
   the header without it.
2. **§6.3 confirmed: no boundary digest, no `writer_version`** (spec §6.3 and arch §4/§6
   amended) — a cohort whose files cut their blocks in different places must never be
   refused; alignment is a property same-block-size files get from the coordinate grid by
   construction, never one a run demands.
3. **The walk's concurrency is ruled** (spec §5.2 and §11 q2 amended):
   `generate-psps` processes its samples **one at a time, in the order given** — no
   samples-in-flight knob, no in-process fan-out. Each sample's generation is independent,
   so a cohort is parallelised by running invocations, typically one sample each. That
   independence is the critical difference from direct mode, which must hold every sample
   open at one shared frontier; the walk stage never does.

## Principles (how the order was chosen)

- **Format before wiring.** The header grows its fields (A) before anything writes a file,
  because **no ng psp exists outside tests yet** — extending `Header` today is free, and
  after `generate-psps` ships every field added is a version negotiation.
- **The gatherer before its command; the source before its caller.** Each run object is
  built and proven against an in-memory oracle before the CLI that drives it (the
  heart-before-plumbing rule).
- **Direct mode is the oracle at every step** (spec §2: "direct mode is then psp mode's
  oracle"). The gatherer's file must read back equal to the walker's own stream; the psp
  route's VCF must equal direct mode's bytes.
- **Reuse over rewrite.** The walker (`AlignmentFilesWalker`), the segmentation
  (`segments_over`), the merge, the calling loop, the VCF writer, the parameters file and
  the census accumulator all exist; this plan writes drivers and one adapter, no new
  algorithmic code.
- **Isolate the silent-failure refactor.** Lifting the calling tail (D2) changes no
  behaviour and could break everything quietly; it lands as its own commit with direct
  mode's VCF byte-identical before and after.
- **Incremental, with pauses.** One milestone, then a checkpoint.

## Preconditions (already in place — verified 2026-09-03)

- The psp store: `src/ng/psp/` complete (`psp_file_format.md` A–H ✅), with
  `PspWriter::create(path, Header)` / `push(&SampleLocusObservations)` /
  `finish(trailer) → WriteStats` (`writer.rs:96,318,388`) and
  `PspReader::open/records/records_from` (`reader.rs:85,257,340`); worker-count byte
  identity pinned (`writer.rs:1394`); repeat-tract records round-trip (`record.rs`,
  `locus-kind` field).
- Direct mode complete (`run_driver_direct_mode.md`, 2026-09-01): `run_call_from_alignments`
  (`src/pop_var_caller_exp/call_from_alignments.rs:687`), `AlignedFilesVariantCaller`
  (`src/ng/run/callers.rs:180`), the per-sample walker seam
  (`AlignedFilesVariantCaller::walkers`, `callers.rs:439`; `AlignmentFilesWalker`,
  `walker.rs:89`).
- The source seam is a trait with a blanket impl for any
  `Iterator<Item = Result<SampleLocusObservations, E>>`
  (`cohort_merge/observation_cache.rs:70,98`).
- `Segmentation` + `SegmentationInputs::first_difference` (`src/ng/run/segments.rs:36,166`)
  — the §6.2 comparison operand.
- The census machinery: accumulator (`parameter_estimation/joint/census.rs`) and file
  (`census_file.rs`: `write_census:200`, `PileupIdentity:81`), driven today only by
  `examples/ng_joint_records_walk.rs:1131-1146` — the pattern Milestone G wires.
- Six of the eight `RunError` refusals built; `AnalysedRegionsDiffer`,
  `SegmentationInputsDiffer` (and `SampleAppearsTwice`) are psp mode's to add
  (arch §6).
- Subcommand names agreed (owner, 2026-08-28): `generate-psps`, `call-from-psps`.

---

## The steps

### Milestone A — the header carries what the run will check (format, no wiring)

**A1. ✅ Typed header fields for the run's identity checks.**
`Header` gains the analysed regions and the segmentation inputs (catalog identity + routing
criteria), as typed TOML sections — they are what §6.2's refusals compare, and
`SegmentationInputs` is the operand type. Round-trip and map-order-independence tests extend
the existing ones. *Depends:* —. *Source:* spec §6.1, §6.2; `psp_file_format.md` §3.1.

**A2. ✅ The sample's read-group table in the header.**
`@RG ID`, library, and the walk-local identifier per group (spec §6.1), so calling can merge
tables into a run-wide numbering (§6.2) without opening blocks. *Depends:* A1 (same encode
path). *Source:* spec §6.1 l.676-732, §6.2 l.759-763.

**A3. ✅ The observation reach ceiling (`max_record_span`).**
The header field `psp_file_format.md` §11 owes and `cohort_merge.md` §13 routed here; the
reader exposes it, nothing consumes it yet (E4 does). *Depends:* A1. *Source:*
`psp_file_format.md` §3.1 l.183-196, §11; spec §6.1.

**A4. ✅ Read filters and the command line into `WriterProvenance`.**
The applied read-filter settings recorded as provenance parameters — recorded, never
compared (goal 4). *Depends:* —. *Source:* spec §6.1.

> **Checkpoint A: the header holds §6.1 minus the record count (routed upstream, item 1).
> Owner confirms the typed-fields shape and the §6.3 stance. Pause for review.**
>
> **Passed (owner, 2026-09-03), three rulings recorded:**
>
> 1. **The typed-fields shape stands** — the catalog is recorded whole in the
>    `[segmentation]` section. Measured: 30,000 digest-carrying scaffolds encode to
>    10,798,518 bytes of the 16,777,187-byte header ceiling.
> 2. **`SegmentationInputs` moves to its own top-level module** (`src/ng/segmentation_inputs.rs`)
>    as Milestone B's first commit, with the `ng::run` re-export kept so no call site churns —
>    psp and run become mutually dependent in B, and the lift resolves the direction.
> 3. **No format-version bump for Milestone A's required header fields.** The reason: no psp
>    outside the test suite predates them, so the only file a missing-field refusal can reach
>    is a pre-A1 scratch file, and refusing it with a missing-field message is accepted.

### Milestone B — the gatherer: one sample's walk into a psp

**B1. ✅ `SampleObservationGatherer` over the existing walker.**
Arch §3.3's object: constructed from one sample's alignment files, the shared
`Arc<Segmentation>`, and a census configuration; internally the direct-mode chain
(`SampleReads` + `GeneratorSet` + `RunSegments` → `SampleLocusObservationsIterator`) that
`AlignedFilesVariantCaller::walkers` builds today, driving `PspWriter::push` per record.
Serial within the sample (spec §5.2; §11 q3 is off the default path). *Depends:* A1-A4 (the
header it writes). *Source:* spec §5.2 l.571-643; arch §3.3 l.432-469.

**B2. ✅ The gatherer's oracle: the file reads back as the walk streamed.**
On the shared run fixtures **and** one real CRAM slice: gather to a psp, read it back, and
compare record-for-record, field-for-field against the same sample walked directly in
memory. This is the plan's north star for the walk stage. *Depends:* B1. *Source:* spec §2
(direct mode is the oracle); §12.9 (an analysed-but-empty ground round-trips).

**B3. ✅ Worker-count/schedule invariance, end to end.**
The format-level byte-identity test (`writer.rs:1394`) re-made through the gatherer on real
reads: the same sample gathered twice gives byte-identical files but the header timestamp.
*Depends:* B1. *Source:* spec §12.1, §6.3.

> **Checkpoint B: a sample's psp is provably the walk, in bytes. Pause for review.**
>
> **Reached 2026-09-03.** Both oracles are green on real reads — one tomato accession
> (`SRR7279481.p1.bench.cram`) over 200 kb of SL4.0 through the real catalog: **183,807
> records, all equal field for field between the file and the walk, header equal; 1,217 of
> them repeat tracts; a second gather byte-identical (948,689 bytes)**
> (`examples/ng_psp_gather_oracle.rs`, which records the run). At fixture scale, 16 tests
> cover the header, the failure paths, the applied-vs-recorded settings, §12.9's
> analysed-but-empty round trip, tract ground, and byte identity.
>
> **Carried into Milestone C:**
> - the run-level shared test-fixture module (deferred at B1, deferred again at B2 — the
>   gatherer is the fifth consumer of fixtures `callers.rs` keeps private);
> - the real slice the oracle runs on lives in `tmp/tomato_slice/` (untracked, 83 MB);
>   a fixed home for it is C1's to decide when `generate-psps` gets its own runs.

### Milestone C — `generate-psps`

**C1. ✅ The subcommand.**
Reference, catalog, one `--alignment` per sample, optional BED, `--output-dir`; assembles
reference/segmentation exactly as `call_from_alignments.rs` does (reuse `segments_over`,
`analysed_regions`, `build_read_groups`), then one gatherer per sample. Refusals shared with
direct mode (the five pre-open checks). *Depends:* B1. *Source:* spec §2, §5.2; the agreed
names (`run_driver_direct_mode.md` l.108).

**C2. ✅ Samples one at a time, and a per-sample report.**
The loop over samples is sequential, in the order given (owner's ruling 2026-09-03, spec
§5.2) — no concurrency knob; a cohort is parallelised by running invocations. A finished
sample prints its `WriteStats`; the run ends with a per-sample report line. *Depends:* C1.
*Source:* spec §5.2; §11 q2 (walk half, answered).

**C3. ✅ An interrupted run leaves nothing that reads as whole.**
Kill-mid-write fixture at the command level: the half-written file is refused as interrupted
(the format already guarantees it; the command's exit and message are what this pins), and a
`--force`-less rerun refuses to overwrite finished files. *Depends:* C1. *Source:* spec §8
(l.1000-1004); `psp_file_format.md` §10.

> **Checkpoint C: a cohort of psps from the command line. Pause for review.**
>
> **Reached 2026-09-03.** `pop_var_caller_exp generate-psps` walks a cohort and writes one psp
> per sample. On one tomato accession over two 100 kb intervals: **193,603 loci stored,
> 914,715 bytes, about 3 s**, and the run says what ground it spoke for —
> *311 of 318 typed regions, 199,672 of 200,000 bases walked, 99.8%*, with the 328 bases it
> could not store named as *clusters of repeats too close together to have clean flanks*.
>
> **What C3 settled about a run that stops**: each sample is walked into
> `<sample>.psp.<pid>.partial` and renamed only once whole, so a stopped walk leaves nothing
> at the sample's own path and a stopped **re-walk** leaves the psp it was replacing intact.
> Without `--force` a run refuses as soon as it finds a psp already there — before walking
> anything, so a cohort is never left half-replaced.
>
> **Two prerequisites landed inside this milestone, both recorded rather than silent:**
> `ng::run`'s shared test fixtures (`70385f5b`, the carry-forward from B1 and B2), and the
> **ground assembly lifted out of direct mode** into `src/pop_var_caller_exp/run_ground.rs`
> (`f00d56e9`) — which is what C1's "reuse `segments_over`/`analysed_regions`" required, and
> what makes §6.2's cohort-agreement check meaningful: both modes now compute the segmentation
> inputs from one copy.
>
> **Carried into Milestone D:**
> - the reference-open block is still duplicated verbatim between the two commands (22 lines);
>   `run_ground`'s own doc already claims to own it;
> - the five repeat-routing `#[arg]` blocks are byte-identical in both commands —
>   `#[command(flatten)]` over the existing `RepeatRouting` would delete both copies and both
>   `ground_request`s (deferred at Medium confidence: no third consumer is coming, and the
>   drift hazard is closed by a whole-struct comparison plus per-flag tests);
> - the on-disk cohort fixture is duplicated between the two commands' test modules;
> - the read filters and the five locus-generator knobs are hard-coded in **both** commands and
>   invisible at either surface — four of the five are recoverable from no field of the psp.
>   Whether they become flags is the owner's call, not a step's.

### Milestone D — the psp source, and one calling tail for both modes

**D1. ✅ The psp-backed source.**
An iterator adapter over `PspReader::records()` mapping `StreamedRecord` to its record and
`PspReadError` into `RunError::SourceFailed` naming the sample (arch §5) — the blanket
`ObservationSource` impl does the rest. A source whose file yields observations out of order
returns the error the merge's contract asks for, not an assertion (arch §8's cohort-merge
item). Unit tests over reader fixtures, including a head-only-skip block and a tract record.
*Depends:* — (parallel to B/C). *Source:* spec §3.1, §3.4; arch §5, §8.

**D2. ✅ Lift the calling tail. Own commit, do not bundle.**
Everything in `call_cohort_handing_each_record_over` from `parameters.view()` on
(`callers.rs:705-…`) becomes a free function over `ObservationCache<S>` + parameters +
segmentation + sink; `AlignedFilesVariantCaller` calls it. **No behaviour change: direct
mode's VCF is byte-identical before and after on the fixtures and the six-accession tomato
slice — that oracle green on both sides of this one commit is the point of isolating it.**
*Depends:* —. *Source:* spec §3.1 l.177-204 ("both callers drive one merge and one calling
loop").

> **Checkpoint D: two sources, one tail, direct mode provably untouched. Pause for review.**
>
> **Reached 2026-09-04.** `PspObservationSource` decodes a stored sample behind the same trait
> direct mode's walker implements, and the calling loop both modes will drive is now a free
> function over any source (`call_cohort_from_sources_handing_each_record_over`), with
> `AlignedFilesVariantCaller` adding only what alignment files add: opening the walkers, and
> turning them, spent, into per-sample tallies.
>
> **The oracle this milestone exists for is green.** Direct mode's VCF over six tomato
> accessions and the first two 100 kb intervals of `benchmarks/tomato1/regions.bed` — **598
> records, sha256 `5f0903cf…`** — is byte-identical either side of the lift, and so are the
> parameters file and the whole run report.
>
> **Two things the reviews caught that the code did not**: after refusing a record the psp
> source went on handing over the *next* one, so a swallowed error lost an observation in
> silence (a refusal now ends the source and says so on every later draw); and the rule that a
> refused record outranks a source failing afterwards had no test — the lift is what made one
> writable, because the loop now takes a `Vec`'s iterator as a source.
>
> **Carried into Milestone E:**
> - the E-milestone interface question is settled by compilation, not by argument: a reviewer
>   wrote `PspVariantCaller`'s method as E will have to write it — `Vec<PspReader>` out of
>   `self`, one source each, into `ObservationCache::over`, through the lifted loop, sources
>   read back afterwards — and it type-checks. Nothing in the signature forces E to change it;
> - arch §3.4's `PspVariantCaller::open` sketch still lists `callers_in_flight: CallersInFlight`,
>   which arch §8 struck on 2026-09-01 and which no type in `src/` implements. The sketch is
>   stale, not a constraint;
> - Milestone C's four carry-forwards are all still open, none touched by D.

### Milestone E — `PspVariantCaller`

**E1. ✅ Open the cohort: every header, every refusal, before any block.**
`PspVariantCaller::open` (arch §3.4): read every header; refuse per §6.2 — analysed regions
equal across the cohort (`AnalysedRegionsDiffer`, new variant), segmentation inputs match
the run's (`SegmentationInputsDiffer`, new variant, via `first_difference`), two files
naming one sample (`SampleAppearsTwice`), contigs agreeing with the run's reference, sample
sets matching the parameters both ways. Analysed regions come from the headers, not a flag.
Every refusal has a provoking test naming what differs (§12.5). *Depends:* A1-A2, D1.
*Source:* spec §5.3, §6.2; arch §3.4, §6.

**E2. ✅ Merge the read-group tables into the run-wide numbering.**
Each sample numbers from zero; build the run table by merging and remap each observation's
walk-local identifier as records are drawn — the remap lives in the psp source, keyed by the
header's table (A2). Fixture: two samples whose local ids clash and whose libraries differ;
pinned by the per-read-group outputs landing under the right sample. *Depends:* E1, A2.
*Source:* spec §6.2 l.759-763, §6.1 l.703-708 (the rejected alternative).

**E3. ✅ Call the cohort through the lifted tail.**
`call_cohort_handing_each_record_over` on `PspVariantCaller`: sources into
`ObservationCache::over`, the D2 tail does the rest — merge, genotyping, records to the
sink. *Depends:* D2, E1, E2. *Source:* spec §3.1, §3.5; arch §3.4.

**E4. ✅ The reach ceiling read where the merge needs it.**
The A3 header field consumed at open per `cohort_merge.md` §13's routing (the deferred
reader at `observation_cache.rs:420`). *Depends:* A3, E1. *Source:* `cohort_merge.md` §13.

> **Checkpoint E: a cohort of psps calls, refusals all fire at construction. Pause for
> review.**
>
> **Reached 2026-09-04.** `OpenPspCohort` reads every header and settles what the cohort is —
> the ground the files agree on, one run-wide read-group numbering, the widest observation reach
> any of them declares — and `PspVariantCaller` checks each file against the run and drives the
> lifted calling loop over one source per open psp. `ng::run` 475 → 499; 24 mutations across the
> milestone, 23 killed.
>
> **Three checks the spec asks for that nothing had**, all found by review and all now made: the
> descriptor budget of §7.1a (before the first `open(2)`, or a cohort of thousands dies at the
> 249th file on a macOS default limit, blaming an innocent one), direct mode's
> catalog-against-reference refusal, and the psp header's own whole-assembly digest — the
> field's one documented consumer.
>
> **Two rulings owed by the owner, both recorded where the code makes the choice:**
>
> 1. **§6.2's duplicate-`@RG ID` refusal is not made.** The spec's reason is that such a table
>    "cannot be renumbered without guessing"; this format guesses nothing, because identity is
>    the walk-local number, which is the entry's own position, and nothing in the merge reads
>    the id. `psp/header.rs`'s own validator declares the case legal — a psp holds one *sample*,
>    not one alignment file, and a sample sequenced across lanes may carry two entries with one
>    id and different libraries — and direct mode calls that cohort. Refusing it would break
>    §1.1's goal 1 for every multi-lane sample whose lanes reuse an id. **Recommendation: §6.2's
>    clause should say what it means, which is the empty table.**
> 2. **§6.2's by-name parameters match is F1's, not E1's.** `RunParameters` is assembled per
>    sample by position and carries no names at all, so the match has to happen where the
>    parameters *file* meets this cohort's sample list. E1 makes the count, exactly as direct
>    mode does. **Recommendation: amend E1's entry above to say so.**
>
> **Carried into the rest of the plan:**
> - **[`psp_head_compared_reads.md`](psp_head_compared_reads.md) goes next, before Milestone F**
>   — see this plan's ordering note at the top;
> - a run over stored files returns the calling tallies alone, where direct mode also returns
>   per-sample walk tallies: a psp source has none to give, so **what such a run says about each
>   sample is F1's to decide**;
> - a psp records the library a walk resolved and not whether the file declared it, so the
>   merged table marks every library `Synthesized` — the claim that cannot be false. A header
>   field would let it say the true one; nothing reads the origin today;
> - Milestone C's four carry-forwards are all still open, none touched by D or E.

### Milestone F — `call-from-psps`, and the oracle that justifies the design

> **⚑ Precondition: [`psp_head_compared_reads.md`](psp_head_compared_reads.md)'s Milestone H is
> committed.** After F, a change to the record head costs a format version
> ([`record.rs:133-137`](../../../../src/ng/psp/record.rs)); before it, nothing. If H has not
> landed, stop and land it first — see this plan's ordering note above.

**F1. ✅ The subcommand.**
One `--psp` per sample (or a directory), `--parameters`/`--defaults` exactly as direct mode
(`run_parameters` reused; the census argument stays `None` until the fit plan), the same
VCF + parameters-file + run-report outputs. *Depends:* E3. *Source:* spec §2, §5.3; the
agreed names.

> **Done 2026-09-04** ([report](../../reports/implementations/ng_psp_mode_f1_2026-09-04.md),
> [review](../../reports/reviews/ng_psp_mode_f1_2026-09-04.md)). **What a run over stored files
> says about each sample is ruled** (owner, 2026-09-04): how many stored loci it read out of each
> file and how many reads went into the comparison at one of them, both measured by the run as it
> decodes; a file holding no loci over this ground named as contributing nothing; and one line,
> printed only where the cohort's psps disagree, naming the files whose walk applied other read
> filters — the one comparison of a field spec §6.1 records and nothing checks.
>
> **The numbers both calling commands score with now live in one module**
> (`src/pop_var_caller_exp/calling_run.rs`), `run_ground`'s sibling: the parameters, the `NG_*`
> measurement switches, the round width, the ploidy, the output refusals, the VCF header metadata
> and the report printer. Two copies of those would be two places for the modes to drift while
> both kept running, which is what F2 compares. Direct mode's VCF, parameters file and run report
> are byte-identical across the lift.

**F2. ✅ Mode equivalence. Own commit, do not bundle.**
Same cohort, same parameters: `call-from-alignments` and `generate-psps` +
`call-from-psps` produce **the same VCF** — bytes, at the default block size — on the run
fixtures and on six tomato accessions over 400 kb through the real catalog. (Comparing
*across block sizes* is the weaker §10.1 tolerance oracle and is not this test.) This is
§12.3, "the oracle that justifies the design", and goal 1's proof. *Depends:* C1, F1.
*Source:* spec §12.3, §1.1 goal 4; `psp_file_format.md` §10.1.

> **Done 2026-09-04** ([report](../../reports/implementations/ng_psp_mode_f2_2026-09-04.md)).
> **599 records, byte-identical apart from `##commandline`** on six tomato accessions over the
> two 100 kb intervals, parameters file identical too —
> `scripts/ng_mode_equivalence_oracle.sh`, in the repository so the run reproduces. In the suite,
> `pop_var_caller_exp::mode_equivalence` compares the two routes' VCFs **whole**, nothing
> filtered: inside one process both routes record the same command line, so the script's one
> exemption is not needed there.
>
> **The fixture had to be built for it**, because the cohort the commands' own tests use has an
> all-`A` reference and writes no record at all. Each of its four discriminating properties
> closes a defect measured surviving without it: two samples varying in different places, a
> repeat tract with a length variant, alternative reads leaning to one strand, and — because
> `--defaults` scores every read group alike — a second comparison under parameters whose three
> read groups carry different multipliers, which is also the only place psp mode's supplied-file
> path meets direct mode's.
>
> **Where it stops, stated rather than assumed**: it compares VCFs, so a stored locus with no
> variant is not compared (578 of 581 a sample); neither route fits, so what only a fit reads can
> be destroyed on write with both comparisons green.

**F3. ✅ The remaining run-level invariances.**
File order does not matter (§12.6, same VCF sample-for-sample); a cohort of
separately-walked samples calls (§12.7 — E2's fixture, end to end); analysed-but-empty
ground round-trips to the VCF's absence of records (§12.9); concurrency invariance of the
psp-route VCF at pools of 1/2/4/8 (§12.2). *Depends:* F1. *Source:* spec §12.

> **Checkpoint F: psp mode exists and equals direct mode. Pause for review.**
>
> **Reached 2026-09-04.** `pop_var_caller_exp call-from-psps` calls a cohort of stored files and
> writes the VCF, the parameters file and the run report direct mode writes. **On six tomato
> accessions over the two 100 kb intervals: 599 records, byte-identical to direct mode's output
> apart from the `##commandline` line, and the parameters file identical too**
> (`scripts/ng_mode_equivalence_oracle.sh`). All four of §12's remaining run-level invariances
> hold (`scripts/ng_psp_concurrency_invariance.sh` for the thread sweep, tests for the other
> three).
>
> **Two rulings are recorded in this milestone:**
>
> 1. **What a run over stored files says about each sample** (owner, 2026-09-04): the loci it
>    read out of each file and how many reads went into the comparison at one of them, both
>    measured by the run; a file holding none named as contributing nothing; and a line, only
>    where the cohort's psps disagree, naming the files whose walk applied other read filters.
>    **The number is not depth and the line does not say depth** — the head's count excludes
>    filtered reads, depth-capped reads and reads that cover a locus without anchoring it.
> 2. **`##commandline` is exempted from the shell-level comparison and from nothing else.** The
>    two routes are two commands, so the line differs by construction; inside one test process
>    both record the same line and the suite compares the files whole.
>
> **Carried into Milestone G:**
> - the oracle compares VCFs, so a stored locus that produces no record is not compared — 578 of
>   the 581 a sample in the fixture. Field-for-field equality with the walk is Milestone B's;
> - neither route fits parameters, so what only a fit reads (the stored sum of squared mapping
>   qualities, the count of reads that covered a locus without an observation) can be destroyed
>   on write with every comparison green. **G writes the census, whose fit does read them**;
> - a psp written against a reference with a different contig count gets the catalog's refusal
>   rather than its own, because the segmentation is built before the cohort is checked against
>   the reference;
> - the per-sample report section has no cap at the thousand-sample end, and neither does direct
>   mode's — it belongs to both or neither;
> - Milestone C's four carry-forwards: two are closed (the on-disk cohort fixture is shared by
>   all three commands; the parameters assembly and the VCF header are lifted into
>   `calling_run`), and two remain — the duplicated reference-open block, and the read filters
>   and locus-generator knobs being invisible at every surface.

### Milestone G — the census beside the psp

**G1. ✅ The gatherer feeds the census accumulator.**
At the gatherer's ordered yield point, each locus into the accumulator the joint-records
walk example already drives (`examples/ng_joint_records_walk.rs:1131-1146`); `finish()`
writes the census file beside the psp via `write_census`, its `PileupIdentity` built from
the psp's own header — the identity's first real construction site. *Depends:* B1, A1-A4.
*Source:* spec §5.2 l.611-633, §1.2 l.84-87; `parameter_prepass_joint_records.md` §6.1.

**G2. ✅ `generate-psps` writes both files; the walk is once.**
The command's report names both; a census that cannot be written fails the sample's walk
(spec §2: alignments read exactly once — a psp without its census would force a re-walk).
*Depends:* G1, C1. *Source:* spec §2 l.152-171.

> **Both done 2026-09-04, in one commit and deliberately**
> ([report](../../reports/implementations/ng_psp_mode_g_2026-09-04.md)): G1 changes
> `write_psp`'s signature, so G2's command has to move in the same breath or the tree does not
> build — the two steps share one loop iteration rather than being split artificially.
>
> **Measured on six tomato accessions over the two 100 kb intervals: 3,592,149 bytes of psp and
> 1,305,915 bytes of census, in 5 s**, both files named per sample in the run's report. The
> census is fed at the walk's yield point, not by the psp writer, so it records what the walk saw
> rather than what was stored.
>
> **The selection's numbers**: about two million positions and five thousand tracts a stratum are
> the design's own figures; **the seed is a compiled-in constant**, which is what lets two
> invocations of one cohort keep the same positions — a seed that differed between them would
> keep disjoint sets and the samples could not be pooled. Whether the three become flags is the
> same open question Milestone C recorded about the read filters.
>
> **A defect the milestone's own test caught first**: `PspWriter::create` amends the header
> before writing it, so a census built from the header the gatherer *holds* names a psp that does
> not exist, and every freshness check would have said *rebuild* for ever. `WriteStats` now
> carries the digest of the header as written.

> **Checkpoint G: the walk stage is spec §2's, whole. The fit stage and the census-equality
> oracle (§7.12) hand to the next plan. Pause for review.**
>
> **Reached 2026-09-04, and the plan is finished.** `generate-psps` reads each sample's
> alignment files once and writes both files spec §2 gives the walk stage; `call-from-psps` calls
> a cohort of the first kind and produces direct mode's VCF. Nothing in this plan is open.
>
> **Carried to the next plan** ([`parameter_prepass_runs.md`](parameter_prepass_runs.md), written
> 2026-09-04): reading a census back, building one *from* a psp, and §7.12's byte-for-byte
> census-equality oracle. **Two of this plan's own gaps close there rather than here** — the
> mode-equivalence oracle cannot see the stored fields only a fit reads (the sum of squared
> mapping qualities, and the count of reads that covered a locus without producing an
> observation), and nothing yet reads a census at all.
>
> **Ruled by the owner, 2026-09-04: they all stay constants until the fit stage exists.** The
> census selection's seed and two counts, the read filters, and the five locus-generator knobs —
> one question, since Milestone C recorded the second half and G added the first. The reason is
> that **nothing can read a census yet**, so a knob added now is one whose effect nobody can
> check; and the seed in particular is a way to break a cohort silently, since two invocations
> that seeded differently walk perfectly and are refused hours later at the fit. Revisit when the
> fit stage lands and there is something to vary them against.

---

## Verification summary

| milestone | proven by |
|---|---|
| A — header fields | round-trip + order-independence + refusal tests at the format level |
| B — gatherer | the file read back equals the direct walk's stream, field for field, fixtures + real CRAM; schedule invariance in bytes |
| C — generate-psps | command-level fixtures; interrupted-write refusal; per-sample `WriteStats` |
| D — source + lifted tail | direct-mode VCF **byte-identical** across the D2 commit; adapter unit tests incl. tract records |
| E — PspVariantCaller | every §6.2 refusal provoked and named; the RG-clash fixture's outputs land under the right samples |
| F — call-from-psps | §12.3 **same VCF as direct mode** (fixtures + tomato slice); §12.2/6/7/9 |
| G — census | census written beside psp, identity tied to the psp header; walk count stays one |

## Out of scope (next plans)

- **[`parameter_prepass_runs.md`](parameter_prepass_runs.md)** (written 2026-09-04): the fit
  stage, census-from-psp, §7.12's byte-for-byte census equality, `generate-census`.
- **psp-mode performance** (after the first measured run): the cheap-numbers read (spec
  §3.3/§10), the shared contig list (§7.2/§10), leasing through `spare`, §11 q7's psp half,
  and q2's remaining callers-in-flight half.
- **The trailer's contents** — opaque bytes until something needs them
  (`psp_file_format.md` §3.4).
