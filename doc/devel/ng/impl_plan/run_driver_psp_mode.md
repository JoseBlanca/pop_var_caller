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
  (`parameter_prepass_runs.md`, unwritten). This plan writes the census beside the psp
  (Milestone G) so no sample ever needs re-walking; it does not read one back.
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

**B2. ☐ The gatherer's oracle: the file reads back as the walk streamed.**
On the shared run fixtures **and** one real CRAM slice: gather to a psp, read it back, and
compare record-for-record, field-for-field against the same sample walked directly in
memory. This is the plan's north star for the walk stage. *Depends:* B1. *Source:* spec §2
(direct mode is the oracle); §12.9 (an analysed-but-empty ground round-trips).

**B3. ☐ Worker-count/schedule invariance, end to end.**
The format-level byte-identity test (`writer.rs:1394`) re-made through the gatherer on real
reads: the same sample gathered twice gives byte-identical files but the header timestamp.
*Depends:* B1. *Source:* spec §12.1, §6.3.

> **Checkpoint B: a sample's psp is provably the walk, in bytes. Pause for review.**

### Milestone C — `generate-psps`

**C1. ☐ The subcommand.**
Reference, catalog, one `--alignment` per sample, optional BED, `--output-dir`; assembles
reference/segmentation exactly as `call_from_alignments.rs` does (reuse `segments_over`,
`analysed_regions`, `build_read_groups`), then one gatherer per sample. Refusals shared with
direct mode (the five pre-open checks). *Depends:* B1. *Source:* spec §2, §5.2; the agreed
names (`run_driver_direct_mode.md` l.108).

**C2. ☐ Samples one at a time, and a per-sample report.**
The loop over samples is sequential, in the order given (owner's ruling 2026-09-03, spec
§5.2) — no concurrency knob; a cohort is parallelised by running invocations. A finished
sample prints its `WriteStats`; the run ends with a per-sample report line. *Depends:* C1.
*Source:* spec §5.2; §11 q2 (walk half, answered).

**C3. ☐ An interrupted run leaves nothing that reads as whole.**
Kill-mid-write fixture at the command level: the half-written file is refused as interrupted
(the format already guarantees it; the command's exit and message are what this pins), and a
`--force`-less rerun refuses to overwrite finished files. *Depends:* C1. *Source:* spec §8
(l.1000-1004); `psp_file_format.md` §10.

> **Checkpoint C: a cohort of psps from the command line. Pause for review.**

### Milestone D — the psp source, and one calling tail for both modes

**D1. ☐ The psp-backed source.**
An iterator adapter over `PspReader::records()` mapping `StreamedRecord` to its record and
`PspReadError` into `RunError::SourceFailed` naming the sample (arch §5) — the blanket
`ObservationSource` impl does the rest. A source whose file yields observations out of order
returns the error the merge's contract asks for, not an assertion (arch §8's cohort-merge
item). Unit tests over reader fixtures, including a head-only-skip block and a tract record.
*Depends:* — (parallel to B/C). *Source:* spec §3.1, §3.4; arch §5, §8.

**D2. ☐ Lift the calling tail. Own commit, do not bundle.**
Everything in `call_cohort_handing_each_record_over` from `parameters.view()` on
(`callers.rs:705-…`) becomes a free function over `ObservationCache<S>` + parameters +
segmentation + sink; `AlignedFilesVariantCaller` calls it. **No behaviour change: direct
mode's VCF is byte-identical before and after on the fixtures and the six-accession tomato
slice — that oracle green on both sides of this one commit is the point of isolating it.**
*Depends:* —. *Source:* spec §3.1 l.177-204 ("both callers drive one merge and one calling
loop").

> **Checkpoint D: two sources, one tail, direct mode provably untouched. Pause for review.**

### Milestone E — `PspVariantCaller`

**E1. ☐ Open the cohort: every header, every refusal, before any block.**
`PspVariantCaller::open` (arch §3.4): read every header; refuse per §6.2 — analysed regions
equal across the cohort (`AnalysedRegionsDiffer`, new variant), segmentation inputs match
the run's (`SegmentationInputsDiffer`, new variant, via `first_difference`), two files
naming one sample (`SampleAppearsTwice`), contigs agreeing with the run's reference, sample
sets matching the parameters both ways. Analysed regions come from the headers, not a flag.
Every refusal has a provoking test naming what differs (§12.5). *Depends:* A1-A2, D1.
*Source:* spec §5.3, §6.2; arch §3.4, §6.

**E2. ☐ Merge the read-group tables into the run-wide numbering.**
Each sample numbers from zero; build the run table by merging and remap each observation's
walk-local identifier as records are drawn — the remap lives in the psp source, keyed by the
header's table (A2). Fixture: two samples whose local ids clash and whose libraries differ;
pinned by the per-read-group outputs landing under the right sample. *Depends:* E1, A2.
*Source:* spec §6.2 l.759-763, §6.1 l.703-708 (the rejected alternative).

**E3. ☐ Call the cohort through the lifted tail.**
`call_cohort_handing_each_record_over` on `PspVariantCaller`: sources into
`ObservationCache::over`, the D2 tail does the rest — merge, genotyping, records to the
sink. *Depends:* D2, E1, E2. *Source:* spec §3.1, §3.5; arch §3.4.

**E4. ☐ The reach ceiling read where the merge needs it.**
The A3 header field consumed at open per `cohort_merge.md` §13's routing (the deferred
reader at `observation_cache.rs:420`). *Depends:* A3, E1. *Source:* `cohort_merge.md` §13.

> **Checkpoint E: a cohort of psps calls, refusals all fire at construction. Pause for
> review.**

### Milestone F — `call-from-psps`, and the oracle that justifies the design

**F1. ☐ The subcommand.**
One `--psp` per sample (or a directory), `--parameters`/`--defaults` exactly as direct mode
(`run_parameters` reused; the census argument stays `None` until the fit plan), the same
VCF + parameters-file + run-report outputs. *Depends:* E3. *Source:* spec §2, §5.3; the
agreed names.

**F2. ☐ Mode equivalence. Own commit, do not bundle.**
Same cohort, same parameters: `call-from-alignments` and `generate-psps` +
`call-from-psps` produce **the same VCF** — bytes, at the default block size — on the run
fixtures and on six tomato accessions over 400 kb through the real catalog. (Comparing
*across block sizes* is the weaker §10.1 tolerance oracle and is not this test.) This is
§12.3, "the oracle that justifies the design", and goal 1's proof. *Depends:* C1, F1.
*Source:* spec §12.3, §1.1 goal 4; `psp_file_format.md` §10.1.

**F3. ☐ The remaining run-level invariances.**
File order does not matter (§12.6, same VCF sample-for-sample); a cohort of
separately-walked samples calls (§12.7 — E2's fixture, end to end); analysed-but-empty
ground round-trips to the VCF's absence of records (§12.9); concurrency invariance of the
psp-route VCF at pools of 1/2/4/8 (§12.2). *Depends:* F1. *Source:* spec §12.

> **Checkpoint F: psp mode exists and equals direct mode. Pause for review.**

### Milestone G — the census beside the psp

**G1. ☐ The gatherer feeds the census accumulator.**
At the gatherer's ordered yield point, each locus into the accumulator the joint-records
walk example already drives (`examples/ng_joint_records_walk.rs:1131-1146`); `finish()`
writes the census file beside the psp via `write_census`, its `PileupIdentity` built from
the psp's own header — the identity's first real construction site. *Depends:* B1, A1-A4.
*Source:* spec §5.2 l.611-633, §1.2 l.84-87; `parameter_prepass_joint_records.md` §6.1.

**G2. ☐ `generate-psps` writes both files; the walk is once.**
The command's report names both; a census that cannot be written fails the sample's walk
(spec §2: alignments read exactly once — a psp without its census would force a re-walk).
*Depends:* G1, C1. *Source:* spec §2 l.152-171.

> **Checkpoint G: the walk stage is spec §2's, whole. The fit stage and the census-equality
> oracle (§7.12) hand to the next plan. Pause for review.**

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

- **`parameter_prepass_runs.md`** (unwritten): the fit stage, census-from-psp, §7.12's
  byte-for-byte census equality, `generate-census`.
- **psp-mode performance** (after the first measured run): the cheap-numbers read (spec
  §3.3/§10), the shared contig list (§7.2/§10), leasing through `spare`, §11 q7's psp half,
  and q2's remaining callers-in-flight half.
- **The trailer's contents** — opaque bytes until something needs them
  (`psp_file_format.md` §3.4).
