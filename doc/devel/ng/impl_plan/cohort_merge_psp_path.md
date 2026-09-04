# The cohort merge's psp path — summaries first, evidence on demand: implementation plan

*Draft, 2026-09-04. Turns [`../spec/cohort_merge_psp_path.md`](../spec/cohort_merge_psp_path.md)
into build order — no new design here; a design question surfacing mid-plan goes back to that
spec. Successor to [`run_driver_psp_mode.md`](run_driver_psp_mode.md), which reserves this
plan's slot under "psp-mode performance". **The standing oracle for every milestone is that
plan's F2: `call-from-alignments` and `generate-psps` + `call-from-psps` produce the same VCF,
byte for byte, on the run fixtures and the tomato slice.** A milestone that moves a byte is
wrong, not almost done.*

## Scope

**In:** the two-state cached observation; the two-cursor psp source with its retained window;
the record-equality oracle against the decode-everything source; the shared contig list; the
spare offer, measured; the end-to-end timing before and after, and the scaling verdict it
gates.

**Out:**

- **anything structural for scaling** (per-sample cover parallelism, overlapping cover with
  building) — gated by Milestone E's shares; if justified, its own plan, its conclusion owed to
  [`run_streaming.md`](../spec/run_streaming.md) §11 q7 either way;
- **dropping stored reference bases** ([`run_streaming.md`](../spec/run_streaming.md) §11 q4) —
  untimed leaning, untouched;
- **any change to a verdict, a format byte, or a CLI flag** (spec §1.2).

## Principles (how the order was chosen)

- **Baseline before the change.** Nobody has timed a psp-mode run end to end; the same
  measurement is the before/after evidence and the gate on the structural work, so it runs
  first and last.
- **The merge-side refactor lands alone, in direct mode, provably inert** — the D2 pattern
  from the predecessor plan: reshape the seam with the VCF byte-identical either side, then
  put the new source behind it.
- **The earlier implementation is the next one's oracle.** The decode-everything source is
  demoted to tests, never deleted: the two-phase source must equal it record for record before
  the VCF oracle is even interesting.
- **Adopt allocator and memory changes on numbers, not shapes** — the spare offer and the
  window's storage are both "measure, then keep or revert".

## Preconditions (verify before step A1)

- **A callable psp path** — `PspVariantCaller::open` and
  `call_cohort_handing_each_record_over` ([`psp_caller.rs:286,425`](../../../../src/ng/run/psp_caller.rs)),
  both public and merged to main with [`run_driver_psp_mode.md`](run_driver_psp_mode.md)
  Milestone E. **Met.** Milestone A drives them from a probe; it does not need the
  `call-from-psps` subcommand, which is that plan's F1 and still open.
- **F2's mode-equivalence oracle, before Milestone B changes anything** —
  [`run_driver_psp_mode.md`](run_driver_psp_mode.md) F2, **open as of 2026-09-04** and owned by
  the conversation running that plan. Milestone A needs no oracle, being measurement only;
  **B onward do**. **Met 2026-09-04**: F1, F2 and F3 are merged, `call-from-psps` exists, and
  the oracle passes on this branch over three tomato accessions — the two modes' VCFs identical
  apart from the recorded command line, 2,253 records, same sha256, parameters file identical
  too (`scripts/ng_mode_equivalence_oracle.sh`).
- **The head carries the keep rule's denominator** —
  [`psp_head_compared_reads.md`](../spec/psp_head_compared_reads.md) Milestone H, **met**
  (main, `47d4c7e1`). Without it phase 1 cannot apply the keep rule at depth. The locus kind
  is *not* among the prerequisites: the owner's ruling of 2026-09-04 (`fe9df2a3`) keeps the
  tag in the body and has a summary look its kind up from the run's segmentation by
  coordinate (spec §2).
- The seam and machinery the spec's reuse map names, at their cited lines: the head-only and
  selective walks (`ng/psp/walk.rs`), the bounded body decode (`ng/psp/record.rs`), the
  verdict's number pair (`run/cohort_merge/close.rs`), the source adapter (`run/psp_source.rs`),
  the cache's cover/evict (`run/cohort_merge/observation_cache.rs`).
- Benchmark inputs: the 63 tomato CRAMs and BEDs (`benchmarks/tomato1/`), HG002, and psps
  generated from them by `generate-psps`.

## The steps

### Milestone A — the baseline: what a psp-mode run costs today

**A1. ✅ Time `call-from-psps` end to end and by stage.** *Done 2026-09-04, spec §5 carries the
numbers. Both depth corners measured against direct mode on the same ground: the calling phase
is 2.05 s against 13.91 s at 63 tomato accessions and 0.05 s against 3.34 s at one HG002 sample,
with reading records back 43% and 91% of it. The skip is worth 2.57× at the cohort's true keep
rate of one record in eight and 2.66× at 280 reads a position — it does not shrink with depth,
which is a prediction of the spec's this refuted. **Checkpoint A's condition is answered: decode
plus merge is not a small share, so B–D are worth building.*** Wall and peak resident at 1, 16 and
63 samples on the tomato benchmark and 1 sample on HG002; the run's time split into decode
(the source), merge, and calling — instrumentation behind a feature or a throwaway probe, not
shipped in the hot path. Report the merge's and decode's shares: these gate Milestone E.
Median of five runs with the spread, per the measurement hygiene in
[`../research/cohort_merge_parallel_cost_plan.md`](../research/cohort_merge_parallel_cost_plan.md) §2.
*Depends:* preconditions. *Source:* spec §5.

> **Checkpoint A: the numbers exist. Pause for review — if decode plus merge is a small share
> of the run, the owner decides how much of B–D is still worth building before anything is
> coded.**

### Milestone B — the two-state cached observation, in direct mode only

**B1. ✅ The item type** — *done 2026-09-04 (`d22b505b`), and **narrower than this step was
drafted**, which matters for C.*

What this asked for was the cache's item becoming spec §3.1's two-state observation, with
direct mode constructing the already-built state everywhere. **What landed is the half of that
which can be proven inert:** `LocusSummary` — the region, the kind's discriminant, and the keep
rule's two counts — with the closing walk routed through it, so no decision reaches past the
summary any more. **The cache still holds a built record for every observation and still builds
every one.** The second state, where the summary comes from a psp record's head and the
evidence is built only where a locus survived, is C1's work.

Splitting it there was not planned and is worth stating: a cached item with one state is not a
type worth having, and giving it a second state that direct mode never constructs would have
been dead code carried through a checkpoint. Routing the *decision* through the summary is the
property C depends on, and it is testable today. *Depends:* —. *Source:* spec §3.1.

**B2. ✅ Prove it inert. Own commit, do not bundle.** Direct mode's VCF byte-identical either
side of B1 on the run fixtures and the tomato slice; the merge's own suite (partition
invariance, k = 1, the failed-locus set) green unchanged. *Depends:* B1. *Source:* spec goal 2.

> **Checkpoint B: the seam exists and direct mode never noticed. Pause for review.**

### Milestone C — the two-phase psp source

**C1. ☐ The cache's second state, and the summary cursor.** Two things, because B1 left the
first of them: the cached item gains the state spec §3.1 describes — a summary that came from a
record's head, with the evidence not yet built — and the psp source walks heads to produce it,
retaining each record's raw bytes; cover advances it, eviction drains it (spec §3.2, §3.3).
Storage per §3.4's leaning (per-record boxes), eviction verified to return memory.
**A summary carries no locus kind** — the width verdict reads one per *locus*, from the opening
observation, and the per-member read is only a release assertion segments already guarantee
(spec §2). What this step must settle is therefore narrower than an earlier draft of it said:
where that one per-locus kind comes from when no body has been decoded. The ruling of
2026-09-04 supplies the route — look it up from the coordinate against the run's segmentation —
and this step builds that lookup. *Depends:* B1. *Source:* spec §2, §3.1, §3.2, §3.3, §3.4.

**C2. ☐ The build cursor.** Evidence on demand: monotonic per-sample builds through the
retained window, chain-id changes replayed, bodies built by the existing bounded decode; the
two refusals (backwards ask, evicted ground) provoked by tests; `ObservationBodyNotBuilt`
retired with its job handed over (spec §3.2). *Depends:* C1. *Source:* spec §3.2, §7.

**C3. ☐ The oracles. Own commit, do not bundle.** The decode-everything source moves to
test-support and the two-phase source is compared against it record for record over fixtures
and a real psp — including a wide-deletion fixture whose record reaches past its region, and a
tract record. Then the standing oracle: F2's mode-equivalence rerun green. *Depends:* C2.
*Source:* spec goal 1.

> **Checkpoint C: the psp path decides from summaries and builds one locus in a hundred, and
> no byte moved. Pause for review.**

### Milestone D — the run-level companions

**D1. ☐ One contig list for the run.** `PspVariantCaller::open` checks the headers' lists
agree and hands every reader the shared one; the open-cost probe
([`examples/ng_psp_open_cost.rs`](../../../../examples/ng_psp_open_cost.rs)) re-measured —
the spec's expectation is 480 kB → 123 kB per open sample on a human reference. *Depends:* —
(parallel to B/C). *Source:* spec §4; [`run_streaming.md`](../spec/run_streaming.md) §10.

**D2. ☐ The spare offer, measured.** The source builds evidence into offered records where
their buffers fit. Adopt on the allocator's share of the A1-versus-E1 profile; revert if it
does not pay, and record the number either way. *Depends:* C2. *Source:* spec §4.

> **Checkpoint D: the run's per-sample residents are the reader's own. Pause for review.**

### Milestone E — the after, and the scaling verdict

**E1. ☐ A1 rerun on the finished path.** Same corpora, same protocol; report before/after
wall, peak resident, and the stage shares; rerun
[`examples/ng_psp_skip_value.rs`](../../../../examples/ng_psp_skip_value.rs) at depth so the
skip's survival at 300× is a number beside the 3× one. *Depends:* C3, D1. *Source:* spec §5, §6.

**E2. ☐ The gated decision, written down.** From E1's shares: what per-sample cover
parallelism could recover of a whole run, and whether it gets a plan. The answer lands in
[`run_streaming.md`](../spec/run_streaming.md) §11 question 7 — *"not worth building"* closes
it just as well — and the psp half of
[`../research/cohort_merge_parallel_cost_plan.md`](../research/cohort_merge_parallel_cost_plan.md)
is folded in or closed by the same entry. *Depends:* E1. *Source:* spec §5, §9.

> **Checkpoint E: the trick is in, sized at both ends of both axes, and the structural
> question has a written answer. Pause for review.**

## Verification summary

| milestone | proven by |
|---|---|
| A — baseline | the numbers, with spread, at 1/16/63 samples; reviewed before any code |
| B — the seam | direct-mode VCF byte-identical across B1; merge suite green unchanged |
| C — the source | record-for-record equality with the decode-everything oracle; F2 mode equivalence green; both refusals provoked; eviction returns memory |
| D — companions | open-cost re-measured (≈480→123 kB); spare adopted or reverted on a profile number |
| E — the after | before/after wall, resident and shares on the same corpora; §11 q7 answered in writing |

## Out of scope (next plans)

- **Structural scaling of the psp path** — its own plan if E2 says so.
- **Reference bases in the record** ([`run_streaming.md`](../spec/run_streaming.md) §11 q4) —
  time it when builds are rare, which is after this plan.
- **Callers-in-flight default** ([`run_streaming.md`](../spec/run_streaming.md) §11 q2) —
  unchanged by this work.
