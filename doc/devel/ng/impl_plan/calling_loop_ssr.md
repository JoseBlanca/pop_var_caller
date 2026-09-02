# ng STR calling loop — implementation plan

**Status:** draft, 2026-09-02. Turns the settled design of
[`../spec/calling_loop_ssr.md`](../spec/calling_loop_ssr.md) into build order: tract
candidate selection (whose own design is
[`../spec/candidate_alleles_ssr.md`](../spec/candidate_alleles_ssr.md) and whose steps this
plan executes from [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md)), the driver's
dispatch, the record's tract inputs, the QUAL experiment, and allele discovery with its
measurement. **Not a place for new design.**

**Runs in parallel with** [`run_ssr_observations.md`](run_ssr_observations.md), on its own
branch and worktree (§Branch). That plan owns the two seams this one consumes: the
`CohortObservation.kind` field (its Milestone A — on `main` before this plan's first step)
and the driver's temporary set-aside guard (its Milestone C — replaced by this plan's
Milestone C).

---

## Scope

**In:** `calling/allele_candidates/ssr.rs` (new — `select_ssr`, the ladder) with the scratch
additions in `allele_candidates/mod.rs`; the driver dispatch in `run/callers.rs` and the
record wiring in `run/records.rs`; the run report's tract-outcome partition;
`calling/inference/discovery.rs` (the mechanism the bake-offs plan designed, built here per
its ownership note); the QUAL-experiment harness and its report; the probes and benchmark
runs that prove each of those.

**Out (with owners):**

- **Routing, the walk's slot, the merge kind, the report's ground wording** —
  [`run_ssr_observations.md`](run_ssr_observations.md).
- **The tract QUAL *decision* and its spec** — `calling_quality_ssr.md`, written from
  Milestone D's report (spec §3.3); this plan produces the report, never the choice.
- **The per-locus slippage re-fit** — stays with
  [`calling_bakeoffs.md`](calling_bakeoffs.md) Milestone D (spec §3.5).
- **Bundles** — no design yet (spec §1.2).
- **The fit-mode command** — deferred future work (owner, 2026-09-02; spec §3.4).

## Principles (how the order was chosen)

- **The algorithmic heart before the plumbing.** Selection — the one unbuilt statistic — is
  built and proven on hand-built loci and against its two-ended differential before the
  driver learns to dispatch to it.
- **Reuse over rewrite, by reference.** The selection steps and the discovery mechanism are
  already planned in reviewed documents; this plan *executes* those milestones rather than
  restating them, so there is one source of truth for each step and one set of checkboxes.
- **Verify against ground truth.** Selection's oracle is two-ended (production reproduced
  with its rules switched in; the measured HG002 improvement with them switched out); the
  end-to-end oracle is the GIAB tract ground where today's recall is exactly 0.000, with
  production's and freebayes' numbers as the stated bar.
- **A report before a decision.** The QUAL experiment and the discovery measurement each end
  in a written report; the decisions they feed (`calling_quality_ssr.md`, discovery's
  default) are taken from the report, not in the code that produced it.
- **Container builds**: all `cargo` via `./scripts/dev.sh`; a native host build at completion.

## Preconditions (already in place, or owed by the parallel plan)

- The shared selection module, merged: `select_generic`, `SelectionVerdict`,
  `UnmatchedSupport`, `AlleleRemap`, `SelectionScratch`
  ([`allele_candidates/`](../../../../src/ng/calling/allele_candidates/)).
- The tract scoring chain, built and tested end to end given candidates:
  `shape_ssr_locus`, `StutterSubstitutionEmission`, `repeat_tract_parameters`, the loop's
  tract arms ([`summarise_condition.rs:7623`](../../../../src/ng/calling/inference/summarise_condition.rs)),
  the encoder's `STR`/`RU`/`PERIOD`/`REPCN` and `FilterVerdict` (spec §2's table).
- **From the parallel plan, before A1:** its Checkpoint A (`CohortObservation.kind`) merged
  to `main`. **Before C1:** its Checkpoint C (slot filled, set-aside guard) merged.
- The baseline and bar: recall 0.000 on repeat-routed GIAB ground; production 0.855–0.990,
  freebayes 0.818–0.874
  ([the loss report](../../reports/ng_str_path_losses_2026-09-02.md) §5).

## Branch, worktree, and the parallel plan

- **Branch** `ng-ssr-calling-loop`, from `main` **after the observations plan's
  Checkpoint A has merged**. **Worktree** `../pop_var_caller-ssr-calling-loop`
  (`git worktree add ../pop_var_caller-ssr-calling-loop ng-ssr-calling-loop`). Use absolute
  paths in commands; remove the worktree when the plan completes.
- The selection milestones executed from
  [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) run **on this branch** — that
  plan's own branch note (`ng-candidate-alleles-ssr`) is superseded by its ownership notes
  and this section.
- **File ownership** is the table in
  [`run_ssr_observations.md`](run_ssr_observations.md) §Branch: this plan does not touch
  `cohort_merge/`, `walker.rs`, `report.rs`'s ground wording, or
  `call_from_alignments.rs`; changes cross only through `main`. `run/callers.rs` and
  `run/records.rs` become this plan's after the observations plan's Checkpoint C merges;
  rebase on `main` at that point, before C1.

---

## The steps

### Milestone A — selection, on hand-built loci

**A1. Execute [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) Milestones B, C and D.**  ✅
The ladder, nomination, and admission/periodicity/outputs, exactly as that plan writes them —
its checkboxes are the live record; this box flips when its Checkpoint D passes. Its
Milestone A is already on `main` (the parallel plan's Checkpoint A — verify before starting).
*Depends:* the parallel plan's Checkpoint A.
*Source:* [`../spec/candidate_alleles_ssr.md`](../spec/candidate_alleles_ssr.md);
[`../spec/calling_loop_ssr.md`](../spec/calling_loop_ssr.md) §3.1.

> **Checkpoint A:** `select_ssr` complete and proven on hand-built loci (that plan's
> Checkpoint D). Pause for review.

### Milestone B — selection's differential, on real data

**B1. Execute [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) Milestone E.**  ✅
Both ends: production's candidate set reproduced on the tomato fixture with its three
replaced rules switched in; the measured HG002 numbers reproduced through the shipped module
with them switched out. The probe's tract slot is filled inside the probe only — the run's
slot belongs to the parallel plan.
*Depends:* A1. *Source:* [`../spec/candidate_alleles_ssr.md`](../spec/candidate_alleles_ssr.md) §10, §13.

> **Checkpoint B:** the differential green at both ends. Merge to `main`. Pause for review.

### Milestone C — the driver calls a tract, and the record says so

**C1. The dispatch.**  ✅
`call_one_generic_locus`'s call site branches on `CohortObservation::kind`
([`callers.rs:813`](../../../../src/ng/run/callers.rs)): `Generic` unchanged; `Ssr` runs
`select_ssr` → `shape_ssr_locus` → the same `genotyper.call_locus`, replacing the
observations plan's set-aside guard. `SsrBundle` stays set aside and counted. The scratch is
already the STR emission's ([`callers.rs:551`](../../../../src/ng/run/callers.rs)).
*Depends:* B1; the parallel plan's Checkpoint C merged, rebase first.
*Source:* spec §3.2.

**C2. The record's two unwired inputs.**  ☐
`records.rs` passes `Some(TractAnnotation::new(motif))` from the observation's kind
(today `None`, [`records.rs:263`](../../../../src/ng/run/records.rs)), the per-allele repeat
counts reach `REPCN`, and the FILTER verdict comes from selection and the loop
(`NotPeriodic`, `TooManyAlleles`, `LowDepth`, `EmDidNotConverge`) rather than from
convergence alone. Round-trip pin: a tract record through `bcftools` unchanged.
*Depends:* C1. *Source:* spec §3.2; [`../spec/vcf_output.md`](../spec/vcf_output.md) §6–§8.

**C3. The run report partitions tract outcomes.**  ☐
Called, refused-by-FILTER, and set-aside-as-unbuilt (bundles) are three counted lines, the
way the generic path's outcomes already partition. Wording only where this plan owns it — the
ground shares stayed with the parallel plan.
*Depends:* C1. *Source:* spec §3.2.

**C4. End to end, measured against the zeros.**  ☐
A GIAB run at 30× and 50× writes tract records. Measure and record: recall on repeat-routed
ground (baseline 0.000; the bar: production 0.855/0.909 indels, 0.990 SNPs; freebayes
0.818–0.874), genotype concordance on the dashboard's panel, and E2 byte-identity at 1–16
threads with tract records in the file. Update the benchmark dashboards' ng rows.
*Depends:* C1–C3. *Source:* spec §1.1, §8.

> **Checkpoint C:** ng calls repeat tracts through its own loop; the numbers are recorded
> beside the stated bar. Merge to `main`. Pause for review.

### Milestone D — the QUAL experiment (a report, not a decision)

**D1. The harness.**  ☐
Per arm — A the inherited fold, C production's caller as external comparator on the same
ground (arm B exists only if A fails, designed later in `calling_quality_ssr.md`) —
calibration (records binned by QUAL against the share truly variant) and the
precision/recall sweep over a QUAL threshold, split by period (homopolymer vs 2+) and by
parameter provenance (fitted vs `Defaulted`). GIAB tract ground at 30×/50×; the STR
simulator at settable slippage.
*Depends:* C4. *Source:* spec §3.3.

**D2. The runs and the written report.**  ☐
Run the harness; write the report with the decision rule's inputs filled in (does arm A
reach the corrected SNP QUAL's standard on the same benchmark?). **Done is the report** —
the decision and any arm-B design belong to `calling_quality_ssr.md`.
*Depends:* D1. *Source:* spec §3.3 (the obligation, owner 2026-09-02).

> **Checkpoint D:** the QUAL report exists with calibration and sweep numbers per arm and
> split. Pause — the owner reads this one before `calling_quality_ssr.md` is written.

### Milestone E — allele discovery, built and measured

**E1. The mechanism.**  ☐
Execute [`calling_bakeoffs.md`](calling_bakeoffs.md) E1 as written there (its ownership note
hands the build here): retrace what the model explains as slippage after convergence, admit
lengths clearing both halves of the bar, append emission columns, recompute the tract's
genotype likelihoods (the STR carry-over rule), stop on a round that adds nothing or at the
cap, prune, re-run.
*Depends:* C1. *Source:* [`../spec/calling_em_loop.md`](../spec/calling_em_loop.md) §4.1;
spec §3.5; [`calling_bakeoffs.md`](calling_bakeoffs.md) E1.

**E2. The three discovery pins.**  ☐
Execute [`calling_bakeoffs.md`](calling_bakeoffs.md) E2: plant/terminate/free-when-off, and
append-only columns — **own commit**; a leaky `Off` quietly changes the default and a
rebuilt table is quietly slow.
*Depends:* E1. *Source:* [`calling_bakeoffs.md`](calling_bakeoffs.md) E2;
[`../spec/calling_em_loop.md`](../spec/calling_em_loop.md) §13 tests 11–12.

**E3. The measurement, and the default it sets.**  ☐
The owner's two questions, as numbers (the F3 report shape): fire rate in loci per ten
thousand first; on firing loci, alleles admitted, surviving the prune, rounds, cost. Then
off-against-on on the same ground and truth: GIAB tract recall and genotype concordance at
30×/50×, simulator exact-truth accuracy, and tomato at three reads with the evidence-bar
sweep. **The shipped default is set to what the report says**, and the report records it.
*Depends:* E2, C4. *Source:* spec §3.5; [`calling_bakeoffs.md`](calling_bakeoffs.md) F3.

> **Checkpoint E:** discovery built, pinned, measured; the default matches the report.
> Merge to `main`. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md)'s own oracles, through its Checkpoint D |
| B | **external, both ends**: production reproduced rules-in on tomato; HG002 numbers reproduced rules-out |
| C | recall moves off 0.000 on GIAB repeat ground, recorded beside production's and freebayes' numbers; `bcftools` round-trip; E2 byte-identity with tract records |
| D | the report itself: calibration + sweep per arm, split by period and provenance — no decision taken |
| E | plant/terminate/free-when-off and append-only pins; the off-vs-on report; default = report |

## Out of scope (next plans)

- **`calling_quality_ssr.md`** — written from Milestone D's report; owns the QUAL decision
  and any arm-B design.
- **The per-locus slippage re-fit** — [`calling_bakeoffs.md`](calling_bakeoffs.md)
  Milestone D and report F2.
- **Bundles** — a spec before a plan (spec §1.2).
- **The fit-mode command** — deferred future work (owner, 2026-09-02).
- **The routing-frontier measurement** — cheap once this plan's Checkpoint C lands;
  [`../spec/typed_regions.md`](../spec/typed_regions.md) §5.2's question.
