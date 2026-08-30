# ng — the VCF module: implementation plan

**Status:** draft, 2026-08-30. The build order for **`src/ng/vcf/`**: the one file a run writes,
for SNPs, indels and repeat tracts together. Design is settled in
[`../spec/vcf_output.md`](../spec/vcf_output.md). This turns that design into build order; it is
**not** a place for new design — the spec's open questions are its §14, and the one that could
have reopened the shape (one interleaved file) is already resolved there.

**Where this sits on the critical path.** The run driver needs a consumer for its ordered record
stream ([`../spec/run_streaming.md`](../spec/run_streaming.md) §10), and "ng emits a VCF" is the
first end-to-end milestone that means anything to a user. This module is that consumer. It does
**not** wait for steps 11a/11b (emission policy): the writer writes what it is handed, and the
policy stages slot in upstream of it later.

---

## Scope

**In:** a new `src/ng/vcf/` module — the written record's own type and the worker-side assembly
that fills it while the evidence is still in hand; the header; the record encoder (fixed columns,
INFO, FORMAT, the anchor rule); the ordered writer with its one legal POS tie; the atomic sink;
and the verification ladder up to a differential against production's generic file on a real
cohort.

**Out (later plans, or already owned):**

- **What is dropped and what is tagged** — steps 11a/11b of
  [`../spec/ng_proposal.md`](../spec/ng_proposal.md). This module encodes any record it is given,
  including a filtered tract locus in its §8 shape; deciding which records exist is not its job.
- **The run wiring** — which stage feeds the writer, where the artifact-correction stage sits on
  the stream — [`run_driver_direct_mode.md`](run_driver_direct_mode.md) and
  [`../spec/calling_quality.md`](../spec/calling_quality.md) §3.4.
- **Real repeat-tract records end to end.** The tract path is not through the merge yet — the
  motif reaches `CohortObservation` in [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md)'s
  Milestone A, and the tract quality (`calling_quality_ssr.md`) is unwritten. This plan builds and
  tests the tract *encoding* (`STR`/`RU`/`PERIOD`, `REPCN`, the tract anchor cases) on constructed
  records, so that path has nothing left to build here when its blockers clear.
- **gVCF, `GP`, `PARALOG_POST`, a tract id in the ID column** — deferred in the spec (§13), each
  with a home.

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **Encoding before I/O, on golden files.** A record's bytes are a pure function of the record
  (spec §11 — byte-determinism is part of the format), so milestones A–B need no file, no writer
  and no run: fixture record in, asserted bytes out. Everything downstream trusts this.
- **The strict external parser is a gate, not a finale.** The tract writer's conventions were
  never pushed through an external parser in anger (spec §15 test 1). `bcftools` is in the dev
  container (`Containerfile:37`); it gates Checkpoint C, before any real-data work.
- **Split data-shaping from encoding.** Assembling a record from a called locus (worker-side,
  evidence in hand) and rendering a record to bytes (anywhere, evidence long gone) are two
  interfaces; neither reaches into the other. Same split, same reason as the quality module's
  nine-count summary ([`../spec/calling_quality.md`](../spec/calling_quality.md) §3.3).
- **Reuse over rewrite.** Production's header builder, ordered writer and atomic sink are ported,
  not redesigned; the spec's reuse map (§12) names what moves and what is left behind.
- **Isolate the silent failures.** Three places produce a quietly-wrong file rather than a crash:
  the anchor rule (a wrong POS parses fine), the `REPCN`-to-`GT` alignment (production's own
  mismatch shipped), and float formatting (a determinism leak is invisible until two runs
  differ). Each is its own commit with its oracle stated.
- **Incremental, with pauses.** One milestone, then stop for review.
- **Container builds.** All `cargo` through `./scripts/dev.sh` by absolute path (`CLAUDE.md`).

## Preconditions (already in place)

- **The called locus** — [`LocusInference`](../../../../src/ng/calling/mod.rs) (region, alleles,
  per-sample calls in the run's sample order, expected allele copies, `converged`, provenance,
  the tract provenance) and [`SampleGenotypeCall`](../../../../src/ng/calling/mod.rs) with its
  `Missing` variant — the no-call the spec's §7 writes as `./.`.
- **The quality outputs** — the corrected site quality and
  [`ArtifactPenalties`](../../../../src/ng/calling/quality/artifact_correction.rs) (`ABPEN`,
  `SPPEN`), plus the one-quality rule the writer must not break
  ([`../spec/calling_quality.md`](../spec/calling_quality.md) §3.5).
- **The evidence, while it lasts** —
  [`CohortObservation` / `SampleSupport`](../../../../src/ng/run/cohort_merge/build.rs): per
  `(allele, read group)` read counts and mapq sums, partial observations, silent-read tallies.
  Everything per-sample the file carries is summed out of this before the locus is released.
- **The genotype-order table** —
  [`genotype_table.rs`](../../../../src/ng/calling/genotype_table.rs), the decode `GT` runs
  through at any ploidy.
- **`ReferenceInfo`** — contig names, lengths, optional md5s
  ([`../spec/reference_info.md`](../spec/reference_info.md)) for the `##contig` lines.
- **Production's writers as frozen reference** — [`src/vcf/`](../../../../src/vcf/) and
  [`src/ssr/cohort/vcf_out.rs`](../../../../src/ssr/cohort/vcf_out.rs) — both for ported code and
  as the differential's other arm.
- **`noodles` as a dependency** (production's header path uses it), and **`bcftools` in the dev
  container**.

**One join this plan builds that no type carries today, named up front:** nothing currently
survives the locus's release carrying per-sample depths, per-allele read counts or pooled mapping
qualities — the calls do not have them and the evidence does not outlive the worker. Milestone A's
record type is that carrier, and Milestone D fills it. This is the same pattern, for the same
reason, as `ArtifactTestCounts`: sum everything cohort-shaped in the worker, carry scalars and
short vectors downstream.

---

## The steps

### Milestone A — what a record carries (types, no I/O)

✅ **A1. The written record's type.** Everything the file needs after the evidence is gone, per
locus: the region; the allele sequences; per sample, the call (or `Missing`), the genotype
quality, the depth tally and per-allele counts of A2; the pooled mapq sums for `MQREF`/`MQALT`/
`MQDIFF`; the corrected site quality and the two penalties; `converged`; and the tract
annotations (`RU`, and `REPCN` derivable per called allele) as an `Option` that is `Some` exactly
when the record is a tract record. Names are the coder's; the content list is the spec's (§5–§8)
and dropping an item is a stop-and-ask.
*Depends:* —. *Source:* spec §3, §5–§8.

✅ **A2. The per-sample tally, pinned.** The one place `DP` and `AD`'s exact composition is
decided, in types and doc comments before any number is summed: `AD[allele]` = complete
observations whose sequence matched that allele, summed over read groups; `DP` = the sample's
complete observations over the merge's **whole** allele table — written alleles or not — plus its
partial observations, so `DP − ΣAD` is every read no written allele explains (spec §7). Record in
the step's report what was included and what was left out (the silent-read tally is *out*: those
reads produced no observation at all, and its own doc says it double-counts across records).
*Depends:* A1. *Source:* spec §7; [`build.rs`](../../../../src/ng/run/cohort_merge/build.rs).

✅ **A3. The header's metadata type and its refusals.** Contigs (name, length, md5 when known),
the sample names in the run's sample order, source, command line, reference path, the parameters
file's name. Ported refusals: duplicate sample names, duplicate contigs, a contig longer than
`i32::MAX` (spec §4; production [`header.rs:296-321`](../../../../src/vcf/header.rs)).
*Depends:* —. *Source:* spec §4.

> **Checkpoint A:** the types can express both record kinds, the no-call, and a filtered tract
> locus; a record the spec forbids (a tract annotation without its sibling, an `AD` wider than
> the alleles) is unrepresentable or refused at construction. Pause for review.

### Milestone B — one record to bytes

✅ **B1. Fixed columns and the anchor rule. Own commit, do not bundle.** CHROM through FILTER for
a record with no empty allele; then the anchor: any empty allele ⇒ every allele prefixed with the
left reference base and POS shifted one left, and at contig position 1 the right base appended
with no shift — the rule that replaces production's invented `N`
(spec §5). **Silent-failure oracle:** the position-1 cases and the shifted-POS cases asserted
byte-for-byte, including the tract-deletion fixture mirroring production's `N` case
([`vcf_out.rs:405-435`](../../../../src/ssr/cohort/vcf_out.rs)) with the corrected output.
*Depends:* A1. *Source:* spec §5.

☐ **B2. INFO.** `AF` from the expected allele copies normalised over `AN` — **and the step
checks, against the loop, that this equals the converged pass's fitted frequency**; if it does
not, stop for a ruling rather than write either. `AC`/`AN` decoded from the calls with `Missing`
samples excluded and the `AC ≤ AN` assertion kept; `DP`; `ABPEN`/`SPPEN`; the MQ family with
production's omission rules (absent key vs `.` entry, spec §11); `STR`+`RU`+`PERIOD` together or
not at all.
*Depends:* B1. *Source:* spec §6; [`record_encode.rs:280-474`](../../../../src/vcf/record_encode.rs).

☐ **B3. FORMAT and the sample columns.** `GT` through the genotype-order table, sorted,
`/`-joined; `GQ` rounded and capped; `DP`/`AD` from A2's tally; the no-call spelling — `./.` with
`GQ` missing but `DP`/`AD` written when the sample had reads (spec §7, and its Q3 note). **Own
commit inside this step: `REPCN` in `GT`'s sorted order** — production computes it in candidate
order while sorting `GT`, so the two fields need not correspond; the fixture that catches it is a
het whose sorted order differs from candidate order.
*Depends:* B2. *Source:* spec §7; [`vcf_out.rs:372-505`](../../../../src/ssr/cohort/vcf_out.rs).

☐ **B4. Number formatting, pinned as a table. Own commit, do not bundle.** QUAL to one decimal;
the stated precision for `AF`, the penalties and the MQ family; integer fields as integers. One
formatting function per type, one test table of adversarial values — the spec makes the rendering
part of the format (§11), so a formatting choice is a format change from here on.
*Depends:* B1–B3. *Source:* spec §11.

> **Checkpoint B:** golden-file tests cover both record kinds, multi-allelic, all four anchor
> cases, the no-call, and a filtered tract locus (`QUAL 0`, `ALT .`, every sample `./.`); every
> golden line also round-trips through noodles' reader. Pause for review.

### Milestone C — the header, the writer, and a strict parser

☐ **C1. The header.** A3's metadata rendered through the ported noodles path, declarations
exactly the spec's §4/§6–§8 set, `##parametersFile` included. Golden-file the whole header.
*Depends:* A3, B4. *Source:* spec §4; [`header.rs:73-247`](../../../../src/vcf/header.rs).

☐ **C2. The ordered writer and the sink.** Port the ordered writer and the atomic sink
(tmp → flush → bgzf EOF → fsync file and directory → rename; bgzf by suffix). **The one change:**
strictly-increasing POS relaxes to exactly the spec's legal tie — same POS admitted once, generic
before tract, anything else refused ([`../spec/vcf_output.md`](../spec/vcf_output.md) §5).
*Depends:* C1. *Source:* spec §5, §11; [`writer.rs`](../../../../src/vcf/writer.rs),
[`sink.rs`](../../../../src/vcf/sink.rs).

☐ **C3. The strict-parser gate.** A fixture file holding every golden case interleaved — generic
and tract records, the POS tie, the filtered locus — parses under `bcftools view` in the dev
container with **zero warnings**, plain and bgzf both. This is spec §15 test 1 and it runs here,
before any real data, because an encoding defect found later costs a re-run of everything after
this point.
*Depends:* C2. *Source:* spec §15 test 1.

> **Checkpoint C:** ng can write a well-formed interleaved VCF from constructed records, and an
> external parser agrees. Pause for review.

### Milestone D — a record from a called locus

☐ **D1. The worker-side assembly.** One function: `LocusInference` + the locus's evidence
(+ `ReferenceInfo` for anchor bases) → A1's record, called where the quality module's summary is
built — the last moment both are in hand. Fills A2's tallies and the MQ pools; copies the calls;
derives nothing the loop already decided.
*Depends:* Checkpoint B, A2. *Source:* spec §3;
[`../spec/calling_quality.md`](../spec/calling_quality.md) §3.3.

☐ **D2. The quality join, and the invariant.** The record's QUAL is the corrected value the
artifact stage wrote, the penalties ride beside it, and **the writer contains no arithmetic on
any of the three** — assert it structurally (the encoder takes them as opaque finished values).
This is the §3.5 one-quality rule; the shipped production defect it prevents is the spec's
whole reason for it.
*Depends:* D1. *Source:* spec §5;
[`../spec/calling_quality.md`](../spec/calling_quality.md) §3.5.

☐ **D3. Assembly fixtures where the joins can lie.** A locus with a `Missing` sample (AN drops,
`./.` written, tally still filled); a locus where candidate selection dropped a supported allele
(`DP − ΣAD` lands in every carrying sample, spec §15 test 6); a sample with two read groups
(tallies summed across groups exactly once); an unconverged locus (`EMNoConv`, still written).
*Depends:* D1, D2. *Source:* spec §7, §8, §15 tests 5–6.

> **Checkpoint D:** called fixtures become correct records with no hand-filled fields left.
> Pause for review.

### Milestone E — real data, and the differentials

☐ **E1. The generic differential against production.** Drive merge → loop → assembly → writer
over real cohort slices (tomato, and HG002) with production's parameters, and compare against
production's generic VCF on the same input: `GT`, `AC`, `AN` agree; `AF` compared with the E-step
caveat B2 recorded; `DP`/`AD` divergences all in the dropped-allele class and counted; anchor
divergences all in the position-1 class. **Every divergence lands in a named class or the run is
red** — the census discipline, not a tolerance.
*Depends:* Checkpoint D. *Source:* spec §15 test 3.

☐ **E2. The file-invariant scanner.** A small check over any produced file, kept as a test
utility: no two records share a reference base; the only POS ties are generic-then-tract; `STR`
⟺ `RU` and `PERIOD`; `AC ≤ AN` per record; diploid `AN + 2·(no-calls) = 2·samples`. Runs over
E1's outputs and every future run's.
*Depends:* E1. *Source:* spec §15 tests 2, 5.

☐ **E3. Determinism, measured not asserted.** The same slice at 1 and at several workers, files
byte-identical — encoding is a pure function plus C2's ordering, so a difference is a bug in one
of exactly those two places. (If the run driver cannot yet vary worker count over this path, the
test is the formatter's purity plus C2's order refusals, and the whole-run form moves to the
driver's plan — say which one shipped.)
*Depends:* E1. *Source:* spec §11, §15 test 4.

☐ **E4. The concordance tooling over an interleaved file** — spec §14 Q1's verification:
`benchmarks/lib/ssr_concordance.py` and one external dumpSTR-class tool over a fixture file with
both kinds interleaved (constructed tract records are fine — the question is about the file
shape, not the caller). A rigid tool gets a documented one-line `bcftools` pre-split, not a
format change.
*Depends:* C3. *Source:* spec §14 Q1.

> **Checkpoint E:** ng's generic VCF stands next to production's with every difference named,
> and the interleaved shape is confirmed against real tooling. Pause for review — and this is
> the plan's end state until the tract path's blockers clear elsewhere.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | the spec-forbidden states are unrepresentable or refused; `DP`/`AD` composition is written down before any sum exists |
| B | golden bytes for every record shape the spec names, round-tripped through noodles; formatting pinned by table |
| C | `bcftools view` accepts the interleaved fixture file with zero warnings, plain and bgzf |
| D | fixtures where each join could silently lie (missing sample, dropped allele, two read groups, unconverged) produce the spec's exact spelling |
| E | the production differential with every divergence in a named class; the file-invariant scanner; byte-identity across worker counts; the concordance tooling runs |

## Out of scope (next plans)

- **Emission policy** — thresholds, hard-drop vs FILTER-tag: steps 11a/11b's documents, unwritten.
- **The stream wiring** — [`run_driver_direct_mode.md`](run_driver_direct_mode.md) hands this
  module its records; nothing here schedules anything.
- **Real tract records** — gated on [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md)
  Milestone A (the motif into the merge) and on `calling_quality_ssr.md` (unwritten). The
  encoding will already be green from B/C.
- **The spec's §13 deferrals** — gVCF, `GP`, `PARALOG_POST`, a tract id in ID.
