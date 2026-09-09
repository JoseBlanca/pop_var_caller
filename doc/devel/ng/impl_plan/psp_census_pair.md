# ng — the psp and its census are a pair: implementation plan

**Status:** plan, 2026-09-09. **No code yet.** It turns the settled design in
[`psp_census_pair.md`](../spec/psp_census_pair.md) into build order and is **not a place for new
design**; where a step meets a question the spec does not answer, it stops at a checkpoint rather
than deciding. It follows [`parameter_prepass_runs.md`](parameter_prepass_runs.md), which built the
four commands this plan reshapes into three and a repair.

---

## 1. What this closes

Today a user names census files at step 2, retypes the walk's repeat criteria there, and can put a
census anywhere. When a census is stale, step 2 reports the first one it meets and stops. After this
plan:

    generate-psps        alignments        ->  <sample>.psp and <sample>.census
    estimate-parameters  psps              ->  cohort.parameters.toml, or a refusal naming every stale pair
    call-from-psps       psps + that file  ->  the VCF
    regenerate-census    psps              ->  <sample>.census beside each, for the pairs step 2 refused

No command takes a census path. No command takes a repeat criterion the psp header already holds.

---

## 2. Scope

**In:**

- one judgement of whether a psp's census is fresh, made from the two files' headers, shared by
  step 2 and the repair command (spec §4.2, §8);
- `estimate-parameters` over `--psp`, with the six duplicated flags removed and every stale pair
  reported before the reference is read (spec §4, §6);
- `regenerate-census` in place of `generate-census`: beside the psp, no `--output-dir`, no
  `--force`, fresh pairs skipped unless `--all` (spec §8);
- an older-format census named as such rather than as damage (spec §4.2);
- a segmentation entry point that takes the criteria the header holds (spec §6, the trap);
- everything that spells the old shape: subcommand docs, scripts, the pipeline line in
  `PROJECT_STATUS.md`;
- the two measuring probes on branch `census-vs-psp-perf` landed under `examples/`, since the spec
  cites their numbers.

**Out:**

- **a record count in the psp footer** — [`psp_file_format.md`](../spec/psp_file_format.md) §3.3;
  the pair check stays header-only, as spec §5 records.
- **a fit-side budget or cap** — [`parameter_prepass_runs.md`](parameter_prepass_runs.md)
  Checkpoint C's knobs question; spec §7 records why it would be cheap.
- **the census's per-position rule under a widened record** — spec §12, owned by
  [`parameter_prepass_census_sites.md`](../spec/parameter_prepass_census_sites.md) §2.
- **the refusal's time estimate** — spec §13's one open question; confirm at Checkpoint B before
  adding the sentence.

---

## 3. Principles (how the order was chosen)

- **The judgement before the commands.** Whether a pair is fresh is a pure function of two headers
  and a version word. It is built and tested on fixtures first, then both commands call it; a
  command that made its own judgement would be the two-spellings-of-one-test the spec's reuse map
  exists to prevent.
- **Isolate the step whose failure is silent.** Reading the criteria from the header instead of the
  flags (B1) changes where a number comes from and nothing about the arithmetic. If it is wrong the
  fit still produces a file — a plausible one. It lands as its own commit, with the fitted file
  byte-identical to the flag-driven one on the tomato fixture before and after.
- **Reuse over rewrite.** `freshness_by_header`, `header_and_its_digest`, `census_from_psp`, the
  `.partial`-then-rename, and today's directory expansion are called as they are. The plan names
  what it reuses and re-derives nothing.
- **Verify against ground truth.** The oracle for the rename is the existing byte-for-byte test
  between the walk's census and the rebuilt one; the oracle for the whole is
  `scripts/ng_fit_stage_end_to_end.sh` producing the same parameters file and the same VCF as
  before.
- **Types first, then implementation** (project rule).
- **Incremental, with pauses.** Four milestones, a checkpoint after each.
- **Builds through `./scripts/dev.sh` where a container runtime exists**, plain `cargo` where it
  does not (CLAUDE.md).

---

## 4. Preconditions (already in place)

Confirm each before step A1.

- `generate-psps` writes `<sample>.psp` and `<sample>.census` from one pass, census renamed into
  place first ([`generate_psps.rs:676-735`](../../../../src/pop_var_caller_exp/generate_psps.rs)).
- The identity a census carries and the header-only check of it: `PileupIdentity` and
  `freshness_by_header` ([`census_file.rs:83`, `:176`](../../../../src/ng/parameter_estimation/joint/census_file.rs)).
- The cheap header read: `psp::header_and_its_digest` ([`psp/mod.rs:192`](../../../../src/ng/psp/mod.rs)).
- The settings in the header: `SegmentationInputs` ([`segmentation_inputs.rs:25-40`](../../../../src/ng/segmentation_inputs.rs)).
- One census rebuilt from one psp: `census_from_psp` ([`census_from_psp.rs:202`](../../../../src/ng/run/census_from_psp.rs)).
- The cohort opener whose loop this plan replaces: `open_census_cohort`
  ([`census_cohort.rs:175-215`](../../../../src/ng/run/census_cohort.rs)).
- The fixtures both command test modules build on: `a_walked_cohort` in
  [`estimate_parameters/tests.rs:37`](../../../../src/pop_var_caller_exp/estimate_parameters/tests.rs)
  and [`generate_census/tests.rs:36`](../../../../src/pop_var_caller_exp/generate_census/tests.rs).
- The end-to-end oracle: [`scripts/ng_fit_stage_end_to_end.sh`](../../../../scripts/ng_fit_stage_end_to_end.sh),
  which runs the four commands on real reads and compares the two routes' censuses.

---

## 5. The steps

### Milestone A — one judgement of a pair

☐ **A1 — the verdict type.** An enum saying what a psp's census is: fresh; absent; of an older
format, naming the version found and the one this build reads; built from another psp, naming the
field `freshness_by_header` names. A noun, with its own `Display` in the words of spec §4.2. No
logic yet.
*Depends:* —. *Source:* spec §4.2.

☐ **A2 — the judgement.** One function: a psp path in, its verdict out. Derives the census path
from the stem (one helper, replacing both `psp_beside` and `census_path_for` — spec §9), reads the
census's version word before its magic-and-version check can call it malformed, reads the psp's
header digest, and calls `freshness_by_header`. Unit tests on fixtures for every verdict, including
a census file whose version word is `VERSION − 1` and one whose psp has been rewritten.
*Depends:* A1. *Source:* spec §4.2, §5, §9 (the *older format is `Malformed`* trap).

☐ **A3 — the cohort agrees on its settings.** Every psp's `SegmentationInputs` compared against the
first's with `first_difference`, which has no caller outside its own tests today; a disagreement is
refused naming the sample and the field. A test with two psps walked under different `min_copies`.
*Depends:* —. *Source:* spec §6, second paragraph.

☐ **A4 — a cohort judged whole.** Every psp of a list judged, no early return; the result is one
verdict a sample, in the order given, plus a directory expansion identical to today's `--psp`
handling in `generate-census`. A test with three stale pairs in five that asserts all three are
named.
*Depends:* A2, A3. *Source:* spec §4.1.

> **Checkpoint A: a pair can be judged, and a cohort's verdicts come back together.** Pause for
> review.

### Milestone B — `estimate-parameters` over psps

☐ **B1 — the criteria from the header, own commit.** Split `segments_over` after
`routing_criteria` so a caller holding a `StrRepeatCriteria` enters there
([`run_ground.rs:242-281`](../../../../src/pop_var_caller_exp/run_ground.rs)); `estimate-parameters`
builds its segmentation from the cohort's `SegmentationInputs` instead of from its flags, **with
the flags still present and ignored**. Oracle: the parameters file written on the tomato fixture is
byte-identical before and after. **Own commit, do not bundle** — a wrong criteria here fits a
plausible file.
*Depends:* —. *Source:* spec §6 and its trap.

☐ **B2 — `--psp` in, `--census` and the five criteria flags out.** The argument struct, the
directory expansion from A4, and the six deletions. `SUBCOMMAND` unchanged. Tests that build
argument lists (`args_over`, `a_shortest_run`) rewritten to the new shape; the byte-identity of B1
re-asserted on the new arguments.
*Depends:* A4, B1. *Source:* spec §6.

☐ **B3 — the refusal, before the reference.** A4's verdicts are taken first; if any sample is not
fresh, the run stops with spec §4.3's report — every stale sample, its census path, its cause, and
one `regenerate-census --psp …` line built from the arguments given — and nothing else has been
read. A3's agreement check runs first. `a_census_without_its_psp_is_refused` becomes *a psp without its census is refused and the
report names it*; a second test has two stale samples and asserts both are in the message and the
reference was never opened.
*Depends:* B2. *Source:* spec §4, §4.1, §4.3.

☐ **B4 — the fit's own refusal says what to run.** `CohortFitError::AnotherSelection`
([`census_fit.rs:139`](../../../../src/ng/run/census_fit.rs)) is rendered at the command as *these
censuses were built under another selection; run `regenerate-census --all`*, since this is the one
staleness the header cannot show and `--all` is what forces it.
*Depends:* B3. *Source:* spec §4.2's fourth row, §8.

> **Checkpoint B: step 2 takes psps, reads nothing it could be told wrongly, and refuses a stale
> cohort whole.** Decide spec §13's open question here — whether the refusal estimates the
> regeneration's duration. Pause for review.

### Milestone C — `regenerate-census`

☐ **C1 — the rename.** Module, `SUBCOMMAND`, `cli.rs` variant and its doc, the fifteen tests, and
the two `SUBCOMMAND`-spelling tests. A3's agreement check replaces the analysed-regions-only one. `--output-dir` and `--force` go: the census is written beside
its psp through the same `.partial`-then-rename as `generate-psps`, replacing what is there. The
byte-for-byte test (`each_census_it_writes_equals_the_one_the_walk_wrote`) is the oracle that the
rename changed nothing about the file.
*Depends:* A3, A4. *Source:* spec §8, §3.2.

☐ **C2 — fresh pairs skipped unless `--all`.** Each psp is judged with A2 before it is opened; a
fresh one is reported as skipped and not read; `--all` regenerates every one. Tests: a cohort with
one stale pair regenerates one file and names the rest as skipped; `--all` rewrites all and every
file is byte-identical to before.
*Depends:* C1. *Source:* spec §8.

☐ **C3 — a run stopped part-way costs only what is left.** Regenerate three, kill the process
between the second and third (a test that makes the third psp unreadable after two succeed, and
then readable again), run again: the third alone is rebuilt, the first two skipped as fresh, and no
`.partial` file is left behind.
*Depends:* C2. *Source:* spec §8, §10 (errors).

> **Checkpoint C: the repair command writes beside the psp and does only what is owed.** Pause for
> review.

### Milestone D — everything that spells the old shape, and the oracles re-run

☐ **D1 — scripts.** `scripts/ng_fit_stage_end_to_end.sh` (calls `generate-census`, passes
`--census`), and whichever of `ng_census_route_cost.sh` and `ng_census_agreement_mutations.sh`
name either. Run the end-to-end script on the six tomato accessions over the two 100 kb intervals:
same parameters file, same VCF as the last recorded run.
*Depends:* B4, C3. *Source:* spec §9 (last item), §11.

☐ **D2 — the probes landed.** `examples/ng_census_read_vs_psp.rs` and
`examples/ng_census_locus_spans.rs` from branch `census-vs-psp-perf`, unchanged; they are what
spec §2 and §12 cite.
*Depends:* —.

☐ **D3 — words.** The subcommand docs in `cli.rs` for the three commands and the repair;
`generate-psps`'s help, which still describes the census as written beside the psp and now says
that is the only place; `PROJECT_STATUS.md`'s pipeline line; the report for this plan under
`doc/devel/reports/implementations/`.
*Depends:* D1.

> **Checkpoint D: three commands and a repair, and the oracles unchanged.** Pause for review.

---

## 6. Verification summary

| milestone | proven by |
|---|---|
| A — the judgement | unit tests on fixtures for every verdict; two psps under different criteria refused by field; a five-sample cohort with three stale pairs names all three |
| B — step 2 over psps | the parameters file byte-identical to the flag-driven one on the tomato fixture (B1, B2); refusal tests that assert the message names every stale sample and the reference was never opened (B3) |
| C — the repair | `each_census_it_writes_equals_the_one_the_walk_wrote` after the rename; skip/`--all` tests; the stopped-and-resumed test leaves no `.partial` |
| D — the whole | `scripts/ng_fit_stage_end_to_end.sh` on six tomato accessions: same parameters file, same VCF, both routes' censuses still byte-identical |

---

## 7. Out of scope (next plans)

- **A record count in the psp footer** — a psp format change; [`psp_file_format.md`](../spec/psp_file_format.md).
- **A fit-side budget or cap, served by filtering a generous census** — the knobs question at
  [`parameter_prepass_runs.md`](parameter_prepass_runs.md) Checkpoint C.
- **Alleles under a widened record in the census** — [`parameter_prepass_census_sites.md`](../spec/parameter_prepass_census_sites.md) §2.
