# ng — the census lives inside the psp: implementation plan

**Status:** plan, 2026-09-09. **No code yet.** It turns the settled design in
[`psp_census_pair.md`](../spec/psp_census_pair.md) into build order and is **not a place for new
design**; where a step meets a question the spec does not answer, it stops at a checkpoint rather
than deciding. It follows [`parameter_prepass_runs.md`](parameter_prepass_runs.md), which built the
four commands this plan reshapes into three and a repair.

---

## 1. What this closes

Today a sample is two files, a user names census files at step 2 and retypes the walk's repeat
criteria there, and a stale census stops step 2 at the first one it meets. After this plan:

    generate-psps        alignments        ->  <sample>.psp, census inside
    estimate-parameters  psps              ->  cohort.parameters.toml, or a refusal naming every stale sample
    call-from-psps       psps + that file  ->  the VCF
    regenerate-census    psps              ->  each psp's trailer replaced, for the samples step 2 refused

No census file exists. No command takes a repeat criterion the psp header already holds.

---

## 2. Scope

**In:**

- the census written into the psp's trailer when the walk seals the file, and the second file
  gone (spec §3, §3.1, §3.2);
- the fit's lazy census reader taking a path and an offset (spec §5);
- one judgement of whether a psp's census is fresh, shared by step 2 and the repair command (spec
  §4.2, §8);
- a cohort whose psps disagree on the catalog or the criteria refused in `OpenPspCohort::open`, so
  every command that opens a cohort refuses it, `call-from-psps` included (spec §6);
- `estimate-parameters` over `--psp`, with the six duplicated flags removed and every stale sample
  reported before the reference is read (spec §4, §6);
- `regenerate-census` in place of `generate-census`: a read pass then `replace_trailer`, no
  `--output-dir`, no `--force`, fresh psps skipped (spec §8);
- an older-format census named as such rather than as damage (spec §4.2);
- a segmentation entry point that takes the criteria the header holds (spec §6, the trap);
- the sidecar's machinery deleted once nothing reads it: `PileupIdentity`, `freshness`,
  `freshness_by_header`, `psp_beside`, `census_path_for`, `CensusInCohort`;
- everything that spells the old shape: subcommand docs, scripts, `PROJECT_STATUS.md`; and the two
  measuring probes on branch `census-vs-psp-perf` landed under `examples/`.

**Out:**

- **a footer rebuilt from an intact index** — spec §12; [`psp_file_format.md`](../spec/psp_file_format.md) §6.5.
- **a fit-side budget or cap** — [`parameter_prepass_runs.md`](parameter_prepass_runs.md)
  Checkpoint C's knobs question; spec §7.
- **the census's per-position rule under a widened record** — spec §12;
  [`parameter_prepass_census_sites.md`](../spec/parameter_prepass_census_sites.md) §2.

---

## 3. Principles (how the order was chosen)

- **The file first, then everything that reads it.** Until the walk writes the census into the
  trailer nothing downstream can be tested against a real psp, so Milestone A is the writer and
  the reader, proven by the existing byte-for-byte oracle before any command changes.
- **The judgement before the commands.** Whether a psp's census is fresh is one function over a
  footer, a version word and a digest; both commands call it, so it is built and tested on
  fixtures first.
- **Isolate the step whose failure is silent.** Reading the criteria from the header instead of
  the flags (C1) changes where a number comes from and nothing about the arithmetic; wrong, it
  fits a plausible file. It lands as its own commit with the fitted file byte-identical before and
  after.
- **Reuse over rewrite.** `PspWriter::finish`, `replace_trailer`, `write_census`/`decode_census`,
  `census_from_psp`, `OpenPspCohort::open`, `first_difference` — called as they are.
- **Verify against ground truth.** The oracle for Milestone A is the existing walk-versus-rebuild
  byte comparison, now over trailer bytes; the oracle for the repair is a psp copied, regenerated
  and identical to the original whole; the oracle for the whole is
  `scripts/ng_fit_stage_end_to_end.sh` producing the same parameters file and VCF as before.
- **Types first, then implementation** (project rule).
- **Incremental, with pauses.** Five milestones, a checkpoint after each.
- **Builds through `./scripts/dev.sh` where a container runtime exists**, plain `cargo` where it
  does not (CLAUDE.md).

---

## 4. Preconditions (already in place)

Confirm each before step A1.

- The trailer as a replaceable payload: `PspWriter::finish(trailer)`
  ([`writer.rs:672`](../../../../src/ng/psp/writer.rs)), `replace_trailer`
  ([`trailer.rs:103`](../../../../src/ng/psp/trailer.rs)), the footer's `trailer_offset` and
  `trailer_bytes`, and `PspReader::trailer` ([`reader.rs:228-247`](../../../../src/ng/psp/reader.rs)).
- The census encoded to bytes and back: `write_census` into any `Write`, `decode_census` from a
  slice ([`census_file.rs:237`, `:456`](../../../../src/ng/parameter_estimation/joint/census_file.rs)).
- The walk builds its census and finishes the writer line: `write_psp` and `write_census_beside`
  ([`gatherer.rs:483-560`](../../../../src/ng/run/gatherer.rs)); the line's `finish`
  ([`psp_writer_line.rs:244`](../../../../src/ng/run/psp_writer_line.rs)).
- The settings in the header: `SegmentationInputs` and `first_difference`
  ([`segmentation_inputs.rs:25-55`](../../../../src/ng/segmentation_inputs.rs)).
- One census rebuilt from one psp: `census_from_psp` ([`census_from_psp.rs:202`](../../../../src/ng/run/census_from_psp.rs)).
- The cohort opener every command shares: `OpenPspCohort::open` ([`psp_caller.rs:111`](../../../../src/ng/run/psp_caller.rs)).
- The fixtures: `a_walked_cohort` in
  [`estimate_parameters/tests.rs:37`](../../../../src/pop_var_caller_exp/estimate_parameters/tests.rs)
  and [`generate_census/tests.rs:36`](../../../../src/pop_var_caller_exp/generate_census/tests.rs);
  the walk fixture in `generate_psps/tests.rs`.
- The end-to-end oracle: [`scripts/ng_fit_stage_end_to_end.sh`](../../../../scripts/ng_fit_stage_end_to_end.sh).

---

## 5. The steps

### Milestone A — the census in the trailer, written and read back

✅ **A1 — the walk hands the census to `finish`.** `PspWriterLine::finish` takes the trailer's bytes
and passes them to `PspWriter::finish`; `write_psp` finishes and encodes its census *before*
finishing the line, and hands the bytes across; the census is written with its pileup identity
absent (`None`, already an `Option`); `write_psp`'s census-path argument goes. A gatherer opened
without a plan writes an empty trailer. Test: a walked psp's trailer decodes as a census equal to
the one `CensusWriter::finish` returned.
*Depends:* —. *Source:* spec §3.1, §9 (first item).

✅ **A2 — `generate-psps` writes one file.** The census path, its `.partial`, its rename and the
ordering comment go ([`generate_psps.rs:676-735`](../../../../src/pop_var_caller_exp/generate_psps.rs));
the per-sample line and the report say the census's size inside the psp. Every test naming
`<sample>.census` (l.479-487, 699, 975, 1029, the block at 1128) rewritten to look in the trailer.
*Depends:* A1. *Source:* spec §3.2.

✅ **A3 — the lazy reader at an offset.** `open_census` and `SampleCensusEvidence::backed` take the
path and the trailer's extent; every section seek adds the offset
([`census.rs:1280`, `1326-1330`](../../../../src/ng/parameter_estimation/joint/census.rs)). The
resident `read_census` over the trailer bytes is the parity oracle: every section read lazily
equals the resident one, and `bytes_read()` shows only the sections asked for were read.
*Depends:* A1. *Source:* spec §5.

✅ **A4 — the two producers agree, byte for byte, over trailer bytes.**
`each_census_it_writes_equals_the_one_the_walk_wrote` rewired: the trailer the walk wrote against
`census_from_psp`'s output encoded the same way, on the fixture carrying a repeat tract and three
read groups. `census_from_psp` loses its identity argument.
*Depends:* A1, A3. *Source:* spec §11 (parity oracle).

> **Checkpoint A: a psp carries its census, and the fit can read it lazily from there.** Pause
> for review.

### Milestone B — one judgement, and a cohort that agrees on its settings

✅ **B1 — the verdict type.** An enum: fresh; no census (empty trailer); an older format, naming
the version found and the one this build reads; another selection. A noun with its own `Display`
in spec §4.2's words. No logic.
*Depends:* —. *Source:* spec §4.2.

✅ **B2 — the cheap half of the judgement.** From an open `PspReader`: the footer's trailer length,
then the version word at the trailer's front, read before `decode_census`'s magic-and-version check
can call it malformed. Fixtures for an empty trailer and a trailer whose version word is
`VERSION − 1`.
*Depends:* B1. *Source:* spec §4.2 rows 1-2, §9 (the `Malformed` trap).

✅ **B3 — the cohort agrees on its settings, in the psp cohort opener.** `OpenPspCohort::open`
compares every psp's `SegmentationInputs` against the first's with `first_difference` beside its
analysed-regions check ([`psp_caller.rs:677-690`](../../../../src/ng/run/psp_caller.rs)); a
disagreement is refused naming the sample and the field. `call-from-psps` gets it for free. Tests:
two psps under different `min_copies` refused by the opener, and refused at `call-from-psps`.
*Depends:* —. *Source:* spec §6, the decision paragraph.

✅ **B4 — a cohort judged whole.** Every psp of an opened cohort judged with B2, no early return;
one verdict a sample, in the order given. A test with three stale psps in five asserts all three
are named.
*Depends:* B2, B3. *Source:* spec §4.1.

> **Checkpoint B: a psp can be judged, a cohort's verdicts come back together, and a cohort that
> disagrees on its criteria is refused everywhere.** Pause for review.

### Milestone C — `estimate-parameters` over psps

✅ **C1 — the criteria from the header, own commit.** Split `segments_over` after
`routing_criteria` so a caller holding a `StrRepeatCriteria` enters there
([`run_ground.rs:242-281`](../../../../src/pop_var_caller_exp/run_ground.rs));
`estimate-parameters` builds its segmentation from the cohort's `SegmentationInputs` instead of
from its flags, **with the flags still present and ignored**. Oracle: the parameters file written
on the tomato fixture is byte-identical before and after. **Own commit, do not bundle.**
*Depends:* —. *Source:* spec §6 and its trap.

✅ **C2 — `--psp` in; `--census` and the five criteria flags out.** The cohort opened with
`OpenPspCohort::open`; the censuses read from the trailers with A3; the six deletions. The tests
that build argument lists (`args_over`, `a_shortest_run`) rewritten; C1's byte-identity re-asserted.
`open_census_cohort` and `CensusInCohort` go.
*Depends:* A3, B4, C1. *Source:* spec §5, §6.

✅ **C3 — the refusal, before the reference.** B4's verdicts are taken first; if any sample is not
fresh, the run stops with spec §4.3's report — every stale sample, its psp path, its cause, and one
`regenerate-census --psp …` line built from the arguments given — and the reference has not been
opened. `a_census_without_its_psp_is_refused` becomes *a psp without a census is refused and the
report names it*; a second test has two stale samples and asserts both are in the message.
*Depends:* C2. *Source:* spec §4, §4.1, §4.3.

✅ **C4 — the fit's own refusal says what to run.** `CohortFitError::AnotherSelection`
([`census_fit.rs:139`](../../../../src/ng/run/census_fit.rs)) is rendered as *these censuses were
built under another selection; run `regenerate-census`*.
*Depends:* C3. *Source:* spec §4.2 row 3.

> **Checkpoint C: step 2 takes psps, reads nothing it could be told wrongly, and refuses a stale
> cohort whole.** Pause for review.

### Milestone D — `regenerate-census`

☐ **D1 — the rename, and the write into the psp.** Module, `SUBCOMMAND`, `cli.rs` variant and its
doc, the fifteen tests and the two `SUBCOMMAND`-spelling tests. `--output-dir` and `--force` go;
per sample, `census_from_psp` then `replace_trailer`. The cohort opened with `OpenPspCohort::open`,
so B3's check applies.
*Depends:* A4, B3. *Source:* spec §8.

☐ **D2 — fresh psps skipped.** Each psp judged with B2 and, since the selection is rebuilt here
anyway, with the digest check as well; a fresh one is reported as skipped and its records not read.
Tests: a cohort with one stale psp rewrites one trailer and names the rest as skipped; a cohort
walked under a selection with another seed regenerates every one.
*Depends:* D1. *Source:* spec §8.

☐ **D3 — a run stopped part-way costs only what is left.** Regenerate three; make the third psp
unreadable after two succeed, then readable again; run again: the third alone is regenerated, the
first two skipped as fresh.
*Depends:* D2. *Source:* spec §8, §10 (errors).

☐ **D4 — the whole-file oracle.** A walked psp copied, `regenerate-census` run on the copy with
the same reference and catalog, and the copy identical to the original **byte for byte, whole** —
header, blocks, index, trailer, footer. On the fixture with a repeat tract and three read groups,
and in `scripts/ng_fit_stage_end_to_end.sh` on real reads (E1).
*Depends:* D1. *Source:* spec §11 (parity oracle).

> **Checkpoint D: the repair command rewrites only the tail, only where owed, and a regenerated
> psp is the walked one to the byte.** Pause for review.

### Milestone E — the sidecar removed, and everything that spells the old shape

☐ **E1 — scripts.** `scripts/ng_fit_stage_end_to_end.sh` (calls `generate-census`, passes
`--census`, diffs census files → the copy-regenerate-`cmp` of D4); whichever of
`ng_census_route_cost.sh` and `ng_census_agreement_mutations.sh` name either;
`examples/ng_census_route_cost.rs`'s `write_psp(path, Some(census))`. Run the end-to-end script on
the six tomato accessions over the two 100 kb intervals: same parameters file, same VCF as the last
recorded run.
*Depends:* C4, D4.

☐ **E2 — the sidecar's machinery deleted.** `PileupIdentity`, `freshness`, `freshness_by_header`
and their tests; `psp_beside`, `census_path_for`, `CENSUS_FILE_EXTENSION`; what is left of
`census_cohort.rs`. `cargo` says what still reads them; nothing should.
*Depends:* C2, D1.

☐ **E3 — the probes landed.** `examples/ng_census_read_vs_psp.rs` and
`examples/ng_census_locus_spans.rs` from branch `census-vs-psp-perf`; the first writes its census
to a scratch file of its own and is unaffected by E2.
*Depends:* —.

☐ **E4 — words.** The subcommand docs in `cli.rs` for the three commands and the repair;
`generate-psps`'s help, which describes a census beside the psp; `PROJECT_STATUS.md`'s pipeline
line; the report for this plan under `doc/devel/reports/implementations/`.
*Depends:* E1, E2.

> **Checkpoint E: one file per sample, three commands and a repair, and the oracles unchanged.**
> Pause for review.

---

## 6. Verification summary

| milestone | proven by |
|---|---|
| A — the census in the trailer | a walked psp's trailer decodes to the census the writer returned; lazy reads at an offset equal resident reads and read only the sections asked for; the walk-versus-rebuild byte comparison over trailer bytes |
| B — the judgement | fixtures for every verdict; two psps under different criteria refused by the opener and by `call-from-psps`; a five-sample cohort with three stale psps names all three |
| C — step 2 over psps | the parameters file byte-identical to the flag-driven one on the tomato fixture (C1, C2); refusal tests asserting every stale sample is named and the reference was never opened (C3) |
| D — the repair | skip tests; the stopped-and-resumed test; a copied psp regenerated and identical whole (D4) |
| E — the whole | `scripts/ng_fit_stage_end_to_end.sh` on six tomato accessions: same parameters file, same VCF; a regenerated psp `cmp`-identical to the walked one |

---

## 7. Out of scope (next plans)

- **A footer rebuilt from an intact index** — [`psp_file_format.md`](../spec/psp_file_format.md) §6.5.
- **A fit-side budget or cap** — [`parameter_prepass_runs.md`](parameter_prepass_runs.md) Checkpoint C.
- **Alleles under a widened record in the census** — [`parameter_prepass_census_sites.md`](../spec/parameter_prepass_census_sites.md) §2.
