# ng — the hidden-duplication filter: implementation plan

*Draft, 2026-09-06. Turns [`../spec/hidden_paralog_filter.md`](../spec/hidden_paralog_filter.md)
into build order — no new design here; a design question surfacing mid-plan goes back to that
spec, and the spec's one open decision (every record scored, non-SNPs on coverage alone, §3.2) is
confirmed with the owner before step C1 is coded. There is no architecture document: the spec's
§3.7 type blocks are the code shape the steps cite. This plan starts after
[`window_coverage.md`](window_coverage.md)'s Checkpoint C. **The standing oracle is that
`--paralog-fdr 0` writes byte for byte what the run writes today**, on the run fixtures and the
tomato slice.*

## Scope

**In:** production's statistics copied into `src/ng/paralog/` with their tests and a
bit-identity differential; the spill's entry, codec and lifecycle; the writer's line entry point
and the two-column patch; the scoring context built from the histograms and the parameters file;
the three passes wired behind the sink both modes share; the two command-line flags, the header
line, the two INFO fields, the FILTER id, the run report's lines; the runs at the ends of the
range and the measurements the spec leaves open.

**Out (later plans):**

- **a tract-aware allele term** and **the coverage-only score's power** beyond the count this
  plan reports — spec §8;
- **scoring in parallel** — gated on D2's pass shares;
- **streaming compression of the spill** — spec §3.4, if it proves easy;
- **any change to which loci the merge builds** — the merge's own documents.

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone.
- **The copy before the wiring, and a differential before anything reads it.** The statistics
  are transcribed from `src/paralog/` and proven bit-identical to production's on shared inputs
  before a single record is scored — the port is correct when it matches what it ports.
- **The algorithmic heart before the plumbing.** Score → spill → passes → flags, in that order;
  each is tested alone before the next depends on it.
- **Isolate the silent steps.** Three failures here produce a plausible file rather than a
  panic: a non-SNP sample skipped at zero reads (every tract scores nothing), a `NaN` turned into
  a zero (every unscored record folds into π), and a patched line whose other columns moved.
  Each lands as its own commit with its oracle named.
- **Off first.** The spill and the passes are wired behind a flag whose default-off run is
  byte-identical before the flag's on-path is finished; the oracle exists before the behaviour.
- **Verify against ground truth.** Production's scorer for the ratio; the run before this work
  for the file; the tag-mode file for the drop-mode file.
- **Container builds**: `cargo` via `./scripts/dev.sh`; a native host build at completion.

## Preconditions (verify before step A1)

- **[`window_coverage.md`](window_coverage.md) Checkpoint C is passed**: every written record
  reaches the sink with its per-sample pairs, and the cache hands back one histogram per sample.
  This plan does not start without it.
- **Production's statistics are where the spec's reuse map says**
  ([`src/paralog/`](../../../../src/paralog/): `mod.rs`, `coverage_model.rs`, `locus_score.rs`,
  `prior.rs`; [`calibrate.rs:64-118`](../../../../src/var_calling/paralog_filter/calibrate.rs) for
  the calibration and its verdict) and production is frozen — nothing under `src/paralog/` or
  `src/var_calling/paralog_filter/` is edited.
- **The parameters file carries one inbreeding coefficient per sample**, exposed as
  [`RunParameters::inbreeding_coefficient_by_sample`, `run_parameters.rs:695`](../../../../src/ng/calling/run_parameters.rs).
- **The writer takes a record and checks its order**
  ([`write_record`, `writer.rs:114`](../../../../src/ng/vcf/writer.rs)); the header writes meta
  lines as `##key=value` ([`header_text`, `header.rs:431`](../../../../src/ng/vcf/header.rs)); the
  declarations are the two constant lists beside it.
- **The varint primitives** ([`psp/varint.rs:46-126`](../../../../src/psp/varint.rs)) are `pub`.
- **Ground to run on**: the six-accession tomato slice
  ([`scripts/ng_fit_stage_end_to_end.sh`](../../../../scripts/ng_fit_stage_end_to_end.sh)), a
  single accession from it, and the HG002 psp.

---

## The steps

### Milestone A — the statistics, copied and proven equal

**A1. ✅ The constants and the coverage model.** `src/ng/paralog/mod.rs` and
`coverage_model.rs`: `ParalogModelParams` and its defaults, `CoverageFitConfig`,
`SingleCopyCoverageModel::fit` and `relative_copy_number`, `CoverageModelError`, with production's
tests, reading ng's histogram type. Green as transcribed. *Depends:* —. *Source:* spec §3.1, §7.

**A2. ☐ The score. Own commit, do not bundle.** `locus_score.rs`: `SampleObservation`,
`LocusObservations`, `ParalogScorePrecompute`, `score_locus_for_paralogy`, with production's
tests; plus a differential that feeds both implementations the same randomised inputs — cohort
sizes 1, 2, 10, 63; zero-read samples included — and asserts the ratio, both log-likelihoods and
the counts **bit-identical**. `N = 1` is in the differential deliberately (spec §4). *Depends:*
A1. *Source:* spec §3.2, §7, §10.

**A3. ☐ The prior, the curve, the calibration.** `prior.rs`: `ParalogLrHistogram`,
`ParalogPrior::estimate`, `ParalogFdrCurve`, `lr_threshold_for_fdr`, with production's tests; the
`ParalogCalibration` type with `flags` and `posterior` and the fallback prior, copied from
`calibrate.rs`. *Depends:* —. *Source:* spec §3.3, §7.

> **Checkpoint A: production's filter statistics live in ng, and the score is bit-identical
> to production's on shared inputs including one sample. Pause for review.**

### Milestone B — the spill

**B1. ☐ The entry and its codec.** `src/ng/run/paralog_filter/spill.rs`: `SpillEntry` and
`SpilledSample` as spec §3.7; a writer and a reader over the varint primitives, the encoder
destructuring the entry exhaustively, floats by bits; round-trip tests over entries with absent
pairs, zero samples, a `FILTER` already set, and an `INFO` of `.`. *Depends:* —. *Source:* spec
§3.4, §3.7.

**B2. ☐ The lifecycle.** The file at `<output>.paralog-spill.tmp`, created on first write,
removed by a guard on every exit path — success, a `RunError`, a panic unwinding through the
run — with a test for each. *Depends:* B1. *Source:* spec §3.4.

**B3. ☐ The writer's line entry point and the patch. Own commit, do not bundle.**
`VcfWriter::write_line(place, &[u8])` running `check_order` from the entry's three head fields;
a patch function that splits a line on its first eight tabs and rewrites `FILTER` (join with
`;`, never replace a non-`PASS`) and `INFO` (append, or replace a `.`). Tests: a patched line
with nothing to add is byte-identical to its input; a record with nine columns and one with a
thousand sample columns patch the same two; a line with a `.` INFO. *Depends:* —. *Source:*
spec §3.5, §6 traps 5–7.

> **Checkpoint B: an entry round-trips with its absences intact, the file cannot outlive its
> run, and a line goes out unchanged unless the verdict touches it. Pause for review.**

### Milestone C — the three passes, behind a flag

**C1. ☐ The scoring context. Own commit, do not bundle.** `ParalogScoringContext::new` from the
histograms, the parameters file's coefficients and the model params: one fit per sample with
rejections kept by reason, the σ₀ slice with `NaN` where absent, the precompute; and
`observation_of(entry, sample) -> Option<SampleObservation>` applying spec §3.2's rule — a
biallelic SNP hands the scorer its AD and is skipped at zero reads; **every other record hands
`0/0` and is not skipped**. The silent failure is the skip applied to a tract; the test is a
tract entry whose samples all score, against a SNP entry whose zero-read sample does not.
*Depends:* A1, A2, B1; **the spec's §3.2 decision confirmed with the owner.** *Source:* spec
§3.1, §3.2, §6 traps 1–3.

**C2. ☐ Pass one behind the flag, off-path proven.** `--paralog-fdr` and `--paralog-filter-tag`
on both subcommands; with the target at zero the sink is the VCF writer as today and no spill
exists. **Byte-identical on the run fixtures and the tomato slice before C3 starts.** With the
target above zero the sink writes spill entries — the line from `record_line`, the pairs from the
sink's slice, the AD from the record's columns, the SNP test on the record's own alleles — and
the writer is not opened until pass three. *Depends:* B1, B2. *Source:* spec §3.4, §3.6, §1.1
goal 7.

**C3. ☐ Pass two. Own commit, do not bundle.** Read the spill; score each entry through C1;
keep the ratio in a `Vec<f64>` and fold the finite ones into the LR histogram; estimate π, build
the curve, resolve the cut; warn and fall back when the EM does not converge. The silent failure
is a `NaN` becoming a zero somewhere between the scorer and the histogram; the test folds a
mixed stream and asserts the histogram's total equals the finite count. *Depends:* A3, C1, C2.
*Source:* spec §3.3, §6 trap 4.

**C4. ☐ Pass three.** Read the spill again in step with the ratios; drop or tag per spec §3.5;
append the two INFO fields on every finite ratio; write through `write_line`; delete the spill;
count. The header gains the `##paralogFilter=` line and the three declarations; the run report
gains its lines. *Depends:* B3, C3. *Source:* spec §3.5, §3.6.

**C5. ☐ The oracles, all four.** (i) off is byte-identical (C2's, rerun); (ii) on at an
unreachable target, the two INFO keys stripped, equals off; (iii) the drop-mode file equals the
tag-mode file minus its `hiddenParalog` lines; (iv) mode equivalence with the filter on, both
modes' files identical. All on the six-accession slice; (i) and (iv) also on the run fixtures.
*Depends:* C4. *Source:* spec §10.

> **Checkpoint C: the filter runs in both modes, off is exactly today, and the three
> relations between its outputs hold on real reads. Pause for review.**

### Milestone D — the range, and what the spec left to measurement

**D1. ☐ One sample.** A single tomato accession over the slice, filter on: the fit's outcome,
π and whether the EM converged, records dropped. The check is that it runs, that the precompute
at `N = 1` produces finite ratios, and that the header says what happened. Spec §9's second OPEN
closed. *Depends:* C5. *Source:* spec §4, §9.

**D2. ☐ Six samples, with the pass shares.** The slice at six, filter on: records dropped and
tagged, samples rejected by reason, π, the cut; wall per pass — calling, fit, scoring, writing —
and peak resident, against the filter-off run. The scoring share is what decides whether
parallel scoring gets a plan (spec §8). *Depends:* C5. *Source:* spec §4, §5.

**D3. ☐ HG002.** Filter on over the HG002 psp: the fit accepts (the scaled depth axis, spec §4)
or the report says why not; the records dropped, with the ten highest ratios listed for a human
eye. *Depends:* C5. *Source:* spec §4.

**D4. ☐ The coverage-only score, counted.** On D2's and D3's runs: how many biallelic-SNP,
indel, multiallelic and tract records were scored, how many of each were flagged, and the
distribution of ratios by kind. Beside it, production's filter over the same six accessions as a
sanity comparison — the flagged SNP sets' overlap, reported, not chased. **This is what closes
spec §9's first OPEN or sends §3.2 back to the spec.** *Depends:* D2, D3. *Source:* spec §3.2,
§8, §9.

> **Checkpoint D: the filter has run at one sample, at six, and at three hundred reads a
> position, with numbers at each; the one decision beyond production has a measurement
> behind it. Pause for review.**

## Verification summary

| milestone | proven by |
|---|---|
| A — the statistics | production's tests green as transcribed; the randomised differential bit-identical, `N = 1` included |
| B — the spill | round trips with absences by bits; the guard on every exit path; the patch a no-op where nothing is added |
| C — the passes | the four oracles: off identical, unflagged identical, drop = tag − tagged, modes identical |
| D — the range | runs at 1, 6 and HG002 reported with dropped counts, rejections and π; the coverage-only counts by kind; the sanity comparison against production |

## Out of scope (next plans)

- **A tract-aware allele term**, if D4 says the coverage-only score is not enough at tracts —
  a document under the spec.
- **Parallel scoring**, if D2's share says pass two is worth it.
- **Spill compression** — one writer, one reader, when convenient.
- **The trailer histogram** — [`../spec/window_coverage.md`](../spec/window_coverage.md) §8.
