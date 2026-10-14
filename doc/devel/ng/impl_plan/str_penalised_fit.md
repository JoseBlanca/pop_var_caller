# The repeat-tract fit by penalised EM — implementation plan

**Draft, 2026-10-14.** Turns the settled design — specs
[`../spec/str_penalised_fit.md`](../spec/str_penalised_fit.md) and
[`../spec/str_fit_to_calling.md`](../spec/str_fit_to_calling.md), architecture
[`../arch/str_penalised_fit.md`](../arch/str_penalised_fit.md) — into build order. **It invents no
design**: a question that comes up here goes back to the specs or the architecture.

## Scope

**In:** the penalised fit as `estimate-parameters`' repeat-tract fit; parameters-file format 3 with
reading of older files; calling's new fallback step; the checks of spec §4; deleting the climb.

**Out, each with its home:**
- stutter per slippage group and library multipliers → `spec/str_library_stutter.md`, rebuilt on this
  fit;
- reading the census one motif length at a time → the census-read step (spec §6);
- the full-scale run on kimura's real data → when kimura is back (§D3 measures it on a simulated
  cohort here meanwhile);
- re-recording main's own stale oracle baseline → done on this branch by the merge (`e0b41ad4`).

## Principles (how the order was chosen)

- **Types first, then behaviour.** Milestone A changes shapes only; every number stays the same.
- **The hand-over before the fit.** Calling must hold an absent F (milestone B) before any fit can
  produce one (milestone C); in the other order calling panics (spec §3).
- **The prototype is the oracle.** The research bench's `D_tract_prior` on `str-fit-compare` gives
  the numbers the production fit must reproduce.
- **Isolate the silent steps.** Two steps can change numbers without failing anything: making the fit
  deterministic (C3) and switching `estimate-parameters` to it (C6). Each is its own commit, with its
  oracle run before and after.
- **Delete last.** The climb goes only after the new fit has passed on real data (D2) and in calling
  (D5).

## Preconditions (already in place)

- Branch `str-penalised-fit` at `1183fc0e` or later: main merged in (`e0b41ad4`), specs approved.
- The research code on the branch: `ssr_fit/penalised_fit.rs`, `composite_fit.rs`,
  `counting_fit_errors.rs`, `bench_fixtures/{fit_bench,cohort_draw}.rs`,
  `examples/ng_str_fit_bench.rs`.
- tomato2 psps in format 2.0: `/Users/jose/devel/pop_var_caller-str-guard/tmp/tomato2/psps`
  (59 accessions, written 2026-10-14 by main's `pvc`).
- HG002 (`benchmarks/ssr_hg002`) and the GIAB trio (`benchmarks/giab`) as in the earlier STR reports.
- Known test failures before step 1: three `cli::cross_platform_digests` tests (stale, spec §3) and
  three debug-only panic tests in release.

## The steps

### Milestone A — shapes only, no number moves

- ☐ **A1** `StratumFit::concentration` becomes `Option<f64>`; the climb always fills `Some`. Every
  reader updated. *Depends:* —. *Source:* arch §1.2.
- ☐ **A2** The climb's own fields move into `climb: Option<ClimbRecord>`; its log summaries read
  them there. *Depends:* A1. *Source:* arch §1.2.
- ☐ **A3** `LevelSource::PenalisedCurve`, `ShareSource::PenalisedCurve`, `PenalisedPoint`,
  `StratumError::TooFewLoci` — added, unused. *Depends:* —. *Source:* arch §1.3–1.4.
- ☐ **A4** `FittedAlleleLengthVariation::concentration` becomes `Option<f64>`;
  `AlleleLengthVariationRung::StratumsOwnSharesStatedSpread` added; `allele_length_variation_at`,
  `checked_allele_length_variation` and `median_concentration` handle `None`. Tests for the new rung
  with a hand-built table. *Depends:* A1. *Source:* arch §1.5; companion spec §3.2.

> **Checkpoint A:** the oracle (`scripts/promote_ng_oracle.sh`) gives the same checksums as at the
> precondition commit; the touched tests pass. Pause for review.

### Milestone B — the parameters file, format 3

- ☐ **B1** `between_locus_share_of` / `concentration_of`; the rows' and fallback's keys renamed and
  holding F; `between_locus_share_standard_error`; `FORMAT_VERSION = 3`; validation of F in (0, 1).
  Golden files regenerated, diff checked. *Depends:* A4. *Source:* arch §1.6, §2.3; companion spec
  §2.1.
- ☐ **B2** `migrate_to_format_3` on the parsed TOML; `version_1_as_written.toml` restored to its old
  keys (from `aa42d629~1`); tests that a format-1 and a format-2 file read back with F = 1/(1+c).
  *Depends:* B1. *Source:* arch §2.3; companion spec §2.2.
- ☐ **B3** `LevelSmoothing::PenalisedCurve` and `ShareSmoothing::PenalisedCurve`, written and read.
  *Depends:* A3, B1. *Source:* arch §1.6; companion spec §2.3.

> **Checkpoint B:** a format-3 file round-trips, F absent included; the oracle's VCFs are unchanged
> and its parameters files differ only in the renamed keys and F's values. Pause for review.

### Milestone C — the fit

- ☐ **C1** `composite_fit.rs` renamed `expected_counts.rs`; no code change. *Depends:* —.
  *Source:* arch, Module home.
- ☐ **C2** `PenalisedFitConfig::production()` and the named constants; the research-only fields
  removed; the bench's `D_tract_prior` built from `production()`. *Depends:* C1. *Source:* arch §1.1;
  spec §2.7.
- ☐ **C3** **Own commit, do not bundle.** Every sum in the penalised fit taken in index order, so
  two runs and runs at 1 and 18 threads give the same bits; a test pins it. *Oracle:* the bench's
  `D_tract_prior` numbers at the cells of `compare_D_prior_2026-10-14.md` agree with the prototype's
  to 1e-9 before and after. *Depends:* C2. *Source:* spec §3; arch §2.1.
- ☐ **C4** Memory bounded: at most one class's distinct rows per worker. *Oracle:* bench numbers
  bit-identical before and after; peak memory at 2,000 samples × 30 reads below the counting fit's
  1.1 GB × 2. *Depends:* C3. *Source:* spec §3; arch §2.1.
- ☐ **C5** `fit_repeat_tracts`: per-period grouping, the `StratumOutcome` conversion (`Fitted`,
  `Derived`, `Refused`), provenance as `PenalisedCurve`, errors with `TooFewLoci` below
  `LOCI_FOR_A_SHARE_ERROR`. Tested on drawn cohorts against the bench. *Depends:* C4, A2, A3.
  *Source:* arch §2.1.
- ☐ **C6** **Own commit, do not bundle.** `fit_a_cohort` calls `fit_repeat_tracts`;
  `--str-param-estimates-at-once` removed with an error. *Oracle:* the oracle script, run before
  and after; the move is explained (expected: repeat-tract rows of the fitted file, and the calls
  made with it). *Depends:* C5, B3. *Source:* arch §2.2, §2.4.

> **Checkpoint C:** `estimate-parameters` runs the new fit end to end on tomato2; the oracle's move
> is written down. Pause for review.

### Milestone D — checks (spec §4)

- ☐ **D1** The bench reproduces `compare_D_prior_2026-10-14.md`'s `D_tract_prior` rows with the
  production fit, within run-to-run noise. *Depends:* C6. *Source:* spec §4.2.
- ☐ **D2** tomato2, held-out by whole locus: the production fit against the counting fit, with the
  bench's scorer on real psps. *Depends:* C6. *Source:* spec §4.4.
- ☐ **D3** Full scale, simulated: 2,169 samples with kimura's loci per class (about 66 million
  pairs with reads), the production fit's wall time and peak memory on this Mac. Within one hour,
  with a memory watchdog. *Depends:* C6. *Source:* spec §4.5 (kimura's real run deferred).
- ☐ **D4** The oracle and `cli::cross_platform_digests` re-recorded, each move explained.
  *Depends:* C6. *Source:* spec §4.6.
- ☐ **D5** Calling on HG002 against its repeat-tract truth and on the GIAB trio: not worse than with
  the climb's parameters. *Depends:* C6. *Source:* spec §4.7.
- ☐ **D6** Does a slip rate 0.5% low change calls? Calling on simulated 500 and 2,000 samples × 3
  reads, fitted vs true stutter. *Depends:* C6. *Source:* spec §5.1.

> **Checkpoint D:** a short report under `doc/devel/reports/implementations/` with every D result.
> Pause for review: the owner decides whether D2 and D5 pass before the climb is deleted.

### Milestone E — delete the climb

- ☐ **E1** The climb's code, tests and examples deleted in one commit (`fit_strata*`, `fit_stratum`,
  the quadrature, the walks, `SampleSubsets`, `draw_and_apply_curves`'s blend, `ClimbRecord`).
  *Depends:* D2, D5 and the owner's go-ahead. *Source:* spec §7.
- ☐ **E2** The specs the climb built on get a dated note that the fit they describe was replaced
  (`parameter_prepass_ssr.md`, `str_slippage_level_curve.md`), pointing at the new spec.
  *Depends:* E1. *Source:* spec-authoring's rule for decisions that were acted on.

> **Checkpoint E:** full suite and code review. Pause for review before any merge into main.

## Verification summary

| milestone | proven by |
|---|---|
| A | oracle checksums unchanged; new-rung tests |
| B | round-trip and migration tests; golden-file diffs; oracle VCFs unchanged |
| C | bench parity with the prototype (C3, C4); determinism test; oracle move explained (C6) |
| D | the bench, tomato2 held-out, simulated full scale, HG002 and GIAB calling |
| E | full suite and code review |

## Out of scope (next plans)

Slippage groups and library multipliers (`str_library_stutter.md`); per-period census reading;
kimura's full-scale run; the read-crossing bias for long alleles (`str_allele_window.md` §5.3).
