# The SNP/indel fit converges — implementation plan

**Draft, 2026-10-16.** Turns the settled design into build order: the stopping rule of
[`../spec/fit_precision.md`](../spec/fit_precision.md) §2 (amendment of 2026-10-16), and §2.1.3,
§3.3.1 and §12 tests 11–12 of [`../spec/parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md),
all decided by the owner on 2026-10-16 (that spec's §11 question 13). **It invents no design**: a
question that comes up here goes back to those sections. There is no architecture document; the
types are in the specs.

**How each step is run (owner's process).** A subagent implements the step; a fresh agent reviews
it; the implementer fixes what the review finds; the lead commits. The lead pauses for the owner
only at the checkpoints marked *owner*, where the outcome is a decision.

## Scope

**In:** the Newton-decrement stopping rule; the Newton step with halving; the not-separately-determined
test for the density's pair and the duplicated class's pair, with the errors of the determined
shares; the run log's line; the parameters file's format 4 with its two informational tables; the
checks — the identity oracle, the cross-platform digests, calls on the 4 oracle psps and the 59
tomato2 psps against the current fit, and timings at 1, 4 and 59 samples.

**Out, each with its home:**

- whether the outer-product information is too sharp along the slowest direction, and Louis's
  method for the cohort's eight → `fit_precision.md` §8 question 3, a measurement of its own;
- reparameterising the density around the determined share → `parameter_prepass_joint_fit.md` §11
  question 13 part 2, only if a consumer needs the density's numbers themselves;
- a cohort above 59 samples, or another species or depth → the next full-scale run (kimura); the
  threshold was measured on tomato only (`fit_precision.md` §2).

## Principles (how the order was chosen)

- **Isolate the silent steps.** Three steps change numbers without failing anything: the stopping
  rule (A2), the Newton step (B1) and the file format (D2). Each is its own commit, with the oracle
  run before and after and every moved line explained in the commit message.
- **The prototype is the oracle for the fit.** The two prototypes already ran the built design;
  their runs at threshold 0.001 (`tmp/snp_stop/runs/{s1,o4}_d0.001.*`) are what B1 must reproduce,
  bit for bit.
- **Compute before writing.** The verdicts are computed and tested (milestone C) before anything
  writes them (milestone D), so C changes no output at all.
- **The heart before the plumbing.** The decrement and the verdict are pure functions with unit
  tests (A1, C1–C3) before `maximise` or the file reads them.
- **Light dev loop** (owner, 2026-10-10): touched tests during a step; the oracle at each silent
  step; the full suite, both platforms and a code review once, at milestone E.

## Preconditions (already in place)

- Branch `snp-fit-convergence` from main `51b7be99`, with the spec amendments of 2026-10-16
  committed.
- The prototypes, as reference code, **not to be applied**: every behaviour in them hangs off an
  environment variable. `tmp/snp_conv/prototype_knobs.diff` (the Newton step, in `maximise`) and
  `tmp/snp_stop/prototype_stopping.diff` (the decrement, `newton_decrement`, `cohort_covariance`, the
  judge table).
- The measurement runs: `tmp/snp_stop/runs/` (the references at 10⁻⁴: `s1_d1e-4`, `o4_d1e-4`,
  `t59_d1e-4.judge.tsv`; the runs at 0.001: `s1_d0.001`, `o4_d0.001`) and
  `tmp/snp_plan/dup_trade.py` (the verdicts from a judge table).
- The current fit's files, for the call comparisons: `tmp/snp_conv/runs/o4_base.toml` (4 oracle
  samples) and the main checkout's `tmp/snp_fit_baseline/t2_59.parameters.toml` (59 tomato2
  samples, 2026-10-15). Both are format 3; the new build reads them.
- Inputs: the oracle psps (main checkout `tmp/snp_fit_baseline/oracle_psps/`, `SRS1839219.psp` for
  one sample), tomato2's 59 psps (main checkout `tmp/tomato2/psps/`), the reference
  `$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa` and the catalog
  `benchmarks/tomato1/crams/S_lycopersicum_chromosomes.4.00.repeats.parquet` in the main checkout,
  as `tmp/snp_stop/fit.sh` passes them.
- Before A1 the executor runs `scripts/promote_ng_oracle.sh` and `cli::cross_platform_digests` on
  the precondition commit and records which pass; the suite's known failures in release, if any,
  are listed in A1's commit.

## The steps

Each step names what it touches, its tests, and its **identity rule**: what must stay bit-identical,
what moves, and how the move is checked.

### Milestone A — stop on the Newton decrement

- ✅ **A1** `newton_step` returns the decrement beside the distances (one `bounded_step`, one
  `information_times`); `DECREMENT_THRESHOLD` = 0.001 and the verdict `has_converged(decrement,
  threshold)` in `fit/settled.rs`; `JointFitConfig::decrement_threshold`, defaulted to it. Callers of
  `newton_step` take the distances; nothing reads the decrement yet.
  *Touches:* `fit/standard_errors.rs`, `fit/settled.rs`, `fit.rs` (config, call sites).
  *Tests:* the decrement equals `½ gᵀ I⁻¹ g` with no bound binding and `g·d − ½ dᵀ I d` with one,
  on the whole matrix and on the arrow; the verdict at −2,959, NaN, 0, just below and at 0.001.
  *Identity:* every output bit-identical; the touched tests and the digests pass.
  *Depends:* —. *Source:* `fit_precision.md` §2 (where it is computed; traps), §3.6 item 4.
  *Done 2026-10-16:* `newton_step` returns a `NewtonStep { distances, decrement }`; the decrement and
  `has_converged` are marked `expect(dead_code)` until A2 reads them. The precondition oracle matched
  the baseline on every line; the fit module's 122 release tests and the four digest tests pass.
- ✅ **A2** **Own commit, do not bundle.** `maximise` stops on the decrement, in the order of a
  judging pass: converged → stop at the judged point, the rest of the cycle not run; else agreement;
  else the SQUAREM cycle as today. The progress line prints the decrement; `PassSummary` gains
  `decrement`. The three texts that state the old meaning of *converged* change (the file's note at
  `to_toml.rs:159` and the golden file, the outcome's doc at `parameters_file/mod.rs:714`,
  `settled.rs`' module doc). `not_settled` and `furthest` stay, for the progress line and a
  non-converged start's report.
  *Touches:* `fit.rs` (`maximise`, `PassSummary`), `fit/settled.rs`, `calling/parameters_file/`.
  *Tests:* on a drawn cohort that converges, the start returns the parameters its last judging pass
  judged, after one final pass; a start whose decrement is negative is not stopped (the verdict's
  test from A1, and a drawn case if one is found); the touched fit tests.
  *Identity:* every pass before a start's first judging pass is unchanged. **Moves:** where a start
  converges — its passes, and the values it returns (no longer a cycle past the judged point) — and
  the note text in every fitted file. On the oracle no start converged under the old rule; if none
  now reaches 0.001 within 200 passes, only the note moves, and the commit shows it by diffing the
  fitted file with the note restored. The digests move with the note and with any start that
  converges on the two-sample fixture; re-recorded in this commit, each moved line explained.
  *Depends:* A1. *Source:* `fit_precision.md` §2 (the rule, its order), §3.5.
  *Done 2026-10-16:* on the oracle only the note moved (starts at 200, 199, 200 passes as before;
  the fitted file with the note restored is the baseline's byte for byte). On the digests' two
  two-sample fixtures the starts **stop converging** (18 → 198 passes; 49, 54, 15 → 199, 198, 198):
  no parameter there has an error, so the old test called them settled, while the bounded solve ends
  at its round limit and returns a decrement near 290. Three slow drawn-cohort tests are ignored
  until B1 (two crawl without the Newton step; one fixture needs redrawing).

> Checkpoint A (lead): A2's oracle moves explained. No owner pause unless a move is not.

### Milestone B — the Newton step with halving

- ☐ **B1** **Own commit, do not bundle.** At a judging pass that does not stop the start, try its
  Newton step as `parameter_prepass_joint_fit.md` §3.3.1 sets out: the target from the same
  distances, the first fraction `min(1, 2·f_last)`, a try kept only if strictly higher, the kept
  try's pass reused as the next cycle's judging pass, up to `NEWTON_TRIES` = 4 (beside `JUMP_SLACK`
  in `fit.rs`), a try past the shares' bound halved without a pass, SQUAREM from the judged point
  when no try rises. The progress line for every cycle says which way it went.
  *Touches:* `fit.rs` (`maximise`).
  *Tests:* spec §12 test 11's unit cases — overshoot halved and the halved try kept; the remembered
  fraction used next; every try falling → SQUAREM from the judged point; the share bound → halved
  without a pass; no try once `max_passes` is spent; the same bits at 1, 4 and 8 threads
  (`a_census_of_many_chunks_fits_to_the_same_bits_at_any_pool_width` and the digests' pool-width
  test).
  *Identity:* **the built fit reproduces the prototype's runs at 0.001** on 1 sample and the 4
  oracle samples (`tmp/snp_stop/runs/s1_d0.001`, `o4_d0.001`): passes by start, how each ended, and
  every fitted number in the parameters file, bit for bit; any difference is traced to its cause
  before commit. **Moves:** the oracle's fitted file and the VCF called with it (expected: fewer
  passes, every start converged or agreed; the call differences of the diagnosis's §5.1 size), and
  the digests, re-recorded here with the moves explained. A start that never reaches a judging pass
  is unchanged.
  *Depends:* A2. *Source:* `parameter_prepass_joint_fit.md` §3.3.1, §12 test 11.

> Checkpoint B (lead): B1 reproduces the prototype; the oracle's moves are explained. No owner pause
> unless the reproduction fails for a reason that is not a bug.

### Milestone C — what the data determine (computed, not yet written)

- ☐ **C1** The regularised incomplete Beta function `I_x(a, b)` beside `genetics::lgamma`
  (`src/genetics.rs`), by continued fraction, through `crate::float`.
  *Tests:* against values pinned from an independent implementation (SciPy's `betainc`, recorded in
  the test with their source) at `x` = 0.01 and shapes 0.02, 1 and 50; the symmetry
  `I_x(a, b) = 1 − I_{1−x}(b, a)`; `x` = 0 and 1.
  *Identity:* nothing calls it; every output bit-identical.
  *Depends:* —. *Source:* `parameter_prepass_joint_fit.md` §2.1.3, traps.
- ☐ **C2** `cohort_covariance` beside `newton_step`: the cohort's 8 × 8 block of the inverse
  information, by unit-slope solves on the whole matrix or the arrow, `None` for a number not
  solved for.
  *Tests:* equal to the dense inverse restricted to the cohort (the module's
  `dense_errors_without` helpers) on the whole matrix, the arrow and an arrow of mixed libraries; its
  diagonal's square roots equal `StandardErrors`' estimated errors to rounding; a dropped number's row
  is `None`.
  *Identity:* nothing calls it; every output bit-identical.
  *Depends:* —. *Source:* §2.1.3, how each number is computed.
- ☐ **C3** `fit/determined.rs`: `ReportedNumber`, `DeterminedShare`, `DensityReport`,
  `DuplicatedReport`; `STATED_FREQUENCY`, `SEPARATION_LOG_LIKELIHOOD`; the determined share, its
  error, the errors of `Hexp` and `E[f]`, the trade cost and the verdict for both pairs; a marked
  number's error withheld. `JointFit` gains `density_report` and `duplicated_report`, computed after
  the winning start's final pass from the information its errors come from.
  *Tests:* on hand-built covariances, `V_trade` at or below zero gives no trade, a number with no
  covariance marks the pair, a trade cost just under 2 marks it and 2 does not; **the measured
  verdicts reproduced** from the reference runs' recorded values and covariances
  (`tmp/snp_stop/runs/{s1,o4,t59}_d1e-4.judge.tsv`, copied into the test as constants): density
  1.35 marked and duplicated 0.12 marked at 59 samples, density no trade and duplicated 0.0055
  marked at 4, both marked at 1 — `dup_trade.py` gives the expected numbers; spec §12 test 12's
  drawn cohorts for both pairs.
  *Identity:* no written output changes — the oracle's checksums and the digests are unchanged.
  *Depends:* C1, C2. *Source:* §2.1.3; §12 test 12.

> Checkpoint C (lead): C changes no output.

### Milestone D — what is reported

- ☐ **D1** The run log: `StandardErrors::described` takes the two reports, prints the verdict in
  place of a marked number's error, and adds each determined share and trade cost (§2.1.3's
  wording).
  *Tests:* the line from a hand-built report, marked and not; the duplicated class's part absent
  when the class is off.
  *Identity:* only the log line changes; checksums and digests unchanged.
  *Depends:* C3. *Source:* §2.1.3, *What is reported*.
- ☐ **D2** **Own commit, do not bundle.** Parameters file format 4: the two tables'
  types (`deny_unknown_fields`), writing with their comments, reading, validation;
  `FORMAT_VERSION` = 4 (the test at `validate.rs:1854` flips to 4 read, 5 refused);
  `parameters_file_of` (`run/census_fit.rs:533`) passes the reports through builders beside
  `with_snp_indel_fit_starts`; a calling run carries both tables back (`cli/calling_run.rs:613`);
  `what_was_fitted` does not count them. Golden files regenerated by `regenerate_the_golden_file`
  and `regenerate_the_written_golden_file`, the diffs read.
  *Tests:* round trip; a version-3 file reads with both tables absent; a version-5 file is refused;
  a marked number with a `standard_error` is refused, naming its key; a non-finite value refused; a
  calling run writes both tables back; `what_was_fitted`'s counts unchanged.
  *Identity:* **every VCF bit-identical** (the oracle's two calling checksums and from_fit). The
  fitted file moves only by its version line and the two tables with their comments: shown by
  deleting them and diffing against B1's file. The defaults run's file moves by its version line
  alone. The digests re-recorded for those lines.
  *Depends:* C3. *Source:* §2.1.3, *What is reported* and its traps; `parameters_file.md` §3.6.

> **Checkpoint D (owner):** the oracle's fitted file (4 samples: the density not marked, the
> duplicated class marked) and a one-sample file (both marked), with their log lines. The owner
> judges whether the tables and the log say what §2.1.3 decided, in words a reader can use.

### Milestone E — the checks

- ☐ **E1** Cross-platform: `cli::cross_platform_digests` and `scripts/promote_ng_oracle.sh` give the
  same checksums natively on macOS and in the Linux container. The oracle baseline is re-recorded
  once, its comment naming what A2, B1 and D2 each moved.
  *Depends:* D2. *Source:* `fit_precision.md` §3.6 item 5; the oracle script's rule.
- ☐ **E2** Calls on the 4 oracle psps: called once with the new fit and once with the current fit
  (`o4_base.toml`), **both by the new build**, so only the parameters file differs; records added and
  lost, sample genotypes changed, QUAL changes (`tmp/snp_conv/compare_vcfs.py`, `qualdiff.py`).
  Expected near the diagnosis's §5.1: 2 records added, 5 genotypes changed of 26,720.
  *Depends:* D2. *Source:* `parameter_prepass_joint_fit.md` §3.3.1, *Measured*.
- ☐ **E3** The same on the 59 tomato2 psps, against `t2_59.parameters.toml`. Not measured before:
  the diagnosis skipped it for cost (two 59-sample calling runs).
  *Depends:* D2. *Source:* as E2.
- ☐ **E4** The stopping rule loses nothing at 59 samples: the new fit's returned values (from
  `PVC_JOINT_FIT_TRACE`'s returned-fit rows) against the reference `t59_d1e-4.judge.tsv`, every
  determined quantity within 0.1 of its reference error (predicted 0.047); B1's reproduction already
  covers 1 and 4 samples. Also the 59-sample verdicts and file tables against §2.1.3's example.
  *Depends:* D2. *Source:* `fit_precision.md` §3.6 item 4.
- ☐ **E5** Timings, on a machine running nothing else: `estimate-parameters` wall time, the
  SNP/indel stage, passes by start and peak resident memory, at 1 sample, the 4 oracle samples and
  the 59 tomato2 samples. Main's binary and the new one at 1 and 4; at 59 the new one against
  main's run of 2026-10-15 (61 min, 18 threads), unless the owner asks for main to be run again.
  *Depends:* D2. *Source:* `parameter_prepass_joint_fit.md` §3.3.1 (the measured table).
- ☐ **E6** The report, `doc/devel/reports/implementations/snp_fit_convergence_E_<date>.md`, with
  E1–E5; the full suite in release; a code review of the branch.
  *Depends:* E1–E5.

> **Checkpoint E (owner):** the report. The owner decides whether E2's and E3's call changes are
> acceptable and whether the branch merges.

## Verification summary

| milestone | proven by |
|---|---|
| A | decrement and verdict unit tests; oracle and digests before and after A2, moves explained |
| B | bit-for-bit reproduction of the prototype's 0.001 runs at 1 and 4 samples; §12 test 11; determinism at 1, 4, 8 threads |
| C | incomplete Beta against SciPy values; covariance against the dense inverse; the measured verdicts reproduced from recorded covariances; §12 test 12; no output moves |
| D | round-trip, version and validation tests; VCFs bit-identical; fitted file differs only by the tables and version |
| E | both platforms agree; calls against the current fit at 4 and 59 samples; values within 0.1 errors of the 10⁻⁴ reference at 59; timings |

## Out of scope (next plans)

The observed information against the outer-product estimator (`fit_precision.md` §8 question 3);
reparameterising the density (joint-fit §11 question 13); a full-scale run above 59 samples.
