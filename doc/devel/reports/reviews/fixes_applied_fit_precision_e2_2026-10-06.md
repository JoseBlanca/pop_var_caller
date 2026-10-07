# Fix Application Report: fit_precision_e2_2026-10-06.md

**Date:** 2026-10-06
**Source review:** `doc/devel/reports/reviews/fit_precision_e2_2026-10-06.md`
**Source state reviewed against:** review object `c2f32f0c` (parent `25af661c`, branch `fit-precision`)
**Execution mode:** non-interactive
**Overall status:** Completed, one Major deferred to the owner

---

## 1. Executive summary

- Review: 5 Major, 11 Minor, nits.
- Applied: M1, M2, M4, M5, Mi1–Mi11 (Mi2 but for a version-1 file carrying version-2 keys).
- Deferred: M3 (the owner's decision), Mi2's version-1 case, most nits.
- Validation: `cargo fmt --check` and `cargo clippy --all-targets --all-features -- -D warnings` clean; the touched
  modules' release tests pass; the full suite in the step's implementation report. No hot path changed.

### Unresolved high-priority findings
- M3 — the portability digest does not pin the new errors; for the owner at checkpoint E.

## 2. Per-finding log

### M1 — the excess fitted across samples no read informs — Applied
`fits_homozygote_excess` now takes the count of samples that inform their excess: at the fit, those with a read at an
ordinary position (`samples_with_reads`); in the error code, those whose own excess diagonal carries information
(`samples_informing_their_excess`, on the blocks and on the whole matrix), at all six places that read it. Tests:
`one_sample_with_reads_beside_a_tracts_only_sample_holds_its_excess` (fit), `one_sample_informing_its_excess_holds_it`
(the error code). `a_diagonal_that_is_not_finite_is_no_information` moved from two samples to three: with the second
sample's excess uninformed, two samples kept its old meaning only at three.

### M2 — empty totals for every read group — Applied
Empty totals only for a read group whose rate is `Defaulted` from no observation; any other rate without totals still
stops assembly. The run-level test gives its four drawn samples real totals (a thousand reads at Q30 each) and now
asserts their multipliers are fitted, with errors, beside the defaulted fifth.

### M3 — the portability digest pins no error — Deferred
Needs a fixture whose rates leave their bounds and whose tracts mismatch, and a re-record on both platforms. Put to the
owner at checkpoint E with a recommendation.

### M4 — the starts a calling run read — Applied with adaptation
`TheRunsNumbers::parameters_file` is now the one door both calling commands write through; it adds the starts. Tests:
`a_run_writes_back_the_starts_its_file_carried` (the reviewer's, through the whole psp command) and
`a_calling_run_writes_back_the_fit_starts_it_read` (the door). `of_run`'s doc says it writes none.

### M5 — spaces in two refusals — Applied
The continuations restored; `refused` asserts no message carries a run of spaces, which covers every refusal in the
module.

### Minors — Applied
- Mi1: `agreed_with_start = 0` refused, tested.
- Mi2: refused, each tested — a start agreeing with one that did not converge, an empty list, `samples_fitted_on`
  above the cohort's size, an error on the outlier weight or the fallback concentration. A version-1 file carrying
  version-2 keys is still accepted (deferred: it reads correctly, and the reviewer ranked it last).
- Mi3: the subset sentence applies where `samples_fitted_on` is below the cohort's count, says equal means every sample,
  and gives checkpoint D's figure in the unit it was measured in. The fixture's `samples_fitted_on` went from 512 to 1
  (of its two samples), which the new bound requires.
- Mi4: the header names the zero-rate case and whose error a `supplied` number carries.
- Mi5: every new refusal says how to fix it; the `samples_fitted_on = 0` fix corrected.
- Mi6: both new tests moved above the old tests' doc comments.
- Mi7: the file's type renamed `StartOutcome`, a noun distinct from the fit's `StartEnding`; added to both spelling
  tests (twenty-five variants, seventeen in the golden file).
- Mi8: `ReadsBehindEachCalibration::no_count_in_reads`, documented for the fit's case, used by `parameters_file_of`.
- Mi9: the three stale docs corrected.
- Mi10: the comments say the fit's steps move an uninformed excess, so the default value is written; two figures the
  reviewer measured were not copied into the code.
- Mi11: the doc says a number supplied from a file keeps that run's error.

### Nits
Applied: the `StratumOutcome` match names its variants; the baseline note says what the comparison removed instead of
citing a scratch script; one step label in code replaced by a date. Not changed: the remaining step labels in test docs,
the `expect` notes, `SETTLED_FRACTION`'s "tenth" in the notes, inline paths in tests.

## 3. Checksums after the fixes

The notes changed again, so the cross-platform fitted-parameters checksum moved from `f03c0ab1…` to
`a5e87936e968f9b855f91f88db828083`: the two files differ in 19 note lines and nothing else (calls unchanged), and the
first was proved against the old `22dda760…`. The oracle was rerun on the fixed code for its baseline (the step's
implementation report).
