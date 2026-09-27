# Fitting to the precision the data support — implementation plan

*Status: plan, 2026-09-27. Nothing here is built yet.*

This plan turns [`doc/devel/ng/spec/fit_precision.md`](../ng/spec/fit_precision.md) into build
order. **It is not a place for new design**: every step cites the spec section it builds, and a
question the spec does not settle goes back to the spec. It follows
[`estimation_memory.md`](estimation_memory.md), whose kimura run
(`tmp/estimation_memory_2169.log`, 2,169 samples) is where both slow steps were measured.

## 1. Scope

**In:**

- **A** — standard errors for the SNP/indel fit (spec §3.2–3.3), first as a diagnostic that changes
  no fitted number;
- **B** — the SNP/indel stopping rule in units of those errors, and stopping a start that agrees with
  an earlier one (spec §2, §3.3–3.5);
- **C** — standard errors for each repeat-tract stratum and the climb's new stopping rule (spec
  §4.2–4.3);
- **D** — each repeat-tract stratum fitted on a subset of samples, grown to a precision target (spec
  §4.4);
- **E** — the standard errors in `Estimate` and in the parameters file, version 2 (spec §5).

**Out**, each with its home in spec §7: weighting the curves by the new errors; errors for blended
slippage numbers; calling reading the errors; running the SNP/indel starts in parallel; the
repeat-tract guard.

## 2. Principles that fixed the order

- **Measure before changing behaviour.** Milestones A and C compute and print errors without using
  them. Only once the errors are shown to mean what they say (their oracle tests) do B and C3 stop by
  them. A wrong error in a stopping rule is a fit that stops too early, and nothing crashes.
- **The error's oracle first.** Every error is checked against something independent: finite
  differences of the total log-likelihood, the full matrix on a small cohort, and the spread of the
  estimates over many drawn cohorts.
- **Isolate the steps that move a fitted number.** B1, B2, C3 and D2 each change results silently;
  each is **its own commit, never bundled**, with the cross-platform checksum's movement explained in
  its commit before the checksum is re-recorded
  ([`cross_platform_digests.rs`](../../../src/cli/cross_platform_digests.rs) module doc: re-record
  only after `scripts/promote_ng_oracle.sh` has measured the change on the real cohort).
- **Types first, then implementation**, within every milestone.
- **The file last.** Its format changes once, after the errors it writes are trusted.

## 3. Preconditions

- Main at or after `7a378ded` (issues 1–3 of the memory plan, the progress-line instrumentation,
  `--skip-contamination`).
- `scripts/promote_ng_oracle.sh` runs in the dev container and its baseline passes.
- For checkpoints B and D: the owner can run `estimate-parameters` on kimura (the 2,169-sample
  cohort), as for the memory plan.
- Spec question 1 is decided (owner, 2026-09-27): a parameter with no information is written as
  `defaulted` with the default value, built in E2.

## 4. The steps

### Milestone A — the SNP/indel fit's standard errors, as a diagnostic

1. ✅ **A1 — the per-position scores.** For each parameter kind of spec §3.2's table, the observed-data
   score at one position (the posterior expectation of the complete-data score). Unit-tested against
   a central finite difference of the total log-likelihood on a small drawn cohort, one test per
   kind. *Depends:* —. *Source:* spec §3.2.
2. ✅ **A2 — the block information.** `SampleInformationBlock` ((k + 8)² per sample) and the cohort's
   8 × 8 block, accumulated in an E-pass when asked for, in the fit's fixed chunk order. The fit's
   iterating passes do not ask; nothing moves. *Depends:* A1. *Source:* spec §3.2 (the block
   approximation), §3.6 item 5.
3. ✅ **A3 — errors from the blocks.** The Schur-complement errors for each sample's parameters and
   the cohort's; `None` for a parameter with no information. Tested against the full outer-product
   matrix on drawn cohorts of 4 and 20 samples (the full matrix is small there), reporting the
   difference per parameter kind. *Depends:* A2. *Source:* spec §3.2, §3.6 item 2.
4. ✅ **A4 — the errors mean what they say.** A test (ignored by default for its length, run at the
   checkpoint) drawing 200 cohorts from known parameters at 3 and at 30 reads a position, reporting
   the share of estimates within one and two errors of the truth, per parameter kind. *Depends:* A3.
   *Source:* spec §3.6 item 3.
5. ☐ **A5 — compute and print them at the end of a fit.** The final pass accumulates the blocks; the
   fit prints, per parameter kind, the median and largest error, and the count with none. No fitted
   number moves; **the checksums pass unchanged.** *Depends:* A3. *Source:* spec §3.3 (the final
   pass), §3.5.

> **Checkpoint A — pause for review.** Report A3's full-matrix comparison and A4's coverage, and the
> errors printed on the 4-accession oracle cohort. **Decide `SETTLED_FRACTION`** (spec open question
> 2) from how far the current fit's last cycles move each parameter in units of its error. If A3 or A4
> shows the outer-product errors are not trustworthy for some kind, that goes back to spec question 3
> before B starts.

### Milestone B — the SNP/indel fit stops by its errors

1. ☐ **B1 — the settled test replaces the relative-move rule.** The log-likelihood trigger, the
   information pass when it first holds and every `ERROR_REFRESH_CYCLES` after, the per-parameter
   projection with `MAX_CONTRACTION`, and convergence when every parameter is settled.
   `largest_relative_move` no longer stops the fit. **Own commit; moves fitted numbers.** Validated by
   spec §3.6 item 4: the oracle cohort fitted under the new rule and run to 1,000 passes agree within
   `SETTLED_FRACTION` of each error; the passes saved are reported. The checksums are re-recorded in
   this commit, with the explanation. *Depends:* checkpoint A. *Source:* spec §2, §3.3.
2. ☐ **B2 — a later start stops when it agrees.** The agreement test on projected endpoints against
   the best converged earlier start. **Own commit; moves passes, and moves numbers only where a start
   that would have won by a hair now stops.** *Depends:* B1. *Source:* spec §3.4.
3. ☐ **B3 — reporting.** Per start: converged, at the limit, or agreed; and when not converged, the
   parameter furthest from settled with its distance in errors. The per-cycle progress line gains the
   count not yet settled. `JointFit` gains the per-start record. *Depends:* B2. *Source:* spec §3.5.

> **Checkpoint B — pause for review.** The owner reruns the SNP/indel half on kimura
> (`estimate-parameters` with the new build; the repeat-tract half may be interrupted once the SNP/indel
> fit reports done). Report passes per start and wall time against 9 h 25 min, and the parameters
> against the previous run's in units of their errors.

### Milestone C — each stratum's errors and the climb's stopping rule

1. ☐ **C1 — a stratum's curvature errors.** The central-difference curvature of the stratum's total
   log-likelihood over all its numbers, at its final answer, on the climb's scales, carried to the
   natural scale; `None` where the curvature is not negative-definite in that direction. Computed and
   carried in `StratumFit`, printed as a summary; **no number moves.** *Depends:* —. *Source:* spec §4.2.
2. ☐ **C2 — the stratum errors mean what they say.** A test (ignored by default) fitting strata drawn
   at known slippage many times, at 3 and 30 reads, comparing spread with reported error. *Depends:*
   C1. *Source:* spec §4.5 item 1.
3. ☐ **C3 — the climb stops on its projected remaining total gain**, and a round that loses is not
   convergence: its moves are undone and the walk stops at its best point. **Own commit; moves fitted
   numbers.** Validated against a climb run to 20 rounds on drawn strata and on the oracle cohort
   (spec §4.5 item 2); checksums re-recorded in this commit with the explanation. *Depends:* C2,
   checkpoint A's `SETTLED_FRACTION`. *Source:* spec §4.3.

> **Checkpoint C — pause for review.** C2's coverage, C3's comparison against 20 rounds, and the
> rounds saved on the oracle cohort.

### Milestone D — a stratum read from a subset of samples

1. ☐ **D1 — the sample order.** A pure function from the cohort's sample names and a fixed seed to
   their order, tested to be the same whatever order the names arrive in. *Depends:* —. *Source:*
   spec §4.4 (the order).
2. ☐ **D2 — the growing subset.** Fit on the first `FIRST_SUBSET` samples, add slippage groups'
   samples to `MIN_SAMPLES_A_GROUP`, double until the level's relative error is below
   `LEVEL_RELATIVE_ERROR_TARGET` or every sample is in, each climb warm-started from the last; the
   refusal floor judged on the whole stratum. Evidence counts from the subset. **Own commit.** At 256
   samples or fewer it must be byte-identical to C3's result: **the checksums pass unchanged.**
   *Depends:* D1, C1, C3. *Source:* spec §4.4.
3. ☐ **D3 — the comparison tool.** An example fitting chosen strata of a cohort both on every sample
   and on the grown subset, printing each stratum's three numbers, their errors, the subset size
   reached and the time, for the owner to run on kimura. *Depends:* D2. *Source:* spec §4.5 item 3.

> **Checkpoint D — pause for review.** The owner runs D3 on kimura for about five strata, chosen to
> include the largest. **Decide `FIRST_SUBSET`, `LEVEL_RELATIVE_ERROR_TARGET` and spec question 5**
> (which number carries the target) from it.

### Milestone E — the errors in the file

1. ☐ **E1 — `Estimate` gains `standard_error: Option<f64>`**, its "no uncertainty interval" note
   rewritten (spec §5.1); every constructor sets it, `None` where no error is computed. No output
   changes. *Depends:* B3, C1. *Source:* spec §5.1.
2. ☐ **E2 — the parameters file, version 2.** `standard_error` in the value tables and
   `own_fit_standard_error` in the slippage origin blocks (spec §5.2's table), each start's outcome
   in `fitted_from`, `FORMAT_VERSION` = 2, the reader accepting versions 1 and 2, the golden files
   updated; a parameter with no information written as `defaulted` with the default value (spec
   question 1). The fitted-parameters checksum is re-recorded (the key and the version move it; no fitted
   number does — the commit shows the file's diff). *Depends:* E1.
   *Source:* spec §5.2–5.3.

> **Checkpoint E — pause for review.** The written file on the oracle cohort, and the owner's full
> kimura run with the finished build: wall time of each half against the kimura log, and the file.

## 5. Verification summary

| milestone | proven by |
|---|---|
| A | scores against finite differences; block errors against the full matrix at 4 and 20 samples; coverage over 200 drawn cohorts at 3 and 30 reads; checksums unchanged |
| B | the oracle cohort under the new rule against 1,000 passes, within `SETTLED_FRACTION` of each error; kimura rerun (owner) |
| C | coverage over drawn strata; the new stop against 20 rounds; checksums re-recorded with explanation |
| D | order independence; byte-identical at ≤ 256 samples; full against subset on kimura strata (owner) |
| E | golden files; version-1 files still read; absent errors stay absent through a round trip |

## 6. Out of scope (next plans)

Spec §7, each with its recommended home: the curves weighted by the new errors and errors for blended
numbers (an amendment to `str_slippage_level_curve.md`); calling reading the errors; parallel
SNP/indel starts (after checkpoint B measures the passes a start takes); the repeat-tract guard (its
own investigation).
