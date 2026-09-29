# Fitting to the precision the data support — implementation plan

*Status: plan, 2026-09-27; Milestone A built. **Amended at checkpoint A (owner, 2026-09-28)**: three
steps added before Milestone B (§4, Milestone A′), and the checkpoint's decisions recorded under it.*

This plan turns [`doc/devel/ng/spec/fit_precision.md`](../ng/spec/fit_precision.md) into build
order. **It is not a place for new design**: every step cites the spec section it builds, and a
question the spec does not settle goes back to the spec. It follows
[`estimation_memory.md`](estimation_memory.md), whose kimura run
(`tmp/estimation_memory_2169.log`, 2,169 samples) is where both slow steps were measured.

## 1. Scope

**In:**

- **A** — standard errors for the SNP/indel fit (spec §3.2–3.3), first as a diagnostic that changes
  no fitted number;
- **A′** — three corrections checkpoint A found, made before the errors are used to stop the fit:
  each library's reads scored with its own error rates, the allele-frequency shapes' update made to
  maximise the likelihood the fit computes, and the full information matrix at small cohorts;
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
- **Isolate the steps that move a fitted number.** A6, A7, B1, B2, C3 and D2 each change results silently;
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
5. ✅ **A5 — compute and print them at the end of a fit.** The final pass accumulates the blocks; the
   fit prints, per parameter kind, the median and largest error, and the count with none. No fitted
   number moves; **the checksums pass unchanged.** *Depends:* A3. *Source:* spec §3.3 (the final
   pass), §3.5.

> **Checkpoint A — pause for review.** Report A3's full-matrix comparison and A4's coverage, and the
> errors printed on the 4-accession oracle cohort. **Decide `SETTLED_FRACTION`** (spec open question
> 2) from how far the current fit's last cycles move each parameter in units of its error. If A3 or A4
> shows the outer-product errors are not trustworthy for some kind, that goes back to spec question 3
> before B starts.

**Decided at checkpoint A (owner, 2026-09-28).** The measurements are in the
[A4](../reports/implementations/fit_precision_a4_2026-09-27.md) and
[A5](../reports/implementations/fit_precision_a5_2026-09-27.md) reports.

1. **Each library's reads are scored with that library's own error rates** (step A6). A sample's
   genotype stays shared across its libraries; each library's reads are weighed by its own rates.
   Today the fit adds a sample's libraries together and scores them all with the first library's
   rates, then gives every library that same fitted number — though two libraries are separate
   experiments, and one's error rate says nothing about another's, any more than a library of another
   sample would. The copied rates are not written with a shared error: this replaces the proposal to
   do so.
2. **The allele-frequency shapes' update is fixed before the stopping rule** (step A7). It solves
   the digamma form, not the slope of the likelihood the fit computes, so the fit stops beside the
   likelihood's maximum: on drawn cohorts the shapes and the invariant share sit 0.4 to 3.1 errors
   above the truth, and the gap, in errors, grows as the square root of the positions.
3. **Spec question 3: the outer-product estimator stands; at small cohorts the full matrix replaces
   the blocks** (step A8), up to 20 samples. At 4 samples and 3 reads the estimates scatter 1.23 to
   1.67 times the blocks' errors and 0.89 to 1.02 times the full matrix's; at 20 samples the two are
   about 10% apart.
4. **`SETTLED_FRACTION` = 0.1**, the spec's starting value. The oracle cohort could not decide it —
   its fit stops at the 200-pass limit — so step B1's comparison against 1,000 passes is its check.

### Milestone A′ — the corrections checkpoint A found

Before the errors are used to stop the fit, the fit must estimate what the errors describe. Each
step here that moves a fitted number is its own commit, with the checksums re-recorded after
`scripts/promote_ng_oracle.sh` has measured the change.

1. ✅ **A6 — each library's reads under its own error rates.** At a position, a sample's likelihood
   given its genotype becomes the product over its libraries of each library's reads under that
   library's two rates; the genotype, and everything above it, stays one a sample. Today every stage
   assumes one set of counts and one rate pair a sample, so the step reaches each of them:
   - the per-position reader keeps each library's counts and depth range apart, where it now adds
     them (`EvidenceCursor::next_position`);
   - the likelihood reads each library's own rates (`one_position`, which now reads
     `group_index[s][0]`'s);
   - each library's read tallies are credited with its own reads only, so the maximisation fits
     each library's rates from its own evidence (today every library of a sample is credited with
     the sample's pooled reads, which is why the fitted rates come out identical);
   - the scores, the information and the errors of steps A1–A3 carry each library's two rates: a
     sample with k libraries has 1 + 2k parameters of its own, as spec §3.2 sizes it, where they now
     carry three whatever k. The test that pins a second library's rates as having no slope
     (`a_second_read_groups_rates_have_no_slope`) is replaced by one that finite-differences them.

   A library with no reads has no information and is written `defaulted` at E2, per spec question 1;
   the count of "read groups after a sample's first" leaves the log. **Own commit; moves fitted
   numbers for samples with more than one library only** (up to 482 of kimura's 2,651 read groups).
   The oracle cohort has one library a sample (4 read groups over 4 samples), so its checksums are
   expected not to move, which is checked; the cross-platform fixture has two read groups in one
   sample (`cli::test_fixtures`), so its checksums are expected to move, and are re-recorded in this
   commit with the explanation. Validated on drawn cohorts whose samples carry two libraries at
   different rates: each library's rate recovered within its errors, and the
   scores against finite differences per library. Work per position grows with the libraries rather
   than the samples (2,651 against 2,169 on kimura). *If the step is larger than one reviewed commit,
   it splits at the scores: the likelihood first, with a multi-library sample's errors reported
   absent until the second part.* *Depends:* A5. *Source:* checkpoint A decision 1; spec §3.2's
   per-library block.
2. ✅ **A7 — the allele-frequency shapes maximise the likelihood the fit computes.** The update of
   the density's two Beta shapes, and of the carrier Beta's, which uses the same form
   (`fit_beta_shapes`), is replaced by one whose fixed point is where the likelihood's own slope in
   the shapes is zero: Newton steps on the slopes the quadrature rule gives (`RuleSlopes`, step A1),
   kept within the shapes' bounds and never lowering the log-likelihood. More quadrature nodes are not
   the fix: the likelihood itself is already exact at 16 nodes, and at 96 the gap was still 0.16 to
   0.88 errors. **Own commit; moves fitted numbers.** Validated by A4's coverage test: the flat-point
   column near zero for the shapes and the invariant share, their coverage near nominal at 20
   samples, and the test's bounds re-recorded from the new run; the oracle's change measured and
   explained before the checksums are re-recorded. *Depends:* A6. *Source:* checkpoint A decision 2.
3. ✅ **A8 — the full matrix at small cohorts.** When a cohort has at most `FULL_MATRIX_SAMPLES` = 20
   samples, the errors come from the full outer-product matrix — every parameter paired with every
   other, two samples' included — and from the blocks above that. No fitted number moves. Validated by
   A4's coverage test, whose 2- and 4-sample regimes at 3 reads then report the full matrix's
   coverage as the fit's, and by the cost of the final pass at 20 samples, measured. *Depends:* A6.
   *Source:* checkpoint A decision 3.

> **Checkpoint A′ — pause for review.** A6's and A7's changes on the oracle cohort, A4's coverage
> re-run, and the kimura run of the A5 build if the owner has made it.

### Milestone B — the SNP/indel fit stops by its errors

1. ☐ **B1 — the settled test replaces the relative-move rule.** The log-likelihood trigger, the
   information pass when it first holds and every `ERROR_REFRESH_CYCLES` after, the per-parameter
   projection with `MAX_CONTRACTION`, and convergence when every parameter is settled.
   `largest_relative_move` no longer stops the fit. **Own commit; moves fitted numbers.** Validated by
   spec §3.6 item 4: the oracle cohort fitted under the new rule and run to 1,000 passes agree within
   `SETTLED_FRACTION` of each error; the passes saved are reported. The checksums are re-recorded in
   this commit, with the explanation. `SETTLED_FRACTION` = 0.1 (checkpoint A). *Depends:* A8,
   checkpoint A′. *Source:* spec §2, §3.3.
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
   question 1) — after A6, a library with no reads; no library's rate is copied from another's.
   The fitted-parameters checksum is re-recorded (the key and the version move it; no fitted
   number does — the commit shows the file's diff). *Depends:* E1.
   *Source:* spec §5.2–5.3.

> **Checkpoint E — pause for review.** The written file on the oracle cohort, and the owner's full
> kimura run with the finished build: wall time of each half against the kimura log, and the file.

## 5. Verification summary

| milestone | proven by |
|---|---|
| A | scores against finite differences; block errors against the full matrix at 4 and 20 samples; coverage over 200 drawn cohorts at 3 and 30 reads; checksums unchanged |
| A′ | each library's rate recovered on drawn two-library cohorts, and its scores against finite differences; the shapes' flat-point column near zero and their coverage re-run; the full matrix's coverage at 2 and 4 samples; each moved checksum measured and explained |
| B | the oracle cohort under the new rule against 1,000 passes, within `SETTLED_FRACTION` of each error; kimura rerun (owner) |
| C | coverage over drawn strata; the new stop against 20 rounds; checksums re-recorded with explanation |
| D | order independence; byte-identical at ≤ 256 samples; full against subset on kimura strata (owner) |
| E | golden files; version-1 files still read; absent errors stay absent through a round trip |

## 6. Out of scope (next plans)

Spec §7, each with its recommended home: the curves weighted by the new errors and errors for blended
numbers (an amendment to `str_slippage_level_curve.md`); calling reading the errors; parallel
SNP/indel starts (after checkpoint B measures the passes a start takes); the repeat-tract guard (its
own investigation).
