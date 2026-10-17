# Correctness review: SNP/indel fit convergence, plan steps B2–B5 (2026-10-17)

*Written by the reviewing agent; filed by the lead (the harness refused the agent's file write).*

**Scope.** Commits `06f57ef4..f2ef6574` on `snp-fit-convergence`, reviewed as committed (an export of
`f2ef6574`).
- **B2:** the box solve rewritten as the primal active-set method (`bounded_step`); a solve that does
  not finish reports a decrement that is not a number (NaN).
- **B3:** the swap of a library's two error rates deleted, plus a run-log note naming libraries whose
  rates end inverted.
- **B4:** the decrement measured only over parameters that have a standard error; the others profiled
  out through the arrow (`determined_decrement`).
- **B5:** the tests, and the pass-limit fix in `maximise`.

Judged against `fit_precision.md` §2 (amendment of 2026-10-17) and `parameter_prepass_joint_fit.md` §3.3.

**Method.** Static reading, then four experiments on the export: a stress test of the box solve (6,000
random cases of 1 to 8 samples; information from fewer positions than parameters, diagonal floors 0,
10⁻¹⁵, 10⁻¹³, 10⁻¹⁰, condition numbers to about 10¹⁰, some exactly singular; slopes scaled 10⁻⁴ to 10⁴;
values on an end, just inside one, or interior); a round-by-round trace of unfinished solves; timing at
59, 500 and 2,000 samples; 13 mutations against the fit module's 135 tests and the 4 digest tests.

**Baseline at `f2ef6574`:** `parameter_estimation::joint::fit` 134 passed, 1 ignored; `cross_platform_digests` 4 passed.

## Findings, ranked

| # | Severity | Finding |
|---|---|---|
| 1 | Major | The box solve can cycle on ill-conditioned information; it then hits its round limit, the decrement is NaN, and the start cannot stop on it. |
| 2 | Major | Rounding can make the determined decrement slightly negative; `has_converged` refuses any negative value, so a converged start keeps running. |
| 3 | Major | No test checks that an unfinished solve reports NaN; a mutation removing it passes every test. |
| 4 | Minor | The solve spends one round on each parameter meeting an end at its own distance; cost grows with the square of the sample count (1.85 s a judging pass at 2,000 samples in the bad case). |
| 5 | Minor | Three parts of the box solve are untested: freeing one parameter at a time, the 10⁻⁹ rounding threshold, clamping a negative room to zero. |
| 6 | Minor | B4 is pinned only by the digest checksums. |
| 7 | Minor | On a library whose two rates end inverted, calling uses the higher of the two; the run-log note does not say so. |
| 8 | Minor | `newton_step` takes the standard errors separately from the sums they must come from, tied only by a length check. |

**What holds:** in all 6,000 stress cases the solve stays inside the box and its step is worth at least
nothing; `has_converged` refuses NaN; the pass-limit fix works and its tests catch both ways of breaking
it; errors and step come from the same pass's sums at the same point; with every moved parameter
determined, the old formula is used bit for bit; the new code is serial and the pool-width digest test
passes.

### 1. `standard_errors.rs:808-827` — [Major] the box solve can cycle

The method relies on a freed parameter moving back into its interval on the next solve. That holds in
exact arithmetic when the reduced information is positive definite; it fails when the solve over the
free parameters is nearly singular, or when `invert_identified` keeps a different set once the freed
one joins. The freed parameter's step then points out; its room is clamped to 0 (line 779), the reach is
0, it is fixed again, and the next whole step frees it again. Traced on stress case 203: parameter 9
freed with an inward slope of 0.031, alternating from round 21 to 199 with the model's value constant to
10⁻¹¹. Cases 18 and 28 loop the same way; case 573 progresses by reaches of 4×10⁻⁸ and 2×10⁻⁹ until the
limit. 28 of 6,000 cases did not finish, all with floor 10⁻¹⁰ (2.3% of that group). A cohort landing
there runs to `max_passes`. **Fix:** block a parameter that will not move inwards once freed (fix it
again, never free it again within that solve; the blocked set only grows, so the solve terminates); add
case 203 as a regression test asserting the solve finishes.

### 2. `standard_errors.rs:723`, `settled.rs:69` — [Major] a rounding-negative determined decrement never converges

`½ mᵀ I m` with `m = (d_K, −I_UU⁻¹ I_UK d_K)` is non-negative in exact arithmetic but a small difference
of large terms in floating point. 125 of 3,546 stress cases using the profiled formula came out negative,
all on exactly singular information (−1.2×10⁻¹⁶ against ½dᵀId 1.0×10⁻²; −1.4×10⁻¹⁰ against 2.3×10⁴;
−3.8×10⁻²⁴ against 8.3×10⁻¹⁰). They appear exactly where the start should stop. **Fix:** in
`determined_decrement` only, return `length.max(0.0)` (a negative there can only be rounding); keep the
NaN for an unfinished solve and the `>= 0.0` guard; test on an exactly singular arrow.

### 3. `standard_errors.rs:608-613` — [Major] the NaN for an unfinished solve is untested

Mutation S4 (compute the decrement even when the solve did not finish) passes all tests. **Fix:** use
case 203 as a fixture asserting `decrement.is_nan()` and `!has_converged`; after finding 1's fix, reach
the limit through a round count settable in tests.

### 4. `standard_errors.rs:755, 793-803` — [Minor] quadratic cost when many parameters sit just inside an end

Each round fixes only the parameters whose room equals the shortest exactly, and the golden section
leaves values beside an end (a gap of 4.4×10⁻¹⁰), so each outbred sample at its excess floor costs a
round. Measured: 59 samples 61 µs → 1.7 ms; 500: 0.30 → 115 ms; 2,000: 1.3 ms → 1.85 s. **Fix:** treat a
value within a relative 10⁻⁹ of the width from an end as on it before the first round, or fix every
parameter whose room is within `reach × (1 + 1e-12)` of the shortest; pin with the timing case.

### 5. `standard_errors.rs:779, 816-817` — [Minor] three parts survive mutation

S1 (free all inward-pointing fixed parameters at once), S3 (`SLOPE_POINTING_IN` = 0), S13 (no
`.max(0.0)` on the room) pass every test. The ill-conditioned fixtures of findings 1 and 2 would catch
them; also a value exactly on an end with a rounding-level slope, asserting fewer rounds than parameters.

### 6. `standard_errors.rs:597-610` — [Minor] B4 pinned only by checksums

S7 (every parameter given an error) passes all fit tests; only the digests fail. **Fix:** a fit-level
test on the two-sample cohort with sequencing errors (12 of 16 parameters without error): every start
ends converged in under 50 passes.

### 7. `fit.rs:2876-2894`, `run/census_fit.rs:358` — [Minor] calling's rate on an inverted library

Calling reads the ordinary-position rate as the library's sequencing error rate
(`JointFit::sequencing_error_rates`, `fit.rs:456`); on an inverted library that is the higher rate. On the
21-sample digest fixture 12 of 22 libraries end inverted. **Fix:** say it in the note; record in §3.3
that it is accepted.

### 8. `standard_errors.rs:587-607` — [Minor] errors from other sums

Only an `assert_eq!` on length ties `errors` to `sums`. **Fix:** a `debug_assert_eq!` against
`StandardErrors::of(sums)`, or pass them together.

## Nits

- The determined decrement can exceed `½ dᵀ I d` when the solve over the profiled set drops a parameter
  the step moved (12 of 3,546 cases, by up to 3.3×, always ≤ 2.5×10⁻⁵): harmless at 0.001, but the doc's
  "the smallest that length can be" overstates.
- Accepted by the owner, worth keeping in view: with no parameter determined, the decrement is 0 and a
  start stops at its first judging pass (the no-errors fixture: 16 passes, −44.871 against −44.699 at 41).
  At one low-depth sample, when a start stops then depends only on `log_likelihood_stillness`.

## Mutation results

| id | mutation | caught by |
|---|---|---|
| S1 | free all inward-pointing fixed parameters at once | nothing |
| S2 | fix only when room < reach | 7 tests |
| S3 | `SLOPE_POINTING_IN` = 0 | nothing |
| S4 | no NaN for an unfinished solve | nothing |
| S5 | hold profiled parameters still | `the_determined_decrement_is_the_step_through_the_determined_covariance` |
| S6 | leave error-carrying, unmoved parameters out of the profiled set | same |
| S7 | B4 switched off | 3 digest tests only |
| S8 | free the first candidate, not the one pointing furthest in | the property test |
| S9 | restore the error-rate swap | `plain_passes_never_lower_the_log_likelihood`, `every_judged_decrement_is_at_least_zero` |
| S10 | drop `accepted_pass.is_some() \|\|` | `an_accepted_trys_pass_is_judged_at_the_pass_limit` |
| S11 | break before SQUAREM only after refused tries | same |
| S12 | `has_converged` without `>= 0.0` | `settled`'s own test |
| S13 | room without `.max(0.0)` | nothing |

Scratch: `tmp/review_B/` (export with probes and `REVIEW_MUT`/`REVIEW_TRACE` switches, logs, traces).
