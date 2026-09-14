# Code Review: portable_float_B5
**Date:** 2026-09-14
**Reviewer:** rust-code-review skill (orchestrator, single-pass — see §1)
**Scope:** uncommitted diff of plan step B5 (`doc/devel/implementation_plans/portable_float.md`, Milestone B): std transcendental calls in `src/parameter_estimation/` replaced by `crate::float`
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** `git diff HEAD -- src/` in the worktree `pop_var_caller-portable-float`, on
  top of `ce0ea2f4` (B4). 8 files, +127 / −107.
- **Reviewed against:** `HEAD` = `ce0ea2f4` plus uncommitted changes.
- **In-scope files:** `src/parameter_estimation/{calibration, depth_bins}.rs`,
  `src/parameter_estimation/joint/{census_moments, contamination, fit, share_curve, slippage_curve,
  ssr_fit}.rs`; the callee `src/float.rs`; libm 0.2.16's `log` (through B4's Python transcription).
- **Deliberately out of scope:** the other modules' conversions (B2–B4) and the points their reviews
  settled: prose that names `.ln()`, exact comparisons of a value against libm's own output, the
  flat-gap constant written as bits.
- **Categories:** as for B2–B4, the fan-out of mutation-testing agents was **not** run. The dispatcher
  forbade `cargo` and source edits, which removes what worktree isolation is for. One pass covered
  `reliability` (equivalence of every hunk, test challenge), `refactor_safety` (completeness,
  imports under every configuration), `extras` (hot path, what the parity oracle can and cannot
  show) and `naming`/doc accuracy. `unsafe_concurrency`, `module_structure`, `tooling`, `defaults`:
  no triggers.

## 2. Verdict

**Approve-with-changes.** All 124 converted calls compute what they computed before, with the
operands in the same order and nothing duplicated. No std transcendental call remains anywhere in
`src/` outside `src/float.rs`'s own tests, which compare against std on purpose.

One site will probably keep the fit from reproducing A3's libm-build parameters file (`3f5e5e7f`)
byte for byte. `one_position` takes `ln 3` (`fit.rs:2183`). The old code's `(3 as f64).ln()` had a
literal argument, so the compiler computed it while building, with the build machine's library;
the libm build could not redirect that. glibc and Apple both return the correctly rounded
`0x3ff193ea7aad030b`. libm returns `…030a`, one unit lower, and the converted code now uses libm's
value (M1). It is the only such site in the fit's shipped code. If the parity check fails, this is
where to look first. A one-line experiment settles it (§4 Q1).

## 3. Execution status

- **Commands run by the reviewer:**
  - Completeness sweeps over all of `src/`: method calls
    `\.(ln|exp|powf|powi|log10|log2|log|ln_1p|exp_m1|sin|cos|tan|sinh|cosh|tanh|asin|acos|atan|atan2|cbrt|hypot|exp2|asinh|acosh|atanh)\s*\(`,
    chains that start a line with one of those methods (`^\s*\.(…)\b`), and `(f64|f32)::(…)` paths.
    **Hits:** `src/float.rs:200-281`, which is the module's own std comparison tests (intended). The
    other hits were `self.gamma(…)` in three random-number fixtures and four doc comments naming
    `f64::ln`/`f64::powi` (`alignment/emission.rs:162-163`, `float.rs:3`, `:49`,
    `calling/likelihood/generic.rs:2103`). **No std maths call in `src/run/`, `src/cli/`, or
    anywhere else reached by the fit.** No dependency that does maths (statrs, rand_distr, …) is in
    `Cargo.toml`.
  - Counts over the diff. Added: `float::ln` 64, `exp` 33, `powi` 14, `powf` 10, `cos` 3, `sin` 1
    = 125. Removed: `.ln` 64, `.exp` 32, `f64::exp` 1, `.powi` 14, `.powf` 10, `.cos` 3, `.sin` 1 = 125.
    One `ln` on each side is the doc comment at `fit.rs:2259`. **124 code calls: CHECKED-CORRECT.**
    The plan's "125 calls" counts that comment.
  - `uv run --no-project --with mpmath python tmp/review_2026-09-14_portable_float_B5/literal_sites.py`.
    It checks the literal-argument sites against libm 0.2.16's `log` (B4's transcription), the
    platform's, and a correctly rounded value, plus the distance of each depth-ladder rung from a
    rounding boundary. Verbatim, abridged:
    ```
    ln(3): libm 3ff193ea7aad030a apple 3ff193ea7aad030b correctly-rounded 3ff193ea7aad030b
    ln(TAU): libm 3ffd67f1c864beb4 apple 3ffd67f1c864beb4 correctly-rounded 3ffd67f1c864beb4
    ratio apple 3ff486fd85bcf246 cr-of-double-exponent 3ff486fd85bcf246
      8.0*r^1 = 10.263652972521538  distance from .5 = 0.236 (1.05e+14 relative ulps)
      [...]
      124.0*r^7 = 709.4233733993999  distance from .5 = 0.0766 (4.91e+11 relative ulps)
      [...]
    ```
    Log: `tmp/review_2026-09-14_portable_float_B5/literal_sites.log`.
- **Gates reported by the dispatcher (not re-run):** fmt clean; clippy `-D warnings` clean with
  `--all-targets --all-features` and with `--lib` only; container suite 4,796 passed, 0 failed. The
  parity check of the fit against A3's libm build and the identity oracle were still running.
- **Commands not run:** all `cargo`, by instruction; the macOS suite.
- **Needs verification:** 1 (Q1: whether `ln 3` actually moves the parameters file).

## 4. Open questions and assumptions

1. **Does the one-unit `ln 3` change the parameters file?** Assumption: LLVM folds
   `llvm.log.f64(3.0)` with the host library, which it does for `log`, `exp` and `pow` intrinsics
   with constant operands. A3's own report says the same (§3.1 item 3). The A3 libm fits were
   identical on macOS and Linux, which fits both hosts folding to `…030b`. **Experiment:** if the
   parity fit mismatches, set `let ln_three = f64::from_bits(0x3ff193ea7aad030b);` temporarily and
   rerun `tmp/parity_oracle.sh fit`. A match proves the whole difference is this site. Affects M1.
2. **Do the fitted parameters keep libm's `ln 3`, or the old bits?** Plan §8's rule after B2 is that a
   constant is written as its old bits only when converting it moves a *recorded answer*. The parity
   file is an oracle, not a recorded test answer. **Recommendation: keep libm's value, hoist it (M1),
   and record the experiment's result in the step's commit message as the explained difference.**
   The owner decides. Affects M1.
3. **`estimate-contamination` has no parity check.** `contamination.rs`'s shipped `ln`/`exp` calls
   move its output the same way, but the parity oracle only covers `estimate-parameters` and
   `call-from-psps`. A3's libm build would have produced the reference; B7 can compare against it if
   that output was kept. Affects Mi1.

## 5. Top 3 priorities

1. **M1** — `ln 3` is the one literal-argument site in the fit's shipped code. It is computed
   again at every position and differs by one unit from what the libm build used. Hoist it, and run
   the Q1 experiment if parity fails.
2. **Mi1** — the drawn evidence of `bench_fixtures` and the contamination tests changes realisation,
   so B7's fit benchmark does not time the same workload A3 timed.
3. **Nit** — `fit.rs:2259`'s edited comment is 105 characters, and the plan's call count (125)
   includes it.

## 6. Findings

### Blocker

None.

### Major

#### M1: src/parameter_estimation/joint/fit.rs:2183 — `ln 3` is computed at every position, and at one unit from the value the libm build used

- **Confidence:** High for the bits (computed, §3); Medium that it changes the parameters file (Q1).
- **Problem:** `let ln_three = float::ln(CANDIDATE_ALTERNATIVES as f64);` sits inside `one_position`,
  which runs once per position per expectation pass. Before B5 it was `(CANDIDATE_ALTERNATIVES as
  f64).ln()`, whose argument is a constant, so the compiler folded it to the host library's
  correctly rounded `0x3ff193ea7aad030b`. That happened in A3's libm build too, because symbol
  substitution only reaches calls that survive to run time. Now it is libm's `0x3ff193ea7aad030a`.
  `ln_three` is subtracted from the fixed-alternative, segregating and duplicated branches at every
  position (`:2193`, `:2203`, `:2227`). A one-unit shift there can move the last digits of `log_likelihood`
  and of every share the EM loop returns. The TOML writer prints the shortest representation that
  round-trips, so any bit change shows.
- **Other literal sites checked and cleared:**
  - `ssr_fit.rs:2359`: `ln(TAU)` is the same bits in libm, Apple and correctly rounded.
  - `depth_bins.rs:234-307`: `widening_ratio(8, 124, 11)` and `powi(ratio, k)` have constant operands
    after inlining, so the old code may have folded them. Every rung is rounded to an integer, though,
    and the nearest rung to a `.5` boundary is `124·r⁷ = 709.42`, 4.9 × 10¹¹ units away.
  - Every other `float::` call with a literal argument in `parameter_estimation` is in `#[cfg(test)]`
    or `bench_fixtures`.
  - `powi` in shipped code (`contamination::binomial`, `Slippage` in `ssr_fit.rs:126,138`) has run-time
    operands, where B1 found std's `powi` and the loop bit-identical.
- **Why it matters:** B5 is the step whose parity check is meant to be byte-exact. Without this
  explanation, a mismatch looks like a conversion error in 124 calls. It also puts a libm call back
  on the hot path that the old build had as a constant. That costs one call per position; the fit's
  other logarithms per position are far more numerous (see §7).
- **Suggested fix:** hoist the constant out of the per-position function. It depends on nothing.
  ```rust
  // fit.rs, beside `const CANDIDATE_ALTERNATIVES: usize = 3;`
  /// `ln 3`, the log of the number of candidate alternatives. **Written as bits**, because the
  /// correctly rounded value is `…030b` and libm's is `…030a`; the old build folded it to `…030b`.
  const LN_CANDIDATE_ALTERNATIVES: f64 = f64::from_bits(0x3ff1_93ea_7aad_030b);
  ```
  That keeps parity **only** if the owner chooses the old bits (Q2). Otherwise compute it once in
  `expectation_pass` and pass it in, or read it from a `LazyLock<f64>` like `GAP_EXTEND_PROB` in
  `alignment/ssr_best_path_flat_gap.rs:220`, and record the one-unit difference in the commit.
  Whichever is chosen, add the test in §8 so the choice is explicit.

### Minor

#### Mi1: src/parameter_estimation/joint/fit.rs:3393-3421, ssr_fit.rs:2714-2738, contamination.rs:1484-1520 — the seeded draws change realisation, including the benchmark's evidence

- **Confidence:** High.
- **Problem:** the gamma, normal and Poisson generators in `fit::bench_fixtures`,
  `ssr_fit::bench_fixtures` and `contamination`'s test RNG now use libm's `ln`, `cos`, `powf` and
  `exp`. Rejection sampling compares `ln u` against a bound (`fit.rs:3407`, `ssr_fit.rs:2728`,
  `contamination.rs:1499`). When one unit flips that comparison, the generator consumes a different
  number of uniforms, and every later draw in the stream changes. The conversion is still correct:
  the draws are now the same on every platform, and before B5 they already differed between macOS
  and Linux while the thresholds held on both. But `benches/ng_joint_fit_perf.rs` draws its evidence
  from these fixtures (`Cargo.toml`, feature `bench-fixtures`). B7 plans to compare the fit
  benchmark with A3's, and it will be timing a different cohort draw.
- **Why it matters:** a B7 benchmark delta would mix libm's cost with a changed workload. The EM
  iteration count depends on the draw.
- **Suggested fix:** in B7, report the bench's iteration counts beside the times, or compare against
  a benchmark run of A3's libm build, which drew through libm for all but `cos`. A3 did not define
  `cos` (A3 report §"The libm build"), so even that build's Box–Muller draws differ from B5's. State
  this in the B7 report. The macOS suite is the remaining check that no seeded test threshold sits at
  the edge. The draws are platform-independent now, so a Linux pass should carry over.

### Nits

- `fit.rs:2259`: the edited doc comment is 105 characters, past the 100 the surrounding prose keeps
  (rustfmt does not wrap comments). Rewrap at "floors".
- `calibration.rs:404`: `mean_error_probability() − float::exp(−1.5)` against `1e-12` now compares
  the same function on the same argument. It still checks that `mean_error_probability` exponentiates
  `mean_log_error`; accepted under plan §5 ("test code converts too"). No change needed.
- The plan's B5 line says 125 calls; there are 124 code calls plus one comment. Correct the plan
  when ticking B5.

## 7. Out of scope observations

- **`one_position` takes up to 9 logarithms per position of values fixed for the whole pass**
  (`fit.rs:2207-2236`, pre-existing): `ln(1 − duplicated_share)`, `ln p_invariant`, `ln p_fixed_alt`,
  `ln p_segregating`, `ln duplicated_share` and `ln share`, per class (two classes). It also retakes
  `ln multiplicity[candidate]` in six loops (`:2190`, `:2196`, `:2218`, `:2328`, `:2365`, `:2443`). With libm's `ln` costing 0.5 to 1.5 ns more a call (A2), these are the cheapest part of
  the fit's 30% slowdown to win back. Suggested follow-up: compute the parameter logs once per
  expectation pass into a small struct, and `ln multiplicity` once per position. The same values in
  the same order leave the output bit-identical.
- **Literal-argument sites elsewhere in shipped code**, for completeness of the parity story and
  already reviewed in B2: `ssr_best_path_flat_gap.rs:220` (`exp(−1)`), `:241` and `:952`
  (`ln GAP_OPEN_PROB`). These feed alignment. `estimate-parameters` reads fixed psps, so they could
  reach the fit only if it realigns. The VCF parity (`1e40bd1e`) already reproduced after B2–B4.

## 8. Missing tests to add now

### `one_position` / the `ln 3` constant (M1)

- **`ln_candidate_alternatives_is_within_one_unit_of_libm`**: covers the constant written as bits, if
  chosen. It would catch someone changing `CANDIDATE_ALTERNATIVES` without updating the bits. It
  follows the B2 rule that each bits constant is tied to `float::ln` within one step.
  ```rust
  #[test]
  fn ln_candidate_alternatives_is_within_one_unit_of_libm() {
      let computed = float::ln(CANDIDATE_ALTERNATIVES as f64);
      assert!(LN_CANDIDATE_ALTERNATIVES.to_bits().abs_diff(computed.to_bits()) <= 1);
  }
  ```

## 9. What's good

- Multi-line receivers were rebuilt with the whole expression inside the call, not left as a
  trailing chain: `fit.rs:1845-1850` (`float::exp(f64::from(…unwrap_or(0.0)))`) and
  `slippage_curve.rs:631-635` (`float::exp(…).clamp(…)`). Both keep the original order of operations.
- Closures over `&f64` dereference exactly where a plain function needs `f64` (`fit.rs:2825`
  `float::ln(*f)`), and rely on `f64 − &f64` where the old code did (`fit.rs:2826`, `:3031`), with no
  added copies.
- The ln-gamma reflection keeps its nesting, `ln(π / sin(πx)) − lnΓ(1 − x)` (`ssr_fit.rs:2350-2351`),
  and the Beta front factor keeps its single `exp` over the summed logs (`ssr_fit.rs:2390`). Neither
  split a hoisted `ln_beta`.
- The test-only import in `census_moments.rs:902` is placed in `mod tests`, where it is the only user.
  So `--lib` builds carry no unused import, and the files whose shipped code uses `float` let their
  tests and `bench_fixtures` pick it up through `use super::*`.

## 10. Commands to re-verify

- `rg -n --pcre2 '\.(ln|exp|powf|powi|log10|log2|ln_1p|exp_m1|sin|cos|tan)\s*\(' src/ | rg -v '^src/float.rs'`
  should return only `self.gamma(` hits.
- `uv run --no-project --with mpmath python tmp/review_2026-09-14_portable_float_B5/literal_sites.py`
- `sh tmp/parity_oracle.sh fit tmp/parity/B5_macos` (and `_linux` in the container); if it mismatches,
  repeat with `ln_three` set to `f64::from_bits(0x3ff193ea7aad030b)` (Q1).
- `cargo test --lib --tests --all-features` on macOS.

### Author response convention
Address each finding by its identifier (e.g., "M1", "Mi1") with one of: `fixed in <commit>` /
`disputed because …` / `deferred to <issue>` / `won't fix because …`. Answer open questions from
section 4 first.

## Fixes applied and results (2026-09-14)

| finding | what was done |
|---|---|
| M1 `ln(3)` in `one_position` | **Confirmed as the whole difference from the libm build.** The converted fit wrote `0fb3e50d` on macOS and on Linux — the same file on both — against the libm build's `3f5e5e7f`, differing in 8 of 574 lines (two read groups' error multipliers, four inbreeding coefficients, two concentrations). With `ln_three` temporarily set to `f64::from_bits(0x3ff1_93ea_7aad_030b)`, the old compile-time value, the macOS fit wrote `3f5e5e7f` exactly; the probe was reverted. libm's value is kept, by the plan's rule: a constant is written out as its old bits only where converting it moves a recorded answer, and no recorded answer moved. Hoisting it out of `one_position` is left to the performance follow-up the review suggests |
| Mi1 seeded draws | macOS suite 4,795 passed, 0 failed: no test threshold sat at the edge. B7 will note that the fit bench's drawn cohort differs from A3's |
| open question: `estimate-contamination` parity | not checked; carried to Checkpoint B as a known gap |
| nit comment width | rewrapped |

Results on the release builds of this step, before the comment rewrap: identity oracle's five checksums unchanged (Linux); parity calls `1e40bd1e` on Linux and macOS; container suite 4,796 passed, 0 failed; macOS suite 4,795 passed, 0 failed; clippy clean with `--all-targets --all-features` and with `--lib` alone.
