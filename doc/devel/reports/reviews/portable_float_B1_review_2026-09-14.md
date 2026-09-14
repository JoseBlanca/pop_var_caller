# Code Review: portable_float_B1
**Date:** 2026-09-14
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** Step B1 of `doc/devel/implementation_plans/portable_float.md`: the new `src/float.rs` and its registration in `src/lib.rs`, as the uncommitted working tree of branch `portable-float`
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** a diff. `git diff HEAD` (the one line added to `src/lib.rs`, and the plan's Milestone B text) plus the untracked `src/float.rs`.
- **Reviewed against:** `HEAD` = `75722336` plus the working tree, worktree `pop_var_caller-portable-float`.
- **In-scope files:**
  - `src/float.rs` (new, 250 lines)
  - `src/lib.rs` (`pub mod float;`)
  - `Cargo.toml` (the `libm` dependency the module now rests on; unchanged, but a direct dependency of the change)
  - `doc/devel/implementation_plans/portable_float.md` (the B1 item, checked only for agreement with the code)
- **Checked against, not reviewed:** `libm-0.2.16` source (`src/math/{log,exp,pow,log10,log1p,expm1,sin,cos,k_sin,k_cos,rem_pio2,rem_pio2_large,scalbn,floor}.rs`, `src/math/arch/mod.rs`), reports A1, A2 and A3 under `doc/devel/reports/implementations/`.
- **Deliberately out of scope:** the call sites (B2–B5) and the clippy ban (B6), except where this module's shape decides how they will go.
- **Categories:** applied by a single reviewer rather than fanned out (see §3): reliability, errors, naming, idiomatic, refactor_safety, smells, tooling, extras (stable output, hot path). `unsafe_concurrency` and `module_structure` do not apply (no `unsafe`, no shared state, one file).

## 2. Verdict

**Approve-with-changes.** Every function delegates to the right libm function with the right argument order; `powi` is correct at every edge examined and is the same algorithm, in the same multiplication order, as the run-time `powi` it replaces. No Blocker or Major. The Minor findings are about what the tests can and cannot detect, one doc claim that is broader than its evidence, and a `Cargo.toml` range that lets a routine `cargo update` move every output.

## 3. Execution status

- **Commands run by the reviewer:** none that build. Read-only `git`, `rg` and `grep` over the worktree and the libm source.
- **Commands not run, and why:** `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`, `cargo test`, `cargo doc` — the dispatch forbade running cargo and reported them already green (fmt, clippy `-D warnings`, doc, 4,792 tests in the Linux container; the float tests natively on macOS). `cargo audit` — not run, no new dependency.
- **Deviation from the skill:** no per-category sub-agents. The skill's fan-out exists so agents can mutation-test in their own worktrees; with cargo ruled out there is nothing for them to run, and the scope is one 250-line file. No mutation results are reported; the test-strength findings below (Mi1, Mi3) are argued from the fixtures and carry verification steps.
- **Findings labelled "Needs verification":** 2 (Mi1, Mi3).

## 4. Open questions and assumptions

1. **Is `libm` meant to float within `0.2.x`?** `Cargo.toml` says `libm = "0.2"`. Portability survives an upgrade — both platforms move together — but every output moves, and the only thing that notices is the pinned-bits test. Affects Mi2.
2. **Assumption: std's run-time `powi` is the square-and-multiply loop on every target the caller ships to.** The compiler-builtins source is not on this machine (no `rust-src` component). From its published source, Rust's `compiler_builtins` `__powidf2`, LLVM compiler-rt's `__powidf2` and libgcc's `__powidf2` all compute `1 · x^(bit 0) · (x²)^(bit 1) · …` in increasing bit order and apply `1/r` at the end for a negative exponent — the order `src/float.rs:53-68` uses. LLVM's expansion of `powi` with a constant exponent (`ExpandPowI`) also multiplies in that order. The unit test `powi_matches_std_at_run_time_…` confirms it on aarch64 macOS and aarch64 Linux; x86_64 (CI, `ubuntu-latest`) is confirmed only when CI runs it. Affects Mi3.

## 5. Top 3 priorities

1. **Mi1** — the std-agreement test's tolerance is zero for subnormal results, so `exp(-745.0)` fails on any platform whose `exp` is one step away; replace the relative tolerance with a count of representable steps.
2. **Mi2** — pin `libm` exactly, or at least say in `Cargo.toml` that upgrading it moves output; the dependency's comment still describes it as the `lgamma` crate.
3. **Mi3** — `powi`'s agreement with std is tested only for bases in (0, 2.46] and exponents −40 to 300; add the edges (`i32::MIN`, `i32::MAX`, negative and signed-zero bases, infinities, NaN to the power 0, overflow before a reciprocal) against std at run time.

## 6. Findings

### Blocker

None.

### Major

None.

### Minor

#### Mi1: [src/float.rs](../../../../src/float.rs#L182-L185) — the "within a few units" tolerance is zero for subnormal and tiny results
- **Confidence:** High (the arithmetic); whether any listed argument actually trips it off aarch64 is **Needs verification**.
- **Problem:** `close` accepts `(ours - platform).abs() <= 4.0 * f64::EPSILON * platform.abs()`. For a normal result that is 4 to 8 units in the last place. For a subnormal result the spacing between representable numbers is fixed at 4.9e-324, while `4 · ε · |p|` underflows to 0, so the test demands identical bits. `exp(-745.0)` is in the list (line 192) and its value is the smallest subnormal, 4.9e-324; a platform library that rounds it to 0 or to 9.9e-324 fails the test although it is one unit away. Both aarch64 platforms passed; x86_64 glibc (the CI runner) was never measured in A2.
- **Why it matters:** the test's name promises a few units in the last place; near underflow it is an exact-equality test, and a CI failure there would look like a wrong delegation.
- **Suggested fix:** measure the distance in representable steps, which is what the name says:
  ```rust
  fn close(ours: f64, platform: f64) -> bool {
      // Both are one sign, so the difference of their bit patterns is the number of
      // representable numbers between them; a sign difference only matters at zero.
      if ours.to_bits() == platform.to_bits() || ours == platform {
          return true; // the second equates +0.0 and −0.0
      }
      ours.is_sign_negative() == platform.is_sign_negative()
          && (ours.to_bits() as i64 - platform.to_bits() as i64).abs() <= 4
  }
  ```

#### Mi2: [Cargo.toml](../../../../Cargo.toml#L133-L137) — `libm = "0.2"` lets `cargo update` move every output, and its comment describes a different use
- **Confidence:** High.
- **Assumptions:** Open question 1.
- **Problem:** after B2–B5 every transcendental result in the caller is libm's. The caret range admits any 0.2.x; libm's recent releases have been replacing musl ports with new generic implementations, and `Cargo.lock` (which does pin 0.2.16) is regenerated by a routine `cargo update`. The comment above the line still says libm is there for `lgamma`.
- **Why it matters:** a dependency bump would silently change genotypes and fitted parameters on every platform at once; the reproducibility this plan buys is then only as stable as a lockfile nobody reads.
- **Suggested fix:**
  ```toml
  # libm: every transcendental function the caller evaluates (`src/float.rs`), and
  # `lgamma` for the genotype prior. One Rust implementation compiled into the
  # binary, so macOS and Linux give the same bits. **Pinned exactly: a new libm
  # release can round differently and move output everywhere.** Bump deliberately
  # and update the pinned bits in `float::tests`.
  libm = "=0.2.16"
  ```

#### Mi3: [src/float.rs](../../../../src/float.rs#L217-L243) — `powi`'s agreement with std is tested on a narrow range; the edges are checked against constants, not std
- **Confidence:** Medium. The code is correct at every edge below by hand evaluation; the gap is that no test would catch a change that breaks bit-agreement there. **Needs verification** on x86_64 (open question 2).
- **Problem:** `powi_matches_std_at_run_time_…` uses 200 bases from 1e-6 to 2.457 and exponents −40 to 300: 68,200 pairs, all with a positive base and no overflow before the reciprocal. `powi_handles_its_edges` checks eight values against literals, without `black_box` and without std. Nothing compares with std at `i32::MIN` (where `unsigned_abs` gives 2³¹ and the reciprocal of an infinite product gives 0), `i32::MAX`, negative bases with odd negative exponents, `-0.0` (sign of the infinity), `±inf`, `NaN` to the power 0 (1.0), or bases above 1 whose negative power overflows (`1.0 / inf`).
  Hand evaluation of the loop, for the record: `powi(2.0, i32::MIN)` squares 31 times to `inf`, then `1/inf = 0.0` — matches std and the test. `powi(-0.0, -1)` = `1/-0.0 = -inf`. `powi(NaN, 0)` never multiplies and returns 1.0. `powi(x, 3)` = `(1·x)·(x·x)`, the order LLVM's constant-exponent expansion uses.
- **Why it matters:** B2–B5 route per-read and per-fit-iteration calls through this loop (A1 §2: `alignment/stutter.rs:783`, `ssr_fit.rs:125`, `contamination.rs:1278`), and the module's promise is bit-equality with what those sites computed before.
- **Suggested fix:** see §8, `powi_matches_std_at_run_time_on_edge_operands`.

#### Mi4: [src/float.rs](../../../../src/float.rs#L10-L13) — "libm's only processor-specific code is `sqrt`, `fma` and `rint`" holds on aarch64 only; "every function here … through libm" does not include `powi`
- **Confidence:** High.
- **Problem:** two claims are broader than the source.
  1. `libm-0.2.16/src/math/arch/mod.rs` selects `sqrt`, `fma` and `rint` on aarch64, `sqrt` and `fma` on x86_64 with SSE2, and on 32-bit x86 without SSE also `floor`, `ceil` and an x87 `exp` (`exp.rs:86-90`). The x87 `exp` is not exactly rounded. A2 states the claim with its targets; the module doc drops them. (The second half is correct: `rg fma` over `log`, `exp`, `pow`, `log10`, `log1p`, `expm1`, `sin`, `cos` and their helpers `k_sin`, `k_cos`, `rem_pio2`, `rem_pio2_large`, `scalbn` finds nothing. `pow` does call `sqrt`, and `rem_pio2_large` calls `floor`, both exact.)
  2. Line 10 says every function computes through libm; `powi` is a Rust loop that never calls libm.
- **Why it matters:** the module doc is where B6's ban message will send readers; a portability guarantee stated without its targets reads as covering 32-bit x86.
- **Suggested fix:**
  ```rust
  //! Every function here computes in Rust instead — through the [`libm`] crate, or for `powi` as
  //! plain multiplication — so the same argument gives the same bits wherever the binary runs. On
  //! aarch64 and x86_64, libm's only processor-specific code is `sqrt`, `fma` and (aarch64) `rint`,
  //! which IEEE 754 requires to be rounded exactly; `pow` calls `sqrt`, and none of the functions
  //! this module uses calls `fma`. (32-bit x86 without SSE also swaps in an x87 `exp`; the caller
  //! does not target it.)
  ```

### Nits

- [src/float.rs](../../../../src/float.rs#L53) — `powi` is the only function without `#[inline]`. Release builds are fat LTO with one codegen unit, so the shipped binary is unaffected, but the `profiling` and `soak` profiles (`lto = false`, 16 codegen units) may leave a per-read call un-inlined and distort profiles. Add `#[inline]`.
- [src/float.rs](../../../../src/float.rs#L181-L215) — the std-agreement test's arguments are array literals without `black_box`, so the compiler may evaluate both sides while building. On a native build that uses the same library it runs on, so the result stands; wrap the arguments in `black_box` anyway, as the other two tests do, so the test says what it measures.
- [src/float.rs](../../../../src/float.rs#L105-L106) — "the same on every platform" is recorded on two aarch64 platforms (line 156 says so). Say "same on every platform checked (aarch64 macOS and Linux; x86_64 by CI)".
- [src/float.rs](../../../../src/float.rs#L21-L22) — "`exp` 0.5 to 3.3 ns more": A2 measured the two level past underflow on Linux. "Up to 3.3 ns more" is accurate.
- [src/float.rs](../../../../src/float.rs#L46-L47) — "agreed on every argument measured, on both platforms": A2's `powi` row covered exponents 0 to 300 only. The negative branch rests on this module's unit test; say so.
- Forward note for B6: the std-agreement and `powi` tests are the one legitimate use of the banned std methods, so `float::tests` will need a scoped `#[allow(clippy::disallowed_methods)]` with a line saying why.
- Pinned-bits test ([src/float.rs](../../../../src/float.rs#L108-L175)): it catches a libm upgrade or a rewired delegation, but a function reverted to std would pass on a platform where std happens to round those arguments like libm. A2 says `sin(0.77)` differs between the two platforms' std, so `sin` is covered; for the others, choosing one pinned argument per function from A2's "libm vs std differs on both platforms" output would make the test also prove "this is libm". Low priority once B6's ban forbids std calls in `src/float.rs`.

## 7. Out of scope observations

- `doc/devel/implementation_plans/portable_float.md` B6 lists `log2`, `log`, `tan` and the `f32` twins in the ban. `rg` finds no such call in `src/`, `tests/` or `benches/` (the one `f32` near a transcendental, `parameter_estimation/joint/fit.rs:2278`, is an `f64` `exp` cast afterwards), so B1 correctly omits them; the ban can still list them so a future call is steered to `crate::float` — at which point the missing function must be added there.
- `PROJECT_STATUS.md` has no block for `portable-float`, and the A1–A3 reviews did not add one; not changed here, since the dispatch scoped this review to a report.

## 8. Missing tests to add now

### `powi`

- **`powi_matches_std_at_run_time_on_edge_operands`** — covers the extremes of both operands; catches a change to the reciprocal, the sign handling or the `unsigned_abs` that breaks agreement at `i32::MIN` or on negative and non-finite bases.
  ```rust
  #[test]
  fn powi_matches_std_at_run_time_on_edge_operands() {
      let bases = [
          0.0, -0.0, 1.0, -1.0, 0.5, -0.5, 2.0, -2.0, 1e-160, 1e160,
          f64::MIN_POSITIVE, 4.9e-324, f64::INFINITY, f64::NEG_INFINITY, f64::NAN,
      ];
      let exponents = [0, 1, -1, 2, -2, 3, -3, 1023, -1024, 1075, i32::MAX, i32::MIN, i32::MIN + 1];
      for base in bases {
          for exponent in exponents {
              let (b, e) = (black_box(base), black_box(exponent));
              let ours = powi(b, e);
              let std = b.powi(e);
              assert!(
                  ours.to_bits() == std.to_bits() || (ours.is_nan() && std.is_nan()),
                  "powi({base:e}, {exponent}): ours {ours:e}, std {std:e}"
              );
          }
      }
  }
  ```

### `exp`, `ln` and the others

- **`each_function_agrees_with_std_within_four_representable_steps_near_underflow`** — the Mi1 rewrite of `close`, plus arguments whose results are subnormal (`exp(-740.0)`, `exp_m1(-1e-310)`, `sin(1e-310)`); catches a tolerance that silently becomes exact equality.

## 9. What's good

- [src/float.rs](../../../../src/float.rs#L26-L98) — free functions rather than an extension trait are the only shape that works: a trait method named `ln` is shadowed by the inherent `f64::ln`, and `clippy::disallowed-methods` needs a path to point at. Path-form sites such as `.map(f64::exp)` (`parameter_estimation/calibration.rs:238`) convert one-for-one to `.map(float::exp)`.
- [src/float.rs](../../../../src/float.rs#L44-L52) — the `powi` doc explains why the loop is written out with the mechanism (build-time evaluation with the build machine's `pow`) and its measured size (15 of 20, up to 26 units), both matching A2.
- [src/float.rs](../../../../src/float.rs#L177-L179) — pairing a pinned-bits test with an against-std test separates "the value changed" from "the function is the wrong one"; either alone would misreport the other failure.
- [src/float.rs](../../../../src/float.rs#L220-L231) — `black_box` on both `powi` operands keeps std's call at run time, which is exactly the property the test claims.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo test --lib float::` (Linux container) and `cargo test --lib float::` natively on macOS — after Mi1 and Mi3.
- CI on `ubuntu-latest` (x86_64) — the first run of `float::tests` off aarch64; watch `outputs_are_pinned_to_the_bit` and the std-agreement test.
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings` after adding the new test.

## Fixes applied (2026-09-14)

| finding | what was done |
|---|---|
| Mi1 no slack at subnormal results | `close` counts representable steps (at most 4, same sign), as proposed |
| Mi2 `libm = "0.2"` floats | pinned to `=0.2.16`; the comment says why an upgrade must be deliberate. Decided by the implementer: an unpinned version would let `cargo update` move every output |
| Mi3 `powi` edges not checked against std | `powi_matches_std_at_run_time_on_its_edges`: 14 bases (signed zeros, ±1, ±2, −0.5, extremes, infinities, NaN) × 11 exponents including `i32::MIN` and `i32::MAX`, bit-equal to std at run time on macOS and Linux |
| Mi4 module doc overstated | targets named (aarch64, x86_64; 32-bit x86 without SSE excluded); `powi` named as plain multiplication; pinned bits stated as recorded on aarch64 only |
| nits | `#[inline]` on `powi`; `black_box` on the std-comparison arguments; "up to 3.3 ns"; pinned-bits doc says "both platforms recorded". The B6 `allow` on `float::tests` is noted for that step |

Validation after the fixes, in the container: fmt clean, clippy `-D warnings` clean, doc clean, `cargo test --lib --tests --all-features` 4,793 passed, 0 failed, 4 ignored; `float::tests` 6 passed on macOS natively too.
