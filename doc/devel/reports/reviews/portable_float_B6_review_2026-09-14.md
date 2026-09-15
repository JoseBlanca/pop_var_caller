# Code Review: portable_float_B6
**Date:** 2026-09-14
**Reviewer:** rust-code-review skill (orchestrator, single-pass — see §1)
**Scope:** uncommitted diff of plan step B6 (`doc/devel/implementation_plans/portable_float.md`, Milestone B): a clippy `disallowed-methods` ban on std's transcendental float methods, with the allows it needs
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** `git diff HEAD` in the worktree `pop_var_caller-portable-float` on top of
  `08ec742a` (B5), plus the untracked `clippy.toml`. 24 tracked files, +117 / −3.
- **Reviewed against:** `HEAD` = `08ec742a` plus uncommitted changes.
- **In-scope files:** `clippy.toml`; `src/float.rs` (module doc paragraph, test-module allow);
  `src/parameter_estimation/joint/fit.rs` (comment rewrap only); the 22 example files that gained a
  file-level allow (20 in `examples/`, plus `examples/shared/stutter_model.rs` and
  `examples/shared/stutter_table.rs`); `.github/workflows/ci.yml` and `rust-toolchain.toml` as the
  context that decides whether the ban is enforced.
- **Deliberately out of scope:** the conversions themselves (B2–B5, reviewed separately); maths
  inside dependencies, which no lint on this crate can reach (stated once under §7).
- **Categories:** the fan-out of mutation-testing agents was **not** run. The dispatcher forbade
  `cargo` and edits, which removes what worktree isolation is for, and the diff is configuration
  plus attributes. One pass covered `tooling` (the clippy configuration and CI), `refactor_safety`
  (can new code get past the ban without an allow), `extras` (does the diff match the plan step) and
  `naming`/doc accuracy. `reliability`, `errors`, `unsafe_concurrency`, `module_structure`,
  `defaults`: no code paths changed. `PROJECT_STATUS.md` was not updated, because the dispatcher
  forbade edits other than this report.

## 2. Verdict

**Approve-with-changes.** The ban does what the plan asks for `f64`: every stable std `f64` method
whose result depends on the platform's maths library is listed, and nothing IEEE 754 rounds exactly
is. Its realistic bypasses are closed. A function-pointer use like `.map(f64::ln)` and the `<f64>::ln`
spelling resolve to the same method, and clippy checks path expressions, not only calls. A direct
`extern "C"` call to the C `log` needs `unsafe`, which `Cargo.toml` forbids. No crate offering a
`Float` trait is a direct dependency, and edition 2024 does not expose transitive ones. The attribute
placement is correct in all 22 examples, and CI runs the same clippy command against the same
`clippy.toml`.

Three things to change. The `f32` list covers 16 of the 26 methods, while `src/float.rs`'s new
paragraph says the ban covers the "`f32` twins" (Mi1). The paragraph also says "everywhere in the
crate" without naming the two exemptions (Mi2). And `clippy.toml` is still untracked: if it is left
out of the commit, CI stays green with no ban at all, because an `allow` for a lint that never fires
is silent (Mi3).

## 3. Execution status

- **Commands run by the reviewer** (host, read-only):
  - `git -C … diff HEAD`, `git status --short`: 24 modified tracked files, `?? clippy.toml`.
  - `git check-ignore -v clippy.toml`: exit 1, so the file is not ignored and `git add` will take it.
  - Sweep for banned calls across `src tests benches examples`. The method-call pattern was
    `\.(ln|log|log2|log10|ln_1p|exp|exp2|exp_m1|powf|powi|sin|cos|tan|sin_cos|asin|acos|atan|atan2|sinh|cosh|tanh|asinh|acosh|atanh|cbrt|hypot)\(`,
    run alongside the `(f64|f32)::(…)` path form. **Every one of the 22 allowed examples has at least
    one hit** (from 1 in `examples/shared/stutter_table.rs` to 58 in `examples/ng_multilib_key_harness.rs`).
    **No example without the allow has a hit**, so the list of 22 is neither short nor padded, as far
    as a text search can tell. In `src/`, `tests/` and `benches/`, the path form appears only in doc
    comments (`src/alignment/emission.rs:162-163`, `src/calling/likelihood/generic.rs:2103`), and the
    rarer methods (`tan`, `asin` … `hypot`, `log2`, `exp2`) do not appear at all.
  - The 10 calls clippy flagged before the allows (dispatcher's figure) map onto `src/float.rs`'s
    tests as `ln`, `log10`, `exp`, `exp_m1`, `ln_1p`, `sin`, `cos`, `powf` and `powi` twice. **So 9 of
    the 26 `f64` entries, and none of the 16 `f32` entries, are proven by a hit to resolve.** I
    checked the other 33 names by eye against std's method names and found no misspelling (see Mi4
    for why that is worth making mechanical).
  - `Cargo.toml`: edition 2024; `unsafe_code = "forbid"`; direct dependencies contain no maths crate
    other than `libm`. `Cargo.lock` has `num-traits 0.2.19` and `num-complex` only as transitive
    dependencies (through the arrow/parquet stack).
  - `grep cfg(target_os…)`: no macOS-only code in `src/`, `tests/` or `benches/`, so a Linux clippy
    run (CI, container) lints everything that ships. The `#[cfg(unix)]` blocks are compiled on Linux.
- **Commands not run:** every `cargo` command, by instruction. The dispatcher's gates (fmt, clippy
  `-D warnings`, doc, container tests 4,848 passed / 3 pre-existing failures, macOS 4,795 passed,
  release `calling::` 1,165 passed) are taken as given. A clean `clippy --all-targets` also proves
  the inner attributes are well placed: `#![allow]` after an item is a hard compile error.
- **Needs verification:** 2 (Mi4, and the function-pointer claim in §2, which rests on clippy's own
  UI tests rather than a run here — the probe in §10 settles both).

## 4. Open questions and assumptions

1. **Does clippy 1.98 report a `disallowed-methods` path that resolves to nothing, and does
   `-D warnings` make that fatal?** Recent clippy prints a warning for such a path, and entries
   accept `allow-invalid = true` to silence it. I believe this is a plain diagnostic rather than a
   lint, which would mean `-D warnings` does **not** make it an error. Affects Mi4.
2. **Are `f64::gamma`, `ln_gamma`, `erf`, `erfc` still unstable on 1.98?** They existed as unstable
   std methods when I last saw std (`float_gamma`, `float_erf`). If any is now stable it is an
   unlisted platform-library call. Affects the nit on future-proof entries.

## 5. Top 3 priorities

1. **Mi3** — commit `clippy.toml`. Without it the step's only deliverable disappears and nothing fails.
2. **Mi1** — list the ten missing `f32` methods, so the doc's "`f32` twins" is true.
3. **Mi4** — run the one-minute probe in §10 once, so a future typo in `clippy.toml` cannot switch an
   entry off silently.

## 6. Findings

### Minor

- **Mi1:** [clippy.toml](../../../../clippy.toml) — **[Minor]** The `f32` list stops at 16 of the 26 methods the `f64` list bans
- **Confidence:** High
- **Problem:** The `f64` block lists 26 methods. The `f32` block omits `sin_cos`, `asin`, `acos`,
  `atan`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`. `src/` mentions `f32` 207 times and
  calls none of these today, so nothing escapes now. But `src/float.rs:29-31` says the ban covers
  "the other transcendental methods (and their `f32` twins)", which is not what the file does.
- **Why it matters:** the first `f32` `tanh` or `atan` added to shipped code computes through the
  platform library with no warning, which is the exact failure B6 exists to prevent.
- **Suggested fix:** append, in the file's style:
  ```toml
  { path = "f32::sin_cos", reason = "use crate::float on an f64: std's rounds differently on each platform" },
  { path = "f32::asin", reason = "use crate::float on an f64: std's rounds differently on each platform" },
  { path = "f32::acos", reason = "use crate::float on an f64: std's rounds differently on each platform" },
  { path = "f32::atan", reason = "use crate::float on an f64: std's rounds differently on each platform" },
  { path = "f32::sinh", reason = "use crate::float on an f64: std's rounds differently on each platform" },
  { path = "f32::cosh", reason = "use crate::float on an f64: std's rounds differently on each platform" },
  { path = "f32::tanh", reason = "use crate::float on an f64: std's rounds differently on each platform" },
  { path = "f32::asinh", reason = "use crate::float on an f64: std's rounds differently on each platform" },
  { path = "f32::acosh", reason = "use crate::float on an f64: std's rounds differently on each platform" },
  { path = "f32::atanh", reason = "use crate::float on an f64: std's rounds differently on each platform" },
  ```

- **Mi2:** [src/float.rs](../../../../src/float.rs#L29-L31) — **[Minor]** "everywhere in the crate" omits the two exemptions
- **Confidence:** High
- **Problem:** The paragraph reads as if no allow exists. There are two kinds: this module's own test
  module, and 22 examples. `clippy.toml`'s header comment names both; the module doc, which is where
  a reader of the code lands, names neither.
- **Why it matters:** a reader who finds `x.ln()` in an example, or `#[allow(clippy::disallowed_methods)]`
  twenty lines below this paragraph, is left to wonder whether the ban is broken.
- **Suggested fix:**
  ```rust
  //! **Nothing else may call std's versions.** `clippy.toml` refuses `f64::ln`, `exp`, `powf`, `powi`
  //! and the other transcendental methods (and their `f32` twins) in the library, the binary, the
  //! tests and the benches, each with a message pointing here. Two places allow them: this module's
  //! tests, which compare against std on purpose, and the examples, which are research tools outside
  //! this guarantee.
  ```

- **Mi3:** [clippy.toml](../../../../clippy.toml) — **[Minor]** The ban is an untracked file, and forgetting it fails nothing
- **Confidence:** High
- **Problem:** `git status` shows `?? clippy.toml`. Every other change in the step is an `allow`, and
  `#[allow]` of a lint that never fires produces no warning. A commit made with `git commit -a`, or
  with the modified files staged by name, passes every gate and ships no ban.
- **Why it matters:** the step's whole effect is lost silently, and B7 would then measure a tree with
  no guard against regressions.
- **Suggested fix:** `git add clippy.toml` in the B6 commit, and check that `git show --stat HEAD`
  lists it. Optionally switch the allows to `#[expect(...)]` (see Nits): an expectation that is never
  met *does* warn, but only in the opposite direction, so it does not replace the `git add`.

- **Mi4:** [clippy.toml](../../../../clippy.toml) — **[Minor]** 33 of the 42 entries are unproven, and a misspelt entry may fail nothing — **Needs verification**
- **Confidence:** Low (depends on §4 Q1)
- **Assumptions:** clippy 1.98 prints a warning for an unresolvable path, not a lint that
  `-D warnings` promotes.
- **Problem:** A hit proves an entry resolves. Only `f64::{ln, log10, exp, exp_m1, ln_1p, sin, cos,
  powf, powi}` had one. The other 17 `f64` entries and all 16 `f32` entries have no callers, so a
  misspelling (`f64::atna`) would switch that entry off. The build would then print at most a warning
  that CI does not fail on. All 42 names are spelt correctly today.
- **Why it matters:** the entries most likely to be edited later — the "add it to crate::float" ones,
  when someone does add `tan` — are exactly the unproven ones.
- **Suggested fix:** run the probe in §10 once. If the unresolved-path warning is not fatal under
  `-D warnings`, add a CI step that fails on it:
  ```yaml
  - name: Clippy config resolves
    run: |
      cargo clippy --all-targets --all-features 2>&1 | tee clippy.log
      ! grep -q 'does not refer to' clippy.log
  ```
  (Take the exact message text from the probe's output; the wording above is from memory.)

### Nits

- `clippy.toml` reasons say "std's rounds differently on each platform". For most arguments the two
  libraries agree (A2 measured `exp` differing on 1 argument in 540 to 1 in 2,700), so "may round
  differently" is the accurate wording. The file's header comment already states it that way.
- Future-proofing (§4 Q2): `f64::gamma`, `f64::ln_gamma`, `f64::erf`, `f64::erfc` exist in std as
  unstable methods, so entries for them resolve today and would catch the day any is stabilised.
  `genetics::lgamma` is the portable route, and the reason can say so. Note that the `self.gamma(…)`
  calls in `contamination.rs:1487-1508` and `fit.rs:3388-3396` are the fixtures' own RNG method, not
  `f64::gamma`, and would not be flagged.
- `#[expect(clippy::disallowed_methods, reason = …)]` instead of `#[allow]`, in `src/float.rs`'s tests
  and the 20 top-level examples: an example that stops calling std's maths would then warn, and the
  list of exemptions stays exact. The repo already uses `#[expect]` (`src/run/callers.rs:1110`). Keep
  `allow` in `examples/shared/*.rs`, whose lint levels are inherited from each including example.
  Those two allows are redundant under today's only includer (`ng_str_stutter_rate.rs`, itself
  allowed), and harmless.
- `examples/ng_joint_contamination_harness.rs:40-49` and `examples/ng_multilib_key_harness.rs:89-97`:
  the new attribute now sits directly above the `//` comment that explains the *next*
  `#![allow(clippy::needless_range_loop)]`. A blank line between them keeps that comment visibly
  attached to its own attribute.
- `src/parameter_estimation/joint/fit.rs:2260-2262`: rewrap only, text unchanged. Nothing to review.

## 7. Out of scope observations

- **Dependencies are not linted.** Any maths std does inside noodles, parquet or arrow is outside
  every clippy configuration of this crate. None of them computes a value the caller writes to its
  output from a transcendental function, as far as the call sites in `src/` show, so nothing needs
  doing now. It is a limit of the mechanism worth one sentence in the plan's C-milestone
  cross-platform digest test, which is what would catch it.
- **Bypasses considered and closed.** Recorded so B7 and later reviews need not redo them:
  - *Function pointers and qualified spellings* (`.map(f64::ln)`, `<f64>::ln(x)`,
    `core::primitive::f64::ln`). Clippy's `disallowed_methods` checks path expressions, not only
    calls: its UI tests flag `let indirect: fn(&str, &str) -> Regex = Regex::new`,
    `Box::new(f32::clamp)` and `.map(Regex::new)`. All three spellings resolve to the same method.
    Needs verification by the probe; high expectation.
  - *A `Float` trait* (`num_traits::Float::ln`) resolves to the trait method, not `f64::ln`, so it
    would not be flagged. `num-traits` is only transitive, and edition 2024 does not put transitive
    crates in scope. A future direct dependency on `num-traits`, `statrs` or similar reopens this;
    the `libm` comment in `Cargo.toml` is the natural place to say so.
  - *`libm::` called directly* outside `src/float.rs` (`genetics::lgamma` does): portable by
    construction, so not a gap in what the ban protects. Banning it would only enforce routing.
  - *`std::intrinsics::logf64`* needs `#![feature(core_intrinsics)]`, unavailable on the pinned
    stable toolchain. *`extern "C" { fn log(…) }`* needs `unsafe`, which `Cargo.toml` forbids.
  - *Doctests* are not linted by clippy, but they do not reach the binary.
- **What correctly stays allowed.** `sqrt` is correctly rounded by IEEE 754. `mul_add` is IEEE's
  fused multiply-add, specified as correctly rounded: the hardware instruction and the C library's
  `fma` fallback on glibc and Apple give the same bits, and `src/` does not call it. `floor`, `ceil`,
  `round`, `trunc`, `%`, `rem_euclid`, `recip`, `to_degrees`/`to_radians` (a single multiplication),
  `min`/`max`/`clamp` and `as` casts are exact or singly rounded arithmetic. Leaving them off the list
  is right.
- **CI's toolchain comment is imprecise, pre-existing.** `.github/workflows/ci.yml:21-24` says
  `dtolnay/rust-toolchain@stable` "picks up" `rust-toolchain.toml`. The action installs stable; it is
  rustup, when `cargo` first runs, that switches to the pinned `1.98.0` (with `clippy`, listed in the
  toml). The effect is what the comment intends. `clippy.toml` at the manifest root is found without
  configuration, and the step's command matches the local gate exactly
  (`cargo clippy --all-targets --all-features -- -D warnings`, with `RUSTFLAGS: -D warnings`).

## 8. Missing tests to add now

None as Rust tests: the ban is enforced by the lint run, not by the suite. The equivalent check is
the probe in §10, run once, and optionally the CI step in Mi4.

## 9. What's good

- `clippy.toml`'s reasons split into "use `crate::float::x`" and "add it to `crate::float`". The
  message tells the next author which of the two jobs they have, not only that they were refused.
- The `powi` entry carries its own reason (compile-time folding with the build machine's `pow`),
  matching `src/float.rs:53-59`, rather than the generic rounding message that would be wrong for it.
- The exemptions are scoped as narrowly as the plan allows: a `#[cfg(test)]` module in `src/float.rs`,
  not the file; per example file, not a `[lints]` table entry that would also cover `src/`.
- `clippy.toml`'s header names what stays on std (`sqrt`, `abs`, basic arithmetic) and why, so the
  list's omissions read as decisions.

## 10. Commands to re-verify

- The dispatcher's gates, unchanged:
  `cargo clippy --all-targets --all-features -- -D warnings`, `cargo fmt --all -- --check`,
  `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --all-features`.
- **New: a one-off probe** (settles Mi4 and the function-pointer claim; delete afterwards).
  Add temporarily to the end of `src/float.rs`, and a deliberately misspelt entry
  `{ path = "f64::atna", reason = "probe" }` to `clippy.toml`:
  ```rust
  #[allow(dead_code)]
  fn disallowed_methods_probe(x: f64, y: f32) -> f64 {
      let by_pointer: Vec<f64> = [x].iter().copied().map(f64::ln).collect(); // expect: flagged
      let qualified = <f64>::tan(x);                                         // expect: flagged
      let twin = f64::from(y.tanh());                                        // flagged only after Mi1
      by_pointer[0] + qualified + twin
  }
  ```
  Then `cargo clippy --lib -- -D warnings`, and record: (a) whether `map(f64::ln)` and `<f64>::tan`
  are reported; (b) the exact text printed for `f64::atna`, and whether the command's exit status is
  non-zero because of it.
- `git show --stat HEAD` after committing: `clippy.toml` must be listed (Mi3).

## Fixes applied (2026-09-14)

| finding | what was done |
|---|---|
| Mi1 `f32` list incomplete | the ten missing `f32` entries added; 52 entries in all |
| Mi2 doc overstated reach | `src/float.rs` now names the two exemptions (its tests, the examples) and the check script |
| Mi3 `clippy.toml` untracked | committed with this step |
| Mi4 a typo would silently disable an entry | **Measured:** with an entry `f64::lnn_probe_typo` added, clippy printed nothing about it; with `.map(f64::exp)` added to library code, clippy refused it (function values are caught). New `scripts/check_float_ban.sh` writes a throwaway example calling every banned method, runs clippy, and fails unless clippy refuses exactly the listed set; it names a misspelled method (checked with `f64::cbrtt`: "no method named `cbrtt` found for type `f64`", exit 1) and passes on the real file ("all 52 entries in clippy.toml refuse their method") |
| nits | reasons say "may round differently"; blank lines added in the two examples. `gamma`/`ln_gamma`/`erf` entries not added: the canary could not call unstable methods on stable Rust, so the check could not prove them. `allow` kept over `expect` in examples |

Validation: fmt, clippy `--all-targets --all-features -D warnings`, doc clean in the container; `float::tests` pass; the check script passes. Suite results for this step are in the review's header (container with examples 4,848 passed, 3 failed that fail identically before the plan; macOS 4,795 passed, 0 failed); the fixes after them change only `clippy.toml`, doc comments, two blank lines and a new script.
