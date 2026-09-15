# Code Review: portable_float_B3
**Date:** 2026-09-14
**Reviewer:** rust-code-review skill (orchestrator, single-pass — see §1)
**Scope:** uncommitted diff of plan step B3 (`doc/devel/implementation_plans/portable_float.md`, Milestone B): std transcendental calls in `src/calling/`, `src/genetics.rs` and `src/types.rs` replaced by `crate::float`
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** the working-tree diff `git diff HEAD -- src/` in the worktree
  `pop_var_caller-portable-float`, on top of `79a22fff` (B2). 20 files, +214 / −160.
- **Reviewed against:** `HEAD` = `79a22fff` plus uncommitted changes.
- **In-scope files:** the 20 files the diff touches — `src/calling/genotype_prior/{dirichlet_multinomial,
  hardy_weinberg, mod, seed_generic}.rs`, `src/calling/genotype_table.rs`,
  `src/calling/inference/summarise_condition.rs`, `src/calling/likelihood/{generic, mod, ssr,
  ssr_emission}.rs`, `src/calling/loop_parity.rs`, `src/calling/parameters_file/{defaults,
  from_run_parameters, to_run_parameters, validate}.rs`, `src/calling/quality/{artifact_correction,
  mod}.rs`, `src/calling/run_parameters.rs`, `src/genetics.rs`, `src/types.rs`; the callee
  `src/float.rs`; the recorded-answer tests `src/calling/genotype_table_parity.rs`,
  `src/calling/loop_parity.rs`, `src/calling/quality_parity.rs`, and
  `dirichlet_multinomial.rs`'s production-TSV test; the script `tmp/b3_convert.py` (existence only).
- **Deliberately out of scope:** `paralog` and `parameter_estimation` (B4, B5); `examples/` (not
  converted, by plan); `src/float.rs` itself (reviewed as B1); points the B2 review settled — prose
  that names `.ln()`/`f64::ln` stays (the ban covers calls), and exact comparisons of a recorded
  value against libm's `ln` are left as they are.
- **Categories:** as for B2, the per-category fan-out of mutation-testing agents was **not** run:
  the dispatcher forbade `cargo` and source edits, which removes what the worktree isolation is for.
  The dispatcher's questions (a)–(e) were answered in one pass covering `reliability` (semantic
  equivalence, test challenge), `refactor_safety` (completeness), `naming`/doc accuracy, and
  `extras` (hot path, stable output). `unsafe_concurrency`, `module_structure`, `tooling`,
  `defaults`: no triggers.

## 2. Verdict

**Approve-with-changes.** All 157 converted calls mean what they meant before, apart from the
library that rounds them, and no call is left in the three scopes. The two hand edits are sound.
What needs attention is test prose and one test that has never measured what it says: the
rising-product accuracy test compares `float::ln` against `float::ln` of the same bits (Mi1), and
its doc calls that reference "correctly rounded" — which libm's `ln 3`, one of its own cells, is not.

## 3. Execution status

- **Commands run by the reviewer:**
  - `rg` completeness sweeps over `src/calling`, `src/genetics.rs`, `src/types.rs`: method calls
    `\.(ln|exp|powf|powi|log10|log2|log|ln_1p|exp_m1|sin|cos|tan|…|mul_add)\s*\(`, line-start
    chains `^\s*\.(…)\b` (multiline), `(f64|f32)::(…)` paths, `libm::`. Result: no calls; the only
    hits are prose (`generic.rs:2100` "`f64::ln`") and `genetics.rs`'s `libm::lgamma` wrapper.
  - Counts over the diff: 158 `float::…(` added, 158 `.ln(`/`.exp(`/… removed; one of each pair
    is the doc comment at `generic.rs:1430`, so **157 code calls — CHECKED-CORRECT**. 19 files gain
    `use crate::float`; `loop_parity.rs` is the 20th, doc only.
  - `uv run --no-project python tmp/review_2026-09-14_portable_float_B3/probe.py` — rounding of
    `ln 3`. Verbatim:
    ```
    literal ..7 0x3ff193ea7aad030b 1.0986122886681098
    literal ..6 0x3ff193ea7aad030a 1.0986122886681096
    apple math.log(3) 0x3ff193ea7aad030b
    exact ln3 1.09861228866810969139524523692252570464749055782274945173469
    0x3ff193ea7aad030b dist 9.071297235001530042216639799198193804826531E-17
    0x3ff193ea7aad030a dist -1.3133163257501600766255993562618212445173469E-16
    ulp 2.220446049250313080847263336181640625E-16 pos of exact from lower in ulps 0.59146509152680098617064719872616563302465790176783052570624
    ```
  - `uv run --no-project python tmp/review_2026-09-14_portable_float_B3/rising_cells.py` — the
    correctly rounded `ln` of each exact rising product in Mi1's grid; prints `23` cells and their
    bits (table in §8). Apple's `log` agrees with exact on all 23.
- **Gates reported by the dispatcher (not re-run):** fmt, clippy `-D warnings`, doc; container suite
  4,795 passed / 0 failed; macOS 4,794 passed / 0 failed; identity oracle's five checksums unchanged.
- **Commands not run:** all `cargo`, by instruction.
- **Needs verification:** 2 (Q2, Q3).

## 4. Open questions and assumptions

1. **Is B3 checked against the parity oracle at all, or only at B5?** The plan's parity oracle (A3's
   libm build: parameters `3f5e5e7f`, VCF `1e40bd1e`) is the one check that a *correct* conversion
   reproduces libm's answer bit for bit. The identity oracle being unchanged says B3 moved no output
   on its inputs; it does not say B3's outputs equal libm's. If calling from the libm parameters
   file on a B3 build still reaches unconverted `paralog`/`parameter_estimation` maths, the check
   has to wait for B5 — the commit message should say which. Affects §6 (d).
2. **Needs verification — how far below the floor the SSR score fell.** The dispatcher reports "3
   units". Not reproducible without `cargo`; the comment says "a few units", which is consistent.
   Affects Mi2.
3. **Needs verification — does the rising-product test still pass once its reference is exact?**
   libm returns `…030a` for `ln 3` (the `LOG_SPREAD` test asserts it and passed), while the
   correctly rounded value is `…030b`, so at cell (α 3, 1 copy) ours is one step from exact. Whether
   `ours_ulp <= theirs_ulp` holds there depends on the bits of `lgamma(4) − lgamma(3)`, which the
   reviewer could not compute. Affects Mi1.

## 5. Top 3 priorities

1. **Mi1** — `the_rising_product_is_at_least_as_close_to_exact_as_the_difference_of_two_lgammas`
   computes "ours" and "exact" as the same `float::ln` of bit-identical arguments, so half of what
   it claims cannot fail, and its doc calls libm's `ln` a correctly rounded reference.
2. **Mi2** — the loosened SSR floor test now asserts nothing its sibling on the identical fixture
   does not already assert more tightly, and its doc still says "at least".
3. **Nit on `LOG_SPREAD`'s doc** — say that the literal is libm's rounding, not the correctly
   rounded `ln 3`, so nobody "fixes" it back.

## 6. Findings

### Blocker

None.

### Major

None.

### Minor

#### Mi1: src/calling/genotype_prior/dirichlet_multinomial.rs:620-697 — The accuracy test's reference is the code's own computation, and not correctly rounded

**Confidence:** High
**Categories:** reliability (test challenge), extras (diff's quantitative claims)
**Pre-existing in shape; B3 touches both lines and changes the library the claim rests on.**

[dirichlet_multinomial.rs](../../../../src/calling/genotype_prior/dirichlet_multinomial.rs#L638-L697):

```rust
let exact = float::ln(exact_product as f64);          // line 656
...
let ours = { let mut rising = alpha; for step in 1..copies { rising *= alpha + f64::from(step); } float::ln(rising) };
```

Every partial product is an integer below the final one, which the loop keeps under 2⁵³, so every
multiplication is exact and `rising` has the same bits as `exact_product as f64`. The same function
on the same bits returns the same bits: **`ours_ulp` is 0 in every cell, on every platform, with any
`ln`.** So `ours_ulp <= theirs_ulp` and `worst_ours_ulp <= 1` cannot fail. What the test really
checks is `cells == 23` (CHECKED-CORRECT, §3) and that the `lgamma` difference is at least 100 steps
from `ln` of the product.

The doc (line 627) says the product's "logarithm is a correctly rounded reference". It never was
guaranteed to be, and after B3 it is demonstrably not on this grid: cell (α 3, 1 copy) is `ln 3`,
whose correctly rounded value is `0x3ff193ea7aad030b`, and libm returns `…030a`, 0.59 of a step away
(§3; the `LOG_SPREAD` test pins libm's bits). The test would report "ours is exact" there.

Which wrong implementations pass: any `ours` equal to `float::ln` of the exact product — including a
rising product that is correct but a `float::ln` that is several steps off.

**Fix:** compare against recorded, correctly rounded bits computed outside `f64` (table in §8), and
state the claim honestly: ours within one step of exact, theirs at least 100. If Q3 shows ours is
one step out where `lgamma` happens to be exact, the per-cell `ours_ulp <= theirs_ulp` must become
`ours_ulp <= theirs_ulp.max(1)`.

#### Mi2: src/calling/likelihood/ssr_emission.rs:1259-1281 — The loosened floor test is now subsumed by its sibling, and its doc still says "at least"

**Confidence:** High
**Categories:** reliability (test challenge), naming/documentation

The change `scored >= floor` → `scored >= floor * (1.0 - 1e-12)` is justified. The comment's
mechanism is **CHECKED-CORRECT**: `substitution_probability` returns `float::exp` of a sum of
`FlatEmission` log scores ([ssr_emission.rs:612-631](../../../../src/calling/likelihood/ssr_emission.rs#L612-L631)),
while `floor` uses `float::powi`'s repeated multiplication. Rounding along those two routes is
about 18 steps of 1e-3-sized logs plus one `exp`: of order 1e-15 relative, so 1e-12 leaves a margin
near a thousand. **It is not too loose for the claim's size**: every defect in the same-length term
the test is aimed at — a lost `(1 − ε)` factor, a wrong share, a placement counted twice — moves the
score by at least 1e-3 relative.

But `the_score_is_the_length_factor_times_the_letter_factor`
([ssr_emission.rs:1433-1463](../../../../src/calling/likelihood/ssr_emission.rs#L1433-L1463)) runs
the identical fixture (`a_model()`, the same 18-base tract, 6 repeats, `CAG`, ε = 1e-3) and asserts
`|scored − expected| ≤ expected · 1e-12` with the same `expected`. That implies both of this test's
assertions (`≥ floor·(1 − 1e-12)` and `≤ 1.05·floor`). Before B3 the strict `>=` was the one thing
the sibling did not say; now nothing is. And the doc, "scores at least
`same_length_share × (1 − ε)^length`", is no longer what is asserted.

**Fix** (either):
- keep the test as spec §12's named anchor, and make the doc say what it checks: "scores
  `same_length_share × (1 − ε)^length` to within rounding (1e-12 relative) and no more than 5% above
  it". Optionally spell the slack in steps, `floor * (1.0 - 64.0 * f64::EPSILON)`, which states the
  mechanism rather than a round number;
- or give it a fixture the sibling does not cover, where other stutter terms do contribute (a read
  one part-repeat away), so "at least the same-length term" is a non-trivial inequality again.

### Nits

- [generic.rs:2081-2089](../../../../src/calling/likelihood/generic.rs#L2081-L2089) `LOG_SPREAD`:
  the new literal is correct for its test — `…109_6` parses to `0x3ff193ea7aad030a`, libm's `ln 3`,
  and `…109_7` parsed to `…030b`, the correctly rounded value glibc and Apple return (§3). All 11
  other uses still make sense: 10 compare at `1e-12` with at most 6 × `LOG_SPREAD`, so the one-step
  change moves them by at most 1.3e-15, and one asserts a gap `> 3.0`. No shipped code carries a
  hand-written `ln 3` (`rg '1\.098|1_098|LN_3'` finds only this constant and two test-side
  `1.0986` bands); the shipped row divides by `ERROR_SPREAD_BASES` before one `float::ln`. Suggest
  one sentence in the doc: "libm's rounding of `ln 3`, one step below the correctly rounded
  `…030b`; the equality below is with `float::ln`, so this literal follows libm, not exact `ln 3`."
- [loop_parity.rs:31](../../../../src/calling/loop_parity.rs#L31): the edited doc line is 115
  characters; reflow to the file's 100.
- [loop_parity.rs:28-31](../../../../src/calling/loop_parity.rs#L28-L31): "a difference of units in
  the last place" between production's table-interpolated `ln`/`exp` and libm's is inherited, and
  unverifiable now that production is deleted. Not introduced by B3.
- Commit message: the plan lists B3 at 162 calls, the diff converts 157. The sweep finds none left,
  so the gap is in the counting (A1's table includes `libm::lgamma` in `genetics.rs`); say so.

### Answers to the dispatcher's questions

**(a) Semantic equivalence — every hunk.** Nothing changes meaning beyond std → libm/`float::powi`.
Checked specifically:

- *Unary minus:* `-heterozygote_weight.ln()` → `-float::ln(heterozygote_weight)` (method binds
  tighter than `-`, so both negate the log). Every negative literal receiver was already
  parenthesised — `(-7.0_f64).exp()` → `float::exp(-7.0)`, `(-6.0_f64 / 2.0).exp()` →
  `float::exp(-6.0 / 2.0)`, `(-f64::from(phred) / 10.0 * LN_10).exp().ln()` →
  `float::ln(float::exp(-f64::from(phred) / 10.0 * LN_10))` — so none flips sign.
- *Scalar multipliers:* `2.0 * x.ln()`, `-10.0 * x.log10()`, `f64::from(copies) * p.ln()`,
  `reads * (…).ln()` (`ssr.rs:733`): the multiplication still applies to the log.
- *Clamps:* `(1.0 - inbreeding).max(FLOOR).ln()`, `(alpha / total).max(FLOOR).ln()`,
  `hom_ref.max(FLOOR).ln()` → `max` stays inside the argument. `(-10.0 * tail.log10()).max(0.0)`
  (`artifact_correction.rs:140`) keeps `max` outside. `(self.scale * mean_log_error.exp()).max(MIN)`
  (`likelihood/mod.rs:511`) unchanged in shape.
- *Chains:* `((a - l).exp() + (b - l).exp()).ln()` → `float::ln(float::exp(a - l) + float::exp(b - l))`
  in `log_sum_exp_2` and the Hardy–Weinberg test; `(charged / 0.5_f64).ln()` → `float::ln(charged / 0.5)`.
- *Multi-line receivers:* the four in `generic.rs` (810, 829-832, 904-906, 918-921), `ssr.rs:733-738`,
  the `marginal_probability(…).get().exp()` chain and the `two_sided_binomial_tail(…)` pair: whole
  receiver moved into the argument, operators inside unchanged.
- *Closures:* `|entry| entry.exp()` over `&f64` → `float::exp(*entry)`; `|error| error.ln()` →
  `float::ln(*error)`; `|value| (value - largest).exp()` → `float::exp(value - largest)` (`&f64 − f64`
  is `f64`, no deref needed); `|p| (p.get() - largest).exp()` unchanged shape.
- *Literal suffixes:* every dropped suffix is a literal passed straight to a `fn(f64…)` parameter or
  in arithmetic whose other operand is `f64` (`float::ln(0.5)`, `float::powi(2.0, -53)`,
  `float::log10(7.0 / 3.0)`, `float::powf(10.0, -1.3)`, `float::ln(0.95 * 0.25 + 0.05 * 0.2)`).
  None can infer `f32` (no `f32` parameter exists in `float`) and no integer literal lost a suffix.
  Where a suffixed binding stays (`let scale = 2.5_f64;`, `(own, other) = (0.97_f64, 0.03_f64)`), it
  was kept.
- *`powi` exponents:* all `i32` already — `index.abs_diff(at) as i32`, `index as i32`,
  `(at as i32 - middle as i32).abs()`, `tract.len() as i32 [- 1]`, `-53`, `-52`, `3`.
  One behavioural nuance: std `powi` on literal operands could be folded at build time with the
  build machine's `pow` (A2: 15 of 20 cases differ); `float::powi` cannot. All such sites are tests
  (`2.0.powi(-53)` is exact either way).

**(b) Completeness.** No transcendental call remains in `src/calling`, `src/genetics.rs`,
`src/types.rs` (§3). No `.map(f64::ln)`-style paths. The 37 shipped calls A1 lists for these files
(excluding `libm::lgamma`) are all in the diff: dirichlet-multinomial 10, artifact correction 6,
quality 6, generic 4, genetics 3, Hardy–Weinberg 3, likelihood `mod` 2, and one each in
`genotype_table`, `summarise_condition`, `ssr`, `ssr_emission`.

**(c) Recorded-answer tests.**

| test | how it compares | sees a one-step change? |
|---|---|---|
| `genotype_table_parity` | log multinomial coefficients **as bits** against production's recorded bits | **yes** — and `log_factorial` now sums libm `ln`s and still matches bit for bit, a real confirmation |
| `loop_parity` | called genotypes, iteration count, convergence (discrete), keyed by a digest of literal inputs | only if it flips a call on its 9 tables |
| `quality_parity` | `1e-4` Phred (site), `1e-3` Phred (correction); measured disagreement 1.3e-5 and 2.5e-5 | no, by design — production used table-interpolated `ln`; a libm step is ~1e-13 Phred |
| dirichlet-multinomial TSV | `1e-9` nats absolute | no, by design (production's `lgamma` cancellation) |
| `types.rs` Phred 3000 | `f32` bits of `−10·ln(1e-300)/ln 10` | no — the narrowing to `f32` absorbs any `f64` step; it tests width, which is unchanged |

None of these tolerances was loosened by B3. None became tautological through the conversion: the
tests that now call `float::` on both sides (`generic.rs` fixtures built from `float::ln(0.5)` etc.,
the `SlipCosts`-style reconstructions) still spell the quantity by a different route from the code.
The one test whose two sides are literally the same computation is Mi1, and it was so before B3.

**(d) Shipped hot paths.**

- *No extra work:* every conversion is one call for one call. `generic.rs` still takes one `ln` per
  copy count (810, 904-906, 918-921) and one per genotype in the two-column walk (829-832); `ssr.rs`
  one per genotype × observation; `summarise_condition.rs:365` one `exp` per genotype; the quality
  fold one `exp` per copy count and one `ln` per rescale. No hoisted value was split.
- *Inlining:* `float::ln` is `#[inline]` but `libm::log` is not; release uses `lto = "fat"`,
  `codegen-units = 1`, so it can inline across the crate boundary. The `profiling` and `soak`
  profiles set `lto = false`, where it cannot — profiles taken there overstate libm's cost relative
  to release. Worth one line in B7's report.
- *Ties:* a genotype decided by an exact tie that symmetry creates (two genotypes whose scores come
  from the same inputs by the same expression) stays tied under any deterministic library. A
  coincidental near-tie between differently computed scores can flip on a one-step change; no unit
  test observes that, `loop_parity` would only on its nine tables, and the identity oracle (reported
  unchanged) covers its own inputs. The parity oracle (Q1) is the check that covers it.
- `score_best_genotype`: `best.min(1 − ε)` keeps `1 − best` ≥ `2⁻⁵²` exactly, so `log10` stays
  finite under libm as under std.

**(e) Doc comments.** `loop_parity.rs`'s edit is accurate ("ng calls libm's through
`crate::float`"). `generic.rs:1429-1433` now reads "Computed from `float::ln(3.0)` instead, this
would be a test of the logarithm" — accurate. Inaccurate after B3: `dirichlet_multinomial.rs:627`
"correctly rounded reference" (Mi1) and `ssr_emission.rs:1260` "at least" (Mi2). `genetics.rs`'s
`lgamma` doc and `float.rs`'s module doc agree with the code. Prose naming `.ln()`/`f64::ln`
(`generic.rs:2100`, `to_run_parameters.rs:1429`, `summarise_condition.rs:3115`) is left per B2.

## 7. Out of scope observations

- `quality_parity.rs:271-277` and `:634-636` quote disagreements "measured" under std. After B3 they
  may differ in the last digits; the tolerances have an order of magnitude of margin, so nothing
  fails, but the figures describe the pre-B3 build.

## 8. Missing tests to add now

**`fill_random_mating_log_priors`' rising product (`dirichlet_multinomial.rs`)**

- `the_rising_product_is_within_one_step_of_the_exact_logarithm` — replaces the tautological half of
  Mi1's test. Input class: the 23 exact-product cells. Catches an `ln` more than one step out, or a
  rising product that is not the exact one. Reference bits computed by 80-digit `Decimal.ln` and
  correctly rounded to `f64` (`tmp/review_2026-09-14_portable_float_B3/rising_cells.py`):

  ```rust
  /// ln of the exact rising product, correctly rounded, computed outside f64.
  const EXACT_LN_OF_RISING_PRODUCT: [(u64, u32, u64); 23] = [
      (3, 1, 0x3ff193ea7aad030b), (3, 2, 0x4003e116bcd39e7d), (3, 3, 0x4010609bdc65328b),
      (3, 4, 0x40178b5edaefba8c), (3, 5, 0x401f53fb867c5122), (3, 6, 0x4023d2aa530d136e),
      (3, 7, 0x402837a4f1b85430), (3, 8, 0x402cd29160a5a976),
      (125, 1, 0x4013503179f229e7), (125, 2, 0x40235445e15a44a0), (125, 3, 0x402d047f2b8a5bc1),
      (125, 4, 0x40335c5e3d8bea8d), (125, 5, 0x4038387ae7f95840), (125, 6, 0x403d1691a4d6346c),
      (125, 7, 0x4040fb4f49af7342),
      (2001, 1, 0x401e67d6037b19ca), (2001, 2, 0x402e6817801facb6), (2001, 3, 0x4036ce42b96396b7),
      (2001, 4, 0x403e689a68ab805f),
      (6001, 1, 0x4021663ca3fdbd25), (6001, 2, 0x403166478f7eb5ff), (6001, 3, 0x403a197bb8084ca2),
      (6001, 4, 0x4041665d65923982),
  ];
  ```

  Assert `|ours.to_bits() − exact_bits| ≤ 1` per cell and keep `worst_theirs_ulp ≥ 100` against the
  same table. At (3, 1) ours is expected to be exactly one step out (libm `…030a`).

**`StutterSubstitutionEmission` (`ssr_emission.rs`)**

- `a_read_one_part_repeat_away_scores_at_least_its_same_length_term` — only if Mi2's second fix is
  chosen: a fixture where a non-same-length term contributes, so `scored > floor` is strict by more
  than rounding and the "dominated" claim is tested where it can fail.

## 9. What's good

- `genotype_table_parity` compares the log multinomial coefficients as bit patterns against
  production's recorded bits, and it still passes: that is direct evidence libm's `ln` agrees with
  std's on the small integers `log_factorial` sums, not an inference from tolerances.
- The new comment at [ssr_emission.rs:1271-1272](../../../../src/calling/likelihood/ssr_emission.rs#L1271-L1272)
  names the mechanism (sum of logs against repeated multiplication) instead of only widening a
  bound, and the mechanism matches the code.
- The `LOG_SPREAD` edit changed the literal rather than the test's equality with `float::ln`, which
  keeps the fixtures' constant tied to what the code computes.

## 10. Commands to re-verify

- Reviewer ran: the `rg` sweeps in §3; `uv run --no-project python tmp/review_2026-09-14_portable_float_B3/probe.py`;
  `uv run --no-project python tmp/review_2026-09-14_portable_float_B3/rising_cells.py`.
- New, for the author:
  - `./scripts/dev.sh cargo test --lib calling::genotype_prior::dirichlet_multinomial` after Mi1's
    table, on macOS and in the container (answers Q3).
  - `./scripts/dev.sh cargo test --lib calling::likelihood::ssr_emission` after Mi2.

**PROJECT_STATUS.md** was not updated: as for B2, the plan schedules its entry for step C3, the
feature has no block yet, and the dispatcher asked for the report only.

## Fixes applied (2026-09-14)

| finding | what was done |
|---|---|
| Mi1 rising-product test compared `float::ln` with itself | The reference is now `EXACT_LN_OF_RISING_PRODUCT`: 23 correctly rounded bit patterns, recomputed by the implementer independently (100-digit decimal logarithm, nearest `f64` chosen by exact comparison with both neighbours) and identical to the review's §8 table. The test asserts the rising product stays exact, each cell `ours ≤ max(theirs, 1)` steps from the reference, worst ours ≤ 1, 23 cells and worst lgamma difference ≥ 100. It passes, so open question 3 is answered: libm is within one step on every cell. The doc no longer calls libm's `ln` correctly rounded and says why (3, 1) needs the one-step allowance |
| Mi2 loosened floor test | Doc now says "at least, to within rounding". The fixture stays; the sibling test's tighter bound is noted in the review, not duplicated |
| nit `LOG_SPREAD` | Doc says the literal is libm's `ln 3`, one unit below the correctly rounded value, and must not be rounded back |
| nit `loop_parity.rs` width, production's table claim | Rewrapped; the unverifiable "units in the last place" size is removed |
| nit commit message 157 vs 162 | Explained in the commit message |
| open question 1 | The parity oracle against the libm build waits for B5: `call-from-psps` still reaches unconverted `paralog` code until B4, and the fit until B5. This step is checked by the identity oracle (unchanged) and the recorded-answer tests (all pass) |

Validation after the fixes: container fmt, clippy `-D warnings`, doc clean; 4,795 passed, 0 failed, 4 ignored. macOS `calling::` 1,178 passed, 0 failed; macOS full suite before the fixes 4,794 passed, 0 failed. The identity oracle ran on the release build of the conversion before the fixes, which changed only tests and comments; all five checksums unchanged.
