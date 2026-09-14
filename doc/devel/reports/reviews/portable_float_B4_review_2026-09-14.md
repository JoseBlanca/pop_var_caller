# Code Review: portable_float_B4
**Date:** 2026-09-14
**Reviewer:** rust-code-review skill (orchestrator, single-pass — see §1)
**Scope:** uncommitted diff of plan step B4 (`doc/devel/implementation_plans/portable_float.md`, Milestone B): std transcendental calls in `src/paralog/` replaced by `crate::float`, and the production-parity tests loosened from bit equality
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** `git diff HEAD -- src/` in the worktree `pop_var_caller-portable-float`, on
  top of `06bb201f` (B3). 6 files, +117 / −76.
- **Reviewed against:** `HEAD` = `06bb201f` plus uncommitted changes.
- **In-scope files:** `src/paralog/{calibrate, calibration, locus_score, mod, prior,
  production_parity}.rs`; the callee `src/float.rs`; the fixture
  `src/paralog/testdata/production_parity_answers.tsv` (read, to re-derive the diff's numbers).
- **Deliberately out of scope:** `parameter_estimation` (B5); `src/float.rs` itself (B1); points the
  B2 and B3 reviews settled — prose naming `.ln()` stays, and exact comparisons of a value against
  libm's own output are left exact.
- **Categories:** as for B2 and B3, the fan-out of mutation-testing agents was **not** run: the
  dispatcher forbade `cargo` and source edits, which removes what worktree isolation is for. Instead
  the reviewer ported the calibration arithmetic and libm 0.2.16's `exp`, `log` and `log1p` to
  Python and re-ran the parity tests' probability comparisons against the frozen answers (§3). One
  pass covered `reliability` (equivalence, test challenge, tolerance design), `refactor_safety`
  (completeness), `naming`/doc accuracy, and `extras` (hot path, diff-matches-intent).
  `unsafe_concurrency`, `module_structure`, `tooling`, `defaults`: no triggers.

## 2. Verdict

**Approve-with-changes.** All 36 converted calls mean what they meant before, and no std maths call
remains in `src/paralog`. The loosening of the parity tests is the right move and the tolerance is
wide enough: libm moves 17 of the 518 probability values those tests compare, by at most 22 units in
the last place, a relative 3.6e-15, which is 277 times inside the 1e-12 allowed. What needs fixing
before commit: the new `within_rounding` accepts an infinity against *any* finite value (M1), the
constant's doc misstates its own measurement (Mi1), and four doc comments still say the gap is zero
or the comparison is by bit pattern (Mi2).

## 3. Execution status

- **Commands run by the reviewer:**
  - Completeness sweeps over `src/paralog`: method calls
    `\.(ln|exp|powf|powi|log10|log2|log|ln_1p|exp_m1|sin|cos|tan|sinh|cosh|tanh|atan|atan2|asin|acos|cbrt|hypot|mul_add)\s*\(`,
    line-start chains `^\s*\.(…)\b` (multiline), and `(f64|f32)::(…)` paths. **No hits.**
  - Counts over the diff: 36 `float::` calls added (22 `ln`, 12 `exp`, 2 `ln_1p`) against 38
    `.ln(`/`.exp(`/`.ln_1p(` removed, two of which are prose (`locus_score.rs:176`, `:249`). **36
    code calls — CHECKED-CORRECT.** The plan's "38 calls" includes those two.
  - `uv run --no-project python tmp/review_2026-09-14_portable_float_B4/calibration_port.py src/paralog/testdata/production_parity_answers.tsv`
    — prior, curve, cut, fallback and posterior ported from `prior.rs`/`calibration.rs`/`calibrate.rs`,
    with libm 0.2.16's `exp` and `log` transcribed from the registry source, run on the tests'
    seeds. The port's `exp` disagrees with Apple's on 19,300 of 200,000 arguments (about one in
    ten), which is A2's measured libm-against-platform rate, and it reproduces the dispatcher's two
    observed failures to the unit. Verbatim, abridged to the differing rows:
    ```
    port sanity: exp differs from Apple on 19300/200000, log on 9333/200000
    keys not found: 0 []
    compared 518: bit-identical 501, differing 17
      prior_and_curve/nothing is duplicated/prior_probability: libm 1.2301354082845576e-12 (3d75a40990d7962a) glibc 1.230135408284562e-12 (3d75a40990d79640) 22 units, relative 3.61e-15, absolute 4.44e-27
      prior_and_curve/nothing is duplicated/q_of_lr/30.0: libm 0.06747724424480839 (3fb146304d4172b0) glibc 0.06747724424480828 (3fb146304d4172a8) 8 units, relative 1.65e-15, absolute 1.11e-16
      verdict/nothing is duplicated/0.01/-100.0/posterior: libm 2.78578079052782e-56 (3465dbb8aeebc0e3) glibc 2.7857807905278194e-56 (3465dbb8aeebc0e2) 1 units, relative 1.63e-16, absolute 4.53e-72
      [… the same posterior at the other four targets …]
      verdict/one locus in ten/0.01/-30.0/posterior: libm 1.1325748281984857e-14 (3d0980d970a46153) glibc 1.132574828198486e-14 (3d0980d970a46155) 2 units, relative 2.79e-16, absolute 3.16e-30
      [… the same posterior at the other four targets …]
      verdict/beyond the histogram's range/0.01/0.0/posterior: libm 0.2622832058276867 (3fd0c93f7fd47a14) glibc 0.26228320582768677 (3fd0c93f7fd47a15) 1 units, relative 2.12e-16, absolute 5.55e-17
      [… the same posterior at the other four targets …]
    fallback/empty: compared 92, differing 0
    within_rounding(inf, 0.5) = True
    within_rounding(-inf, inf) = True
    within_rounding(1e308*10, 1e308) = True
    within_rounding(0.0, -0.0) = True
    everything-duplicated pi = 1 - 1.96e-13; relative 1e-12 admits 1 - pi anywhere in [0, 1.2e-12]; within(1.0, pi) = True
    ```
    Every cut (`lr_threshold_for_fdr`, `lr_threshold`: 85 values) is bit-identical.
  - `uv run --no-project --with mpmath python tmp/review_2026-09-14_portable_float_B4/series_sweep.py`
    — `the_series_shortcut_lands_on_log1ps_double_or_its_neighbour` with libm's `log1p`, Apple's,
    and a correctly rounded one. Verbatim:
    ```
    libm     threshold 0.00014: checked 8080, differing 3 (1 in 2693), widest 1 ulp, pairs where x moves hi: 293, differing among those 1 in 97.7
    libm     threshold 0.0014: checked 8104, differing 21 (1 in 386), widest 8641 ulp, pairs where x moves hi: 317, differing among those 1 in 15.1
    apple    threshold 0.00014: checked 8080, differing 3 (1 in 2693), widest 1 ulp, pairs where x moves hi: 293, differing among those 1 in 97.7
    apple    threshold 0.0014: checked 8104, differing 21 (1 in 386), widest 8641 ulp, pairs where x moves hi: 317, differing among those 1 in 15.1
    correct  threshold 0.00014: checked 8080, differing 3 (1 in 2693), widest 1 ulp, pairs where x moves hi: 293, differing among those 1 in 97.7
    correct  threshold 0.0014: checked 8104, differing 21 (1 in 386), widest 8641 ulp, pairs where x moves hi: 317, differing among those 1 in 15.1
    libm vs apple log1p differ on 0 of 1010 swept x
    ```
- **Gates reported by the dispatcher (not re-run):** fmt, clippy `-D warnings`, doc; container suite
  4,795 passed / 0 failed; macOS 4,794 passed / 0 failed; identity oracle's five checksums unchanged;
  parity oracle (call-from-psps with A3's libm-build parameters) reproduces A3's VCF records, md5
  `1e40bd1e`.
- **Commands not run:** all `cargo`, by instruction.
- **Needs verification:** 1 (Q2).

## 4. Open questions and assumptions

1. **Should the tolerance be in units in the last place rather than relative?** The measured
   differences are 1 to 22 units. A units bound handles zero, values near 1 and infinities
   uniformly, which a relative bound does not (M1, Mi4). The recommendation is yes; the author
   decides. Affects M1, Mi4.
2. **Needs verification — the scorer's 3.6e-12 nats.** Not re-derived: that needs the whole scorer
   ported, not only the calibration. It is consistent with the calibration result (tens of units at
   log-likelihood magnitudes in the hundreds, where a unit is 5.7e-14). Affects the 1e-10 pin's
   comment (`production_parity.rs:539-550`).

## 5. Top 3 priorities

1. **M1** — `within_rounding` passes `(inf, 0.5)` and `(-inf, inf)`, so a defect that sends π, a
   q-value, a posterior or a cut to infinity is invisible to four tests.
2. **Mi2** — `locus_score.rs` says the differential's gap is "exactly zero" twice, one paragraph
   after saying it is 3.6e-12; `production_parity.rs`'s header still says "by bit pattern".
3. **Mi1** — the constant's doc gives π's gap as a relative 1.8e-15; it is 3.6e-15, and 15 further
   values moved that the doc does not mention.

## 6. Findings

### Blocker

None.

### Major

#### M1: src/paralog/production_parity.rs:922-929 — `within_rounding` accepts an infinity against any finite value

- **Confidence:** High (reproduced in the Python transcription, §3).
- **Problem:** the relative branch is `(ours - theirs).abs() <= 1e-12 * ours.abs().max(theirs.abs())`.
  When either side is infinite the right-hand side is `inf`, and the left is `inf`, so the test is
  `inf <= inf`, which is true. `within_rounding(f64::INFINITY, 0.5)`,
  `within_rounding(f64::NEG_INFINITY, f64::INFINITY)` and `within_rounding(f64::INFINITY, f64::MAX)`
  all return `true`. The doc says infinities are covered by "equal bits", which is the intent; the
  arithmetic branch then accepts them anyway.
- **Why it matters:** these tests are the only check that the copied calibration computes what
  production did. A port defect that produces an infinite π (dividing by a zero total), q-value,
  posterior or cut now passes all four renamed tests. Before B4 the `to_bits` comparison caught it.
- **Suggested fix:** make non-finite values match exactly, and measure the tolerance in units in the
  last place, as the series test in `locus_score.rs` already does (Q1). The widest measured
  difference is 22 units.
  ```rust
  /// How far apart, in units in the last place, ng's probability may sit from production's. The
  /// widest measured on 2026-09-14 was 22 (π, stream "nothing is duplicated").
  const PROBABILITIES_MAY_DIFFER_BY_UNITS: u64 = 64;

  /// Both `NaN`; or both the same infinity; or finite, of one sign, and within
  /// [`PROBABILITIES_MAY_DIFFER_BY_UNITS`] of each other (`+0.0` and `−0.0` are equal).
  fn within_rounding(ours: f64, theirs: f64) -> bool {
      if ours.is_nan() || theirs.is_nan() {
          return ours.is_nan() && theirs.is_nan();
      }
      if !ours.is_finite() || !theirs.is_finite() {
          return ours == theirs;
      }
      ours == theirs
          || (ours.is_sign_negative() == theirs.is_sign_negative()
              && ours.to_bits().abs_diff(theirs.to_bits()) <= PROBABILITIES_MAY_DIFFER_BY_UNITS)
  }
  ```
  If the relative form is kept, add `ours.is_finite() && theirs.is_finite() &&` in front of the
  arithmetic branch. Either way add the test in §8.

### Minor

#### Mi1: src/paralog/production_parity.rs:913-920 — the constant's doc misstates π's gap and lists 2 of the 17 values that moved

- **Confidence:** High (computed, §3).
- **Problem:** "π, 22 units in the last place away (a relative 1.8e-15)". The units are right; the
  relative figure is not. π is `1.2301354082845576e-12` against `1.230135408284562e-12`, which is
  4.44e-27 apart, a relative **3.61e-15**. The 1.8e-15 is that difference divided by the *sum* of the
  two values. 22 units cannot be 1.8e-15 of any `f64`: one unit is between 1.1e-16 and 2.2e-16 of the
  value, so 22 units is at least 2.4e-15. Separately, "the first values the old comparison refused"
  is literally true, because each test stops at its first failure. But a reader takes it as the
  whole list. The full list is 17 values out of 518: π (22 units); the q-value at a ratio of 30 in the
  same stream (8 units); and three posteriors, each at all five targets (1, 2 and 1 units).
- **Why it matters:** these figures are what the 1e-12 tolerance is justified from. A reader
  checking the headroom finds it halved, and does not know that a q-value moved.
- **Suggested fix:**
  ```rust
  /// … The values it moves, measured 2026-09-14 by re-running these tests' arithmetic with libm's
  /// `exp` and `log`: 17 of the 518 probabilities compared, none by more than 22 units in the last
  /// place — π on the "nothing is duplicated" stream (22 units, a relative 3.6e-15), that stream's
  /// q-value at a ratio of 30 (8 units), and three posteriors at every target (1, 2 and 1 units).
  /// No cut moved. …
  ```

#### Mi2: src/paralog/locus_score.rs:58-60, :515-519; src/paralog/production_parity.rs:40-41, :83-84, :433-435, :458-460 — doc comments still describe the gap as zero, the comparison as bit-exact, and the series as the only departure

- **Confidence:** High.
- **Problem:** the diff added a paragraph at `locus_score.rs:51-55` saying the gap "went from exactly
  zero to 3.6e-12 nats". Five lines later, `:58-60` still says **"The gap it measures is exactly zero,
  so on everything that differential draws the two trees still return the same `f64`."** And
  `log_add_exp`'s doc at `:515-519` says the widest gap "is **exactly zero**". In
  `production_parity.rs`, the module header at `:40-41` says the prior, curve, cut and verdict are
  compared "**by bit pattern**", and `:83-84` says "again by bit pattern". The in-loop comment at
  `:433-435` calls the series shortcut "this port's one departure from production's arithmetic",
  and the constant's doc at `:458-460` calls it "the one place the two trees now round differently".
  Since B4 there are two departures.
- **Why it matters:** the series shortcut was released from the copy guard on the strength of "the
  measured gap is zero". A reader who takes that sentence at face value will read the new 3.6e-12
  as a regression of the shortcut, when it comes from libm.
- **Suggested fix:**
  - `locus_score.rs:58-60`: "Until 2026-09-14 the gap it measured was exactly zero; since the maths
    moved to libm it is 3.6e-12 nats, all of it from `exp`, `ln` and `log1p` rounding differently
    from glibc's."
  - `locus_score.rs:515-519`: "…against production's frozen answers, and before the maths moved to
    libm the widest gap on any field was **exactly zero** — so the shortcut itself moves nothing
    there. (Since 2026-09-14 the gap is 3.6e-12 nats, from libm.)"
  - `production_parity.rs:40-41` and `:83-84`: "…the prior, the curve and the verdict's posterior
    within [`PROBABILITIES_MAY_DIFFER_RELATIVELY_BY`]; the flags, the convergence, the counts and
    the cut exactly" (adjusted to Mi3's outcome).
  - `production_parity.rs:433-435` and `:458-460`: "one of two departures … the other, since
    2026-09-14, is that `exp`, `ln` and `log1p` are libm's".

#### Mi3: src/paralog/production_parity.rs:1000-1008, :1103-1110, :1285-1291 — the cut is loosened although it did not move and is a decision

- **Confidence:** High.
- **Problem:** `lr_threshold_for_fdr` returns a bin centre (`prior.rs:292-298`): `lo + (i + 0.5) ·
  width / n`, computed with no transcendental call. The maths library can change *which* bin is
  chosen, but never the bits of a given bin's centre. All 85 cuts compared are bit-identical under
  libm (§3). Comparing them with `optional_within_rounding` passes exactly the cases `to_bits`
  would pass — two different centres are never 1e-12 apart; the closest centre to zero is ±0.05,
  and neighbours are 0.1 or 0.6 apart. But it contradicts the constant's own doc, which says every
  decision is still compared exactly. It also means a future change that interpolates the cut would
  get a tolerance nobody chose.
- **Why it matters:** it blurs the rule the diff sets out: probabilities within rounding, decisions
  exactly. The cut is the decision the VCF header records.
- **Suggested fix:** restore the bit comparison at the three sites and name the cut among the
  decisions in the doc.
  ```rust
  assert_eq!(
      ours.map(f64::to_bits),
      theirs.map(f64::to_bits),
      "{case}: the cut for a target of {target} differs — ng {ours:?}, production {theirs:?}. \
       The cut is a bin centre, so no maths library moves it; only a different bin does",
  );
  ```

#### Mi4: src/paralog/production_parity.rs:924-929 — a relative tolerance is blind to π near 1

- **Confidence:** High (computed, §3).
- **Assumptions:** that the author keeps the relative form; the units form in M1 removes this.
- **Problem:** in the "everything is duplicated" stream, π is `1 − 1.96e-13`. A relative 1e-12 is an
  absolute 1e-12 there, so any π from `1 − 1.2e-12` up to **exactly 1.0** passes. The information in
  a probability near 1 is in its distance from 1, and that distance is 6 times smaller than the
  tolerance. The verdict test catches π = 1.0 indirectly, because `posterior` returns `None` at
  π ≥ 1. A π of `1 − 1e-12` would pass everything.
- **Why it matters:** low. The values at stake are below any cut. But the blind spot is not stated,
  and the units bound in M1 closes it for free: 64 units at 1 is 7e-15.
- **Suggested fix:** take M1's units form. Otherwise state the blind spot in the constant's doc.

### Nits

- `src/paralog/prior.rs:29-30` — the inserted clause makes line 30 117 characters, where the file
  wraps at 100; rewrap. The sentence also says ng "changes nothing else, except" the maths calls,
  but ng also added the `use crate::float;` line.
- `src/paralog/mod.rs:20-23` — "Three of the five files … byte for byte apart from their maths
  calls, which go through `crate::float`" suggests all three had calls converted. `coverage_model`
  and `model_params` have none; only `locus_score` changed. (And `locus_score` stopped being byte for
  byte on 2026-09-09, which predates this diff.) Suggest: "…`coverage_model` and `model_params` byte
  for byte; `locus_score` apart from the series shortcut and, since 2026-09-14, its maths calls
  through `crate::float`."
- `src/paralog/locus_score.rs:53-54` — "The reason and the measurements are
  `doc/devel/implementation_plans/portable_float.md`": the measurements are in the A2 and A3 reports,
  which the plan links. Name those.
- `src/paralog/production_parity.rs:545-550` — the 1e-10 pin is 28 times the measured 3.6e-12. Every
  input is fixed and libm is pinned `=0.2.16`, so the gap is the same on every run and machine. A pin
  at 1e-11 would be equally stable and would report a change 10 times smaller. Optional.
- Plan `portable_float.md` B4 says "38 calls"; the code has 36, plus two mentions in prose. Say so
  when ticking the step, as B3's commit did.

### Answers to the dispatcher's questions

**1. The tolerances.**
- *Is a relative 1e-12 tight enough to catch a real port defect?* In the ranges these tests probe,
  yes, except for infinities (M1). A defect worth catching — a wrong bin centre, a wrong sign in the
  logit, a mis-clamped π — moves a posterior in the transition region (ratios −1 to 18) by a relative
  1e-3 or more. All six probes there are compared.
  - **Near 0:** posteriors reach 2.8e-56. A relative bound is the right measure there, and the
    libm-against-glibc gap is 1.6e-16 relative. A posterior that underflows to 0.0 on both sides
    matches by bits (`-1e6`).
  - **Near 1:** a relative bound is loose (Mi4).
  - **Near-zero q-values:** these are where a relative bound could be *too tight*. A q-value like
    7.7e-13 is a mean of `1 − posterior` terms that each round in steps of 1.1e-16, so one unit's
    change in one `exp` could move it by a relative 1e-4. None did (§3), and because the inputs are
    fixed and libm is pinned the result cannot change between runs. So this is brittleness on a libm
    upgrade, not flakiness.
- *Thresholds as ratios near 0?* The cuts are bin centres, never closer to 0 than ±0.05 and never
  moved (Mi3).
- *Opposite signs and zero?* Correct: opposite signs make the difference at least the larger
  magnitude, so they fail unless both are zero, and `+0.0`/`−0.0` pass.
- *NaN?* Correct: both NaN pass, one NaN fails.
- *Is 1e-10 nats justified?* Yes, as a stable pin (see Nits on tightening). The measured 3.6e-12 is
  not re-derived here (Q2).
- *A bit-exact comparison left that could flake?* No. Four checks stay exact:
  - `mod.rs:256-263` compares `posterior` with the same expression through the same `float`
    functions, and Rust does not fuse multiply-adds, so this is the same computation twice;
  - `production_parity.rs:738-741` compares the neutral `0.0`;
  - the flags, the convergence flag and the counts are integers or booleans;
  - `scoring_is_deterministic` compares a function with itself.

**2. Doc comments.** `calibration.rs:19-21` and `calibrate.rs:6,16` are accurate. See Mi1, Mi2 and
the Nits for the rest.

**(a) Semantic equivalence.** All 36 hunks are equivalent:
- a method call binds tighter than unary minus, so `-(p * (1.0 - p)).ln()` is
  `-float::ln(p * (1.0 - p))`, and `ln_normal`'s `- sigma.ln() -` is unchanged in precedence;
- `(-log_odds).exp()` becomes `float::exp(-log_odds)`, and `0.5f64.ln()` becomes `float::ln(0.5)`,
  which infers `f64`;
- in `map(|x| (x - max).exp())`, `x` is `&f64`, and `&f64 - f64` is an `f64` through std's
  reference `Sub` impl, so `float::exp(x - max)` type-checks to the same value;
- `.sum::<f64>().ln()` becomes `float::ln(….sum::<f64>())`.

**(b) Completeness.** No std transcendental call, chained or by path, remains in `src/paralog` (§3).

**(c) The series test.** Its numbers did not move: libm's `log1p` equals Apple's on all 1,010 swept
`x`, and the sweep gives 3 differing pairs of 8,080, widest 1 unit, under libm, Apple and a
correctly rounded `log1p` alike. Its doc was not changed by B4, but two of its claims were already
inaccurate (see §7).

**(d) Hot path.** The number of calls is unchanged. `score_locus_for_paralogy`, `h2_log_likelihood`,
`log_add_exp` and `LogSumExp` make the same `exp`, `ln` and `log1p` calls on the same branches, and
the diff adds no allocation or branch. The cost per call is A3's, where libm's `exp` is up to 3.3 ns
slower. It may come out lower here, because the release profile's fat LTO can inline libm into
`log_add_exp`, while a call into glibc or libSystem is never inlined. The series threshold was chosen
from `log1p`'s cost; libm's `log1p` for `x < 0.41` is a division and a seven-term polynomial, so the
shortcut should still pay, but by a different margin. B7's benches are where that shows.

## 7. Out of scope observations

- `src/paralog/locus_score.rs:785-790, :819-828` (pre-existing, not touched by B4). Measured with the
  sweep in §3:
  - "Below about 1 in 500 of these pairs land on the neighbouring double" holds but understates the
    margin: 3 of 8,080, 1 in 2,693. Only 293 of the 8,080 pairs have an `x` large enough to change
    `hi` at all, and among those, 1 in 98 lands on the neighbour.
  - "Raising the threshold one decade takes the disagreement rate past 1 in 20" is **wrong** by the
    test's own count: 21 of 8,104, 1 in 386. The figure is only near 1 in 20 (1 in 15) among the pairs
    where `x` moves `hi`.
  - "the size assertion starts failing" is right: the widest gap becomes 8,641 units.

  Follow-up: restate both rates with their denominator, and consider pinning the count of pairs
  where `x` moves `hi`, so the thousands of trivially equal pairs cannot dilute the rate.
- No `PROJECT_STATUS.md` block exists for the portable-float plan, and the B2 and B3 reviews did not
  add one; left for the author to decide at the milestone's end.

## 8. Missing tests to add now

`within_rounding` (in `production_parity.rs`):

```rust
#[test]
fn within_rounding_refuses_an_infinity_against_anything_but_itself() {
    assert!(!within_rounding(f64::INFINITY, 0.5));
    assert!(!within_rounding(0.5, f64::NEG_INFINITY));
    assert!(!within_rounding(f64::INFINITY, f64::NEG_INFINITY));
    assert!(!within_rounding(f64::INFINITY, f64::MAX));
    assert!(within_rounding(f64::INFINITY, f64::INFINITY));
    assert!(within_rounding(f64::NAN, f64::NAN));
    assert!(!within_rounding(f64::NAN, 0.0));
    assert!(within_rounding(0.0, -0.0));
    assert!(!within_rounding(1e-300, -1e-300));
}

#[test]
fn within_rounding_accepts_the_measured_gap_and_refuses_a_real_one() {
    // π on the "nothing is duplicated" stream: libm against glibc, 22 units apart.
    let glibc = f64::from_bits(0x3d75_a409_90d7_9640);
    let libm = f64::from_bits(0x3d75_a409_90d7_962a);
    assert!(within_rounding(libm, glibc));
    // A probability near 1 that has lost its distance from 1 (Mi4) must be refused.
    let near_one = f64::from_bits(0x3fef_ffff_ffff_f917);
    assert!(!within_rounding(1.0, near_one));
}
```
The second test's last assertion fails on the relative form and passes on M1's units form, which is
the point of it.

## 9. What's good

- `production_parity.rs:545-550` keeps a *measured* pin beside the derived wall instead of relaxing
  the pin to the wall, so the test still reports how far the trees actually are apart.
- The flags, the convergence flag and the counts stayed exact while only the probabilities were
  loosened (`production_parity.rs:950-954`, `:1164-1168`, `:1293-1297`) — the split between numbers
  and decisions is the right one (Mi3 only asks that the cut join the decisions).
- Renaming the three tests from `_bit_for_bit` to `_to_within_rounding`, and updating the
  cross-reference in `calibrate.rs:16`, keeps names from asserting what the body no longer checks.

## 10. Commands to re-verify

- `uv run --no-project python tmp/review_2026-09-14_portable_float_B4/calibration_port.py src/paralog/testdata/production_parity_answers.tsv`
- `uv run --no-project --with mpmath python tmp/review_2026-09-14_portable_float_B4/series_sweep.py`
- After the fixes, in the container: `cargo test --lib --all-features paralog::` (the two new tests
  in §8), then the gates in plan §7.

### Author response convention

Address each finding by its identifier (M1, Mi1, …) with `fixed in <commit>` / `disputed because …` /
`deferred to <issue>` / `won't fix because …`. Answer Q1 first.

## Fixes applied (2026-09-14)

| finding | what was done |
|---|---|
| M1 infinity accepted against a finite value | `within_rounding` now counts representable numbers: NaN only against NaN, an infinity only against the same infinity, otherwise the same sign and at most `PROBABILITIES_MAY_DIFFER_BY_UNITS` = 64 apart. New test `within_rounding_refuses_infinities_sign_changes_and_distant_values` covers the cases the review listed, plus 64/65 steps and π next to 1 |
| Mi1 figures | Measured by the implementer, not taken from the review's replication: a temporary print in `within_rounding` counted, over these tests, 17 probabilities of 465 whose bits differ from production's, identically on macOS and Linux — 10 by one step, 5 by two, 1 by eight, π by twenty-two. The doc states those and the relative 3.6e-15 |
| Mi2 stale docs | `locus_score.rs` "exactly zero" passages, `production_parity.rs` "by bit pattern" passages, and both "the one departure" passages updated |
| Mi3 cut needs no tolerance | All three cut comparisons are bit-for-bit again; they pass |
| Mi4 π near 1 | Covered by the step count (M1) and asserted in the new test |
| nits | `prior.rs` rewrapped; `mod.rs` names only `locus_score` as converted; `locus_score.rs` points at the A2 and A3 reports. The 1e-10 pin stays |

Validation after the fixes: container fmt, clippy `-D warnings`, doc clean; container suite 4,795 passed, 0 failed, 4 ignored before the new comparison test, `paralog` 239 passed on macOS and Linux after it. Identity oracle unchanged and the parity oracle's calls (`1e40bd1e`) reproduced, both on the conversion before these test-only fixes.
