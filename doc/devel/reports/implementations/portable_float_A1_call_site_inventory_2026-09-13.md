# Where the caller calls the platform's maths library

*Step A1 of `doc/devel/implementation_plans/portable_float.md`. Branch `portable-float`, tree
identical under `src/` to commit `75722336` (`git diff --stat 75722336 HEAD -- src Cargo.toml` is
empty). Nothing under `src/` was changed. 2026-09-13.*

Rust's `f64::ln`, `exp`, `powf`, `log10`, `ln_1p` and the trigonometric functions hand their
work to the operating system's maths library, and macOS and Linux round those results
differently in the last binary place. The plan is to route every such call through the
pure-Rust `libm` crate so both platforms produce the same bits, but `libm` may be slower. This
report lists every call in `src/`, separates the calls the tests make from the calls the caller
ships with, and says for each shipped call how often it runs and what range of numbers it is
given. It is the input for step A2, which times `libm` against the standard library over those
ranges.

**Words used below.**

- **Transcendental function** — one whose result a processor cannot produce with a single
  exact instruction (`ln`, `exp`, `powf`, `sin`, …), so a maths library computes it by series or
  tables and chooses how to round. `+ − × ÷` and `sqrt` are not in this group: IEEE 754 requires
  them to be correctly rounded, and every platform gives the same bits.
- **Shipped call** — a call outside test code (the test rule is in §1.2).
- **Reachable** — a shipped call that some path from the `pop_var_caller` binary can reach. 31
  shipped calls are compiled into the library but called only from tests or from `examples/`;
  they are listed and marked *not reachable*.
- **Observation** — at a locus, the reads of one sample that showed the same allele in the same
  read group, counted together. Likelihood rows loop over observations, not over single reads.
- **Pass** — one round of the calling loop that re-estimates a locus's allele frequencies and
  rescores every sample (`run_frequency_loop`, `src/calling/inference/summarise_condition.rs:793`).
- **The parameter fit** — `fit_a_cohort` (`src/run/census_fit.rs:150`), run by
  `estimate-parameters`. It holds three optimisers: the SNP/indel joint fit
  (`parameter_estimation/joint/fit.rs`), the repeat-tract fit (`joint/ssr_fit.rs`) and the
  contamination fit (`joint/contamination.rs`).
- **Log-sum-exp** — computing `ln(eᵃ + eᵇ + …)` by subtracting the largest term first, so every
  `exp` sees a number at or below zero.

## Summary

**202 shipped calls and 232 test calls**, found by the commands in §1. Every shipped call is
`f64`; none is `f32`. Of the shipped calls, 170 are reachable from the binary, 31 are not, and 1
runs only when the environment variable `PVC_MINTED_ERROR_CENSUS=1` is set.

**Shipped calls by function and by how often they run** (all 202; the number reachable from the
binary in brackets where it differs):

| function | per-base | per-read | per-genotype | per-fit-iteration | per-locus | per-table | per-run | unclear | total |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `ln` | – | 30 (9) | 16 (13) | 40 | 8 | 7 | 24 | 1 (0) | 126 (101) |
| `exp` | 1 (0) | 3 (2) | 8 | 23 (22) | 2 | – | 6 | 1 (0) | 44 (40) |
| `sqrt` | – | – | – | 3 | – | – | 11 | – | 14 |
| `powi` | – | 2 (1) | – | 4 | 1 | – | 2 (1) | – | 9 (7) |
| `powf` | – | – | – | – | – | – | 3 | – | 3 |
| `ln_1p` | 1 (0) | – | 1 | – | – | – | – | – | 2 (1) |
| `log10` | – | – | – | – | 2 | – | – | – | 2 |
| `sin` | – | – | – | 1 | – | – | – | – | 1 |
| `libm::lgamma` | – | – | – | – | 1 | – | – | – | 1 |
| **total** | **2 (0)** | **35 (12)** | **25 (22)** | **71 (70)** | **14** | **7** | **46 (45)** | **2 (0)** | **202 (170)** |

The per-read column's "(12)" counts the one census call gated by the environment variable as not
reachable. The class names are defined at the head of §2. One per-table call,
`calling/genotype_table.rs:535`, also runs per locus at repeat tracts with more than 16 candidate
alleles (§2.3); it is counted once, as per-table.

**No reachable call runs per base.** The two per-base calls sit in `ssr_marginal_whole_read.rs`,
an aligner nothing in the binary constructs. The per-base arithmetic that does run — the
repeat delimiter's match and mismatch scores and each read's error logarithm — reads tables
built from written-out bits or from `LN_10` multiplication (`static PER_QUALITY_LN`, `src/alignment/emission.rs:170`,
`src/locus_generation/pileup/open_record.rs:3233`), so it calls no maths library at all.

**Where the hot reachable calls are.**

- **Per read (12 calls).** The repeat-tract path: the delimiter's slip costs, 4 `ln` per
  delimited read (`alignment/ssr_unit_robust.rs:278-288`); the STR emission, 2 `ln` and 1 `exp`
  per observation, candidate allele and slip placement, and for a read that ran off the tract
  also per reachable length change (`alignment/emission.rs:558-559`,
  `calling/likelihood/ssr_emission.rs:629`) and a `powi` for the stutter probability
  (`alignment/stutter.rs:783`). The SNP/indel likelihood row: 3 `ln` per observation per copy
  count (`calling/likelihood/generic.rs:810, 903, 912`) and 1 `exp` per observation
  (`calling/likelihood/mod.rs:510`).
- **Per genotype (22 calls).** The genotype prior, evaluated per sample per genotype per pass
  (`calling/genotype_prior/dirichlet_multinomial.rs`, 10 calls); turning scores into posteriors
  (`calling/inference/summarise_condition.rs:364`); the two likelihood rows' per-genotype
  logarithms (`generic.rs:830`, `calling/likelihood/ssr.rs:736`); the site-quality fold over
  samples (`calling/quality/mod.rs:562-612`, 4 calls); and the paralog scorer's inner loops over
  frequency grid points and samples (`paralog/locus_score.rs:490-573`, 5 calls).
- **Per fit iteration (70 calls).** Nearly all in the parameter fit. The loops over samples and
  over quadrature points only multiply and add (`joint/fit.rs:2115-2127`,
  `joint/ssr_fit.rs:2031-2069`); the maths-library calls sit on either side of them, so each
  call's count is its own product, not the product of every loop:
  - `exp` per census position × 2 classes × candidate (≤ 3) × **sample** × 3 genotypes, per pass
    (`joint/fit.rs:2098`);
  - `ln` per census position × 2 classes × candidate × **16 quadrature nodes**, per pass
    (`joint/fit.rs:2132`);
  - `ln` per repeat tract × **sample × 91 genotypes** × read group × non-empty length bucket (≤ 17),
    each time the slippage numbers change (`joint/ssr_fit.rs:1948`);
  - `ln` and `exp` per repeat tract × **256 quadrature points**, per score evaluation
    (`joint/ssr_fit.rs:2073, 2332`);
  - `exp` and 2 `ln` per bisection step (≤ 60) × 12 sticks × 256 points, each time the length
    spectrum or concentration changes (`joint/ssr_fit.rs:2388`).

  The paralog filter's prior fit adds 3 (`paralog/prior.rs`).

**Many hot calls take an argument that does not change inside their loop.** Converting such a
call to a slower `libm` function costs speed only because it is recomputed; the value could be
computed once outside the loop instead. This bears on how step A2 should weight the per-function
slowdowns. The cases, from the code:

- `alignment/emission.rs:558-559` compute `ln(1 − ε)` and `ln(ε/3)` from the substitution rate,
  which is fixed per (read group, candidate) at a locus; they rerun because
  `substitution_probability` builds a new `FlatEmission` on every call
  (`calling/likelihood/ssr_emission.rs:617`).
- `calling/likelihood/generic.rs:810, 912` take `(1 − c)·k/P + c·q`. Where a read group has no
  contamination fraction (*c* = 0), that is `k/P`: one of at most 16 values for the whole run.
- `alignment/stutter.rs:783` raises a share fixed per (read group, candidate) to one of at most 10
  exponents.
- `alignment/ssr_unit_robust.rs:278-288` take the same four literal values on every read (listed
  as constants below).
- `joint/ssr_fit.rs:1948` takes a probability that depends only on (read group, genotype pair,
  bucket) (`ssr_fit.rs:1946-1947`), so one refresh needs at most read groups × 91 × 17 distinct
  logarithms rather than one per tract × sample.
- `paralog/locus_score.rs:490` takes σ = a sample's single-copy depth SD × 1 or √(T/2), every factor
  fixed once per run; per sample that is at most 5 distinct values.
- The multiplicity logarithms and `ln 3` in `joint/fit.rs` (table below).

**Sites whose argument is computable before the run starts** — they could be written out as
bits, as commit `75722336` did for the emission table, at no speed cost:

| site | what is computed | why it is constant |
|---|---|---|
| `alignment/ssr_best_path_flat_gap.rs:201` | `e⁻¹` | `(-1.0f64).exp()` |
| `alignment/ssr_best_path_flat_gap.rs:221-225` (5) | `ln` of the gap-open, gap-extend and match-to-match probabilities | built from the literals at `:183`, `:196` and `:201` |
| `alignment/ssr_unit_robust.rs:278-288` (4) | `ln` of the stutter shares that price a whole-unit slip | the shipped generator always hands the delimiter `StutterModel::hipstr_shipped()` (`locus_generation/ssr.rs:2331`), whose shares are literals (`alignment/stutter.rs:308-317`) |
| `calling/genotype_table.rs:535` | `ln 2 + … + ln n` | `n` is at most the ploidy: a small integer, so a table of `ln k` covers it |
| `paralog/locus_score.rs:260, 261, 265, 266` (8) | `ln ε`, `ln(1−ε)`, `ln 0.5` | `ε` is `pseudocount_vaf`, and the run passes `ParalogModelParams::default()` (`run/paralog_filter/finish.rs:211`; `DEFAULT_PSEUDOCOUNT_VAF = 0.01`, `paralog/model_params.rs:14`) |
| `paralog/locus_score.rs:454` | `ln(cells)` | cells = 40 grid points × 7 carrier configurations under the same defaults |
| `paralog/locus_score.rs:466, 479, 480` | `√(T/2)`, `ln(m/T)`, `ln(1 − m/T)` | `T` from `DEFAULT_CARRIER_COPY_NUMBERS = [3, 4, 6, 8]` (`paralog/model_params.rs:27`) |
| `parameter_estimation/depth_bins.rs:233, 241, 299` | the depth ladder's widening ratio `(124/8)^(1/11)` and its powers (`:241` is not reachable, §2.7) | literals at `depth_bins.rs:56, 77, 85, 97` |
| `parameter_estimation/joint/fit.rs:2183` | `ln 3` | `CANDIDATE_ALTERNATIVES = 3` (`fit.rs:941`), recomputed per position |
| `parameter_estimation/joint/fit.rs:2190, 2196, 2218, 2327, 2364, 2442` | `ln(multiplicity)` | multiplicity is 1, 2 or 3 (`fit.rs:2047, 2055`), so three values cover every call |
| `parameter_estimation/joint/share_curve.rs:558` | logit of the fallback share | the fallback is `DEFAULT_SHORTER_SHARE` or `DEFAULT_FALL_OFF` (`joint/ssr_fit.rs:1799, 1805`) |
| `parameter_estimation/joint/ssr_fit.rs:2305` | `−ln(points)` | `QUADRATURE_POINTS = 256` (`ssr_fit.rs:529`) under the default configuration (`ssr_fit.rs:653`) |
| `parameter_estimation/joint/ssr_fit.rs:2357`, one of three | `½ ln 2π` | `TAU.ln()`, recomputed on every `ln_gamma` call |

Two of these are in hot loops: the four slip-cost logarithms run once or twice per delimited read,
and `fit.rs:2183` plus the six multiplicity logarithms run per census position per pass.

**Other checks.**

- **`wide` SIMD types:** used at one place, `joint/ssr_fit.rs:2042-2046`, with `ZERO`, `new`,
  `*`, `+=` and `to_array` only. In `wide` 1.4.0 those map to the processor's plain add and
  multiply (`wide-1.4.0/src/f64x4_.rs:54-69, 88-103`; `+=` forwards to `+`, `wide-1.4.0/src/lib.rs:369`), which are
  exact.
- **`mul_add`:** no call anywhere in `src/`. Rust does not fuse a separate `*` and `+` into one
  instruction on its own, including under the `target-cpu=x86-64-v3` flag in
  `.cargo/config.toml:16`.
- **`f32` transcendental calls:** none. The only `f32` on a call line converts an `exp` result
  after the fact (`joint/fit.rs:2278`).
- **`libm::lgamma`** (`genetics.rs:54`) is already the portable library. 15 shipped calls reach it
  through the `genetics::lgamma` wrapper (`calling/quality/artifact_correction.rs:308, 327`;
  `calling/quality/mod.rs:460-473`), and none of them changes under the plan.
- **`sqrt` (14 shipped calls)** is exact on every platform, as defined above. It is listed so the
  inventory is complete, and needs no conversion.

## 1. How the calls were found

### 1.1 The search

One regular expression matches both call forms — the method form `x.ln()` and the path form
`f64::ln(x)` or `.map(f64::exp)` — for every function the plan names, plus any `libm::` path:

```sh
F='ln|exp|powf|powi|log10|log2|log|ln_1p|exp_m1|sqrt|cbrt|sin|cos|tan|tanh|sinh|cosh|asin|acos|atan|atan2|hypot|mul_add|exp2'
# every occurrence, one line per call
rg -n -o --no-heading -e "(\.|f(32|64)::)($F)\b(\(|\))|libm::[a-z_0-9]+" src      # 450 occurrences
# drop occurrences on comment lines (16: doc comments that name f64::ln and libm::lgamma)
rg -n --no-heading -e "(\.|f(32|64)::)($F)\b(\(|\))|libm::" src | rg -v '^[^:]+:[0-9]+:\s*//' \
  | rg -o -e "(\.|f(32|64)::)($F)\b(\(|\))|libm::[a-z_0-9]+" | wc -l             # 434 remain
# (without the last `rg -o` stage the same pipeline counts 381 matching lines)
```

A second search without the trailing parenthesis, `rg -n -e "f(32|64)::($F)\b" src`, found one
path-form call in code (`parameter_estimation/calibration.rs:238`, `.map(f64::exp)`), which the
first search already counts. No project type defines a method named `ln`, `exp`, `powf`, `powi`,
`log10`, `log2`, `ln_1p`, `exp_m1`, `sqrt`, `sin` or `cos`
(`rg -n -e 'fn (ln|exp|powf|powi|log10|log2|ln_1p|exp_m1|sqrt|sin|cos)\b' src` returns nothing), so
every match is a standard-library float call. None of `log2`, `log` with a base, `exp_m1`, `cbrt`,
`tan`, the hyperbolic or inverse trigonometric functions, `hypot`, `mul_add` or `exp2` occurs.

The working files are in `tmp/` of this worktree: `portable_float_A1_occ.txt` (the 434
occurrences), `portable_float_A1_classified.tsv` (each marked shipped or test),
`portable_float_A1_sites.tsv` (each shipped site with its frequency class, reachability and
whether its argument is constant).

### 1.2 Telling test code from shipped code

In each of the 44 files with a match I listed every `#[cfg(test)]` and
`#[cfg(any(test, feature = "bench-fixtures"))]` attribute
(`rg -n -A2 '#\[cfg\((test|any\(test|feature = "bench-fixtures")' <file>`), found the line where
test-only code starts, and then checked that no shipped item follows the test module's closing
brace. Three files have a test-only item *before* shipped code, so their boundary is the later
`mod tests`: `calling/likelihood/ssr_emission.rs` (test helpers at 656 and 735, shipped code to
821, tests from 823), `parameter_estimation/joint/census_moments.rs` (a test helper at 537,
tests from 899) and `paralog/mod.rs` (`mod production_parity` at 56, tests from 84). Calls in the
`bench_fixtures` modules (`joint/fit.rs:3331-3683`, `joint/ssr_fit.rs:2692-2841`; 18 calls, a
random-number generator for synthetic cohorts) count as test code, as the brief asks. No file
whose name contains `parity` or `test_fixtures` holds a match. `paralog/calibration.rs` has no test
module; both its calls are shipped.

### 1.3 Totals

| function | shipped | test | all |
|---|---:|---:|---:|
| `ln` | 126 | 126 | 252 |
| `exp` | 44 | 51 | 95 |
| `powi` | 9 | 27 | 36 |
| `sqrt` | 14 | 10 | 24 |
| `powf` | 3 | 11 | 14 |
| `log10` | 2 | 3 | 5 |
| `cos` | 0 | 3 | 3 |
| `ln_1p` | 2 | 1 | 3 |
| `sin` | 1 | 0 | 1 |
| `libm::lgamma` | 1 | 0 | 1 |
| **total** | **202** | **232** | **434** |

Shipped calls per file. *Reachable* is the part a path from the binary reaches (§2 says why each
unreachable call is unreachable):

| file | shipped | reachable | functions |
|---|---:|---:|---|
| `parameter_estimation/joint/fit.rs` | 38 | 37 | ln 23, exp 12, sqrt 3 |
| `parameter_estimation/joint/ssr_fit.rs` | 24 | 24 | ln 13, exp 8, powi 2, sin 1 |
| `paralog/locus_score.rs` | 21 | 21 | ln 16, exp 3, ln_1p 1, sqrt 1 |
| `parameter_estimation/joint/slippage_curve.rs` | 11 | 11 | ln 5, exp 2, powf 2, sqrt 2 |
| `calling/genotype_prior/dirichlet_multinomial.rs` | 10 | 10 | ln 8, exp 2 |
| `parameter_estimation/joint/contamination.rs` | 10 | 10 | sqrt 4, ln 3, powi 2, exp 1 |
| `parameter_estimation/joint/share_curve.rs` | 9 | 9 | ln 4, sqrt 3, exp 2 |
| `alignment/ssr_best_path_flat_gap.rs` | 6 | 6 | ln 5, exp 1 |
| `calling/quality/artifact_correction.rs` | 6 | 6 | ln 4, exp 1, log10 1 |
| `calling/quality/mod.rs` | 6 | 6 | ln 3, exp 2, log10 1 |
| `alignment/ssr_unit_robust.rs` | 4 | 4 | ln 4 |
| `calling/likelihood/generic.rs` | 4 | 4 | ln 4 |
| `genetics.rs` | 4 | 4 | ln 3, libm::lgamma 1 |
| `alignment/ssr_anchor_firm.rs` | 4 | 0 | ln 4 |
| `alignment/ssr_anchor_robust.rs` | 4 | 0 | ln 4 |
| `alignment/ssr_best_path_unit_slip.rs` | 4 | 0 | ln 4 |
| `alignment/ssr_noise_robust.rs` | 4 | 0 | ln 4 |
| `alignment/ssr_robust_indel.rs` | 4 | 0 | ln 4 |
| `calling/genotype_prior/hardy_weinberg.rs` | 3 | 0 | ln 3 |
| `paralog/prior.rs` | 3 | 3 | exp 2, ln 1 |
| `parameter_estimation/depth_bins.rs` | 3 | 2 | powi 2, powf 1 |
| `alignment/emission.rs` | 2 | 2 | ln 2 |
| `alignment/stutter.rs` | 2 | 2 | powi 2 |
| `paralog/calibration.rs` | 2 | 2 | ln 1, exp 1 |
| `calling/likelihood/mod.rs` | 2 | 1 | exp 1, ln 1 |
| `alignment/ssr_marginal_sequence.rs` | 2 | 0 | ln 1, powi 1 |
| `alignment/ssr_marginal_whole_read.rs` | 2 | 0 | exp 1, ln_1p 1 |
| `locus_generation/pileup/minted_error_census.rs` | 2 | 0 | exp 2 (one gated by an environment variable) |
| `calling/genotype_table.rs` | 1 | 1 | ln 1 |
| `calling/inference/summarise_condition.rs` | 1 | 1 | exp 1 |
| `calling/likelihood/ssr.rs` | 1 | 1 | ln 1 |
| `calling/likelihood/ssr_emission.rs` | 1 | 1 | exp 1 |
| `parameter_estimation/calibration.rs` | 1 | 1 | exp 1 |
| `parameter_estimation/joint/census_moments.rs` | 1 | 1 | sqrt 1 |
| **34 files** | **202** | **170** | |

For scale only: the same search over `examples/` finds 454 occurrences, comment lines included
(`rg -o --no-heading -e "(\.|f(32|64)::)($F)\b(\(|\))" examples | wc -l`). They are outside this
inventory.

## 2. Every shipped call site

One row per source line and function; *×2* marks two calls of the same function on one line.
Paths are under `src/`.

**Frequency classes.**

- **per-base** — inside an alignment matrix cell or once per base of a read.
- **per-read** — once or a few times per read or per observation, often per candidate allele.
- **per-genotype** — inside a locus's loops over samples, genotypes or allele counts, usually
  once per pass.
- **per-fit-iteration** — inside the parameter fit's (or the paralog prior's) optimisation loops:
  once per run overall, but each loop iterates, and several run per census position.
- **per-locus** — a fixed number of times per locus, or per sample per locus.
- **per-table** — while building a lookup table that later loops read.
- **per-run** — during setup, or a small fixed number of times per run.
- **unclear** — no reachable caller was found.

The **Constant** column says whether the argument can be computed before the run: *literal* (from
literals in the code), *default* (from default parameters the run always passes), *small integer*
(one of a handful of integers), or *no*.

### 2.1 Alignment

**How the reachable delimiter is reached.** The calling run builds one repeat-tract generator per
walker with `SsrGenerator::with_default_aligner` (`run/walker.rs:1638`), which constructs
`SsrUnitRobustAligner::new(PerQualityEmission::new())` (`locus_generation/ssr.rs:2297`). The
generator classifies each kept read (`locus_generation/ssr.rs:2519-2526`), and classifying calls
the delimiter once, or twice when it retries on a wider window (`locus_generation/ssr.rs:932, 942`).

| site | fn | enclosing fn | the value | how often | constant |
|---|---|---|---|---|---|
| `alignment/emission.rs:558` | ln | `FlatEmission::from_unchecked_rate` | log probability that a base is read correctly, `1 − ε`, floored at the smallest positive double (`emission.rs:139`) | **per-read.** `FlatEmission::try_new` (`emission.rs:539`) calls it; `substitution_probability` (`calling/likelihood/ssr_emission.rs:617`) builds a fresh `FlatEmission` on every call, and that runs per observation per candidate allele per slip placement (`ssr_emission.rs:418, 507, 533, 544`), and at an interrupted tract per reachable length change for a read that ran off it (`:437-448`), inside `fill_emissions` (`calling/likelihood/ssr.rs:762-778`) | no — `ε` is the fitted per-stratum substitution rate |
| `alignment/emission.rs:559` | ln | same | log probability of one particular wrong base, `ε / 3` | **per-read**, same chain | no |
| `alignment/ssr_unit_robust.rs:278` | ln | `SlipCosts::from_model` | log share of reads with no length change | **per-read.** `delimit` calls it on every read (`ssr_unit_robust.rs:429`), reached through `align` (`:1023-1030`) from the generator's per-read classification above | **default** — `StutterModel::hipstr_shipped()`, set in `SsrGenerator::new` (`locus_generation/ssr.rs:2331`); same-length share 0.88 |
| `alignment/ssr_unit_robust.rs:282` | ln | same | log cost of opening a run of gained repeat units, `ln(longer share × one-step share)` | **per-read**, same | **default** — `ln(0.05 × 0.95)` (`stutter.rs:310, 312`) |
| `alignment/ssr_unit_robust.rs:286` | ln | same | log cost of opening a run of lost units | **per-read**, same | **default** — `ln(0.05 × 0.95)` |
| `alignment/ssr_unit_robust.rs:288` | ln | same | log cost of each further unit, `ln(1 − one-step share)` | **per-read**, same | **default** — `ln(0.05)` |
| `alignment/ssr_best_path_flat_gap.rs:201` | exp | `static GAP_EXTEND_PROB` (a `LazyLock`) | the gap-extension probability, `e⁻¹` | **per-run** — evaluated once on first use | **literal** |
| `alignment/ssr_best_path_flat_gap.rs:221` | ln | `TransitionCosts::default` | log match-to-match transition, `ln(1 − 2 × 2.9e-5)` | **per-run.** `TransitionCosts::new` runs in each aligner's constructor (`ssr_unit_robust.rs:381`), once per generator | **literal** (`:183`) |
| `alignment/ssr_best_path_flat_gap.rs:222` | ln | same | log gap-open in the flank, `ln 2.9e-5` | **per-run**, same | **literal** |
| `alignment/ssr_best_path_flat_gap.rs:223` | ln | same | log gap-open in the tract, `ln 1e-2` | **per-run**, same | **literal** (`:196`) |
| `alignment/ssr_best_path_flat_gap.rs:224` | ln | same | log gap-close, `ln(1 − e⁻¹)` | **per-run**, same | **literal** |
| `alignment/ssr_best_path_flat_gap.rs:225` | ln | same | log gap-extend, `ln e⁻¹` | **per-run**, same | **literal** |
| `alignment/stutter.rs:486` | powi | `StutterModel::unreachable_mass` | chance a geometric slip stays within `steps` repeats, `1 − (1 − one-step share)^steps` | **per-locus** — per candidate allele per read group, from `SsrScoringContext::new` (`calling/likelihood/ssr_emission.rs:138`) | no |
| `alignment/stutter.rs:783` | powi | `Regime::probability` | probability of a slip of a given size, `share × one-step × (1 − one-step)^(size−1)` | **per-read** — `StutterModel::probability` (`stutter.rs:369-397`) calls it from `emission` per observation per candidate (`ssr_emission.rs:355`), and for a read that ran off the tract per reachable length change, through `reachable_length_changes` (`stutter.rs:594`) from `censored_emission` (`ssr_emission.rs:437-439`) and `probability_at_least_this_much_longer` (`stutter.rs:687`, `ssr_emission.rs:400-401`) | no — the shares are fitted |

**Compiled but not reachable from the binary.** Four sibling delimiters are referenced only from
`examples/` (`rg -l` for their type names finds `ng_ssr_anchor_firm_validate.rs`,
`ng_ssr_synthetic_bakeoff.rs` and, for unit-slip, the bake-off examples and
`locus_generation/ssr.rs`'s tests from line 2715 on). The two marginal aligners are used only by
tests.

| site | fn | enclosing fn | the value | how often, if called | constant |
|---|---|---|---|---|---|
| `alignment/ssr_anchor_firm.rs:171, 175, 179, 181` | ln ×4 lines | `SlipCosts::from_model` | the four slip costs, as in `ssr_unit_robust.rs` | per-read (`delimit`, `:301`) | depends on the model the caller passes |
| `alignment/ssr_anchor_robust.rs:212, 216, 220, 222` | ln ×4 lines | same | same | per-read (`:343`) | same |
| `alignment/ssr_best_path_unit_slip.rs:141, 145, 149, 151` | ln ×4 lines | same | same | per-read (`:269`) | same |
| `alignment/ssr_noise_robust.rs:209, 213, 218, 221` | ln ×4 lines | `SlipCosts::from_model` (with a margin) | same | per-read (`:343`) | same |
| `alignment/ssr_robust_indel.rs:196, 200, 205, 208` | ln ×4 lines | same | same | per-read (`:330`) | same |
| `alignment/ssr_marginal_sequence.rs:145` | powi | `SsrSequenceMarginal::equal_length_probability` | probability a read identical to the reference was read without error, `(1 − ε)^length` | per-read | no |
| `alignment/ssr_marginal_sequence.rs:294` | ln | `marginal_probability` | log of the summed alignment probability | per-read | no |
| `alignment/ssr_marginal_whole_read.rs:239` | exp, ln_1p | `ln_sum2` | log-sum-exp of two path scores | **per-base** — matrix cells (`:192, 201`) | no |

### 2.2 Calling: likelihoods

| site | fn | enclosing fn | the value | how often | constant |
|---|---|---|---|---|---|
| `calling/likelihood/generic.rs:810` | ln | `assemble_genotype_log_likelihood_row` | log probability that a read showing this allele came from a genotype carrying *k* copies of it, `ln((1 − c)·k/P + c·q)` (*c* contamination fraction, *P* ploidy, *q* contaminant's allele frequency) | **per-read** — per observation, per copy count up to the ploidy. The row is built once per sample per locus (`calling/inference/summarise_condition.rs:2481`) and again every pass at a contaminated locus (`:2480`, `:912-918`) | no |
| `calling/likelihood/generic.rs:830` | ln | same | log probability of that read under a genotype carrying none of the allele, `ln((1 − c)·ε̄/spread + c·q)` | **per-genotype** — per observation per genotype with zero copies | no |
| `calling/likelihood/generic.rs:903` | ln | `score_partials` | log probability of a partial read under a genotype carrying no compatible allele | **per-read** — per partial observation, same chain | no |
| `calling/likelihood/generic.rs:912` | ln | same | same, for 1 to *P* compatible copies | **per-read** — per partial observation per copy count | no |
| `calling/likelihood/mod.rs:510` | exp | `ReadGroupCalibration::charged_error` | the error probability a group of reads is charged: the calibration scale times the geometric mean of their minted error probabilities | **per-read** — per observation and per partial, from `fill_generic_emissions` (`generic.rs:626, 644`), once per sample per locus | no |
| `calling/likelihood/ssr.rs:736` | ln | `genotype_log_likelihood_row` (repeat tracts) | log probability of one observed tract length under a genotype, mixing this individual, junk and contamination | **per-genotype** — per observation per genotype (`ssr.rs:719-736`); the row runs once per sample per tract (`summarise_condition.rs:1859-1862`) | no |
| `calling/likelihood/ssr_emission.rs:629` | exp | `substitution_probability` | probability that a read's bases differ from a candidate's in the observed way, from the summed per-base log scores | **per-read** — per observation per candidate per placement (chain in §2.1), and for a read that ran off an interrupted tract also per reachable length change (`letters_over` inside the loop at `ssr_emission.rs:437-448`) | no |
| `calling/likelihood/mod.rs:438` | ln | `ReadGroupCalibration::log_scale` | log of the calibration scale | **unclear** — every caller is a test (`generic.rs:2615`, `parameters_file/defaults.rs:563`); not reachable | no |

### 2.3 Calling: genotype priors and posteriors

The shipped prior is `MarginalizedDirichletPrior` (`cli/call_from_alignments.rs:605`,
`cli/call_from_psps.rs:547`). `score_one_sample` fills it for every sample on every pass
(`summarise_condition.rs:337`, called from the pass loop at `:839-843`).

| site | fn | enclosing fn | the value | how often | constant |
|---|---|---|---|---|---|
| `calling/genotype_prior/dirichlet_multinomial.rs:78` | ln | `fill_random_mating_log_priors` | log of each allele's concentration (its pseudo-count in the prior) | **per-genotype** — per allele per sample per pass | no |
| `…/dirichlet_multinomial.rs:136` | ln | `one_genotypes_log_prior` | log of the rising product `α(α+1)…(α+k−1)` for an allele a genotype carries 2+ copies of | **per-genotype** — per such genotype per sample per pass | no |
| `…/dirichlet_multinomial.rs:282` | ln | `fill_inbreeding_mixture_log_priors` | log of the total concentration | **per-genotype** — per sample per pass | no |
| `…/dirichlet_multinomial.rs:292` | ln | same | log rising product of the total concentration over the ploidy | **per-genotype** — per sample per pass | no |
| `…/dirichlet_multinomial.rs:307` | ln | same | log of the inbreeding coefficient *F* | **per-genotype** — per sample per pass | no |
| `…/dirichlet_multinomial.rs:308` | ln | same | log of `1 − F`, floored at 1e-300 | **per-genotype** — per sample per pass | no |
| `…/dirichlet_multinomial.rs:322` | ln | same | log concentration of the allele a homozygote carries | **per-genotype** — per homozygous genotype per sample per pass | no |
| `…/dirichlet_multinomial.rs:372` | exp ×2, ln | `log_sum_exp_2` | log-sum-exp of the random-mating and identical-by-descent branches | **per-genotype** — per homozygote per sample per pass (`:325`); also per sample per genotype in site quality (`calling/quality/mod.rs:418`) | no |
| `calling/inference/summarise_condition.rs:364` | exp | `score_one_sample` | a genotype's unnormalised posterior weight, `exp(score − best score)` | **per-genotype** — per genotype per sample per pass | no |
| `calling/genotype_table.rs:535` | ln | `log_factorial` | `ln n!`, summed term by term | **per-table** — the multinomial coefficients of a genotype table (`:543-545`), built once per ploidy and allele count and cached per thread (`:173-211`). **Per locus × genotypes at repeat tracts with more than 16 candidate alleles**: the cache holds only tables up to 16 alleles and ploidy 8 (`:54-57`, `:175-176`), every locus asks for its table (`calling/inference/summarise_condition.rs:2315`), and the repeat-tract candidate cap ships at 32 (`calling/allele_candidates/ssr.rs:477-478`). Such a tract rebuilds its table, up to 528 diploid genotypes × at most 2 `ln` calls | **small integer** — `n` ≤ ploidy; at diploid always `ln 2` |

Not reachable: `PlugInWrightPrior` is referenced only by `genotype_prior/mod.rs`'s tests, so
`calling/genotype_prior/hardy_weinberg.rs:145, 146, 150` (`ln` of *F*, of `1 − F`, and of each
allele's frequency, in `fill_plug_in_mixture_log_priors`) would be per-genotype if called.

### 2.4 Calling: quality

| site | fn | enclosing fn | the value | how often | constant |
|---|---|---|---|---|---|
| `calling/quality/mod.rs:191` | log10 | `score_best_genotype` | the genotype quality, `−10 log₁₀(1 − best posterior)` | **per-locus** — once per sample per locus (`summarise_condition.rs:1294`) | no |
| `calling/quality/mod.rs:562` | exp | `fold_samples_into_allele_counts` | each copy count's likelihood relative to the sample's best | **per-genotype** — per sample per copy count (0 to ploidy), once per locus (`quality/mod.rs:436`, via `score_uncorrected_site_quality`, `summarise_condition.rs:1343`) | no |
| `calling/quality/mod.rs:588` | ln | same | log of the running rescaling factor | **per-genotype** — per sample, once per locus | no |
| `calling/quality/mod.rs:597` | ln | same | log probability of each cohort allele count | **per-genotype** — per allele count (up to samples × ploidy + 1), once per locus | no |
| `calling/quality/mod.rs:612` | exp | `log_sum_exp_over` | each allele count's weight relative to the largest | **per-genotype** — per allele count, once per locus (`quality/mod.rs:488`) | no |
| `calling/quality/mod.rs:613` | ln | same | log of the normaliser | **per-locus** | no |
| `calling/quality/artifact_correction.rs:139` | log10 | `two_sided_binomial_tail_phred` | a bias test's p-value as Phred, `−10 log₁₀ p` | **per-locus** — up to three tests per record (`artifact_correction.rs:445, 509, 514`), from `correct_site_quality` (`run/records.rs:276`) | no |
| `…/artifact_correction.rs:309` | ln ×2 | `log_binomial_probability` | log of the share and of `1 − share` in a binomial probability | **per-locus** — a bisection makes about log₂(depth) calls per test (`:210, 232`) | no |
| `…/artifact_correction.rs:327` | ln ×2 | `regularised_incomplete_beta` | `ln x` and `ln(1 − x)` in the incomplete-beta prefactor | **per-locus** — up to two per test (`:256, 284`) | no |
| `…/artifact_correction.rs:328` | exp | same | the prefactor itself | **per-locus** | no |

### 2.5 Shared helpers and the error census

| site | fn | enclosing fn | the value | how often | constant |
|---|---|---|---|---|---|
| `genetics.rs:54` | libm::lgamma | `lgamma` | log-gamma, already through `libm` | **per-locus** — the site-quality count prior (`quality/mod.rs:460-473`, over every allele count, rebuilt only when its shape changes) and the bias tests (`artifact_correction.rs:308, 327`) | no |
| `genetics.rs:86, 87, 88` | ln (3 lines) | `wright_genotype_log_priors` | log prior of hom-ref, het and hom-alt at frequency *p* and inbreeding *F*, each floored at 1e-300 | **per-table** — the paralog scorer's precompute, 200 frequency points × samples (`paralog/locus_score.rs:228`), built once per run (`run/paralog_filter/scoring_context.rs:209`) | no |
| `locus_generation/pileup/minted_error_census.rs:151` | exp | `record_read` | one read's minted error probability, summed for a diagnostic | **per-read**, and **only when `PVC_MINTED_ERROR_CENSUS=1`** (`fast_column.rs:412-419`, `open_record.rs:913`) | no |
| `locus_generation/pileup/minted_error_census.rs:85` | exp | `MintedErrorTotals::geometric_mean` | geometric mean of the census | **unclear** — only `examples/ng_minted_error_means.rs` and tests call it; not reachable | no |

### 2.6 Paralog filter

The filter runs after calling in `call-from-alignments` (`cli/call_from_alignments.rs:637`). It
builds `ParalogScoringContext` once with default parameters (`run/paralog_filter/finish.rs:208-213`),
scores each record (`scoring_context.rs:354`) and fits the paralog rate from a histogram
(`paralog/calibrate.rs:39`).

| site | fn | enclosing fn | the value | how often | constant |
|---|---|---|---|---|---|
| `paralog/locus_score.rs:226` | ln | `ParalogScorePrecompute::new` | log weight of each frequency grid point, `−ln(p(1−p))` | **per-table** — 200 points (`model_params.rs:32`), once per run | no — the grid depends on cohort size |
| `paralog/locus_score.rs:246, 247` | ln (2 lines) | same | log probability a sample is or is not a carrier, per carrier-frequency point and sample | **per-table** — 40 points × samples, once per run | no — depends on each sample's *F* |
| `paralog/locus_score.rs:260, 261` | ln ×3 each | same | `ln ε`, `ln 0.5`, `ln(1 − ε)` for the real-variant hypothesis | **per-run** | **default** (ε = 0.01) |
| `paralog/locus_score.rs:265, 266` | ln | same | `ln ε`, `ln(1 − ε)` | **per-run** | **default** |
| `paralog/locus_score.rs:454` | ln | `h2_log_likelihood` | log of the number of grid cells averaged over | **per-locus** | **default** (280) |
| `paralog/locus_score.rs:466` | sqrt | `enumerate_carrier_configs` | how much wider the depth spread is at *T* copies, `√(T/2)` | **per-run** | **default** |
| `paralog/locus_score.rs:479, 480` | ln (2 lines) | same | `ln(m/T)` and `ln(1 − m/T)`, the expected alt share at a carrier | **per-run** — 7 configurations | **default** |
| `paralog/locus_score.rs:490` | ln | `ln_normal` | `ln σ` in a normal log density of a sample's relative copy number | **per-genotype** — per usable sample in `h1_log_likelihood` (`:365`), and per sample, once plus once per carrier configuration, in `h2_log_likelihood` (`:419-421`) | no |
| `paralog/locus_score.rs:528` | exp | `log_add_exp` | the smaller term relative to the larger in a two-term log-sum-exp | **per-genotype** — twice per sample per frequency point in `h1` (`log_sum_exp3`, `:390`), once per sample per cell in `h2` (`:447`); both per locus | no |
| `paralog/locus_score.rs:532` | ln_1p | same | `ln(1 + x)` for that ratio when `x ≥ 1.4e-4` | **per-genotype**, same | no |
| `paralog/locus_score.rs:570, 573` | exp (2 lines) | `LogSumExp::push` | a streamed log-sum-exp's running ratio | **per-genotype** — per frequency point (×2) and per cell, per locus | no |
| `paralog/locus_score.rs:581` | ln | `LogSumExp::value` | log of the streamed sum | **per-locus** (3 per locus) | no |
| `paralog/calibration.rs:109` | ln | `ParalogCalibration::posterior` | logit of the fitted paralog rate | **per-locus** — per scored record (`run/paralog_filter/pass_three.rs:158`) | no — fitted once per run, recomputed per record |
| `paralog/calibration.rs:110` | exp | same | the posterior through the logistic function | **per-locus** | no |
| `paralog/prior.rs:301, 303` | exp (2 lines) | `sigmoid` | the logistic function | **per-fit-iteration** — per non-empty histogram bin (up to 2,000) per EM iteration (up to 500) in `ParalogPrior::estimate` (`prior.rs:206-213`; limits at `:45, 56`) | no |
| `paralog/prior.rs:311` | ln | `logit` | `ln(p/(1−p))` | **per-fit-iteration** — once per EM iteration (`:207`) | no |

### 2.7 Parameter estimation: calibration, depth ladder, curves

| site | fn | enclosing fn | the value | how often | constant |
|---|---|---|---|---|---|
| `parameter_estimation/calibration.rs:238` | exp | `MintedReadErrors::mean_error_probability` | mean minted error probability of a read group | **per-run** — per read group, from `ReadGroupCalibration::from_fitted_rate` (`calling/likelihood/mod.rs:390`, called at `calling/run_parameters.rs:219`) | no |
| `parameter_estimation/depth_bins.rs:233` | powf | `widening_ratio` | ratio between successive widening depth bins | **per-run** | **literal** |
| `parameter_estimation/depth_bins.rs:241` | powi | `ladder_tops` | top of each widening bin | **per-run**, and **not reachable**: its only caller is `DepthBinEdges::new` (`:270`), whose only shipped caller is `impl Default for DepthBinEdges` (`:201-204`); every call of `new` or `default` is in the test module (from `:443`) or `examples/`, and the shipped run builds ladders only through `for_census` | **literal** |
| `parameter_estimation/depth_bins.rs:299` | powi | `DepthBinEdges::for_census` | top of each extra census bin | **per-run** (`run/gatherer.rs:311`) | **literal** |
| `parameter_estimation/joint/census_moments.rs:410` | sqrt | `standard_error_of_the_mean` | standard error of the census frequency and heterozygosity | **per-run** (`:369, 374`) | no |
| `joint/share_curve.rs:133` | ln | `FittedShare::logit` | logit of a stratum's fitted share | **per-run** — curve fits over a handful of strata | no |
| `joint/share_curve.rs:147, 178` | sqrt (2 lines) | `logit_standard_error`, `relative_standard_error` | a share's standard errors | **per-run** | no |
| `joint/share_curve.rs:290` | exp | `ShareCurve::share_on_the_curve` | a curve's share through the logistic function | **per-run** — `share_at` per stratum (`joint/ssr_fit.rs:1556-1557`) and leave-one-out scoring (`share_curve.rs:494`) | no |
| `joint/share_curve.rs:495` | ln | `held_out_error_of` | logit of a held-out prediction | **per-run** | no |
| `joint/share_curve.rs:558` | ln | `share_curve_for_a_period` | logit of the fallback share | **per-run** | **literal** |
| `joint/share_curve.rs:666, 669, 683` | ln, sqrt, exp | `blend_share` | blending a stratum's share with the curve's on the logit scale | **per-run** — per stratum (`ssr_fit.rs:1849`) | no |
| `joint/slippage_curve.rs:132` | sqrt | `FittedCell::relative_standard_error` | a cell's relative standard error | **per-run** | no |
| `joint/slippage_curve.rs:207` | exp | `SlippageCurve::level_on_the_line` | slippage level on a multiplying curve | **per-run** — `level_at` per stratum (`ssr_fit.rs:1555`) and leave-one-out (`slippage_curve.rs:491`) | no |
| `joint/slippage_curve.rs:214` | powf | same | slippage level on a power curve, `line^(1/s)` | **per-run** | no |
| `joint/slippage_curve.rs:303, 305` | ln, powf | `fit_line` | the level transformed onto the line's scale | **per-run** | no |
| `joint/slippage_curve.rs:608, 609, 622, 623, 624` | ln ×2, sqrt, ln, ln, exp | `blend_level` | blending a cell's level with the curve's on the log scale | **per-run** — per stratum (`ssr_fit.rs:1684`) | no |

### 2.8 Parameter estimation: the joint SNP/indel fit (`joint/fit.rs`)

**The loop structure.** `fit_jointly` alternates an expectation pass and a maximisation step until
the parameters stop moving. Each `expectation_pass` (`fit.rs:1763`) rebuilds its quadrature
(16 nodes by default, `fit.rs:606`) and its per-read-group logarithm tables, then calls
`one_position` for every census position (`fit.rs:1852`), in parallel chunks. Inside, the work is
per sample × 2 classes (clean or noisy) × up to 3 candidate alleles × 3 genotypes, and per
quadrature node.

| site | fn | enclosing fn | the value | how often | constant |
|---|---|---|---|---|---|
| `fit.rs:1051` | ln | `ReadLogs::of` | log probabilities that a read shows the candidate, the reference or neither, per carried copy count | **per-fit-iteration** — per pass × 2 classes × read groups (`fit.rs:1789`), 7 each | no |
| `fit.rs:1167` | ln | `ln_reference_reads` | log of a sum over the depths a binned sample could have had | **per-fit-iteration** — per position × sample × class × 3, from `reference_terms` (`fit.rs:2075`), when a sample's depth spans more than one value | no |
| `fit.rs:1176` | ln | `count_times_ln` | a count times the log of a probability | **per-fit-iteration** — inside the maximisation step's golden-section searches (`fit.rs:2675-2677, 2702`) | no |
| `fit.rs:1199` | exp, ln | `ln_sum_exp` | log-sum-exp over a slice | **per-fit-iteration** — about ten times per position per pass (`fit.rs:2193-2450`) | no |
| `fit.rs:2098` | exp | `one_position` | each genotype's likelihood relative to the sample's best | **per-fit-iteration** — per position × sample × class × candidate × 3 | no |
| `fit.rs:2132` | ln | same | log of the product over samples at one quadrature node | **per-fit-iteration** — per position × class × candidate × node | no |
| `fit.rs:2175` | ln | same | same, for the duplicated-locus branch | **per-fit-iteration** — same, when the duplicated branch is on | no |
| `fit.rs:2183` | ln | same | `ln 3` | **per-fit-iteration** — per position | **literal** |
| `fit.rs:2190, 2196, 2218` | ln (3 lines) | same | log multiplicity of a candidate allele | **per-fit-iteration** — per position × class × candidate | **small integer** (1, 2 or 3) |
| `fit.rs:2207, 2209, 2212, 2214, 2226, 2235` | ln (6 lines) | same | logs of the class and branch shares (duplicated share, invariant, fixed, segregating, noisy), floored at the smallest positive double | **per-fit-iteration** — per position × class; the arguments change only between passes | no |
| `fit.rs:2278` | exp | same | the noisy-class posterior, stored as `f32` | **per-fit-iteration** — per position, final pass only | no |
| `fit.rs:2283` | exp | same | a class's posterior | **per-fit-iteration** — per position × class | no |
| `fit.rs:2293, 2294, 2295, 2296` | exp (4 lines) | same | each branch's posterior within the class | **per-fit-iteration** — per position × class | no |
| `fit.rs:2327` | ln | same | log multiplicity | **per-fit-iteration** — per position × class × candidate, when the fixed branch has weight | **small integer** |
| `fit.rs:2332` | exp | same | a candidate's share of the fixed branch | **per-fit-iteration** — same | no |
| `fit.rs:2364` | ln | same | log multiplicity | **per-fit-iteration** — per position × class × candidate | **small integer** |
| `fit.rs:2377` | exp | same | a (candidate, node) share of the segregating branch | **per-fit-iteration** — per position × class × candidate × node | no |
| `fit.rs:2442` | ln | same | log multiplicity | **per-fit-iteration**, duplicated branch | **small integer** |
| `fit.rs:2456` | exp | same | a (candidate, node) share of the duplicated branch | **per-fit-iteration** — per position × class × candidate × node | no |
| `fit.rs:2824, 2825, 2828` | ln (3 lines) | `BetaQuadrature::new` | logs of the quadrature nodes, of one minus each node, and of the weights | **per-fit-iteration** — per node, each time a pass rebuilds a quadrature (`fit.rs:1779`) and in the maximisation step (`fit.rs:2638`) | no |
| `fit.rs:2865` | sqrt | `gauss_jacobi` | off-diagonal entries of the Jacobi matrix | **per-fit-iteration** — per node per rebuild | no |
| `fit.rs:2902, 2903` | sqrt (2 lines) | `symmetric_eigen` | the rotation in a Jacobi eigenvalue sweep | **per-fit-iteration** — up to 100 sweeps × 120 pairs per rebuild | no |
| `fit.rs:2950` | ln | `digamma` | the asymptotic series' `ln x` | **per-fit-iteration** — up to 50 Newton steps × 4 calls (`digamma(a)`, `digamma(b)`, `digamma(a + b)` twice) in `fit_beta_shapes` (`fit.rs:2745-2748`) | no |

Not reachable: `fit.rs:1850` (`exp` of each sample's coverage log-odds at a position) runs only
when `JointFitConfig::coverage_odds` is filled. It is empty by default (`fit.rs:617`), and only
`examples/ng_joint_duplicated_in_fit.rs` fills it (`examples/ng_joint_records_walk.rs` names the
field but passes an empty vector, `:302`).

### 2.9 Parameter estimation: the repeat-tract fit (`joint/ssr_fit.rs`)

**The loop structure.** Per stratum, `climb_from` runs up to `max_rounds` rounds
(`ssr_fit.rs:786-807`). Each round climbs every slippage number, every length-spectrum class and the
concentration by golden section, 18 evaluations each (`climb_scalar`, `ssr_fit.rs:2460-2480`).
Each evaluation calls `Scorer::score`, which rebuilds the per-tract likelihoods when slippage
changed (`refresh`, `:2180-2215`), rebuilds the Dirichlet quadrature of 256 points when the
spectrum or concentration changed (`:2143-2152`), and then scores every tract at every point
(`ln_tract`, `:2169, 2174`).

| site | fn | enclosing fn | the value | how often | constant |
|---|---|---|---|---|---|
| `ssr_fit.rs:125` | powi | `Slippage::read_probabilities` | probability that a slip is exactly `step` repeats, `(1 − fall_off) × fall_off^(step−1)` | **per-fit-iteration** — per read group × 13 allele offsets × up to 14 steps, per `refresh` (`:2192`) | no |
| `ssr_fit.rs:130` | powi | same | probability of slipping past the recorded range, `fall_off^steps` | **per-fit-iteration** — per group × offset × direction per `refresh` | no |
| `ssr_fit.rs:984, 1005` | ln (2 lines) | `climb_one_round` | the starting point of a climb on a length-spectrum class and on the concentration | **per-fit-iteration** — per class / once, per round | no |
| `ssr_fit.rs:988, 1002` | exp (2 lines) | same | the trial value at each golden-section evaluation | **per-fit-iteration** — per evaluation | no |
| `ssr_fit.rs:995, 1008` | exp (2 lines) | same | the accepted value | **per-fit-iteration** — per class / once, per round | no |
| `ssr_fit.rs:1948` | ln | `TractLikelihoods::of` | log probability of a read-length bucket under a genotype, floored at 1e-300 | **per-fit-iteration** — per tract × sample × genotype × read group × non-empty bucket, per `refresh` | no |
| `ssr_fit.rs:1957` | exp | same | each genotype's likelihood relative to the sample's best | **per-fit-iteration** — per tract × sample × genotype, per `refresh` | no |
| `ssr_fit.rs:2073` | ln | `ln_tract` | log of the product over samples at one quadrature point | **per-fit-iteration** — per tract × 256 points, per evaluation | no |
| `ssr_fit.rs:2305` | ln | `dirichlet_points` | quadrature weight, `−ln 256` | **per-fit-iteration** — per quadrature rebuild | **default** |
| `ssr_fit.rs:2332` | exp, ln | `ln_sum_exp` | log-sum-exp over the 256 point terms | **per-fit-iteration** — exp per tract × 256, ln per tract, per evaluation (`:2076`) | no |
| `ssr_fit.rs:2349` | sin, ln | `ln_gamma` | the reflection formula for arguments below one half | **per-fit-iteration** — `ln_beta` (`:2366`) per stick of the stick-breaking draw, per point per rebuild (`:2255-2260`), when a concentration is below 0.5 | no |
| `ssr_fit.rs:2357` | ln ×3 | same | `½ ln 2π`, `ln t` and `ln(series)` of the Lanczos approximation | **per-fit-iteration** — same chain | one of three **literal** |
| `ssr_fit.rs:2388` | exp, ln ×2 | `regularised_incomplete_beta_with` | the incomplete beta function's prefactor | **per-fit-iteration** — up to 60 bisection steps (`beta_quantile_with`, `:2436-2446`) × 12 sticks × 256 points per quadrature rebuild | no |
| `ssr_fit.rs:2452` | ln | `logit` | logit of a slippage number, clamped to [1e-9, 1 − 1e-9] | **per-fit-iteration** — per slippage number per round (`:975`) | no |
| `ssr_fit.rs:2456` | exp | `expit` | the logistic function | **per-fit-iteration** — per golden-section evaluation (`:972`) and per accepted value (`:978`) | no |

### 2.10 Parameter estimation: contamination (`joint/contamination.rs`)

`fit_contamination_over` runs inside `fit_jointly` (`fit.rs:1526`), which both
`estimate-parameters` and `estimate-contamination` (`cli/estimate_contamination.rs:480`) call. It
builds ancestry coordinates once, then for each library runs golden-section searches of 30 steps
(`contamination.rs:1180-1201`) over the contamination fraction and each coordinate
(`:1132-1161`); every score sums `ln_marker` over the library's covered markers (`:1088-1110`).

| site | fn | enclosing fn | the value | how often | constant |
|---|---|---|---|---|---|
| `contamination.rs:1250` | ln ×2 | `ln_marker` | log prior weight of an (own genotype, contaminant genotype) pair, and log of the averaged binomial probability over the depth bin | **per-fit-iteration** — per marker × up to 9 pairs × score evaluation | no |
| `contamination.rs:1278` | powi ×2 | `binomial` | `p^k (1 − p)^(n−k)` | **per-fit-iteration** — per marker × pair × depth in the bin × evaluation; also per marker × sample × 3 once when dosages are built (`:805`) | no |
| `contamination.rs:1286` | exp, ln | `ln_sum_exp` | log-sum-exp over the 9 pairs | **per-fit-iteration** — per marker × evaluation | no |
| `contamination.rs:613` | sqrt | `fit_contamination_over` | spread of the coordinates along each axis | **per-run** | no |
| `contamination.rs:880` | sqrt | `ancestry_coordinates` | a marker's scaling, `√(p(1−p))` | **per-run** — per marker, once | no |
| `contamination.rs:1311, 1312` | sqrt (2 lines) | `leading_eigenvectors` | the rotation in a Jacobi eigenvalue sweep over the sample similarity matrix | **per-run** — up to 100 sweeps × samples²/2 pairs | no |

## 3. Input ranges the hot calls see

For each function with calls in the per-base, per-read, per-genotype or per-fit-iteration classes
(reachable calls only), the arguments those calls receive, and the lines that bound them. "Not
bounded" means I found no clamp or floor in the code; the practical range is then a property of
the data.

### 3.1 `ln`

| call sites | argument | bounded by |
|---|---|---|
| per-read: `emission.rs:558` | `1 − ε`, ε ∈ [0, 1]: so [2.2e-308, 1] | `FlatEmission::try_new` refuses ε outside [0, 1] (`emission.rs:535-537`); floor `f64::MIN_POSITIVE` (`emission.rs:139`) |
| per-read: `emission.rs:559` | `ε / 3`: [2.2e-308, 1/3] | same |
| per-read: `ssr_unit_robust.rs:278-288` | four fixed values in the shipped run: 0.88, 0.0475, 0.0475, 0.05 | `stutter.rs:308-317`; for a fitted model, shares are clamped to [0, 1] and one-step shares to [0.01, 0.99] (`stutter.rs:104-106, 857, 868`), and the same-length share to [0.01, 1] (`stutter.rs:287`) |
| per-read: `generic.rs:810, 912` | `(1 − c)·k/P + c·q`, with *k/P* from 1/16 to 1 | *P* ≤ 16 (`calling/likelihood/mod.rs:154, 173`); *c* and *q* are fractions, so the argument is in about [(1 − c)/16, 1]. *c* itself is not clamped at this site |
| per-read and per-genotype: `generic.rs:830, 903` | `(1 − c)·ε̄/spread + c·q` | ε̄ ≥ 1e-12 (`MIN_BASE_ERROR`, `calling/likelihood/mod.rs:207, 510`) and **not capped above** (the scale can push it past 1; `mod.rs:2917-2919` is a test of exactly that); spread is 1 or 3 (`generic.rs:26, 31, 157-161`). Lower end about 3e-13 when *c* = 0 |
| per-genotype: `ssr.rs:736` | the mixture `(1 − λ − c)·explained + λ·allowed/reachable + c·contaminant` | at least λ / (number of reachable lengths) because junk is always added, λ = 0.20 by default (`ssr.rs:155, 696-697`); at most about 1. `1 − λ − c` is floored at 1e-12 (`ssr.rs:254`) |
| per-genotype: `dirichlet_multinomial.rs:78, 322` | an allele's concentration α | ≥ 1e-12 (`MIN_ALT_CONCENTRATION`, `genetics.rs:107`, applied where each seed is built: `calling/genotype_prior/seed_generic.rs:417, 425`, `calling/genotype_prior/seed_ssr.rs:183`, `genetics.rs:147`); α = seed + the other samples' expected copies (`genotype_prior/mod.rs:781`), so up to about seed + (samples − 1) × ploidy. Not bounded above by a clamp |
| per-genotype: `dirichlet_multinomial.rs:136, 292` | a rising product `α(α+1)…(α+k−1)`, *k* ≤ ploidy | same lower bound on α; grows as α^ploidy |
| per-genotype: `dirichlet_multinomial.rs:282` | total concentration | same |
| per-genotype: `dirichlet_multinomial.rs:307, 308` | *F* and `1 − F` | *F* ∈ [0, 1] (asserted in debug, `:268-273`); `ln 0 = −∞` is possible at `:307`; `:308` floors at 1e-300 |
| per-genotype: `dirichlet_multinomial.rs:372`, `quality/mod.rs:613` (per-locus), `joint/fit.rs:1199`, `joint/ssr_fit.rs:2332`, `joint/contamination.rs:1286` | the sum in a log-sum-exp | [1, number of terms]: every term is `exp` of a value ≤ 0 and the largest is exactly 1 |
| per-genotype: `quality/mod.rs:588` | the peak of the count vector, taken before it is rescaled | (0, ploidy + 1]: each entry sums up to `ploidy + 1` products of a previous entry ≤ 1 and a weight ≤ 1 (`quality/mod.rs:560-576`); not clamped |
| per-genotype: `quality/mod.rs:597` | a rescaled count probability | (0, 1] (`value > 0` checked, `:596`) |
| per-genotype: `paralog/locus_score.rs:490` | σ = the sample's single-copy depth SD × √(T/2) | σ₀ finite and > 0 (`locus_score.rs:311`); √(T/2) ∈ [1.22, 2]. Not bounded above |
| per-fit-iteration: `paralog/prior.rs:311` | `p/(1−p)` | p clamped to [1e-12, 1 − 1e-12] (`prior.rs:316`): [1e-12, 1e12] |
| per-fit-iteration: `joint/fit.rs:1051, 1176, 2207-2235` | probabilities | floored at `f64::MIN_POSITIVE`; ≤ 1 |
| per-fit-iteration: `joint/fit.rs:1167` | a weighted sum whose first term is 1 | ≥ 1 (the first depth weight and the first power are both 1, `fit.rs:1095, 1162-1164`); upper end not bounded by a clamp |
| per-fit-iteration: `joint/fit.rs:2132, 2175`, `joint/ssr_fit.rs:2073` | a product of per-sample terms ≤ 1, rescaled by 1e150 whenever it falls below 1e-150 | [about 1e-300, 1] (`fit.rs:2000-2001, 2123-2126`; `ssr_fit.rs:2065-2068`) |
| per-fit-iteration: `joint/fit.rs:2183, 2190, 2196, 2218, 2327, 2364, 2442` | 3, or 1, 2, 3 | `fit.rs:941, 2047, 2055` |
| per-fit-iteration: `joint/fit.rs:2824, 2825` | quadrature nodes and one minus them | [1e-12, 1 − 1e-12] (`fit.rs:2820`) |
| per-fit-iteration: `joint/fit.rs:2828` | normalised quadrature weights | floored at `f64::MIN_POSITIVE`, ≤ 1 |
| per-fit-iteration: `joint/fit.rs:2950` | *x* after the recurrence | ≥ 8 (`fit.rs:2944`); the shapes are clamped to [0.02, 50] after each Newton step (`fit.rs:2758-2759`), so `a + b` ≤ 100 and, since the recurrence only raises arguments below 8, *x* is at most 100 from the second step on. The first step's starting shapes are not clamped |
| per-fit-iteration: `joint/ssr_fit.rs:984` | a length-spectrum share | floored at 1e-9 (`:984`), ≤ 1 |
| per-fit-iteration: `joint/ssr_fit.rs:1005` | the concentration | not bounded by a clamp at this site |
| per-fit-iteration: `joint/ssr_fit.rs:1948` | a bucket probability | floored at 1e-300 (`:1948`), ≤ 1 |
| per-fit-iteration: `joint/ssr_fit.rs:2349` | `π / sin(πx)`, *x* ∈ [1e-3, 0.5) | the concentrations are floored at 1e-3 (`:2241`); argument in (π, about 1,000] |
| per-fit-iteration: `joint/ssr_fit.rs:2357` | 2π (literal); *t* = *x* + 6.5 ≥ 7; the Lanczos series | `:2350-2357`; *x* ≥ 0.5 on this branch |
| per-fit-iteration: `joint/ssr_fit.rs:2388` | *x* and `1 − x`, *x* a bisection midpoint in (0, 1) | `:2377-2386` returns early at 0 and 1 and mirrors *x* above the mode |
| per-fit-iteration: `joint/ssr_fit.rs:2452` | `p/(1−p)` | p clamped to [1e-9, 1 − 1e-9] (`:2451`): [1e-9, 1e9] |
| per-fit-iteration: `joint/contamination.rs:1250` | a prior weight in (0, 1], and an averaged binomial probability | the second is floored at `f64::MIN_POSITIVE` (`:1250`) |

### 3.2 `exp`

| call sites | argument | bounded by |
|---|---|---|
| per-read: `calling/likelihood/mod.rs:510` | a group's mean log error | [−58.72, 0]: each read's log error is `max(ln 10^(−q/10) at the base quality, the mapping-quality term)` with *q* a `u8` (`open_record.rs:3212, 3233`; `read/prepared_read.rs:154`), so no read goes below −25.5 × ln 10 ≈ −58.72; the sum is asserted ≤ 0 (`mod.rs:498-499`). No cap at Phred 93 exists in this path |
| per-read: `ssr_emission.rs:629` | the sum over an observation's bases of `ln(1 − ε)` or `ln(ε/3)` | ≤ 0. **Not bounded below** by the code: length × ln(ε/3), so several hundred to a few thousand below zero for a long read at ε near 1e-3, where the result underflows to 0 |
| per-genotype: `summarise_condition.rs:364`, `quality/mod.rs:562, 612`, `dirichlet_multinomial.rs:372`, `paralog/locus_score.rs:528, 570, 573` | a log score minus the largest (or the larger) score | ≤ 0. Not bounded below; a genotype the reads rule out gives −∞, and the three `quality`/`summarise` sites skip or tolerate it (`quality/mod.rs:561`; `summarise_condition.rs:355` asserts only the largest is finite) |
| per-fit-iteration: `paralog/prior.rs:301, 303` | `−|x|`, *x* = bin centre + logit(π) | bin centres in [−100, 100] (`prior.rs:50, 52`), logit(π) in [−27.6, 27.6] (`prior.rs:316`): argument in [−127.7, 0] |
| per-fit-iteration: `joint/fit.rs:1199, 2098, 2278-2296, 2332, 2377, 2456`, `joint/ssr_fit.rs:1957, 2332`, `joint/contamination.rs:1286` | differences from the largest term | ≤ 0, not bounded below (−∞ for impossible branches) |
| per-fit-iteration: `joint/ssr_fit.rs:988, 995` | log of a length-spectrum share ± 2 | the start is ≥ ln 1e-9 ≈ −20.7 (`:984`), the climb spans ±2 (`:993`) — so about [−22.7, 2] |
| per-fit-iteration: `joint/ssr_fit.rs:1002, 1008` | log concentration ± 2.5 | not bounded by a clamp |
| per-fit-iteration: `joint/ssr_fit.rs:2388` | `a ln x + b ln(1−x) + ln B(a, b)` | ≤ about 0 in practice; not clamped. *a*, *b* ≥ 1e-3 (`:2241`) |
| per-fit-iteration: `joint/ssr_fit.rs:2456` | `−x`, *x* = logit of a slippage number ± 3 | logit clamped to about ±20.7 (`:2451`), the climb spans ±3 (`:976`): argument in about [−23.7, 23.7] |

### 3.3 `ln_1p`

| call site | argument | bounded by |
|---|---|---|
| per-genotype: `paralog/locus_score.rs:532` | `x = exp(lo − hi)` | [1.4e-4, 1]: below 1.4e-4 a series replaces the call (`locus_score.rs:529-531, 542`); *x* ≤ 1 because `lo ≤ hi` (`:527`) |

### 3.4 `powi`

| call site | base | exponent | bounded by |
|---|---|---|---|
| per-read: `alignment/stutter.rs:783` | `1 − one-step share`, in [0.01, 0.99] | `size − 1`, 0 to 9 | `stutter.rs:868` (share clamp), `:772` (sizes above 10 return early), `:65, 92` |
| per-fit-iteration: `joint/ssr_fit.rs:125, 130` | `fall_off`, in (0, 1) | 0 to 14 | `fall_off` is `expit` of a clamped logit (`:2451-2456`); steps ≤ `read_span − allele_offset` with `read_span = 8` (`joint/census.rs:420`, `ssr_fit.rs:2603`) and |offset| ≤ 6 (`ssr_fit.rs:536`) |
| per-fit-iteration: `joint/contamination.rs:1278` | *p* and `1 − p`, in [1e-12, 1 − 1e-12] | *k* and `n − k`: a marker's alternative count and its depth in the bin, 0 up to the summed depth of the library's read groups | `:1277`; each read group's depth is capped at 255 (`joint/census.rs:841`) and a library's groups are summed (`contamination.rs:849-851`). How many groups a library holds is not bounded in this file |

### 3.5 `sin`

| call site | argument | bounded by |
|---|---|---|
| per-fit-iteration: `joint/ssr_fit.rs:2349` | π*x*, *x* ∈ [1e-3, 0.5) | the branch is `x < 0.5` (`:2348`); concentrations ≥ 1e-3 (`:2241`); argument in [3.1e-3, π/2) |

### 3.6 `sqrt`

The three per-fit-iteration calls (`joint/fit.rs:2865, 2902, 2903`) take non-negative numbers
(`:2865` applies `max(0.0)`; `:2902-2903` add 1 to a square). `sqrt` is correctly rounded on every
platform, so these need no micro-benchmark.

## 4. `powi` in shipped code

`powi` is not a maths-library call in Rust: it compiles to an LLVM intrinsic that becomes a short
multiplication loop, the same source on every platform. Whether its bits agree across platforms is
for step A2 to measure (the plan's step B1 depends on it). The 9 shipped calls:

| site | base | exponent | exponent at run time | class | reachable |
|---|---|---|---|---|---|
| `alignment/stutter.rs:486` | `1 − one-step share` | `steps`, capped at 10 | **variable** | per-locus | yes |
| `alignment/stutter.rs:783` | `1 − one-step share` | `size − 1`, 0–9 | **variable** | per-read | yes |
| `alignment/ssr_marginal_sequence.rs:145` | `1 − ε` | read length | **variable** | per-read | no |
| `parameter_estimation/depth_bins.rs:241` | the widening ratio | `step`, 1–11 | **variable** (loop counter; base and range from literals) | per-run | no |
| `parameter_estimation/depth_bins.rs:299` | the widening ratio | `step`, 1–10 | **variable** (loop counter; base and range from literals) | per-run | yes |
| `parameter_estimation/joint/ssr_fit.rs:125` | `fall_off` | `step − 1` | **variable** | per-fit-iteration | yes |
| `parameter_estimation/joint/ssr_fit.rs:130` | `fall_off` | `inside_steps` | **variable** | per-fit-iteration | yes |
| `parameter_estimation/joint/contamination.rs:1278` (first) | `p` | `k` | **variable** | per-fit-iteration | yes |
| `parameter_estimation/joint/contamination.rs:1278` (second) | `1 − p` | `n − k` | **variable** | per-fit-iteration | yes |

No shipped `powi` has a literal exponent. The 27 test calls include literal exponents (for example
`fit.rs:3401`, `powi(3)`, in the bench fixtures).

## 5. What this inventory does not establish

- **How many times each class runs on real data.** The classes come from reading the loops, not
  from a profile. The plan allows a sampling profile where the reading is not decisive; I did not
  need one to assign a class, but the relative cost of, say, the joint fit's per-position loop
  against the calling loop's per-genotype loop is not measured here.
- **The number of passes and optimiser rounds.** Per-genotype calls repeat every pass and
  per-fit-iteration calls every round; neither count is known from the code alone.
- **x86_64.** Nothing in this inventory depends on the processor, but step A2's timings will.
- **Calls in `examples/`**, which the plan leaves out of scope.

## Review applied

The review is `doc/devel/reports/reviews/portable_float_A1_review_2026-09-13.md`. Each finding was
checked against the code before the text changed; I agreed with all of them. The working file
`tmp/portable_float_A1_sites.tsv` was updated with the report. Recounted from it: 170 reachable,
31 not reachable, 1 gated.

| finding | what I checked | what changed |
|---|---|---|
| 1 (Fix) `depth_bins.rs:241` not reachable | `ladder_tops` ← `DepthBinEdges::new` (`:270`) ← `Default` (`:201-204`); every `new()`/`default()` call is in the test module or `examples/`; no struct holding one derives `Default`; `CensusWriter.edges` is filled from `for_census` at both shipped `CensusWriter::new` calls (`run/gatherer.rs:301-311`, `cli/estimate_contamination.rs:785-791`) and `JointFitConfig`'s own `Default` uses `for_census` (`fit.rs:614`) | site marked not reachable in §2.7, §4 and the working file; reachable 171 → 170, not reachable 30 → 31 in the Words list, Summary and per-file table (`depth_bins.rs` 3 shipped, 2 reachable); `powi` per-run "2 (1)", `powi` total "9 (7)", per-run total "46 (45)", grand total "202 (170)" |
| 2 (Fix) §1.1 command printed 381, not 434 | ran it: 381 lines without `-o`, 434 with the extra stage | added the `rg -o … \| wc -l` stage and a note on the 381 |
| 3 (Fix) Summary multiplied loop dimensions of different calls | the sample and node loops at `fit.rs:2115-2127` and `ssr_fit.rs:2031-2069` only multiply and add | per-fit-iteration bullet rewritten call by call |
| 4 (Fix) `genotype_table.rs:535` per locus past 16 candidates | cache bound at `:57`, `:175-176`; table built per locus at `summarise_condition.rs:2315`; repeat-tract cap 32 at `allele_candidates/ssr.rs:477-478` | §2.3 row gives both frequencies; a note under the Summary table says the call is counted once, as per-table; working file's constant column corrected to `small-integer` |
| 5 (Nit) `digamma` × 4 calls, *x* ≤ 100 | `fit.rs:2747-2748` makes four calls; clamp at `:2758-2759` bounds `a + b` by 100 | §2.8 and §3.1 corrected |
| 6 (Nit) arguments fixed inside their loops | each case read at the lines the review gives; for `locus_score.rs:490` the factors are 1 and √(T/2), 5 distinct per sample | new Summary paragraph listing the cases |
| 7 (Nit) partial reads multiply the STR emission calls | `censored_emission` loops over `reachable_length_changes` (`ssr_emission.rs:437-448`); `probability_at_least_this_much_longer` sums the same iterator (`stutter.rs:687`) | added to the §2.1 rows for `emission.rs:558` and `stutter.rs:783`, the §2.2 row for `ssr_emission.rs:629`, and the Summary's per-read bullet |
| 8 (Nit) `ng_joint_records_walk.rs` does not fill `coverage_odds` | `examples/ng_joint_records_walk.rs:302` is `Vec::new()` | §2.8 corrected; the not-reachable conclusion stands |
| 9 (Nit) concentration floor cited the constant only | floors at `seed_generic.rs:417, 425`, `seed_ssr.rs:183`, `genetics.rs:147` | §3.1 cites them |
| 10 (Nit) `quality/mod.rs:588` peak can exceed 1 | the `ln` takes the peak before the rescale (`:579-588`) | range changed to (0, ploidy + 1] |
| 11 (Nit) table citation | `static PER_QUALITY_LN` at `emission.rs:170` | Summary cites `:170` |
