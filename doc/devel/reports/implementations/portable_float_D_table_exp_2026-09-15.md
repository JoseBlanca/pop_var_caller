# A table-driven `exp`: what it wins back, and what it moves

*Milestone D of [`portable_float.md`](../../implementation_plans/portable_float.md). Branch
`portable-float`, on top of commit `caca7b62` (Milestone C). Measured 2026-09-15.*

**The answer.** Computing `exp` with the table-driven algorithm glibc and musl use, ported to Rust,
wins back most of what the portable maths cost the parameter fit. Measured directly against the
unchanged caller (the one that called the platform's maths library), interleaved:

| | unchanged build | portable maths, `libm`'s `exp` (Milestones B–C) | portable maths, table `exp` (this milestone) |
|---|---|---|---|
| `estimate-parameters`, macOS | baseline | +27% (B7) | **+3.4%** |
| `estimate-parameters`, Linux | baseline | +29% (B7) | **+8.0%** |

The same bits on both platforms, as before; the fitted parameters moved again, in their last digits,
and calls with default parameters did not.

**Terms.** *The unchanged build* is commit `c6a4394b`, whose source is identical to `75722336`:
std's maths methods, the platform's library. *The C build* is this branch at Milestone C
(`libm`'s `exp`, the fit's fixed chunking). *The D build* is this milestone. *A step* is the gap
between one `f64` and the next.

**Machines and data** as in the earlier reports: macOS is an Apple M5 Pro host (18 cores); Linux is
the arm64 glibc container VM on it (8 CPUs); four tomato accessions at about 3× over 20 regions for
the fit, 160 regions for calling. One low-depth corner of the caller's range. x86_64 is not covered.

## 1. The function

`src/float/table_exp.rs` is Arm's `exp` (optimized-routines, 2018, MIT) ported line for line from
musl's `src/math/exp.c`: `eˣ = 2^(k/128) · eʳ`, the first factor from a 128-entry table, the second a
degree-5 polynomial in the small remainder `r`. The constants are written as bits; the table is
generated at 200 bits by `scripts/float_exp_table.py`, which checks six entries against musl's
`exp_data.c` and, with `--check`, that the Rust table is its output. Arm's copyright and permission
notice heads the file and is appended to `LICENSE`.

**Correctness, checked three ways.**

- The review (`doc/devel/reports/reviews/portable_float_D1_review_2026-09-15.md`) recomputed all 10
  constants and all 256 table words independently, and a Python transliteration of musl's C gave the
  same bits as the Rust port on all 5,242,880 arguments of §2.
- `the_grid_has_musls_bits` pins a checksum of the port's output over 200,001 arguments spanning
  every finite, non-zero result, computed from that transliteration; it passes on macOS and Linux.
  Unlike a within-one-step test, it fails for an error smaller than a step: the review showed a
  Taylor coefficient, one zeroed table tail and the removed subnormal rounding each change it.
- Ten boundary arguments with bits checked against correct rounding at 300 bits by the review, plus
  `±f64::MAX`, a negative NaN and the edges (zeros, NaN, infinities, overflow and underflow).

## 2. Per call: speed and accuracy

`examples/float_exp_check.rs` times the platform's `exp`, `libm`'s and the port over five ranges of
524,288 arguments each, the median of 21 passes, rotating which goes first, and writes every output
(`tmp/exp_D2/{macos,linux}/`, timings in `timings.tsv`). Nanoseconds per call:

| range | macOS: Apple / libm / table | Linux: glibc / libm / table |
|---|---|---|
| log-probabilities, [−700, 0] | 1.46 / 1.80 / **1.49** | 2.10 / 1.80 / **1.48** |
| log-ratios, [−20, 20] | 1.42 / 1.96 / **1.39** | 1.39 / 1.98 / **1.39** |
| past underflow, [−745, −700] | 1.80 / 2.15 / **2.13** | 2.54 / 2.15 / **2.17** |
| near overflow, [0, 709.7] | 1.39 / 1.79 / **1.38** | 2.07 / 1.80 / **1.39** |
| near zero, [−1e-3, 1e-3] | 1.41 / 1.39 / **1.39** | 1.39 / 1.39 / **1.39** |

The port takes 0.70 to 0.83 of `libm`'s time on the three wide ranges, the same past underflow and
near zero, and about the platform library's time or less everywhere except past underflow on macOS
(1.18 times Apple's). The macOS run started at a host load of 10.5 as an earlier job finished. A
first run at a load of 6.8 (`tmp/exp_D1/`) gave the port's and `libm`'s columns to within 0.04 ns of
these; the platform libraries' columns moved by up to 0.20 ns (Apple, log-ratios) and 0.16 ns (glibc,
near overflow).

**Accuracy** (`scripts/float_exp_accuracy.py`, 32,768 arguments a range spread over its values,
against exp at 200 bits rounded to the nearest double): how many are *not* the correctly rounded
double. Every miss by every implementation is one step.

| range | table | libm | glibc (Linux) | Apple (macOS) |
|---|---:|---:|---:|---:|
| [−700, 0] | 32 | 3,187 | 23 | 42 |
| [−20, 20] | 40 | 3,152 | 23 | 67 |
| [−745, −700] | 3 | 862 | 3 | 15 |
| [0, 709.7] | 26 | 3,220 | 22 | 37 |
| [−1e-3, 1e-3] | 11 | 10 | 11 | 69 |

So the port rounds the wrong way about one argument in a thousand, `libm` about one in ten. The
port's and `libm`'s outputs are identical on macOS and Linux for every argument; the two platforms'
own libraries disagree on 196 to 1,067 arguments a range. glibc, the same algorithm, differs from the
port on 274 to 404 of the 524,288 arguments in the three wide ranges and on none past underflow or
near zero (`tmp/exp_D2/linux/timings.tsv`).

## 3. The caller

### 3.1 Against the C build

Interleaved on the same psps and regions, both builds alternating first (`tmp/d_measure/`):

| command | platform | rounds | C build | D build | change in mean |
|---|---|---:|---|---|---:|
| `estimate-parameters` | macOS | 3 | 443.0–448.1 s | 369.3–374.8 s | **−16.8%** |
| `estimate-parameters` | Linux | 2 | 682.4–685.7 s | 539.2–544.8 s | **−20.8%** |
| `call-from-alignments`, 160 regions | Linux | 5 | 19.26–19.57 s | 19.19–19.26 s | −1.1% |
| `call-from-psps`, 160 regions | Linux | 5 | 8.19–8.24 s | 8.06–8.20 s | −1.1% |

The slowest D fit is faster than the fastest C fit on both platforms. **The macOS calling rounds are
not quoted**: both builds' times spread by up to a quarter across rounds (for example
`call-from-psps` 6.48 to 8.22 s) with no process other than the run using more than about a core at
the logged samples, so a difference of a few per cent cannot be read from them. The host load during
the fits was 8 to 27, median 14, mostly the fit itself on every core.

**Benches**, same session, passes in the order C, D, D, C (`tmp/d_measure/benches_table.txt`): every
case is faster with the D build.

| | macOS | Linux |
|---|---|---|
| joint fit, repeat-tract cases (6) | −7.6 to −10.5% | −3.5 to −7.3% |
| joint fit, SNP/indel cases (5) | −9.3 to −19.4% | −10.8 to −16.6% |
| site quality, by sample count (4) | −1.2 to −17.8% | −1.2 to −11.2% |
| site quality, by allele count (3) | −3.2 to −9.6% | −2.4 to −8.3% |

Each build's two passes agreed to within 3% except in four cases, spread 3.0 to 4.9%. In one of them
the change is inside the spread: Linux, 32 tracts, whose C passes were 4.9% apart against a change of
−4.4%.

### 3.2 Against the unchanged build

The question the owner's request turns on: how much of the conversion's cost does this win back?
Measured directly, interleaved, the D build against the unchanged build's original binaries
(`tmp/ab/bin/std_*`, the same ones B7 timed; `tmp/d_vs_std/`):

| command | platform | rounds | unchanged | D build | change in mean |
|---|---|---:|---|---|---:|
| `estimate-parameters` | macOS | 3 | 347.3–372.2 s | 366.5–370.6 s | **+3.4%** (per round −0.4, +6.1, +4.8%) |
| `estimate-parameters` | Linux | 2 | 481.1–483.3 s | 519.8–522.1 s | **+8.0%** (+7.6, +8.5%) |
| `call-from-psps`, 160 regions | Linux | 5 | 8.00–8.12 s | 8.01–8.10 s | +0.2% |
| `call-from-alignments`, 160 regions | Linux | 5 | 19.17–19.52 s | 19.08–22.22 s | +4.9%; three of five rounds within ±1.2% |

- **The fit.** B7 measured the conversion with `libm`'s `exp` at +27% on macOS and +29% on Linux, in
  another session. Now +3.4% and +8.0%. The macOS figure is the less certain: the unchanged build's
  first round (372.2 s) was 7% slower than its other two, during a burst of two background services
  at about a core each; without it the D build is +5.4% (both later rounds).
- **Calling.** `call-from-psps` is level on Linux. `call-from-alignments` had two D rounds at 21.1 and
  22.2 s among three at 19.1–19.5 s, so its mean is not a measurement of the build. The macOS calling
  rounds were as noisy as in §3.1 and are not quoted.
- The host load over this run was 3.5 to 33, median 8.9 (`tmp/d_vs_std_hostload.tsv`).

## 4. Output

- **Calls with default parameters did not move**: over 160 regions, `call-from-alignments` and
  `call-from-psps` wrote the same records with the C and D builds on both platforms (120,538
  records; `tmp/d_measure/calls_{macos,linux}/`).
- **The identity oracle**, both platforms identical (`tmp/d_measure/oracle_{linux,macos}/digests.txt`;
  the macOS run compared itself against Linux's and matched): the five default-parameter checksums
  unchanged; the fitted parameters `502c3d5e…` → `be31e78a…` and the calls with that fit
  `b7cbf433…` → `18f4fc11…`. `scripts/promote_ng_oracle.baseline` is re-recorded.
- **The fitted parameters**, on B7's psps (Linux): 55 of 574 lines differ from the C build's fit and
  32 from the unchanged build's. With `libm`'s `exp` the conversion had moved 56 from the unchanged
  build's, so the table `exp`, glibc's algorithm, brings the fit closer to what the unchanged caller
  wrote on Linux.
- **Calls made with each build's own fit**, B7's psps, Linux (`tmp/d_measure/calls_with_fit_d.records`):

  | against | records differing, of 6,735 | fields | genotypes, FILTER or sites |
  |---|---:|---|---:|
  | C build | 19 | AF 5, PARALOG_POST 12, PARALOG_LR 1, QUAL 1 | 0 |
  | unchanged build | 18 | AF 8, PARALOG_POST 9, PARALOG_LR 1 | 0 |

- **The cross-platform test's fit checksum** moved from `28722984…` to `3edab375…`, the same on both
  platforms (`tmp/digests_D1/{macos,linux}.log`: the test passes against the new value); its calls
  checksum did not move.
- **One pinned value in `float`'s tests moved**: `exp(13.5)` is now the correctly rounded
  `0x41264290bd5cad8b`, 0.44 of a step from the true value, where `libm` gave the next double up,
  0.56 away.

## 5. Tests

`cargo test --lib --tests --all-features`: **Linux container 4,804 passed, 0 failed, 4 ignored;
macOS 4,803 passed, 0 failed, 4 ignored.** In the container `cargo fmt --check`,
`cargo clippy --all-targets --all-features -- -D warnings` and `cargo doc --no-deps --all-features`
with warnings as errors passed. Five tests were added in `float::table_exp`.

## 6. What was not covered

- **x86_64.** The port has no processor-specific code, but no x86_64 machine ran it.
- **macOS calling commands**, whose rounds were too noisy in both sessions (§3).
- **Other `libm` functions.** `ln` is the next-largest (about 10% of the unchanged fit's CPU) and
  stays on `libm`; its speed against glibc was 1.2 to 1.4 times in A2.
- **GIAB and larger cohorts**, as in B7.

## 7. Deviations

- **§3.2 was added.** The plan asked for the D build against the Milestone B build; the owner's
  question at Checkpoint B was whether this recovers the conversion's slowdown, which only a direct
  comparison with the unchanged build answers. The baseline for §3.1 is the C build, not B, because
  Milestone C changed the fit's chunking (plan §8).
- **The accuracy figures were corrected by review**: the first scoring sampled only the first
  sixteenth of each range.
