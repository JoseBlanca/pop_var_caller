# The platform maths library against libm, function by function

*Step A2 of [`portable_float.md`](../../implementation_plans/portable_float.md). Branch
`portable-float`. Measured 2026-09-14. Nothing under `src/` changed.*

Rust's `f64::ln`, `exp`, `powf` and the rest call the operating system's maths library: glibc on
Linux, Apple's libm on macOS. The two round differently in the last binary place, which is what
moved a repeat tract by a base on macOS on 2026-09-13. The candidate fix routes those calls through
the `libm` crate, one Rust implementation compiled into the binary. For each function the caller
uses, this report answers:

1. how often macOS and Linux disagree today;
2. whether libm gives the same bits on both;
3. how much slower libm is, and how often it rounds differently from each platform's own library —
   which is how much Linux's own output would move, once.

It measures isolated calls. What the change costs the caller's commands, and whether it moves a
VCF, is step A3.

## How it was measured

The program is [`examples/float_libm_vs_std.rs`](../../../../examples/float_libm_vs_std.rs). For each
function and each range of arguments it:

- evaluates 524,288 arguments spread evenly over the range, visited in a scrambled order, through
  std and through libm;
- runs 20 untimed passes of each, then 21 timed passes, alternating which side goes first, and
  reports the median;
- counts the arguments whose two outputs are not the same bits, and the largest distance between
  them in *units in the last place* — steps between adjacent representable numbers;
- writes every output to a file. Its `compare` mode reads a macOS directory and a Linux directory and
  counts, argument by argument, where they differ.

The arguments are generated through libm, so both platforms evaluate the same numbers. The ranges
are those A1's inventory (§3) found at the caller's hot calls, widened to its extremes where §3 says
the code does not bound them (`ln` up to 1e12 and down to the smallest normal number, `exp` below
−708 where results become subnormal). `powf`, `log10` and `exp_m1` have no hot call in A1, so their
ranges are plausible arguments, not measured ones.

A first set of runs, reviewed in
[`portable_float_A2_review_2026-09-14.md`](../reviews/portable_float_A2_review_2026-09-14.md), is
superseded: it visited the arguments in a regular sawtooth, which the branch predictor learned, and
warmed up for one pass only. Every figure below is from the four runs after those fixes.

**The machines.** macOS is the host, an Apple M5 Pro, built natively. Linux is the Apple `container`
VM on the same host, arm64 with glibc, 8 CPUs. The two builds do not share compiler flags (the host's
cargo configuration targets `apple-m1`), so compare std with libm *within* a platform, not macOS
with Linux. **Linux on x86_64 — the CI runner and the Linux dev box — is not measured.** Its glibc
picks between two builds of `exp` and `log` when a program starts, by whether the processor has
fused multiply-add, and neither the speeds nor the disagreement counts here transfer to it.

## Speed

Nanoseconds a call, median of 21 passes in the first of the four runs; *libm ÷ std* is that run's
ratio with the range over all four runs in brackets. Above 1, libm is slower.

**The first row is the loop with no function in it**, about 1 ns. Every other row includes that
nanosecond on both sides, so the ratios understate the functions' own ratios; the nanosecond
differences, which are what a hot loop pays, are exact. Absolute times moved between runs by up to
half (the Linux VM's fourth run was about twice as fast as its first three on every row); the ratios
did not.

| function | arguments | macOS std | macOS libm | macOS libm ÷ std | Linux std | Linux libm | Linux libm ÷ std |
|---|---|---:|---:|---|---:|---:|---|
| — | loop only | 0.95 | 0.95 | 1.00 (1.00) | 1.06 | 1.04 | 0.98 (0.98–0.99) |
| `ln` | probabilities and the fit's rescaled products, [1e-300, 1] | 2.46 | 3.37 | 1.37 (1.37–1.43) | 2.56 | 3.54 | 1.38 (1.38) |
| `ln` | near 1: 1 − d, d in [1e-12, 0.5] | 2.75 | 3.37 | 1.23 (1.22–1.23) | 3.01 | 3.54 | 1.18 (1.16–1.18) |
| `ln` | counts and concentrations, [1, 1e4] | 2.87 | 3.49 | 1.22 (1.21–1.22) | 2.65 | 3.53 | 1.33 (1.33) |
| `ln` | large concentrations and odds, [1e4, 1e12] | 2.45 | 3.49 | 1.42 (1.39–1.43) | 2.55 | 3.54 | 1.39 (1.38–1.39) |
| `ln` | at the floor, [2.2e-308, 1e-300] | 2.45 | 3.49 | 1.42 (1.38–1.43) | 2.55 | 3.53 | 1.38 (1.38–1.39) |
| `exp` | log-probabilities, [−700, 0] | 2.11 | 3.70 | 1.75 (1.40–1.85) | 2.94 | 3.69 | 1.26 (1.11–1.26) |
| `exp` | log-ratios, [−20, 20] | 2.11 | 4.11 | 1.95 (1.76–2.05) | 2.02 | 4.28 | 2.12 (1.85–2.13) |
| `exp` | past underflow, [−745, −700] | 2.87 | 4.83 | 1.68 (1.61–2.00) | 5.03 | 4.97 | 0.99 (0.98–1.08) |
| `powf` | 10^(−q/10), q in [0, 93] | 7.63 | 28.14 | 3.69 (3.68–3.77) | 5.58 | 27.30 | 4.89 (4.74–4.89) |
| `powf` | base in [1e-30, 1], exponent in [0.1, 10] | 7.15 | 27.63 | 3.86 (3.84–3.89) | 6.77 | 26.14 | 3.86 (3.82–3.86) |
| `log10` | probabilities, [1e-300, 1] | 3.12 | 4.60 | 1.47 (1.47) | 6.98 | 5.22 | 0.75 (0.75) |
| `ln_1p` | the small term of a two-term log-sum-exp, [1.4e-4, 1] | 4.73 | 3.34 | 0.71 (0.70–0.71) | 4.61 | 3.42 | 0.74 (0.74–0.77) |
| `exp_m1` | −d, d in [1e-12, 50] (no shipped call) | 2.79 | 3.38 | 1.21 (1.09–1.22) | 4.56 | 3.46 | 0.76 (0.69–0.76) |
| `sin` | πx, x in [1e-3, 0.5) | 2.89 | 3.66 | 1.27 (1.18–1.33) | 4.75 | 2.81 | 0.59 (0.56–0.60) |
| `powi` | base in (0, 1), exponent 0 to 300: std against a Rust loop | 14.91 | 18.06 | 1.21 (1.21–1.22) | 11.98 | 15.32 | 1.28 (1.28–1.30) |

What a hot loop pays, per call, over the four runs:

- **`ln`**: 0.6 to 1.5 ns more on macOS, 0.5 to 1.3 ns more on Linux, in every range.
- **`exp`**: 1.4 to 3.3 ns more on macOS; on Linux 0.5 to 2.9 ns more, except past underflow, where
  glibc slows down and the two are level.
- **`powf`**: 20 to 33 ns more on macOS, 11 to 27 ns more on Linux. A1 puts all three shipped
  calls in the class that runs a few times per run, so this does not reach a hot loop.
- **`log10`, `ln_1p`, `exp_m1` and `sin`** are as fast or faster through libm on Linux; on macOS
  `log10` costs 1.5 to 1.9 ns more, `exp_m1` and `sin` up to 1 ns more, and `ln_1p` is faster.

## Bits

Arguments, out of 524,288, whose outputs are not the same bits. Every differing pair, in every
column, was one unit in the last place apart (from the first run's `measure` and `compare` output).

| function | arguments | macOS std vs Linux std | libm macOS vs libm Linux | macOS std vs libm | Linux std vs libm |
|---|---|---:|---:|---:|---:|
| `ln` | [1e-300, 1] | 0 | 0 | 62 | 62 |
| `ln` | 1 − d | 152 | 0 | 2,765 | 2,791 |
| `ln` | [1, 1e4] | 131 | 0 | 7,304 | 7,369 |
| `ln` | [1e4, 1e12] | 0 | 0 | 0 | 0 |
| `ln` | [2.2e-308, 1e-300] | 0 | 0 | 0 | 0 |
| `exp` | [−700, 0] | 971 | 0 | 51,334 | 51,181 |
| `exp` | [−20, 20] | 979 | 0 | 51,136 | 51,165 |
| `exp` | [−745, −700] | 196 | 0 | 13,423 | 13,361 |
| `powf` | 10^(−q/10) | 861 | 0 | 51,461 | 51,494 |
| `powf` | general | 814 | 0 | 51,097 | 51,105 |
| `log10` | [1e-300, 1] | 1,232 | 0 | 42 | 1,222 |
| `ln_1p` | [1.4e-4, 1] | 12,119 | 0 | 12,119 | 238 |
| `exp_m1` | [−d] | 9,283 | 0 | 9,366 | 143 |
| `sin` | πx | 17,739 | 0 | 20,061 | 16,398 |
| `powi` | std against the loop | 0 | 0 | 0 | 0 |

**libm gives the same bits on macOS and Linux on every argument of every range.** The source agrees.
In libm 0.2.16, the newest release, the processor-specific code on aarch64 is `sqrt`, `fma` and
`rint`, and on x86_64 `sqrt` and `fma` — operations IEEE 754 requires to be rounded exactly
(`src/math/arch/mod.rs`). `exp` has one further variant, for x86 without SSE. `log`, `exp`, `pow`,
`log1p`, `log10` and `sin` are ports of FreeBSD's maths library, and neither they nor the helpers
they call use `fma`.

**Today macOS and Linux disagree on `exp` and `powf` about once in 540 to 2,700 arguments** (971,
979, 861 and 814 of 524,288, and 196 past underflow). `ln` disagrees at most once in 3,400; `log10`
once in 430; `ln_1p`, `exp_m1` and `sin` far more often — `ln_1p` once in 43.

**Linux's own output would move where libm and glibc round differently**: `exp` and `powf` on about
one argument in ten (one in 40 past underflow), `ln` on up to one in 70 depending on the range,
`sin` on one in 32, `log10` on one in 430, `ln_1p` and `exp_m1` on one in 2,200 and one in 3,700.

## Constants the compiler computes while building

**`powi` with literal operands is computed by the compiler, and differently from at run time.** Of
20 `powi` calls written with a literal base and exponent, 15 came out of the build different from
the same call made at run time with the operands hidden from the optimiser, by up to 26 units in
the last place. At run time std's `powi` and the Rust square-and-multiply loop agreed on every
argument, on both platforms, so the run-time value is portable; the folded one is not the loop's.
Whether a given call is folded depends on inlining and on what the optimiser can prove constant, so
a portable `powi` means writing the loop out, not calling std's. The loop costs 2.0 to 3.8 ns a
call more than std's (1.21 to 1.30 times).

**For `ln`, `exp`, `powf` and `sin` this check cannot tell.** Twenty literal calls of each gave, in
every case, the same bits as std at run time, on both builds. That is what folding with the build
machine's own library would give on a build that runs on the machine it was built on, and also what
no folding at all would give. Two of the twenty `ln`, `exp` and `powf` values differ from libm (and
one `sin` value on macOS), so writing a constant out as bits and switching its function to libm are
not the same thing either. The caller's constant sites are listed in A1 §2.

The macOS-built and Linux-built binaries produced the same bits for 99 of the 100 literal calls. The
one that differs, `sin(0.77)`, equals each platform's own std at run time — the platform difference
this report is about, whether or not the compiler folded it.

## Limits

- **aarch64 only**, as said above.
- **Isolated calls.** The caller may inline libm, which std's calls cannot be, or pay for call
  stubs libm avoids. Step A3 measures the commands.
- **Uniform and log-uniform samples** of the ranges; real arguments cluster, and the disagreement
  rate in a cluster can differ.

## Deviations from the plan

- The ranges were first written before A1 finished and revised when it reported; the review then
  asked for the extremes A1 §3 leaves unbounded. The first runs' sawtooth order and single warm-up
  pass were replaced (review F5, F6), and the loop-only row added (F8). All figures are from the
  runs after that.
- The build-time check was added, first for `powi` (plan step B1 depends on whether `powi` is
  bit-stable) and then for `ln`, `exp`, `powf` and `sin` (review F3).

Sources: `tmp/float_A2v2/{macos,linux}_{1,2,3,4}.tsv` and `tmp/float_A2v2/compare.tsv` for the speed and bit tables; `tmp/float_A2v2/{macos,linux}_final.tsv` and `compare_final.tsv` for the literal calls, rerun after three literals were changed to satisfy clippy (the bit counts of that rerun equal the first run's).

## Review applied

| finding | what was done |
|---|---|
| F1 `sin` and `log10` do not move rarely on Linux | stated with their rates |
| F2 the `depth_bins.rs` example cannot move output | dropped from the report |
| F3 folding may apply to `ln`/`exp`/`pow`/`sin` | checked for all four; the check cannot tell, and says so |
| F4 ranges not all from A1 | sources stated per function; extremes added |
| F5 sawtooth argument order | replaced by a permutation; all timings rerun |
| F6 first row timed cold | 20 warm-up passes; rerun |
| F7 `powi` loop cost compared two copies of one algorithm | with the scrambled order the gap fell from 2.4–2.6× to 1.21–1.30×; the earlier gap was the order, not the loop |
| F8 no loop-only row | added; ratios described as understating |
| N1 distance measure | measured in the smaller value's gap; one-finite pairs counted infinite; `compare` now reports distances |
| N2, N3, N4 figures | superseded by the rerun |
| N5 `rint` and the x87 condition | corrected |
| N6 "all twelve rows" included `powi` | removed |
| N7 per-run class wording | "a few times per run" |
| N8 compiler flags differ | stated under the machines |
| N9 folded-value distance 26, not 18 | corrected |
| N10 unexplained terms | replaced or explained |
| N11 stale files | the new runs are in their own directory |
