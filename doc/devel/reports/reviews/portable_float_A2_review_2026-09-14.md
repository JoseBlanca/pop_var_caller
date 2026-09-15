# Review: the portable-float A2 comparison of std against libm

*Reviews `examples/float_libm_vs_std.rs` and
`doc/devel/reports/implementations/portable_float_A2_libm_vs_std_2026-09-14.md` (step A2 of
`doc/devel/implementation_plans/portable_float.md`), both uncommitted in the worktree
`pop_var_caller-portable-float`. 2026-09-14. No build was run and neither file was edited. The
reviewer's one scratch file is `tmp/float_A2_review/bitdist.awk`.*

## Verdict

**The report can carry the Checkpoint A decision. No finding is a Blocker.** Every number in its
two tables matches the raw logs. The central result, that libm gives the same bits on macOS and
Linux, holds. So does the claim that every difference is exactly one step between adjacent
floating-point numbers: I checked that independently of the program.

There are eight Fix findings. Four correct statements in the prose. One of those is wrong: `sin`
does not move "rarely" on Linux, but on one argument in 32. The build-time `powi` hazard is argued
from two `depth_bins.rs` calls whose outputs cannot move. One call is not reachable, and both
round to integers far from a rounding boundary. The folding finding is also not specific to
`powi`. The ranges are not all from the A1 inventory, and the ranges do not cover all of the
inventory. The other four Fix findings are in the program. The argument "shuffle" is a sawtooth.
The first timed row on macOS was measured before the machine settled. The `powi` timing compares
two copies of one algorithm that differ 2.5-fold in speed, and nothing explains why. No empty
pass was timed, so the slowdown percentages understate the cost. There are eleven Nits.

## What was checked, and held

1. **Table 1, all 72 cells**: the 24 std and libm times and 12 ratios per platform against
   `tmp/float_A2_{macos,linux}.tsv`, and each "(4 runs)" range as the minimum and maximum over
   the final run and `run2`–`run4`. Every one matches. Computed with
   `paste <(awk -F'\t' 'NR>2 && NF==7{print $1"\t"$5}' tmp/float_A2_$p.tsv) <(… run2) <(… run3) <(… run4)`.
2. **Table 2, all 48 cells**, against the "arguments differing" column of the two platform logs
   and `tmp/float_A2_compare.tsv`. All match.
3. **The repeats used the same arguments and code.** Each of the 25 output files in
   `tmp/float_A2/{macos,linux}` is byte-identical to its counterpart in all six
   `tmp/float_A2_noise/*` directories (`cmp -s`, no difference printed).
4. **"Every such difference was one unit in the last place" is true, for all three pairings.**
   The program checks only std against libm on each platform, so I checked all three pairings
   from the output files. I dumped each file as 64-bit words (`od -An -v -t x8`) and counted
   the integer distance between the bit patterns of each differing pair
   (`tmp/float_A2_review/bitdist.awk`). Every differing pair, in all 12 ranges, for macOS std
   against Linux std and for std against libm on each platform, has the same sign and is exactly
   1 step apart. Example output: `sin.0 mac-vs-lin std: differing=17739 signdiff=0 steps1:17739`.
5. **The prose frequencies**: 524,288 / 12,119 = 43.3; / 152 = 3,449; / 51,181 = 10.2;
   / 51,494 = 10.2; / 62 = 8,456; / 7,369 = 71. Speed: `ln` 1.00–1.16 macOS and 1.11–1.26 Linux
   over four runs; `exp` 1.53–1.62 and 1.18–2.10; `powf` 3.29–4.10, 9.58–10.87 ns extra; `powi`
   2.43–2.59, 3.88–3.96 ns extra; `log10` +47% and `exp_m1` +20% on macOS. The folded `powi`
   figures: 15 of 20 differ, `0x…1e23 − 0x…1e11` = 18, and 0 of 20 differ between the builds.
6. **The arguments are the same on both platforms.** The index permutation is integer
   arithmetic. Every argument is built from `+`, `-`, `*`, `/`, which IEEE 754 rounds exactly,
   plus `trunc` and `floor`, which are exact, and `libm::exp`/`libm::log`. Rust never fuses a
   multiply and an add into one instruction. The libm outputs agreeing on all 6.3 million
   arguments confirms it: different arguments would have given different libm outputs.
7. **The timing loop.** Every pass is warmed once untimed, so page faults on the output vector
   land outside the timing. The order alternates. With 21 trials, `timings[10]` is the true
   median (`float_libm_vs_std.rs:256-262`). The std calls take `black_box`ed arguments, so they
   cannot be folded. `x.ln()` and the others compile to calls into the platform library.
8. **The run-time `powi` comparison is not folded** (`black_box` on both operands,
   `float_libm_vs_std.rs:406-407`). **The folded side did fold**: 15 of its values differ from
   both the run-time value and the Rust loop, by up to 26 steps.
9. **`compare`** reads the same file names for both directories. It asserts equal lengths and
   compares raw `u64` bits, so it counts what its header says.
10. **libm source** (`~/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/libm-0.2.16`, the
    version in `Cargo.lock:905-906`). `exp.rs`, `log.rs`, `pow.rs`, `log1p.rs`, `log10.rs`,
    `sin.rs` and `expm1.rs` each open with `origin: FreeBSD /usr/src/lib/msun/…`. Neither they
    nor the helpers `sin` uses (`k_sin`, `k_cos`, `rem_pio2`, `rem_pio2_large`, `scalbn`)
    contain `fma` or `mul_add`. `pow` calls `sqrt` (`pow.rs:184`), which is exact. The only
    `select_implementation!` among them is `exp`'s x87 branch (`exp.rs:86-90`).
11. **Inventory ties.** `powf` has 3 shipped calls, all per-run (A1 Summary table, `powf` row).
    No shipped call uses `exp_m1` (A1 §1.1, "None of `log2`, …, `exp_m1` … occurs").

## Fix

### F1. `sin` and `log10` do not move "rarely" on Linux

Report, "What the numbers say", third paragraph: "`exp` and `powf` on about one argument in ten,
`ln` on one in 8,500 to one in 70 depending on the range, and the rest rarely."

Linux std ≠ libm, from `tmp/float_A2_linux.tsv`:

| function | differing | one in |
|---|---:|---:|
| `sin` | 16,398 | 32 |
| `log10` | 1,222 | 429 |
| `ln_1p` | 238 | 2,203 |
| `exp_m1` | 143 | 3,666 |

`sin` moves more often than `ln` does in any range. List `sin` with its number, and say "the rest"
only of `ln_1p` and `exp_m1`.

### F2. The `depth_bins.rs` hazard is argued from calls whose output cannot move

Report, `powi` section, third paragraph: "A1 found two shipped calls (`depth_bins.rs:241, 299`)
whose base is a literal and whose exponent is a loop counter over a literal range … A change in
optimisation can then move their bits."

Three things are wrong or missing:

- **`:241` is not reachable from the binary** (A1 §2.7 and §4: its only caller is
  `DepthBinEdges::new`, used only by `Default`, tests and `examples/`). "Two shipped calls" is
  A1's word for it, but a hazard to the binary is one call, `:299`, reached through `for_census`
  (`run/gatherer.rs:311`).
- **The base is not a literal.** It is `widening_ratio(8, 124, 11)`, itself a `powf` of literals
  (`src/parameter_estimation/depth_bins.rs:232-233, 238, 293-297`). Folding the `powi` requires folding
  that `powf` first.
- **The value is rounded to an integer at once** (`… * ratio.powi(step as i32)).round() as u32`,
  `depth_bins.rs:241, 299`), and no rung is near a rounding boundary. For the census ladder
  124·r^k, k = 1–10, the closest fraction to .5 is k = 7, 709.423, which is 0.077 away. For the
  unreachable ladder 8·r^k the closest is 0.152 away (awk,
  `r=exp(log(124/8)/11); v=124*r^k`). A 26-step move at 709 is about 3e-12. So neither call's
  output can change.

The recommendation to write the loop out may still be right for a general `float` module, but
these two sites are not evidence for it. Say that one reachable site has this shape and that its
rounding makes it immune. Keep the hazard as a statement about any future `powi` whose operands
the compiler can see.

### F3. Build-time folding is not a `powi` property

The section heading "`powi` is not safe as it stands, for a different reason" and its body treat
folding as specific to `powi`. The report's own explanation is that LLVM folds `powi` with the
compiler's host `pow`. LLVM's constant folder (`ConstantFolding.cpp`) folds calls to
`llvm.pow`, `llvm.exp`, `llvm.log`, `llvm.sin` and the other intrinsics the same way, whenever
their operands are constant: it evaluates them with the host C library. I did not verify this on
this toolchain, because cargo was off limits.

If it holds, it matters for Checkpoint A in two ways:

- **Every std transcendental call with literal operands carries the build machine's rounding, not
  the run machine's.** A1 marks such calls "literal" (for example `depth_bins.rs:233`, the
  `powf` above, and the `LazyLock` constants in `ssr_best_path_flat_gap.rs:201-225`). A Linux
  binary built on macOS would carry macOS's bits in those constants. The plan's "pin the
  constants" option depends on this.
- **libm is immune.** It is Rust code, so folding it means evaluating IEEE basic operations,
  which LLVM does exactly.

One extra line in `check_folded_powi` per function would settle it: `f64::exp(-3.7)` folded
against `black_box(-3.7).exp()`.

### F4. The ranges are not all A1's, and do not cover all of A1's

The report says twice that its ranges come from A1 §3 ("The ranges come from the step A1
inventory's §3", and in Limits, "The ranges are the hot calls', from A1 §3"). But §3 has sections
only for `ln`, `exp`, `ln_1p`, `powi`, `sin` and `sqrt`. The `powf` ranges (10^(−q/10), and base
in [1e-30, 1]), the `log10` range and the `exp_m1` range are the author's own choice. Those
functions have no hot call, so the choice is reasonable, but the text should say so.

The ranges also leave out parts of §3, and the disagreement counts vary with the range: `ln`
disagrees on 0, 131 and 152 arguments in its three ranges. The ranges miss:

- **`ln` above 1e4.** `paralog/prior.rs:311` takes `p/(1−p)` up to 1e12. Rising products at
  `dirichlet_multinomial.rs:136, 292` grow as α^ploidy, and α reaches about ploidy × samples (A1 §3.1), so
  about 1e8 at 5,000 diploid samples. This is the large-cohort end that `CLAUDE.md` requires.
- **`ln` in [2.2e-308, 1e-300]** (`emission.rs:558-559`, floored at `f64::MIN_POSITIVE`).
- **`exp` below −708**, where results become subnormal, meaning smaller than any normal `f64`,
  with fewer significant bits. A1 §3.2 says `ssr_emission.rs:629` and the log-sum-exp sites are
  "not bounded below". Libraries treat that band specially, and the sampled [−700, 0] stops short
  of it.

Either widen the ranges, or state per function that the counts describe the sampled ranges only.

### F5. The argument order is a sawtooth, not a shuffle

`float_libm_vs_std.rs:237-241`: "shuffled by a multiplicative step so a branch predictor sees no
monotone sweep", using `index.wrapping_mul(0x9E37_79B9) % 2^19`.

`0x9E3779B9 mod 2^19` = 489,913. Successive arguments therefore step by 489,913/524,288 of the
range, which is −34,375/524,288, or −6.6%. They wrap about every 15 calls: index 0, 1, 2, 3 give
0, 0.934, 0.869, 0.803. The branch predictor sees a regular descending pattern with a period of
about 15. The golden-ratio multiplier spreads points only when you take the top bits, not the low
ones: `((index as u32).wrapping_mul(0x9E37_79B9) >> 13) as f64 / 2^19` keeps a permutation. A
modulus with 324,019 (odd, about 2^19/φ) also works.

Two related effects. `powf`'s second range and `powi` derive the exponent from the same `x` as
the base (`(x * 7919.0).fract()`). So the exponent is a fixed function of the base, and the
exponents follow the same period-15 pattern.

What this does to the timings is not known. An implementation with branches that depend on the
argument, such as `exp`'s range thresholds or `powi`'s bit loop, gains from a predictable
pattern, and the two sides do not branch alike. The counts of differing arguments are unaffected,
because the argument set is still every point on the grid. Fix the comment at least, and
preferably the permutation. The A2 timings then need a rerun.

### F6. The first macOS row was timed before the machine settled

`ln` over [1e-300, 1] is the first thing `measure` times. In the final macOS run its medians are
2.23 ns (std) and 2.58 ns (libm), with passes ranging over 66% and 55% of the median. In the three
repeats the same row gives 1.63–1.84 and 1.89–2.06, with a spread of at most 25%. Table 1 quotes
the final run, so the most-called function's macOS times are about 37% above what the other runs
show. The ratio survived (1.12–1.16), as the report notes. The single warm-up pass, about 1 ms,
is too short for an Apple Silicon core to reach full speed or for the scheduler to settle the
thread.

Quote a repeat's times for this row, or warm up for a fixed time (say 200 ms) before the first
range.

### F7. The `powi` speed gap is between two copies of one algorithm

Report: "The loop written in Rust is 2.4 to 2.6 times slower than std's, about 4 nanoseconds a
call. That cost would only be paid if `powi` were replaced." The same section then recommends
replacing it.

The report itself shows that std's run-time `powi` is the same square-and-multiply loop, bit for
bit on 524,288 arguments. Two implementations of one algorithm that differ 2.5-fold in time
point to how this loop was compiled into the pass, not to the algorithm. The likely cause is a
branch on each exponent bit where the library version is branch-free, or the reverse, together
with F5's regular exponent pattern. So the figure is not what B1's loop will cost. State that the
gap is unexplained, or time a branch-free variant (`product *= if remaining & 1 != 0 { factor }
else { 1.0 }`) before quoting a cost for B1.

### F8. No empty pass was timed, so the slowdown percentages understate the cost

Every iteration of a pass does two `black_box` stack round-trips, loop bookkeeping and a write,
identical on both sides (`float_libm_vs_std.rs:47-50`). The program's own comment
(`:38-40`) says a cost common to both sides "pull[s] their ratio towards one". That applies to
this overhead as much as to the indirect call it avoided. At 1.5–2 ns a call it can be a large
share. The nanosecond differences in the report are unaffected. The percentages ("`ln` costs 0
to 16% more", "`exp` 53 to 62%") are lower bounds.

Time `side!(|x, _| x)` beside each range and subtract it, or say in the report that the
percentages include this overhead.

## Nits

**N1. The program does not check the "one unit" claim for the platform comparison, and its
distance measure is loose.** `compare` counts differing values but never measures their size.
Point 4 above shows the claim holds, but the report cites sources that do not show it.
`ulps_apart` (`:266-270`) divides by the gap above the *larger* magnitude. So a one-step
difference across a power of two reads 0.5, and a two-step difference there reads 1. It also
silently skips pairs where one value is finite and the other is not (`:345`), and it reports
+0 against −0 as 0. Counting the steps between the bit patterns of same-sign values, as
`bitdist.awk` does, is exact. Put that in both `measure` and `compare`.

**N2. "about one argument in 540 to 650".** 524,288 / 979 = 536 and 524,288 / 814 = 644. Write
"540 to 640", or "about one in 600".

**N3. "0.4 to 1 extra nanosecond"** for `exp` is from the final run alone. The "110%" beside it
comes from Linux runs 3 and 4 of [−20, 20], where the extra time is 1.32 and 1.30 ns (2.52 − 1.20
and 2.50 − 1.20). Say "0.4 to 1.3".

**N4. The spread sentence.** "Under 10% for most rows" holds for 18 of 24 sides on macOS and 19
of 24 on Linux in the final run. By row it is 8 of 12 and 7 of 12. "The worst is the macOS `ln`
row" is true of the final run. Across all four runs the worst is macOS run 2, `powf`
10^(−q/10) on libm, where passes ran from 14.68 to 24.92 ns (69%). Say "in the final run".

**N5. libm's processor-specific code.** Report: "the only processor-specific code on aarch64 and
x86_64 is `sqrt`, `fma` and `rint`". In `src/math/arch/mod.rs`, x86_64 takes the `sse2` branch:
`sqrt`, `sqrtf`, `fma`, `fmaf`. `rint` is on aarch64 only. "`exp` has one other variant, for
32-bit x86 without SSE2": the condition is no SSE at all (`configure.rs:116`,
`!target_features.any(== "sse")`). "0.2.16, the newest release" cannot be checked offline; it is
the version the project locks. None of this changes the conclusion.

**N6. "for all twelve rows".** The `powi` row compares std against a Rust loop, not against
libm. Say "all eleven libm rows, and the loop".

**N7. "Its three shipped calls run once per run."** A1's per-run class is "during setup, or a
small fixed number of times per run" (A1 §2, class definitions). Say "a few times per run".

**N8. The two builds used different compiler flags.** `.cargo/config.toml` sets
`target-cpu=apple-m1` for macOS aarch64 and nothing for Linux aarch64. So the cross-platform
speed comparison (macOS std `ln` 2.23 ns against Linux 1.49 ns) mixes libraries and code
generation. It does not affect the bits. Add it to "The machines".

**N9. `check_folded_powi` cannot tell a build that did not fold.** In an unoptimised build the
"folded" side runs at run time too, and the check prints "0 of 20 differ", which reads as "folding
is harmless". Print a warning when every pair agrees. Also, the largest folded difference is 26
steps (`1.000001^250`: `0x3ff001062d3840cd` against `…40e7`), not the 18 the report quotes as its
example. Quoting the maximum is more useful.

**N10. Clarity** (`ai/skills/clear-technical-writing/SKILL.md`, Rules 3 and 8).

- "msun" is not explained. Say "FreeBSD's maths library (msun)".
- "The second half says how much Linux's output would move once": "once" is ambiguous. Say
  "…how much Linux's output would change, a one-time change, on switching."
- "call stubs", "unroll" and "constant propagation" are compiler jargon the reader, a
  geneticist, needs explained or replaced. For example: "the compiler can see that the base and
  every exponent are known while building, and compute the results then".
- The report defines "units in the last place" in the Results preamble. The `powi` section
  reuses it well, but write "steps between adjacent floating-point values" once for a reader
  who does not know the idiom.

**N11. Stale files beside the final outputs.** `tmp/float_A2/{macos,linux}` still hold
`cos.0.*`, `ln_1p.1.*` and `powf_as_exp_y_ln_x_.*` from the 21:17 run, before A1 changed the
ranges. `compare` ignores them, but anyone reusing the directory may not. Also,
`target-container/release/examples/float_libm_vs_std` is dated 00:34, after the Linux repeats
(00:27). The binary behind those runs is gone, although identical outputs (point 3) show it
computed the same thing.
