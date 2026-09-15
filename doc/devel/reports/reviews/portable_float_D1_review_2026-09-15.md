# Code Review: portable_float_D1
**Date:** 2026-09-15
**Reviewer:** rust-code-review skill (orchestrator, single pass — see §1)
**Scope:** uncommitted diff of step D1 (`doc/devel/implementation_plans/portable_float.md`, Milestone D): `crate::float::exp` computed by a Rust port of Arm's table-driven `exp` (via musl's `src/math/exp.c`), its table generator, an accuracy scorer, a timing and output example, the re-pinned `exp(13.5)` bits and the re-recorded fit checksum
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** `git diff` in the worktree `pop_var_caller-portable-float` on top of
  `caca7b62` (`src/float.rs`, `src/cli/cross_platform_digests.rs`), plus the untracked
  `src/float/table_exp.rs`, `scripts/float_exp_table.py`, `scripts/float_exp_accuracy.py` and
  `examples/float_exp_check.rs`, and the measurement outputs under `tmp/exp_D1/`.
- **Deliberately out of scope:** `target-bench/`, `target-container-bench/`, `tmp/inline_exp/`
  (not part of this step); step D2's measurements.
- **Categories:** one pass, no sub-agent fan-out, because `cargo` and builds were forbidden while a
  measurement runs on this machine. Covered: the port's correctness against musl's C, the
  `float_portability` checklist, licensing, test strength, the measurement tools, and doc accuracy
  under `clear-technical-writing`. Nothing was compiled; every check of the port's arithmetic was
  done by an independent Python transliteration of musl's C and by the Rust port's own outputs
  already on disk (§3).

## 2. Verdict

**Approve-with-changes.** The port is correct. Every constant has the bits of musl's hex float;
all 256 table words match a recomputation at 300 bits with a different rounding method; and a
Python transliteration of musl's `exp.c`, written for this review from the C and run on IEEE
doubles, gives the same bits as the Rust port on every one of the 5,242,880 arguments the example
wrote (524,288 per range, five ranges, both platforms). The Rust port's outputs are byte-identical
on macOS and Linux, and they are within one step of the correctly rounded value everywhere scored.
The re-pinned `exp(13.5)` is the correctly rounded double (0.44 of a step below the true value;
libm's was 0.56 above).

Three things to change before committing:

- **The accuracy figures describe only the first sixteenth of each range** (M1). The example
  writes arguments in a scrambled order in which every 16th position falls in the first 1/16 of its
  range, and the scorer takes every 16th position. The "past underflow" row scored only
  [−745, −742.2], where all three implementations happen to be exact. Corrected figures are in M1;
  the conclusion survives, the numbers change.
- **The Arm copyright notice is missing** (M2). The MIT licence the code was taken under requires
  it to travel with the port.
- **The unit tests cannot see a port error smaller than one step** (M3). Replacing Arm's `C4` by
  the Taylor coefficient 1/24, zeroing a table tail, or deleting the subnormal rounding step each
  moves 46 to 711 outputs on the test grid, and every test still passes. The test's doc says such
  an error "fails here".

## 3. Execution status

- **Commands run by the reviewer** (host, read-only; Python via `uv run --no-project --with mpmath
  python`; scripts and logs in `tmp/review_D/`):
  - `git status`, `git diff`, `git log`.
  - `scripts/float_exp_table.py`, output diffed against the `TABLE` in `table_exp.rs`: identical
    (`table_script_out.txt`).
  - `tmp/review_D/verify_port.py` (logs `verify_port.log`, `verify_port_specific.log`,
    `verify_port_digest.log`):
    - the ten constants from `float.fromhex` of musl's hex floats against the Rust `from_bits`: all
      equal. `top12` of 2⁻⁵⁴, 512, 1024 and ∞ are `0x3c9`, `0x408`, `0x409`, `0x7ff`.
    - the table recomputed at 300 bits with Python's correctly rounded `int / int`: 0 of 256 words
      differ. The closest `H[k]` to a rounding tie is 0.0035 of a step from it, so 200 bits is ample.
    - a transliteration of musl's `exp.c` and `specialcase` bit-compared with the Rust outputs in
      `tmp/exp_D1/{macos,linux}/exp.{0..4}.table.f64`: 0 differences in 5,242,880.
    - edge values (see §9) and eight mutations of the algorithm run over the unit test's grid (M3).
  - `scripts/float_exp_accuracy.py tmp/exp_D1/linux 256 tmp/exp_D1/macos`
    (`accuracy_stride256_linux.txt`); `cmp` of every `arguments`, `libm` and `table` file between
    the two platforms: all identical.
  - `tmp/review_D/check_range2.py`, `debug_range2.py`: found M1. `rescore_uniform.py` re-scores
    with positions taken in argument order (`rescore_uniform_{linux,macos}_16.txt`).
  - WebFetch of musl's `src/math/exp.c`, `src/math/exp_data.c` and `COPYRIGHT`, for the licence
    headers and the error-bound comments.
- **Commands not run:** every `cargo` command (tests, clippy, fmt, doc, the example), by
  instruction. `tmp/clippy_D.log` and `tmp/doc_D.log` show clippy and rustdoc finishing without
  errors on this tree; there is no test log for D1 under `tmp/`.
- **Needs verification:**
  - the new `FITTED_PARAMETERS_MD5` `3edab375…` "the same on macOS and in the Linux container": no
    log under `tmp/` contains it (Mi4).
  - the digest `0xb1d16a18c6a42d12` proposed in M3 was computed from the Python transliteration,
    not from Rust; the test that pins it must be seen to pass on both platforms.

## 4. Open questions and assumptions

1. **Are release binaries distributed?** If so, the Arm notice (M2) belongs in a notices file
   shipped with them as well as in the source.
2. **Is the D1 timing run kept anywhere?** The example prints timings to stdout only, and no copy
   is under `tmp/exp_D1/` (Mi6). D2 is where speed is judged, but D1's example is the only
   per-call comparison of the three implementations.

## 5. Top 3 priorities

1. **M1** — fix the scorer's sampling and replace the figures before they reach the D report.
2. **M2** — add the Arm copyright and permission notice.
3. **M3** — pin the port's bits over a grid, so a sub-step error in a coefficient, a tail or the
   subnormal path fails a test.

## 6. Findings

### Major

- **M1:** [scripts/float_exp_accuracy.py:61](../../../../scripts/float_exp_accuracy.py#L61) — **[Major]** Scoring every `stride`-th position samples only the first 1/`stride` of each range
- **Confidence:** High (measured)
- **Problem:** `float_exp_check.rs`'s `arguments()` places index `i` at the value
  `reverse_bits((i · 0x9E3779B9) mod 2¹⁹)` / 2¹⁹. When `i` is a multiple of 16, the product's low
  four bits are zero, so after reversal the top four bits are zero and the argument lies in the
  first sixteenth of the range. The scorer walks `range(0, len, stride)` over file positions, so at
  stride 16 it scores only that sixteenth, and at stride 256 only the first 1/256.
- **Evidence:**
  - The stride-16 positions of range 2 ("past underflow") cover arguments −745.0 to −742.19
    (`rescore_uniform_linux_16.txt`, last line). Every result there is a deep subnormal, which is
    why `accuracy.txt` reports 0 misrounds for all three implementations in that row, while the
    same files show libm and glibc differing on 13,361 of its 524,288 arguments.
  - Of the 13,361 positions where libm and the table port differ in range 2, none is a multiple of
    16 (`debug_range2.py`).
  - Re-scored with positions taken in argument order, stride 16 (32,768 a range), Linux directory;
    "wrong" means not the correctly rounded double, and the largest miss is one step in every cell:

    | range | table port wrong | libm wrong | glibc (arm64) wrong | Apple (macOS) wrong |
    |---|---:|---:|---:|---:|
    | log-probabilities, [−700, 0] | 32 | 3,187 | 23 | 42 |
    | log-ratios, [−20, 20] | 40 | 3,152 | 23 | 67 |
    | past underflow, [−745, −700] | 3 | 862 | 3 | 15 |
    | near overflow, [0, 709.7] | 26 | 3,220 | 22 | 37 |
    | near zero, [−1e-3, 1e-3] | 11 | 10 | 11 | 69 |

    The claim to be corrected: "table misrounds 25–33 of 32,768 per range, libm 3,249–3,366, glibc
    20–33". Over the three wide ranges it is 26–40 for the port (about 1 argument in 1,000), 3,152
    to 3,220 for libm (about 1 in 10) and 22–23 for glibc. Past underflow and near zero are lower
    for everyone and should be quoted separately. The correct-rounding computation itself is sound:
    it agreed with an independent integer-arithmetic rounding on all 163,840 arguments re-scored.
- **Fix:** sort positions by argument before striding, as `rescore_uniform.py` does:
  ```python
  order = sorted(range(len(arguments)), key=lambda i: as_double(arguments[i]))
  for position in order[::stride]:
  ```
  and say in the docstring that positions are spread over the argument values. Re-run on both
  directories and quote those figures.

- **M2:** [src/float/table_exp.rs:1-21](../../../../src/float/table_exp.rs#L1-L21) — **[Major]** The port carries no copyright notice or permission notice for Arm's MIT-licensed code
- **Confidence:** High
- **Problem:** The module says the algorithm is "Arm's … (optimized-routines, 2018, released under
  the MIT licence)" and "ported line for line from musl's `src/math/exp.c`". The MIT licence's one
  condition is that "the above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software". A line-for-line port with the same constants and
  table layout is a substantial portion. musl's `exp.c` and `exp_data.c` both carry
  `Copyright (c) 2018, Arm Limited.` / `SPDX-License-Identifier: MIT`, and musl's `COPYRIGHT` lists
  "Copyright © 2017-2018 Arm Limited" for these files. The project's own `LICENSE` names only its
  authors, and nothing in the repository reproduces Arm's notice.
- **Fix:** add, at the top of `table_exp.rs`, before the module doc:
  ```rust
  // Ported from musl's src/math/exp.c and src/math/exp_data.c, which are
  // Copyright (c) 2018, Arm Limited. SPDX-License-Identifier: MIT
  //
  // Permission is hereby granted, free of charge, to any person obtaining a copy of this software
  // and associated documentation files (the "Software"), to deal in the Software without
  // restriction, including without limitation the rights to use, copy, modify, merge, publish,
  // distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
  // Software is furnished to do so, subject to the following conditions:
  //
  // The above copyright notice and this permission notice shall be included in all copies or
  // substantial portions of the Software.
  //
  // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
  // BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  // NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
  // DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  // OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
  ```
  and the same notice under a "Third-party code" heading appended to `LICENSE` (or a
  `THIRD_PARTY_NOTICES` file), so it also travels with a binary. `license = "MIT"` in `Cargo.toml`
  stays correct. Put `// Copyright (c) 2018, Arm Limited. SPDX-License-Identifier: MIT` in the
  generator script's header too, since it reproduces the table's definition and six entries.

- **M3:** [src/float/table_exp.rs:296-312](../../../../src/float/table_exp.rs#L296-L312) — **[Major]** The tests pass with a wrong coefficient, a wrong tail or the subnormal rounding removed; the doc says they fail
- **Confidence:** High (mutations run in the transliteration over the test's exact 200,001-point grid)
- **Problem:** `agrees_with_libm_to_within_one_step_everywhere_it_is_finite` says "a wrong polynomial
  coefficient or reduction constant, which errs by far more, fails here". The port and libm each
  land within one step of the truth, so the test tolerates any error that keeps the port within one
  step, and most plausible port errors do. `every_table_entry_is_two_to_the_j_over_128` checks only
  that a tail is smaller than half a step, not its value. Outside this file, `float.rs` pins three
  `exp` values (−0.7, −700.25, 13.5), none in the subnormal path.
- **Evidence** (`verify_port_specific.log`), outputs moved on the grid / moved by more than a step:

  | mutation | moved | > 1 step | caught today |
  |---|---:|---:|---|
  | `C4` = 1/24 (Taylor) | 177 | 0 | no |
  | `C4` off by 2²⁰ of its own steps | 0 | 0 | no (invisible to any output test on this grid) |
  | tails of entries 37 and 38 swapped | 1,169 | 0 | no |
  | tail of entry 37 zeroed | 711 | 0 | no |
  | subnormal `hi`/`lo` rounding removed | 46 | 0 | no |
  | high parts of entries 37 and 38 swapped | 3,117 | 3,112 | yes |
  | threshold `0x408` → `0x409` | 5,100 | 5,050 | yes |
  | threshold `0x408` → `0x407` | 0 | 0 | harmless: 256 ≤ \|x\| < 512 then takes the exact slow path |

  All 46 outputs moved by removing the subnormal rounding are wrong (not correctly rounded) where
  the port is right, for example `exp(-711.25242065)`: the port gives `0x0000eb833419deff`, the
  mutation `…df00`.
- **Fix:** keep the two tests (they catch gross errors with a readable message) and correct the
  doc sentence to "a mistyped table entry or a threshold error fails here; an error smaller than a
  step does not — the digest below catches those". Add:
  ```rust
  /// The bits over the grid, folded into one number. Computed from a transliteration of musl's
  /// exp.c run on IEEE doubles (review portable_float_D1, tmp/review_D/verify_port.py), so a
  /// coefficient, tail or subnormal-path error that stays within one step still fails here.
  #[test]
  fn the_grid_has_musls_bits() {
      let (low, high) = (-745.13, 709.78);
      let points = 200_000;
      let mut digest: u64 = 0xcbf2_9ce4_8422_2325;
      for i in 0..=points {
          let x = black_box(low + (high - low) * f64::from(i) / f64::from(points));
          digest = (digest ^ exp(x).to_bits()).wrapping_mul(0x0100_0000_01b3);
      }
      assert_eq!(digest, 0xb1d1_6a18_c6a4_2d12);
  }
  ```
  With `C4` = 1/24 the digest would be `0x9a73d4d8007071ed`, with entry 37's tail zeroed
  `0xe59df08c2913e003`, without the subnormal rounding `0x7757c8921cbcd118`. If the assertion fails
  on first run, the port and the transliteration differ somewhere, which must be found, not
  re-recorded. Also add to `keeps_its_edges` (all bits checked against correct rounding at 300
  bits):

  | argument | expected bits | covers |
  |---|---|---|
  | `-711.25242065` | `0x0000_eb83_3419_deff` | subnormal rounding step |
  | `0.20036285688060917` (≈ 37·ln2/128) | `0x3ff3_8cae_6d05_d865` | entry 37's tail (zeroed: `…866`) |
  | `f64::from_bits(1)`, `-f64::from_bits(1)` | `1.0` | subnormal arguments |
  | `512.0`, `-512.0` | `0x6e19_4765_04ba_852e`, `0x11c4_4109_edb2_0931` | first arguments of the slow path |
  | `f64::from_bits(0x4086_2e42_fefa_39ef)` | `0x7fef_ffff_ffff_ff2a` | largest argument with a finite result |
  | `f64::from_bits(0x4086_2e42_fefa_39f0)` | `f64::INFINITY` | first that overflows |
  | `f64::from_bits(0xc087_4910_d52d_3051)` | `0x0000_0000_0000_0001` | smallest argument with a non-zero result |
  | `f64::from_bits(0xc087_4910_d52d_3052)` | `0.0` (bits) | first that underflows to zero |
  | `f64::MAX`, `-f64::MAX` | `∞`, `+0.0` | the \|x\| ≥ 1024 branch with finite x |
  | a negative NaN (`-f64::NAN`) | `is_nan()` | sign of NaN does not reach the `-∞` check |

### Minor

- **Mi1:** [src/float/table_exp.rs:9-11](../../../../src/float/table_exp.rs#L9-L11) — **[Minor]** "glibc's own build does not [give the same bits]: its C source lets the compiler fuse a multiply and an add on some processors" states a mechanism nobody measured
- **Confidence:** Medium
- **Problem:** The difference is real and measured here: on arm64 Linux, glibc's `exp` differs from
  this port on 274, 404 and 375 of 524,288 arguments in the three wide ranges, and on none past
  underflow or near zero (`verify_port.log`, "table vs platform", linux). That fused multiply-add is
  the cause is plausible (Arm's comments give separate error bounds with and without it) but was
  not shown. The checklist asks for findings with sizes and without unverified mechanisms.
- **Fix:** "On arm64 Linux, glibc's `exp`, the same algorithm, gave different bits from this port
  on 274 to 404 of 524,288 arguments per range (`tmp/exp_D1/linux`). Arm's source is written to be
  built with or without fused multiply-add, which Rust never uses; that is the likely reason, not a
  measured one."

- **Mi2:** [src/float/table_exp.rs:3-7](../../../../src/float/table_exp.rs#L3-L7), [src/float.rs:47-50](../../../../src/float.rs#L47-L50) — **[Minor]** Two speed claims are wider than their source
- **Confidence:** High
- **Problem:**
  - "1.3 to 2.1 times as long a call as glibc's on the arguments the parameter fit passes it": A2's
    Linux ratios are 1.26 ([−700, 0]), 2.12 ([−20, 20]) and 0.99 past underflow. Past underflow
    they are the same speed.
  - "`exp` is the parameter fit's most-called function": A3 measured CPU samples, not calls — `exp`
    65,484 of 326,383 working samples (20%), `log` 31,534.
- **Fix:** "1.3 to 2.1 times as long a call as glibc's, except past underflow, where they are
  equal (A2)"; "`exp` took about a fifth of the parameter fit's CPU time, the most of any maths
  function (A3)".

- **Mi3:** [src/float/table_exp.rs:91](../../../../src/float/table_exp.rs#L91) — **[Minor]** "A valid scale only while −1023·N < k < 1024·N, which the check above guarantees" is false for 512 ≤ |x| < 1024
- **Confidence:** High
- **Problem:** For |x| up to 1024, `k` reaches about ±1,477·128, outside the valid range. The
  check does not exclude those arguments; it sends them to `near_the_limits`, which moves the
  exponent back into range before using the scale.
- **Fix:** "A valid scale only while −1023·N < k < 1024·N, which holds for |x| < 512; larger |x|
  goes to `near_the_limits`, which shifts the exponent first."

- **Mi4:** [src/cli/cross_platform_digests.rs:62-78](../../../../src/cli/cross_platform_digests.rs#L62-L78) — **[Minor]** The new checksum has no log, and the calls checksum's doc is now stale
- **Confidence:** High on the text; the claim itself needs verification
- **Problem:**
  - "Recorded 2026-09-15 with `float::exp` computed by table lookup, the same on macOS … and in the
    Linux container": no file under `tmp/` contains `3edab375`. The old checksum and the teeth
    measurement are cited with logs; the new one is not.
  - `CALLS_MD5`'s doc says "Recorded with the checksum above". It was recorded with `28722984…` and
    did not move with the table `exp`; that is worth saying, since it shows the change moved the
    fit's numbers and not a call.
  - Line 64 is 131 characters; the file's other doc lines keep to 100.
  - "plan Milestone D" is a plan-internal label; name the file.
- **Fix:** keep the two platforms' test logs as `tmp/digests_D1/{macos,linux}.log` and cite them;
  "Unchanged when `exp` moved to the table-driven version, which moved only the parameters file";
  rewrap; "(`doc/devel/implementation_plans/portable_float.md`, Milestone D)".

- **Mi5:** [src/float.rs:15](../../../../src/float.rs#L15), [src/float.rs:180-181](../../../../src/float.rs#L180-L181) — **[Minor]** `float.rs` docs not updated for the new `exp`
- **Confidence:** High
- **Problem:**
  - "(On 32-bit x86 without SSE libm swaps in an x87 `exp`; …)" no longer concerns this module:
    `exp` does not call libm.
  - `PINNED`'s doc says "Recorded on 2026-09-14 … both gave these bits", and the test's doc says a
    change "means libm … changed". `exp(13.5)` was re-recorded on 2026-09-15 from Rust code that is
    not libm.
- **Fix:** drop the parenthetical — A2 (line 118) names `exp` as libm's only function with an x87
  variant, and this module's `exp` no longer reaches it; "a change here means libm, the
  table-driven `exp`, or the function a name delegates to, changed"; date the `exp(13.5)` entry.

- **Mi6:** [examples/float_exp_check.rs:93-161](../../../../examples/float_exp_check.rs#L93-L161) — **[Minor]** The timings are printed but not written with the outputs
- **Confidence:** High
- **Problem:** `tmp/exp_D1/{macos,linux}/` holds the argument and output files and no timings. The
  table port's speed against libm and the platform is therefore not traceable to a log, which the
  plan requires of every figure a report quotes.
- **Fix:** write the printed table to `<output-dir>/timings.tsv` as well (or `tee` it in the
  command the report records).

- **Mi7:** [scripts/float_exp_table.py:1](../../../../scripts/float_exp_table.py#L1), [src/float/table_exp.rs:20-21](../../../../src/float/table_exp.rs#L20-L21) — **[Minor]** The generator names the wrong file, and nothing checks that the Rust table is its output
- **Confidence:** High
- **Problem:** the script's docstring says "The table `src/float/exp.rs` reads"; the file is
  `src/float/table_exp.rs`. The module doc says the script "checks it against the entries musl
  publishes": it checks six of 128. The script prints a block to paste; no step compares the
  pasted table with it (this review did, and they match).
- **Fix:** correct the path; "checks six of its 128 entries against musl's `exp_data.c`"; give the
  script a `--check` mode that reads `src/float/table_exp.rs` and exits non-zero when its `TABLE`
  differs.

### Nits

- **N1:** [src/float/table_exp.rs:16-17](../../../../src/float/table_exp.rs#L16-L17) — "so the lookup
  loses nothing": the tail is itself rounded, so the lookup is exact to about 2⁻¹⁰⁶ relative, far
  below a step. Say that.
- **N2:** [src/float/table_exp.rs:34-35](../../../../src/float/table_exp.rs#L34-L35) — `SHIFT`'s
  bits carry 2⁵¹ + round(z) in their low 52 bits, not round(z); only the low seven bits
  (round(z) mod 128) and the shift into `top` are used.
- **N3:** [src/float/table_exp.rs:132-134](../../../../src/float/table_exp.rs#L132-L134) — "for any
  integer k" holds only within −1023·128 < k < 1024·128.
- **N4:** [src/float/table_exp.rs:17-18](../../../../src/float/table_exp.rs#L17-L18) — cite where the
  0.511 bound is: musl's `exp_data.c` ("ulp error: 0.509 (0.511 without fma)") and `exp.c`
  ("Without fma the worst case error is 0.25/N ulp larger"). Add the measured rate beside it (about
  1 argument in 1,000 not correctly rounded, against about 1 in 10 for libm; M1's table).
- **N5:** [src/float.rs:43](../../../../src/float.rs#L43) — `mod table_exp;` sits between `ln` and
  `exp`; the crate declares modules at the top of a file.

## 7. Out of scope observations

- **The port versus the platform libraries on accuracy.** On the three wide ranges the port misses
  correct rounding on 26 to 40 arguments of 32,768 and glibc on 22 to 23; the two differ on 274 to
  404 of 524,288. So the port is not glibc's bits on arm64 Linux, and is about as accurate. Apple's
  library, on macOS, misses on 37 to 67.
- **A2's disagreement counts are unaffected by M1.** A2 compared outputs over all 524,288 arguments,
  not a stride.

## 8. Missing tests to add now

- The grid digest and the edge rows in M3.

## 9. What's good

- **The port matches musl operation for operation.** Checked line by line:
  - `r = x + kd * NEG_LN2_HI_N + kd * NEG_LN2_LO_N` and `tmp = tail + r + r2 * (C2 + r * C3) + r2 *
    r2 * (C4 + r * C5)` parse in Rust to the same trees as in C (left to right, `*` before `+`), and
    Rust does not fuse them.
  - `lo += 1.0 - hi + y` is `lo + ((1.0 - hi) + y)`; musl's `((1.0 - hi) + y) + lo` differs only in
    the order of the last addition's operands, and IEEE addition is commutative bit for bit (no NaN
    can occur there).
  - `abstop.wrapping_sub(0x3c9)` reproduces C's `uint32_t` wraparound; `ki & 0x8000_0000 == 0` and
    `x.to_bits() >> 63 != 0` parse as `(…) == 0` in Rust, where `&` and `>>` bind tighter than `==`.
  - `ki % N` and `ki << 45` recover `k mod 128` and `k·2⁴⁵ mod 2⁶⁴` because `z + SHIFT` lies in
    [2⁵², 2⁵³), so the double's low 52 bits are 2⁵¹ + round(z); `wrapping_add`/`wrapping_sub` with
    `1009 << 52` and `1022 << 52` (typed `u64` by the call) match C's `uint64_t` arithmetic.
  - The special cases return the same bits as musl: `1 + x` for |x| < 2⁻⁵⁴ (±0 and subnormal
    arguments give 1.0), `+0.0` for −∞, NaN or +∞ from `1 + x`, and `+0.0`/`+∞` where musl calls
    `__math_uflow(0)`/`__math_oflow(0)` (which compute `2⁻⁷⁶⁷·2⁻⁷⁶⁷` and `2⁷⁶⁹·2⁷⁶⁹`). The only
    omission is `fp_force_eval`, which raises a floating-point flag Rust does not read.
  - Twenty edge arguments, among them ±2⁻⁵⁴, ±512, the last finite and first infinite result, the
    last non-zero and first zero result, and ±`f64::MAX`, give the correctly rounded (or IEEE
    special) result in the transliteration (`verify_port.log`).
- **Every constant and table word verified independently**, as §3 lists: 10 constants, 256 table
  words at 300 bits with a different rounding method, and four pairs quoted from musl's current
  `exp_data.c` beyond the script's six found in the table.
- **The re-pinned `exp(13.5)` is right**: `…ad8b` is 0.442 of a step below e¹³·⁵ and `…ad8c` 0.558
  above; the comment's "0.44 of a step" holds.
- **Bit-identity across platforms holds on the evidence**: the `table` and `libm` output files are
  byte-identical between `tmp/exp_D1/macos` and `tmp/exp_D1/linux` for all five ranges, and the
  Rust outputs equal the transliteration's, which is an arithmetic statement independent of either
  machine.
- **The accuracy script's correct rounding is sound**: Python's `int / int` is correctly rounded
  into the subnormals, and it agreed with an independent integer rounding on every argument
  re-scored. The step count by bit difference is valid for same-sign results, subnormals included.
- **The example's timing design is sound**: the scrambled order defeats the branch predictor as in
  A2, the start side rotates each trial, and warm-up passes precede timing. Only the scorer's use of
  the order is wrong (M1).

## 10. Commands to re-verify

- `uv run --no-project --with mpmath python tmp/review_D/verify_port.py --all` (about 5 s): constants,
  table, and the Rust outputs against the transliteration.
- After M1: `uv run --no-project --with mpmath python scripts/float_exp_accuracy.py tmp/exp_D1/linux
  16 tmp/exp_D1/macos`, and the same on `macos`; the rows should match
  `tmp/review_D/rescore_uniform_{linux,macos}_16.txt`.
- After M3: `cargo test --lib float::table_exp` on both platforms; then each mutation in M3's table
  applied by hand must fail `the_grid_has_musls_bits`.
- After Mi7: `uv run --no-project --with mpmath python scripts/float_exp_table.py --check`.
- The step's gates as in plan §7, and `cargo test --lib cross_platform_digests` on both platforms
  with the logs kept (Mi4).

## Fixes applied

| finding | resolution |
|---|---|
| M1 scorer sampled the first sixteenth | `scripts/float_exp_accuracy.py` now takes positions in argument order; re-scored on both platforms' outputs (figures in the D report) |
| M2 no Arm notice | the full MIT notice with `Copyright (c) 2018, Arm Limited.` heads `src/float/table_exp.rs`, is appended to `LICENSE` under "Third-party code", and the copyright line is in the generator's docstring |
| M3 sub-step errors pass | `the_grid_has_musls_bits` pins the grid digest `0xb1d16a18c6a42d12` computed from the review's transliteration; `keeps_the_bits_at_each_paths_boundary` adds the ten bit-exact boundary cases plus `±f64::MAX` and a negative NaN; the one-step test's doc says what it does not catch |
| Mi1 fused multiply-add as cause | reworded as the likely reason, with the measured size (274 to 404 of 524,288 arguments) |
| Mi2 speed claims | "except past underflow, where the two are equal"; "about a fifth of the fit's CPU time, the most of any maths function (A3)" |
| Mi3 valid-scale comment | reworded for 512 ≤ \|x\| < 1024 |
| Mi4 digest docs | log paths `tmp/digests_D1/`, the calls checksum's history, rewrapped, the plan named by file |
| Mi5 float.rs docs | x87 remark removed, `PINNED` doc names the table-driven `exp` and the re-recording date, `mod table_exp;` moved to the top |
| Mi6 timings not saved | the example writes `timings.tsv` beside its outputs |
| Mi7 generator | path corrected, "six of 128 entries", and a `--check` mode confirming the Rust table is its output (passes) |
| N1–N4 | "exact to about 2⁻¹⁰⁶"; `SHIFT`'s bits hold `2⁵¹ + round(z)`; the `k` range in the `TABLE` doc; the 0.511 bound cited to `exp_data.c` with the measured rate beside it |
