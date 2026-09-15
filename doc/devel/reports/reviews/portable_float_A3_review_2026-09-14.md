# Review: the portable-float A3 baseline, and the build that routes the maths through libm

*Reviews `doc/devel/reports/implementations/portable_float_A3_caller_baseline_2026-09-14.md` and
`scripts/portable_float_baseline.sh` (step A3 of `doc/devel/implementation_plans/portable_float.md`),
both uncommitted in the worktree `pop_var_caller-portable-float`. 2026-09-14. No build or run was
made and neither file was edited. The raw data read are under this worktree's `tmp/` and the
throwaway worktree `pop_var_caller-portable-float-probe/tmp/`. The reviewer's one scratch file is
`tmp/review_A3_bench_parsed.tsv`, every criterion median re-read from the 20 baseline logs.*

In this review, **the libm build** means the throwaway build whose source defines the C symbols
`exp`, `log`, `pow`, `log10`, `log1p` and `sin` itself, each calling the `libm` crate. **A psp** is
the per-sample stored pileup file that `generate-psps` writes and `estimate-parameters` and
`call-from-psps` read.

## Verdict

**The report can carry Checkpoint A's decision. No finding is a Blocker.** Every timing, memory
figure, bench percentage, digest, line count and record count in it matches the raw data, with two
exceptions: one profile total, and three statements about which significant digit differs. Its
two main results hold:

- The fit slows by about a third on both platforms.
- The fitted parameters become the same file on macOS and Linux.

A cross-check the report does not make supports the fit's slowdown. In the macOS profile, `exp` took
17.4% of working samples, `log` 8.4% and the call stubs 3.9%. A2's per-call costs, with the empty
loop taken out, put libm's `exp` at about 2.5 times the system library's and `log` at 1.6 times. The
libm build has no stubs. That predicts 17.4 × 2.5 + 8.4 × 1.6 − 29.7 ≈ +27% for the whole fit,
against the +31% measured from the means.

There are ten Fix findings:

- Two concern which binary was measured.
- One is a noisy bench pass behind the report's largest number, +56%.
- Two make the baseline impossible for Milestone B to repeat from what will be committed.
- Five correct or complete statements.

There are twelve Nits.

## What was checked, and held

1. **§2.1, all 42 ranges** (7 commands × time and memory × 2 platforms) against the seven
   `runs.tsv` / `fit.tsv` files. All match, rounding half up. "Repeats spread by at most 5.1%":
   the largest is Linux `estimate-parameters`, 469.60 against 446.72 s. The spread is measured as
   the largest run over the smallest.
2. **§3.2's calling-command table, all 12 libm ranges**, against the probe's two `runs.tsv`. All
   match, and every libm range does overlap its baseline range.
3. **§3.2's fit figures.** macOS 367.06 / 374.72 / 377.90 s and Linux 585.67 / 591.02 / 586.98 s
   match. Peak memory 237.2–239.9 and 193.5–203.7 MB match. The §1 ranges of +28 to +35% and +25 to
   +36% are the smallest and largest ratios of any libm run to any baseline run: 367.06 / 287.64 =
   1.276, and 377.90 / 280.12 = 1.349. On Linux they are 585.67 / 469.60 = 1.247 and
   591.02 / 434.57 = 1.360.
4. **The benches.** I re-read every criterion median from the 20 baseline logs, skipping the
   `change:` lines, which also start with `time:`. All 180 values equal `benches_wide.tsv` up to
   unit conversion. That is 45 cases × 2 passes × 2 platforms, and the four apparent mismatches
   are float formatting. The 36 values in `probe_benches.tsv` equal the four probe logs.
5. **§3.2's bench table, all 72 cells.** Each libm change is taken from the mean of the two
   baseline passes. Each run-to-run figure is |pass 1 − pass 2| ÷ their mean. All 36 percentages
   and all 36 spreads round to the printed values. So do §1's summary ranges: +13 to +30 / +9 to
   +44 for the joint fit, and +1 to +28 / +2 to +56 for site quality.
6. **"Within 3.9% except 5.0% and 5.5%" over all 45 cases** holds under the same spread measure.
   The largest three are 5.49% (macOS `walk_one_in_100/deep_280_reads`), 4.99% (macOS
   `walk_heads/deep_280_reads`) and 3.85% (Linux `region_grain/100000`).
7. **The fit profile.** Sum every line of "Sort by top of stack". Drop everything in
   `libsystem_kernel.dylib` (thread waits, spin yields and locks) and rayon's three idle functions
   (`steal`, `wait_until_cold`, `try_advance`). That leaves **326,383**, exactly the report's
   figure. Maths: 56,814 + 27,517 + 8,670 + 4,017 = **97,018**. `one_position`: **224,566**. All
   three hold.
8. **The parameter files.** Each platform's three fits are byte-identical. The macOS and Linux
   baseline files differ in 7 of 574 lines. The Linux fit of Linux's psps is byte-identical to its
   fit of macOS's psps. All six libm fits have md5 `3f5e5e7fb3892f7619471f159657ecac`. The Linux
   baseline file and the Linux libm file differ in 57 lines.
9. **The calls.** In `tmp/calls_by_fit/`, both VCFs have 6,735 records at the same sites. 46
   differ: AF 19, PARALOG_POST 23, PARALOG_LR 3, QUAL 1, with no genotype, GQ, DP or AD
   difference. Every difference is the last printed digit moving by one, for example
   PARALOG_POST 0.161475 against 0.161476, and QUAL 104.1 against 104.0. In the probe's
   `calls_libm/`, `std.records` against `libm.records` gives 39: AF 16, PARALOG_POST 20,
   PARALOG_LR 2, QUAL 1, and no genotype. `std.records` against `libm_callstd_fit.records` is
   byte-identical. The binary `calls_libm/pop_var_caller_std` is byte-identical to the unchanged
   container binary.
10. **The oracle digests.** `tmp/oracle_A0`, `tmp/oracle_A_macos_std` and the probe's
    `oracle_probe_{macos,linux}` all list the brief's five digests, over 7,787 records. The
    `.comparable` files are the whole VCF minus `##commandline`
    (`scripts/promote_ng_oracle.sh:53`), so the comparison covers every printed field.
11. **The substitution is in force in the binaries on disk.** macOS `nm -m` on the probe's
    `pop_var_caller` and both bench binaries shows `_exp`, `_log`, `_log10` and `_pow` defined in
    `__TEXT`, and none imported. `nm` on the Linux ELF binaries shows `T exp`, `T log`, `T log10`,
    `T log1p`, `T pow` and `T sin`. On Linux, `exp` is at 0x6a0878 and `log` at 0x6a087c: the `exp`
    wrapper is one 4-byte instruction, a jump into libm's code. The remaining imports are not the
    caller's maths:
    - `_log2` / `log2@GLIBC`, called only from `parquet::arrow` writer code;
    - `___exp10` in the macOS bench binaries, called only from `plotters`, criterion's plotting
      library;
    - `cos` in the joint-fit bench.
12. **The psps carry the same content, and this can be shown directly.** The macOS and Linux psps
    differ in their headers: `command-line` (`target/` against `target-container/`, 10 bytes) and
    `created`. After shifting by those 10 bytes, `cmp -l -i 20000:20010` finds 66–68 differing
    bytes a file. They sit 13 bytes apart and each is exactly 10 larger on Linux, which makes
    them file offsets. No payload byte differs.
13. **The source.** `git diff 75722336 HEAD -- src Cargo.toml benches` is empty in both worktrees,
    and the probe's diff is the `libm_interpose` module, `unsafe_code = "deny"` and the
    `PROBE_CHECK` bit check.

## Fix

**F1. Some "unchanged binary" and libm runs used a bench build, not the release build the
baseline used.** Cargo writes `target*/release/pop_var_caller` again when it builds the benches,
using the `bench` profile and, for the joint-fit bench, `--features bench-fixtures`.
`find_binary` then takes it because it is newest, and the chain scripts `touch` it to make sure.
Evidence from `.fingerprint/*/bin-pop_var_caller.json`, file sizes and hard links:

- Base worktree, Linux. `target-container/release/pop_var_caller` is hard-linked to
  `deps/pop_var_caller-1acd7d4e…`, whose fingerprint lists `bench-fixtures` under a profile
  hash different from the 21:10 release build `a290c6d4…`. Its mtime is 00:18:24, the moment
  `fit_chain.sh` touched it. So the **434.57 s "same psps" Linux fit**, the `calls_by_fit` VCFs
  and `calls_libm/pop_var_caller_std` ran on the bench build. The 446.7–469.6 s runs used the
  release build.
- Base worktree, macOS. `target/release/pop_var_caller` is 10,711,472 bytes. That is the size of
  `deps/pop_var_caller-ac982f72…`, built 22:24:47 during the macOS bench pass with
  `bench-fixtures`, not the 21:12 release build (10,713,264 bytes) that the macOS baseline runs
  used.
- Probe worktree, macOS. `ac982f72…` (bench profile, `bench-fixtures`) was built at 00:34:16, 39 s
  before `probe_speed_chain.sh` started the macOS calling runs. Those runs most likely timed that
  build, against the baseline's release build.

The codegen difference is probably small: debug info plus a feature that adds code behind
`cfg`. But the report treats all of these as one build. §1 and §3.2 use 434.6 s as a baseline,
and it is the only run below 446 s. Fix:

- Say which build each comparison used.
- For Milestone B, have the script build with `cargo build --release` itself, or refuse a binary
  older than the source.
- Write the binary's checksum and `git rev-parse HEAD` into `runs.tsv` and `fit.tsv`.

**F2. The Linux site-quality libm pass is disturbed, and it holds the report's largest number.**
In `probe_linux_benches/ng_site_quality_perf.1.log`:

| case | criterion interval | outliers | baseline interval |
|---|---|---:|---|
| `alleles/2` | not listed | 18 in 100 | not listed |
| `alleles/4` | [1.2202 1.2976 1.4517] ms, −6% / +12% around the median | 14 in 100 | 0.05% wide |
| `alleles/6` | [1.9930 2.1114 2.2187] ms, ±5% | 12 in 100 | 0.05% wide |

§3.2's "run-to-run" column describes only the baseline passes, so a reader has no way to see that
the one libm pass behind "+13%" and "+56%" was this noisy. Fix: rerun the Linux site-quality
bench for the libm build, or print the libm pass's interval beside those rows and treat +56% as
"+47 to +64%".

**F3. The baseline Milestone B must repeat lives only in gitignored `tmp/`.** `tmp/` is ignored
(`.gitignore:14`). The report prints 18 of the 45 bench cases. It also omits the exact
invocations, which exist only in `tmp/run_A3_chain.sh`:

- which psps each platform's fit read (`$OUT/{macos,linux}_runs/work/psps/*.psp`, the last
  repeat's);
- that `regions20.bed` is `head -n 20 benchmarks/tomato2/regions_n160_200kb.bed`;
- the order runs → fit → benches;
- the separate `cargo bench --no-run -j 2` prebuild in the container.

Once `tmp/` is cleaned, Milestone B has neither the numbers nor the commands. Fix: put all 45
cases in an appendix table, or commit the tsv, and put the invocations in the report or in the
script's usage.

**F4. The script does not hold everything Milestone B needs to repeat it.**

- `benches` builds inside `cargo bench` with default parallelism. The plan (§7) says the
  container build takes `-j 2`, and the chain did that as a separate step the script lacks.
- `runs` deletes `work/` at the start of each repeat. The psps a later `fit` reads are therefore
  whatever the last `runs` left, and nothing records that.
- Nothing records the machine, the commit or the toolchain.
- Running `benches` before `runs` or `fit` silently changes which binary is timed (F1).

Fix:

- Take `-j` for the build.
- Record `git rev-parse HEAD`, `rustc -V` and the binary checksum.
- Say in the usage text that `runs` and `fit` must follow a release build.

**F5. §1's fit percentages are envelopes that mix in a warm-up run and a different build.** On
Linux, the "+25%" lower end comes from the first baseline run, 469.60 s. The other two runs took
446.72 and 448.21 s. The "+36%" upper end comes from the 434.57 s run on the bench build (F1). The
means say it plainly:

| platform | baseline mean | libm mean | change |
|---|---:|---:|---:|
| macOS | 284.4 s | 373.2 s | **+31%** |
| Linux | 454.8 s | 587.9 s | **+29%** |

Against Linux's second and third baseline runs alone, the change is +31%. Fix: lead with the
means, and give the ranges as noise.

**F6. The deviation that justifies the 160-region runs is not applied to the libm build.** §5
says the 160 regions were added because "fixed start-up cost would hide a small slowdown at 20
regions". The libm build's calling commands were timed only at 20 regions, and §1 concludes "no
change beyond repeat-to-repeat spread" from them. The conclusion is still very likely right,
because the `call-from-alignments` profile puts 55 of about 34,500 working samples in the
maths library. Fix: either run the 160-region set for the libm build, or rest the "no change"
row on the profile (maths is about 1 sample in 630, so no measurable change is expected) rather
than on overlapping 20-region ranges.

**F7. Three statements about significant digits are wrong.** §3.3: "The four samples'
inbreeding coefficients differ from the eighth significant digit on … and two concentrations
from the ninth and tenth."

| value | macOS | Linux | first differing significant digit |
|---|---|---|---:|
| inbreeding coefficient, SRS1839219 | 0.941173069… | 0.941173079… | 8th |
| inbreeding coefficient, SRS1839231 | 0.9283287869… | 0.9283288017… | **7th** |
| inbreeding coefficient, SRS1839234 | 0.9863939635… | 0.9863939723… | 8th |
| inbreeding coefficient, SRS1839262 | 0.99501638282… | 0.99501638822… | **9th** |
| `reference_concentration` | 1.892485070… | 1.892485063… | 9th |
| `alternative_concentration_total` | 0.000826137472… | 0.000826137466… | **8th** |

Fix: "from the seventh to the ninth significant digit", and "the eighth and ninth".

**F8. "57 lines move" gives no subject and no size.** The 57 lines hold 259 differing numbers:

- 36 slippage rows (`share_of_reads_that_slip`, `shorter_share`, `fall_off` and the curves
  behind them);
- 14 per-length concentration rows;
- 3 inbreeding coefficients, 2 read groups' error multipliers, and the 2 concentration totals.

The two platforms agree on every slippage and per-length concentration row today; their 7
differing lines are all in the last group. The largest relative move among values
above 1e-4 is 0.3% (0.000107729 against 0.000108046). Most moves are near 1 part in 10 million.
The single largest relative move is 9%, on a share of 2.7e-9. Fix: say this in §3.3. It is the
information that says whether 57 lines matter. Also fill §1's macOS cell, which reads "—": the
macOS baseline file and the libm file differ in 58 lines.

**F9. §3.3's statement of which psps each comparison read is wrong for two of the three call
rows.** "All fits below read the same four psp files, written on macOS" is literally false for the
Linux fit used in `calls_by_fit`. `call_with_fits.sh` reads `tmp/baseline_A3/linux_fit20/`, which
was fitted from Linux's psps, and its output is byte-identical to the fit of the macOS psps. The
calls table then says "On Linux, calling the same psps", but:

- row 1 (`calls_by_fit`) calls **Linux's** psps (`P=tmp/baseline_A3/linux_runs/work/psps`);
- rows 2–3 (`calls_libm/run.sh`) call **macOS's** psps (`tmp/probe_fit`, byte-identical to
  `macos_runs/work/psps`).

The results are unaffected, because check 12 shows the psp payloads are identical. Fix: state
which psps each row read, and replace the indirect argument ("the same fit from both, so the same
content") with the direct one from check 12.

**F10. The substitution check leaves no record.** §3.1's "all of 200,000 arguments" and the `nm`
listing are not in any file under either `tmp/`. `grep -rl differing` finds nothing, and the plan
(§6) requires every figure to trace to a log. The report also says "`nm` shows the six symbols
defined in the binary". The macOS `pop_var_caller` now on disk exports five: `_log1p` is neither
exported nor imported, so fat LTO resolved it internally. Fix: save the `PROBE_CHECK` output and
the `nm` lines for each binary that was timed, and say "none imported" rather than "six defined".

## Nit

1. **"Overstates one cost … so its slowdowns are, if anything, larger than a real conversion's."**
   The direction is right, but the size is one jump a call: on Linux the wrapper is a single
   4-byte instruction (check 11). The main thing the libm build removes is the dynamic-library
   stub, 3.9% of the fit's samples on macOS, and a real conversion removes it too. A real
   conversion could also differ the other way. Calling `libm::exp` as ordinary Rust makes it
   visible to the inliner, which can grow or shrink hot loops. Suggest: "should match a real
   conversion to within one jump a call; inlining could move it either way". Also cut "and
   understates none that matter", which has no evidence behind it.
2. **Two things the libm build does that a real conversion would not.** Neither changes a
   conclusion, but the method paragraph should say both:
   - It also redirects `exp` and `log` calls from dependencies and C code linked into the binary.
   - Calls whose arguments the compiler proves constant are still folded at build time with the
     build machine's own library, which A2's "Constants the compiler computes" covers. A real
     conversion through `libm::` Rust code would not do that.
3. **The joint-fit bench fixtures are drawn with the substituted functions.** `fit.rs:3392–3416`
   and `ssr_fit.rs:2712–2736` use `ln`, `exp` and `powf` in gamma, normal and Poisson draws, so the
   libm build fits slightly different synthetic data. Integer draws flip only when a uniform
   falls within one unit in the last place of a threshold, so the data are almost surely the same
   counts. Say so, or hash the fixture once.
4. **The call-from-alignments profile total is not reproducible from a stated rule.** The rule
   that gives 326,383 for the fit (drop `libsystem_kernel` and rayon's idle functions) gives
   34,494 for `cfa.sample.txt`, not 34,316. No natural subset of wait functions gives 34,316.
   The maths count also changes rule between the two runs:
   - the fit's 97,018 includes the `DYLD-STUB$$` samples;
   - the call run's 46 leaves out `DYLD-STUB$$log`, 9 samples.

   Counted with stubs, it is 55 samples, about 1 in 630 rather than 1 in 750. State the
   exclusion rule and apply it to both.
5. **The profile durations do not match the files.** "`call-from-alignments`, 160 regions, 19 s"
   is Linux's time. The macOS run took 21.6 s, and the profile's main thread holds 13,574 samples.
   "30 s of the run" for the fit does not match its main thread's 19,425 samples. Say how long
   `sample` ran.
6. **Two spread measures under one word.** §2.1's 5.1% is largest ÷ smallest − 1. §2.2 and §3.2 use
   |difference| ÷ mean. Name the measure once.
7. **§4 "the baseline passes agreed within 3.6%"** is true of the 18 cases in §3.2's table, but §2.2
   says 3.9% with two exceptions at 5.5%. Say "within 3.6% on the 18 cases above".
8. **§3.2 Linux memory "against 193–196"** includes the one 192.9 MB run on the bench build.
   §2.1 says 194–196.
9. **"MB"** is mebibytes: the script divides bytes by 1,048,576 and kilobytes by 1,024. Say MiB, or
   say megabytes means 2²⁰ bytes.
10. **Script mechanics:**
    - On macOS the command's stderr goes to `$log.time`, so the parse `/ real /` could match a
      program line. The usage comment says the command's output goes to the log.
    - `[ a -nt b ]` is not in POSIX `test`, although dash, bash and zsh accept it.
    - The rest is sound. `result=$(measure …)` with `exit 1` in the function aborts under `set -e`,
      because an assignment takes its command substitution's status. `ru_maxrss` is bytes from
      `/usr/bin/time -l` on macOS and kB from `getrusage` on Linux, and both are converted.
      `RUSAGE_CHILDREN` sees only the one waited child.
11. **Clarity.** Define on first use:
    - *psp*;
    - *identity oracle* (the five-digest check);
    - *working samples*;
    - *call stubs*, the jump table a program uses to reach a function in a shared library;
    - *repeat-tract stratum* against *all strata* (`strata/stands_alone`).

    "The substitution is in force" should say what was compared. The §1 table's "—" cells read as
    "no change" rather than "not measured".
12. **§4, GIAB.** "Milestone B's plan runs it if the output changes there": the plan's condition is
    that the oracle VCF moves (B(n+2)). The fit moves here while the oracle does not, so say which
    output decides.
