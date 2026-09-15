# Keeping the caller's output the same on every machine

*Milestone C of [`portable_float.md`](../../implementation_plans/portable_float.md). Branch
`portable-float`. 2026-09-15.*

Milestone B sent every transcendental maths call through `crate::float`, so macOS and Linux compute
the same bits. This milestone makes that hard to lose, and in doing so found and fixed a second
way the output depended on the machine. It adds:

1. **a cross-platform test** that walks, fits and calls a small cohort and asserts the checksums of
   the parameters file and the VCF (`src/cli/cross_platform_digests.rs`);
2. **a fix to the SNP/indel fit**, whose output depended on how many threads ran it (commit
   `f04f7c43`);
3. **a macOS job in CI** running the test suite;
4. **an identity oracle that also fits**, with a committed baseline it compares against;
5. **a floating-point portability checklist** in the code-review and performance-review skills.

The owner amended items 1 and 4 at Checkpoint B: both must cover the parameter fit, not only
calling, because calls made with default parameters were already identical across platforms before
the conversion (B7 report §3.1), so a calling-only check would pass with the defect back in place.

**The machines** are those of the earlier steps: macOS is an Apple M5 Pro host (18 cores), running
native builds; Linux is the Apple `container` VM on it (arm64, glibc, 8 CPUs). **x86_64 is not
covered here**; CI's Linux runner will be the first x86_64 evidence once the branch is pushed.

## 1. The cross-platform test

`cli::cross_platform_digests` builds the two-sample, 600-base varying cohort the CLI tests share (a
SNP in each sample, a repeat tract one sample shortens by a copy, two read groups in one sample),
then in-process:

- walks it into psps (`generate-psps`);
- fits its parameters (`estimate-parameters`);
- calls it with that fit and the shipped hidden-duplication filter (`call-from-psps`);
- asserts the MD5 of the parameters file and of the VCF without its `##commandline` and
  `##reference` lines, both of which name this run rather than its calls.

A second test runs the same pipeline in pools of one, four and seven threads and requires the same
bytes. The whole module takes under two seconds.

**The checksums, and the evidence that the test would catch the defect it guards.** "Before the
conversion" is the tree at `c6a4394b`, whose maths still called the platform library, with the fit
change of §2 applied so the only difference is the maths:

| tree | platform | parameters file | VCF |
|---|---|---|---|
| the conversion (this branch) | macOS | `28722984…` | `3432b421…` |
| the conversion (this branch) | Linux | `28722984…` | `3432b421…` |
| before the conversion | macOS | `f3f66a48…` | `3432b421…` |
| before the conversion | Linux | `cbcb75af…` | `3432b421…` |

The parameters file tells the platforms apart before the conversion and not after. **The VCF does
not**: each platform's calls, made with its own differing fit, came to the same checksum. That is
the finding the owner's amendment anticipated, now shown on the fixture too. The VCF's checksum is
pinned all the same, because the calls are what the caller is for. Logs:
`tmp/digests_C/{macos,linux}_c6a4394b_final.log`.

**What else moves these checksums** is written in the module's documentation: a change to the
arithmetic of the walk, the fit or calling (re-record after measuring it on the real cohort), and a
new crate version, which is written into the parameters file's census digests and the VCF's
`##source` line.

## 2. The fit's output depended on the thread count

**What was found.** The pool-width test failed when first written. The SNP/indel fit's expectation
pass — one sweep over every census position per iteration, and most of the fit's time — split the
census into chunks whose size was `positions ÷ threads` when that was below 16,384, and joined the
chunks' floating-point totals with rayon's `reduce`, in a tree the pool chose at run time. A
floating-point sum's last bits depend on where it is split and in what order the parts are added.
Measured on the fixture with the code at `cd056a3f` (`tmp/digests_C/pool_width_cd056a3f.log`):

| threads | parameters file |
|---:|---|
| 1 | `f7e3248e…` |
| 4 | `02f469ad…` |
| 8 | `eff5fabb…` |

The real four-sample census has about two million positions, so there the chunk size was already
16,384 at any thread count, and B7's macOS and Linux fits matched while running on 18 and 8
threads. That match did not show the fit was independent of the width; only the fixture at several
widths exposed it. A small census — a targeted panel, a short genome, a test — was exposed.

**The fix** (commit `f04f7c43`). The chunk size now depends only on the census's length: 128 chunks,
none smaller than 256 positions. The iterating passes join their chunks in a tree of `rayon::join`
split at fixed points, adding the right half's totals to the left's; the pass whose per-position
lists are reported collects its chunks in position order and adds them once each, as before. The
additions happen in the same order at any width. The fit's own unit test,
`a_census_of_many_chunks_fits_to_the_same_bits_at_any_pool_width`, fits a cohort drawn at 3,172
positions at one and four threads and requires the same bits; with the old chunking and `reduce`
put back, it fails.

**Three shapes were measured, because the first two were slow.** Against the build at `cd056a3f`:

| shape | joint-fit bench, SNP/indel cases (macOS, 4 threads, same session) | tomato fit |
|---|---|---|
| fixed 16,384-position chunks, collected in order | not benched; every SNP/indel case is one or two chunks, so one or two cores | +1.6, +1.2, −0.3% (macOS), −1.0, −0.5% (Linux), interleaved |
| fixed 2,048-position chunks, computed 64 at a time and added in order | **+47 to +56%** at 5,000 positions and by sample count; +15% at 20,000 (both platforms) | +2.0, +4.7, +7.0% (macOS), −2.9, −0.6% (Linux), interleaved |
| fixed 256-position chunks, joined in a tree | −1.3 to +2.0% | one run each: **+28%**, peak memory 339 against 239 MiB (macOS) |
| **chosen**: 128 chunks of at least 256 positions, joined in a tree; the reported pass collected in order | **+0.5 to +1.7%** | **+2.8, +2.2, +2.3%** (macOS), **−1.8, −5.7%** (Linux), interleaved; peak memory 252.5 against 244.3 MiB, one run each (macOS) |

Why the 2,048-position chunks were slow on the bench's small cases is inferred, not counted. The
expectation pass measures a census by its first section's length, and at one thread the old and new
builds took the same time on the 5,000-position case (179 against 180 ms), while at four threads the
new one reached only about the speed-up that three chunks allow (70 against 50 ms). The mechanisms
for the other two slow shapes are not measured either: 7,813 chunks' setup at two million positions
for the 256-position tree, and threads waiting for each batch's slowest chunk for the batches, are
guesses.

**Output moved once, the same on both platforms.** The four tomato samples' 20-region fit, on the
psps B7 used, went from `0fb3e50d` to `9e0e7a4a` on macOS and on Linux: 7 of 574 lines, the same
kind B5's `ln(3)` constant moved — two read groups' error multipliers, inbreeding coefficients, the
two genome-wide concentrations — in about the seventh significant digit. B7's recorded fit
checksum is superseded. Timings are in `tmp/c_measure3/fit_{macos,linux}/ab.tsv`, with the host
load in `tmp/c_measure3_hostload.tsv`: 7 to 23 during the fits, which use every core, with a
PDF-preview service holding about one core throughout and falling on both builds alike.

## 3. The identity oracle now fits

`scripts/promote_ng_oracle.sh` calls four tomato samples over 20 regions from their CRAMs and from
psps with default parameters, as before, and now also fits the stored psps and calls them with the
fit. It writes seven checksums instead of five and compares them, with the slice of genome walked,
against `scripts/promote_ng_oracle.baseline`; it exits 1 on any difference or an unreadable
baseline. `PROMOTE_NG_FIT=0` skips the fit, which takes about 7 minutes on the macOS host and 10 in
the Linux container. The `PROMOTE_NG_*` settings have to be set inside the container command,
because `dev.sh` does not forward them.

**The new baseline, identical on both platforms** (`tmp/oracle_C2/{linux,macos}/digests.txt`,
built from the tree of commit `f04f7c43`; the macOS run compared itself against the Linux run's
digests and matched):

| artefact | checksum |
|---|---|
| VCF from the CRAMs, default parameters | `b74f3edf…` (unchanged since before the work) |
| VCF from the psps, default parameters | `b74f3edf…` (unchanged) |
| parameters file beside those VCFs | `5e6d874d…` (unchanged) |
| window rows | `25d388c8…` (unchanged) |
| per-sample coverage histograms | `259baf04…` (unchanged) |
| **fitted parameters file** | **`502c3d5e…`** (new) |
| **VCF called with that fit** | **`b7cbf433…`** (new) |

The oracle fits the psps it has just written, not B7's, so its fitted file is a different file
from §2's `9e0e7a4a` for the same cohort. **What the fit fix did to calls**, on B7's psps in the
Linux container, each build calling with its own fit: against the Milestone B build, 13 of 6,735
records differ — AF in 4, PARALOG_POST in 8 and QUAL in 1 (104.0 against 104.1); no genotype,
FILTER or site (`tmp/c_calls/`).

## 4. CI and the skills

**CI.** A `macos` job on `macos-15` (arm64) runs `cargo test --lib --tests --all-features`, the
suite that holds the cross-platform test; the Linux job already runs it on x86_64. `--examples` is
left out on macOS: three example tests fail on `main` on every platform measured (B7 report §4), so
the Linux job's `--examples` step is likely failing on `main` too; that is not changed here. **The
job has not run**: nothing on this branch has been pushed.

**Skills.** `ai/skills/rust-code-review/code_review/float_portability.md` and
`ai/skills/rust-performance-review/performance_review/float_portability.md` hold this work's
lessons as rules with the measurements behind them: std's transcendental methods and constants
folded from them, parallel float totals whose split or join depends on the pool, hash-map order
into float sums, reordering licences such as `f64::algebraic_add`, tests against recorded answers,
near-ties, portability tests that must digest full-precision numbers and be shown to fail without
the fix, and — for performance proposals — stating the effect on output bits and measuring small
differences interleaved. Both skills' triage tables make the category required when floating-point
code reaches output, and `hot_loops.md`'s advice to reach for `algebraic_*` first is restricted to
values that never reach output.

## 5. Tests

`cargo test --lib --tests --all-features` on the tree committed with this report: **Linux container
4,799 passed, 0 failed, 4 ignored; macOS 4,798 passed, 0 failed, 4 ignored.** In the container,
`cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings` and
`cargo doc --no-deps --all-features` with warnings as errors passed. Three tests were added: the two
in `cli::cross_platform_digests` and the fit's pool-width test.

## 6. What was not covered

- **x86_64.** No x86_64 machine was run; CI's Linux job will be the first once pushed.
- **The macOS CI job itself**, for the same reason.
- **Thread counts above eight on Linux and eighteen on macOS**, and censuses between the fixture and
  the tomato cohort, except as the fit's unit test and the bench sweep them.
- **Other parallel float totals.** The review of this milestone checked every parallel site in
  `parameter_estimation`, `calling`, `paralog` and `paralog_filter` and found all collect in order
  (review §7); nothing outside those modules was searched.

## 7. Deviations

- **The fit fix was not in the plan.** The pool-width test found it; it is its own commit before
  this milestone's, and plan §8 records it.
- **The oracle script is extended**, where plan principle 5 said it would be reused unchanged.
- **The review of this milestone** (`doc/devel/reports/reviews/portable_float_C_review_2026-09-15.md`)
  found the first fix's 16,384-position chunks left small censuses on one core, which led to the
  shapes measured in §2.
