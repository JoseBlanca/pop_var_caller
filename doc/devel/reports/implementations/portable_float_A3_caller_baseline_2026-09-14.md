# The caller's speed before the change, and what routing its maths through libm costs it

*Step A3 of [`portable_float.md`](../../implementation_plans/portable_float.md). Branch
`portable-float`; the caller's source is identical to `75722336`. Measured 2026-09-13 and 2026-09-14.*

This report records two things:

1. **the caller's speed and memory before any change**, taken with commands Milestone B can repeat;
2. **what the change would cost and what it would move**, measured on the caller itself, not
   estimated from step A2's per-function timings. This part goes beyond the plan.

Part 2 used a throwaway build, the **libm build**, which was never committed. Its source differs
from this tree in one respect: it defines the C functions `exp`, `log`, `pow`, `log10`, `log1p` and
`sin` itself, each calling the libm crate. std's `f64::exp` and the rest call those C functions by
name. The linker takes a function defined inside the program before one in the operating system's
library, so every such call in the binary runs libm's code, and not one call site is edited.

**Two terms.** A *psp* is one sample's stored evidence, written by `generate-psps`. The *identity
oracle* is `scripts/promote_ng_oracle.sh`: it calls four samples over 20 regions both ways, straight
from the CRAMs and through psps, and prints checksums of the two VCFs, the parameters file and the
window-coverage files.

**The machines.** macOS is the host, an Apple M5 Pro with 18 cores and 64 GB, running native builds.
Linux is the Apple `container` VM on the same host: arm64 with glibc, 8 CPUs and 16 GB. They differ
in core count and compiler flags, so compare each platform only with itself. **Linux on x86_64 is not
measured.**

**The data** is the identity oracle's four tomato accessions (SRR5079906, SRR5079864, SRR5079878,
SRR5079876), at about three reads a position. They are called over the oracle's 20 regions (4 Mb) and
over all 160 regions of `benchmarks/tomato2/regions_n160_200kb.bed` (32 Mb). That is one corner of
the caller's range: few samples, low depth.

## 1. What it costs and what it moves

Times are the mean of three repeats. The ranges are in §2 and §3.

| | macOS | Linux |
|---|---|---|
| `estimate-parameters`, the 20-region psps | 284 s → 373 s, **+31%** | 454 s → 588 s, **+30%** |
| joint-fit bench, 11 cases | **+13 to +30%** | **+9 to +44%** |
| site-quality bench, 7 cases | **+1 to +28%** | **+2 to +15%**, and +13% and +56% in two cases from a noisy pass (§3.2) |
| `call-from-psps`, 160 regions | +0.2 to +2% | +3% |
| `call-from-alignments`, 160 regions | ranges overlap | 0 to +3% |
| `generate-psps`, 160 regions | ranges overlap | ranges overlap |
| identity oracle checksums | unchanged | unchanged |
| fitted parameters, 574 lines | 58 lines move | 57 lines move; the largest change to a value above 1e-4 is 0.29% |
| VCF called with each build's own fit | — | 39 of 6,735 records differ in a printed decimal; no site, no genotype |

**Today, with no change**, the macOS and Linux fits of the same psps differ in 7 of 574 lines, and
calling one set of psps with each changes 46 of 6,735 records, also in printed decimals only. The
libm build's fits on macOS and on Linux are the same file, byte for byte.

## 2. The baseline

The commands are in [`scripts/portable_float_baseline.sh`](../../../../scripts/portable_float_baseline.sh):
`runs` for the three calling commands, `fit` for `estimate-parameters`, and `benches` for the five
criterion benches. Appendix A lists every invocation. `call-from-alignments` and `call-from-psps` take
`--threads 0`, which means every core. `estimate-parameters` has no thread flag and uses every core;
`generate-psps` has none and walks its samples one after another. Peak memory is the operating
system's high-water mark of resident memory for the finished process, in MiB. Logs are under
`tmp/baseline_A3/`.

### 2.1 The commands

Range over the repeats. The 160-region rows are from the second session (`tmp/baseline_A3/redo/`),
run beside the libm build's; the first session's 160-region ranges overlapped them on every command.

| command | data | repeats | macOS s | macOS MiB | Linux s | Linux MiB |
|---|---|---:|---|---|---|---|
| `call-from-alignments` | 20 regions | 5 | 4.64–4.71 | 939–1,030 | 4.78–4.88 | 875–942 |
| `generate-psps` | 20 regions | 5 | 7.42–7.60 | 487–540 | 8.51–8.65 | 428–447 |
| `call-from-psps` | 20 regions | 5 | 2.58–2.68 | 135–142 | 3.32–3.37 | 96–109 |
| `call-from-alignments` | 160 regions | 3 | 21.55–21.79 | 1,339–1,447 | 19.10–19.31 | 890–941 |
| `generate-psps` | 160 regions | 3 | 34.94–35.24 | 531–592 | 39.52–40.47 | 406–439 |
| `call-from-psps` | 160 regions | 3 | 6.43–6.51 | 162–165 | 7.96–7.97 | 123–132 |
| `estimate-parameters` | the 20-region psps | 3 | 280.1–287.6 | 248–249 | 447.9–459.2 | 192–195 |

The memory ranges of the 160-region rows are from the first session. Repeats of one command spread
by at most 3.8% in time, the macOS `call-from-psps` over 20 regions. The fit takes about a hundred
times as long as calling the same psps (284 s against 2.6 s on macOS), so its cost dominates any
psp-mode run that fits its own parameters.

The macOS fit read the psps written on macOS; the Linux fit read the same files. Linux's fit of its
own psps (446.7–469.6 s, first session) wrote the same parameters file, so the two sets of psps
carry the same content. The review confirmed that directly: they differ only in the header, which
records the command line and the time of writing.

### 2.2 The benches

Criterion's median of each case, from two passes back to back, in
`tmp/baseline_A3/benches_wide.tsv`. The two passes agreed to within 3.9% on 43 of the 45 cases. The
exceptions are `ng_psp_perf`'s `walk_heads/deep_280_reads` and `walk_one_in_100/deep_280_reads` on
macOS, at 5.0% and 5.5% apart. The joint-fit and site-quality cases are in §3.2.

### 2.3 Where the maths time is

The macOS `sample` profiler took one sample a millisecond of every thread. *Working samples* leave
out threads blocked in a lock or condition wait or spinning for work: rows naming `psynch`,
`swtch_pri`, `wait_until_cold`, `Stealer`, `try_advance`, `mutex`, `cvsignal` or `cvbroad`. Calls to
`exp` and `log` go through a *call stub*, a jump the dynamic linker fills in, and the stub's samples
are counted with the function.

- **`call-from-alignments`, 160 regions** (19 s requested): the maths library took 55 of 34,316
  working samples (`log` 30 plus 9 in its stub, `exp` 16), about 1 in 620.
- **`estimate-parameters`, 20-region psps** (30 s requested, from 2 s into the run): `exp` and `log`
  took 97,018 of 326,383 working samples (`exp` 56,814 plus 8,670, `log` 27,517 plus 4,017), about
  30%. Of the rest, 224,566 were in one function, the joint fit's per-position loop
  (`parameter_estimation::joint::fit::one_position`).

The source files are `tmp/profile_A3/cfa.sample.txt` and `fit20.sample.txt`.

## 3. The libm build

### 3.1 Checking that it measures what it claims

- **In force, and logged** (`tmp/baseline_A3/redo/probe_{macos,linux}_check.txt`). Inside the libm
  build, std's `exp` and `ln` returned libm's bits on all 200,000 arguments tried, on both platforms.
  A2 found the two libraries disagree on about one `exp` argument in ten.
- **Not imported.** `nm` shows none of the six functions imported from the system library. They are
  defined in the binary. On macOS `log1p` has no symbol of its own, having been resolved inside the
  binary.
- **The bench binary.** It links the same library; `exp` and `log` are defined in it.

**How far it stands for the real change.** Three differences, none large:

1. The real change calls libm directly, and the compiler can inline it. The libm build reaches it
   through one extra jump per call. Inlining can move speed either way, so the build's costs are
   close to the real change's but not a bound on them.
2. It also redirects maths calls made by the dependencies, which the real change would not touch.
3. Constants the compiler computes while building still use the build machine's library (A2).

The fit timings used a version with the substitution in `src/main.rs`. The command timings and
benches used one with it in the library, which the bench binaries also link.

### 3.2 Speed

**`estimate-parameters`** read the 20-region psps written on macOS, three repeats each.

| | macOS unchanged | macOS libm | Linux unchanged | Linux libm |
|---|---|---|---|---|
| seconds | 280.1, 285.3, 287.6 | 367.1, 374.7, 377.9 | 447.9, 459.2, 454.6 | 585.7, 591.0, 587.0 |
| peak MiB | 248–249 | 237–240 | 192–195 | 194–204 |

The means are +31% on macOS (284.4 to 373.2 s) and +30% on Linux (453.9 to 587.9 s). The smallest
libm time is 28% above the largest unchanged one on both platforms.

A2 predicts a similar figure. The profile's `exp` share of 20% at 1.4 to 3.3 ns more a call, and its
`log` share of 10% at 0.6 to 1.5 ns more, give about +27% on macOS. That estimate is the review's.

**The calling commands.** libm build against unchanged binary, same session, seconds:

| command | regions | repeats | macOS libm | macOS unchanged | Linux libm | Linux unchanged |
|---|---:|---:|---|---|---|---|
| `call-from-alignments` | 20 | 5 | 4.60–4.72 | 4.64–4.71 | 4.78–4.87 | 4.78–4.88 |
| `generate-psps` | 20 | 5 | 7.50–7.75 | 7.42–7.60 | 8.44–8.53 | 8.51–8.65 |
| `call-from-psps` | 20 | 5 | 2.62–2.69 | 2.58–2.68 | 3.33–3.37 | 3.32–3.37 |
| `call-from-alignments` | 160 | 3 | 21.67–21.94 | 21.55–21.79 | 19.41–19.60 | 19.10–19.31 |
| `generate-psps` | 160 | 3 | 34.84–35.34 | 34.94–35.24 | 39.36–40.29 | 39.52–40.47 |
| `call-from-psps` | 160 | 3 | 6.52–6.57 | 6.43–6.51 | 8.20–8.22 | 7.96–7.97 |

At 20 regions every libm range overlaps the unchanged one. At 160 regions, `call-from-psps` is
slower on both platforms: by up to 2% on macOS, where the ranges miss each other by 0.01 s, and by 3% on Linux,
by 0.23 s. `call-from-alignments` on Linux is slower by up to 3%. The other ranges overlap.
The 20-region rows for the unchanged binary are from the first session. The libm rows are in the
libm build's worktree, `tmp/redo_{macos,linux}_runs{,160}/`.

**Benches.** One libm pass against the unchanged binary's two. The change is measured from the mean
of the two unchanged passes; *run-to-run* is how far apart those two passes were.

| case | macOS | run-to-run | Linux | run-to-run |
|---|---:|---:|---:|---:|
| joint fit, repeat-tract stratum, 32 tracts | +18% | 0.6% | +12% | 2.9% |
| joint fit, repeat-tract stratum, 128 tracts | +15% | 2.7% | +11% | 1.0% |
| joint fit, repeat-tract stratum, 1 sample | +21% | 0.1% | +15% | 1.6% |
| joint fit, repeat-tract stratum, 8 samples | +19% | 0.1% | +13% | 2.6% |
| joint fit, repeat-tract stratum, 32 samples | +13% | 0.2% | +9% | 1.8% |
| joint fit, all strata | +21% | 0.2% | +14% | 0.6% |
| joint fit, SNP/indel, 5,000 positions | +23% | 1.8% | +28% | 1.1% |
| joint fit, SNP/indel, 20,000 positions | +22% | 1.8% | +25% | 1.9% |
| joint fit, SNP/indel, 1 sample | +30% | 3.6% | +44% | 0.4% |
| joint fit, SNP/indel, 8 samples | +23% | 1.6% | +27% | 3.1% |
| joint fit, SNP/indel, 32 samples | +15% | 1.8% | +15% | 0.6% |
| site quality, 63 samples | +21% | 1.1% | +15% | 1.2% |
| site quality, 200 samples | +10% | 0.7% | +7% | 1.0% |
| site quality, 1,000 samples | +3% | 0.2% | +2% | 0.4% |
| site quality, 3,000 samples | +1% | 0.0% | +2% | 0.9% |
| site quality, 2 alleles | +3% | 0.2% | +3% | 0.6% |
| site quality, 4 alleles | +8% | 1.2% | +13% | 0.7% |
| site quality, 6 alleles | +28% | 1.3% | +56% | 0.8% |

The data are `tmp/baseline_A3/probe_benches.tsv` against `benches_wide.tsv`.

**Two Linux site-quality numbers are unreliable.** The libm pass's own confidence intervals for 4 and
6 alleles are 1.22–1.45 ms and 1.99–2.22 ms. Criterion flagged 12 to 18 of 100 samples as outliers.
The unchanged passes' intervals were a few hundredths of a percent wide. The +13% and +56% need a
second pass before they are quoted. The joint-fit bench generates its test data with the same
functions, so the libm build benchmarked slightly different data. The bench draws integer counts,
which the review judged almost certainly unchanged.

### 3.3 Output

**Calling with default parameters does not move.** On both platforms the libm build reproduced the
identity oracle's five baseline checksums, and so did the unchanged binary run natively on macOS. The
oracle fits nothing, so it cannot see a change in the fit.

**The fit moves.** All four fits below read the psps written on macOS.

- **Today, unchanged binary, macOS against Linux: 7 of 574 lines differ.**
  - one read group's error-probability multiplier: 4.910629719888638 against 4.910632131808087;
  - the four samples' inbreeding coefficients, from the 7th to 9th significant digit — for example
    0.941173069553551 against 0.9411730799630529;
  - the two genome-wide concentrations, from the 9th and 8th significant digit.

  Each platform reproduces its own file exactly across repeats.
- **libm build, macOS against Linux:** the same file (checksum `3f5e5e7f`) on all six runs.
- **Linux, unchanged binary against libm build: 57 lines differ.** Of those, 36 are rows of the
  repeat-tract slippage model and 14 are per-repeat-length concentrations, where the two platforms
  agree today. The other 7 are the lines above. Among the numbers above 1e-4 that moved, the largest
  change is 0.29%. On macOS, 58 lines differ.

**What that does to calls.** On Linux, with the unchanged release binary (SHA-256 `07a44a40…`) or
the libm build (`b0b151ac…`):

| comparison | psps | records | differ | fields that differ | genotypes that differ |
|---|---|---:|---:|---|---:|
| today: unchanged binary with the macOS fit against it with the Linux fit | written on Linux | 6,735 | 46 | AF 19, PARALOG_POST 23, PARALOG_LR 3, QUAL 1 | 0 |
| libm build with its fit against unchanged binary with its fit | written on macOS | 6,735 | 39 | AF 16, PARALOG_POST 20, PARALOG_LR 2, QUAL 1 | 0 |
| libm build against unchanged binary, both with the unchanged fit | written on macOS | 6,735 | 0 | — | 0 |

AF is the allele frequency written to the VCF. PARALOG_LR and PARALOG_POST are the hidden-duplication
filter's likelihood ratio and posterior. Every differing field differs in its last printed digit, and
no site appeared or disappeared. The first row's fits came from each platform's own psps, which carry
the same content (§2.1).

The data are `tmp/baseline_A3/redo/calls_by_fit/` and, in the libm build's worktree,
`tmp/calls_libm/`.

## 4. Limits

- **One corner of the range.** Four samples at about 3×, over 20 regions. A decimal that moves here
  could move a genotype or a filter decision near a threshold in a larger or deeper cohort. The GIAB
  trio was not run.
- **aarch64 only.** CI and the Linux dev box are x86_64.
- **One libm pass per bench.** See the two Linux site-quality cases above.

## 5. Deviations from the plan, and a correction

- **`estimate-parameters` was added to the baseline.** The profile put about 30% of its CPU in `exp`
  and `log`.
- **160-region runs were added**, because fixed start-up time hides a small slowdown at 20 regions.
- **Everything in §3 is beyond the plan's Milestone A.** It changes nothing shipped.
- **Correction, found by review.** `cargo bench` rebuilds the package's binary, with the bench
  profile, into the same directory the release build uses, and the script's binary finder then chose
  it. That affected four sets of measurements, all taken again from fresh release builds:
  - a Linux fit;
  - both output comparisons;
  - the libm build's calling runs.

  The outputs came back byte-identical; every figure above is from the release builds. The script now
  builds benches into their own directory and writes each run's binary checksum to `build.txt`.

## Appendix A. The invocations

`B=/Users/jose/devel/pop_var_caller/benchmarks`, `REF=/Users/jose/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa`,
`CAT=$B/tomato1/crams/S_lycopersicum_chromosomes.4.00.repeats.parquet`, `CRAMS` the four CRAMs above
in the order given, `regions20.bed` the first 20 lines of `regions_n160_200kb.bed`.

```sh
# release builds, before any measurement
./scripts/dev.sh cargo build --release --bin pop_var_caller          # Linux
cargo build --release --bin pop_var_caller                           # macOS

# calling commands
DEV_EXTRA_MOUNT=$B ./scripts/dev.sh scripts/portable_float_baseline.sh runs <out> 5 $REF $CAT tmp/regions20.bed $CRAMS
DEV_EXTRA_MOUNT=$B ./scripts/dev.sh scripts/portable_float_baseline.sh runs <out> 3 $REF $CAT benchmarks/tomato2/regions_n160_200kb.bed $CRAMS
./scripts/portable_float_baseline.sh runs <out> 5 $REF $CAT tmp/regions20.bed $CRAMS                # macOS, same arguments

# the fit, on the psps the macOS 20-region runs left in <out>/work/psps
./scripts/portable_float_baseline.sh fit <out> 3 $REF $CAT tmp/baseline_A3/macos_runs/work/psps/*.psp
DEV_EXTRA_MOUNT=$B ./scripts/dev.sh scripts/portable_float_baseline.sh fit <out> 3 $REF $CAT tmp/baseline_A3/macos_runs/work/psps/*.psp

# benches (into target-bench / target-container-bench)
./scripts/portable_float_baseline.sh benches <out>
./scripts/dev.sh scripts/portable_float_baseline.sh benches <out>
```

One command at a time: the Linux VM shares the host's cores, so no two measurements overlapped. The
first session's benches were built into the release directories, before the fix above. They ran
after every release-binary measurement of that session, so its figures are unaffected.
