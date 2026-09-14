# What converting the caller's maths to libm cost, and what it moved

*Step B7 of [`portable_float.md`](../../implementation_plans/portable_float.md). Branch
`portable-float`, commit `cd056a3f`. Measured 2026-09-14 and 2026-09-15.*

Steps B1–B6 sent every transcendental call in `src/` — `ln`, `exp`, `powf`, `powi`, `log10`,
`ln_1p`, `exp_m1`, `sin`, `cos` — through one module, `crate::float`, which calls the `libm` crate
instead of the operating system's maths library, and a clippy rule now refuses the std methods. This
report measures the result against the unchanged caller:

1. how much slower each command and each bench is, on macOS and on Linux;
2. what output moved;
3. whether the test suite passes on both.

**The two builds.** *Unchanged* is the caller at commit `c6a4394b`, whose source is identical to
`75722336`, in a second worktree. *The change* is `cd056a3f`. Both were built in release mode on each
platform; the binaries timed are listed with their checksums in Appendix A.

**Two terms.** A *psp* is one sample's stored evidence, written by `generate-psps`. The *identity
oracle* is `scripts/promote_ng_oracle.sh`: it calls four samples over 20 regions both ways, straight
from the CRAMs and through psps, with default parameters, and prints checksums of the two VCFs, the
parameters file and the window-coverage files.

**The machines.** macOS is the host, an Apple M5 Pro with 18 cores and 64 GB, running native builds.
Linux is the Apple `container` VM on the same host: arm64 with glibc, 8 CPUs and 16 GB. Compare each
platform only with itself. **Linux on x86_64 is not measured** (§5).

**The data** is the identity oracle's four tomato accessions (SRR5079906, SRR5079864, SRR5079878,
SRR5079876), at about three reads a position, over the oracle's 20 regions (4 Mb) or all 160 regions
of `benchmarks/tomato2/regions_n160_200kb.bed` (32 Mb). **That is one corner of the caller's range:
four samples, low depth.** The benches cover the other axes in miniature — the joint-fit bench from
one sample to 32, the site-quality bench from 63 samples to 3,000.

## 1. Summary

| | macOS | Linux |
|---|---|---|
| `estimate-parameters`, 20-region psps | 298 s → 379 s, **+27%** | 452 s → 584 s, **+29%** |
| `call-from-psps`, 160 regions | 6.44 s → 6.54 s, **+1.7%** | 8.08 s → 8.24 s, **+2.0%** |
| `call-from-alignments`, 160 regions | 21.62 s → 21.62 s, no change | 19.25 s → 19.36 s, +0.6%, ranges overlap |
| `generate-psps`, 160 regions | 34.48 s → 34.75 s, +0.8% | 40.29 s → 39.87 s, −1.0% |
| joint-fit bench, 11 cases | **+8 to +21%** (about +10% at the low end, §2.3) | **+7 to +46%** |
| site-quality bench, 7 cases | **0 to +26%** | **0 to +23%** |
| repeat-tract delimiter, pileup and psp benches, 27 cases | all within 2.6% of A3's unchanged figures | delimiter, same session: −2.3% to +3.5%; pileup and psp within 2.0% of A3's |
| peak memory, every command | no consistent direction; within 12% of A3's range, from a different session (§2.4) | the same |
| identity oracle, five checksums | unchanged | unchanged |
| calling with default parameters, 160 regions, 120,538 records | identical | identical |
| fitted parameters, 574 lines | 57 lines move; **the same file on macOS and Linux** | 56 lines move; the same file |
| VCF called with each build's own fit, 6,735 records | — | 30 differ, in a printed decimal only; no genotype, QUAL, FILTER or site |
| test suite, `--lib --tests` | 4,795 passed, 0 failed | 4,796 passed, 0 failed |

The owner accepted at Checkpoint A a fit about 30% slower and up to 3% on the calling commands. Both
limits hold on both platforms.

## 2. Speed

### 2.1 How it was measured, and why this way

An earlier interleaved attempt in this step ran while another application loaded the host. Its one
unchanged macOS fit took 496 s, against 280–288 s in A3's quiet session, which swamps a 30% effect
(`tmp/ab/{linux,macos}/ab.tsv`). It is not used here. The figures below come from one session,
2026-09-14 22:04 to 2026-09-15 01:07, with no other application doing work: the busiest other
processes were macOS's own background services, the editor (at most 114% of one core, in one
sample) and the Claude Code session writing this report.

- **The builds were interleaved.** Each round ran the unchanged binary and the change back to back,
  alternating which went first, so a change in machine speed during the session falls on both.
- **The load was logged every 20 seconds on the host** with each run's start and end time
  (`tmp/ab2/hostload.tsv`), because the Linux VM cannot see the host's load.
- **Nothing measured overlapped.** The commands and the fit ran on Linux first, then on macOS; the
  benches on macOS first, then on Linux; then `generate-psps` on Linux and macOS, and finally the
  Linux delimiter bench again (§2.3).

Loads, and what made them:

- **Linux commands and fit:** host load 2.4 to 11 once the runs had started, median 8.4. Almost all
  of it is the VM running the command on its 8 CPUs. The busiest other process at any sample took
  at most about one core of the 18: a media-analysis service and a PDF-preview service. The longest
  burst of those overlapped the second round, in which the change's fit took 579.1 s, against
  578.6 s in the first round outside it.
- **macOS commands and fit:** host load 6.9 to 22, median 18.7 — the fit itself, which uses every
  core (one sample shows the change's fit at 1,671% CPU). The same two services took up to one core
  in bursts.
- **Benches and `generate-psps`:** logged before and after each pass in
  `tmp/ab2/benches_{macos2,linux}/load.txt`, `tmp/ab2/delim_linux/load.txt` and `tmp/ab2/gen.log`;
  macOS load 1.2 to 5.4, Linux 1.6 to 6.6, the upper values being the joint-fit bench's own four
  threads and the VM. The PDF-preview service held about one core almost continuously through the
  Linux benches and `generate-psps`.

### 2.2 The commands

Seconds, range over the rounds, and the change in the mean. `call-from-alignments` and
`call-from-psps` take `--threads 0`; `estimate-parameters` uses every core.

| command | platform | rounds | unchanged | change | change in mean |
|---|---|---:|---|---|---:|
| `call-from-alignments`, 160 regions | macOS | 5 | 21.51–21.75 | 21.51–21.76 | 0.0% |
| `call-from-psps`, 160 regions | macOS | 5 | 6.42–6.46 | 6.52–6.58 | +1.7% |
| `estimate-parameters`, 20-region psps | macOS | 4 | 290.4–303.6 | 375.3–387.8 | +27.1% |
| `call-from-alignments`, 160 regions | Linux | 5 | 19.08–19.61 | 19.27–19.54 | +0.6% |
| `call-from-psps`, 160 regions | Linux | 5 | 7.98–8.20 | 8.20–8.30 | +2.0% |
| `estimate-parameters`, 20-region psps | Linux | 3 | 446.3–458.3 | 578.6–592.9 | +29.1% |
| `generate-psps`, 160 regions | macOS | 3 | 34.39–34.60 | 34.69–34.82 | +0.8% |
| `generate-psps`, 160 regions | Linux | 3 | 40.03–40.54 | 39.60–40.29 | −1.0% |

**The fit's slowest unchanged run and fastest change-build run are 24% apart on macOS and 26% on
Linux**, so the fit's cost is not noise.

**The calling commands move by 2% or less, and at that size the binary alone can move either way.**
`call-from-psps` is slower on both platforms: on macOS the fastest change-build run is 0.9% above
the slowest unchanged one, and on Linux the ranges touch at 8.20 s. `generate-psps` shows the same
pattern on macOS — all three change-build runs slower than all three unchanged, by 0.8% in the mean —
but on Linux the change was the *faster* build in each of the three rounds, by 1.0% in the mean. A build
that is slower cannot be faster in every round, so a difference of about 1% in either direction can
come from the two binaries differing at all, and not from the maths. `call-from-alignments` is 0.0% on macOS; on
Linux it is +0.6% in the mean and slower in 4 of 5 rounds, with the ranges overlapping. The
`generate-psps` rounds ran after the benches, at 00:50–00:57, at a host load of 2.3 to 5.0.

**Against A3's libm build** — the throwaway build that sent every call through libm without editing
a call site — the fit costs 27.1% here against A3's 31% on macOS, and 29.1% against 30% on Linux.
The unchanged macOS fit here (mean 298.4 s) is 5% above A3's (284.3 s), while the change's (379.1 s)
is 1.6% above the libm build's (373.2 s); the sessions differ, so the two percentages cannot be
compared closely. The plan kept this step partly to see whether the real conversion, which lets the
compiler inline libm, costs less than that build. The benches hint that it does on macOS (§2.3), but
nothing here settles it.

### 2.3 The benches

Criterion's median of each case. **The joint-fit and site-quality benches were run in this session on
both trees**, four passes in the order unchanged, change, change, unchanged. *Run-to-run* is how far
apart a tree's two passes were, as a share of their mean.

**Joint fit** (four threads; the SNP/indel cases run eight passes of the fit, the repeat-tract cases
one round from one starting point):

| case | macOS change | macOS run-to-run, unchanged / change | Linux change | Linux run-to-run, unchanged / change |
|---|---:|---|---:|---|
| repeat tracts, 32 tracts | +18.0% | 0.4% / 0.0% | +13.4% | 4.0% / 0.4% |
| repeat tracts, 128 tracts | +13.1% | 0.7% / 0.3% | +8.8% | 0.8% / 0.4% |
| repeat tracts, 1 sample | +20.8% | 0.0% / 0.2% | +12.3% | 2.6% / 0.4% |
| repeat tracts, 8 samples | +17.2% | 0.1% / 0.2% | +10.6% | 1.5% / 0.8% |
| repeat tracts, 32 samples | +8.7% | 2.5% / 0.2% | +7.2% | 1.9% / 0.1% |
| six repeat-tract strata | +20.4% | 0.1% / 0.5% | +14.4% | 0.3% / 0.2% |
| SNP/indel, 5,000 positions | +19.4% | 1.4% / 0.1% | +26.2% | 0.6% / 0.6% |
| SNP/indel, 20,000 positions | +17.3% | 3.6% / 0.2% | +25.9% | 1.0% / 0.5% |
| SNP/indel, 1 sample | +20.9% | 2.6% / 1.2% | +45.9% | 0.4% / 1.7% |
| SNP/indel, 8 samples | +13.7% | 7.5% / 0.2% | +26.6% | 0.5% / 1.6% |
| SNP/indel, 32 samples | +8.3% | 6.2% / 0.8% | +17.6% | 0.9% / 0.3% |

**Site quality** (one site's quality score, over sample count, and over allele count at 1,000
samples):

| case | macOS change | macOS run-to-run, unchanged / change | Linux change | Linux run-to-run, unchanged / change |
|---|---:|---|---:|---|
| 63 samples | +10.8% | 10.4% / 0.1% | +12.3% | 1.2% / 0.2% |
| 200 samples | +5.1% | 4.1% / 0.2% | +6.7% | 1.9% / 0.1% |
| 1,000 samples | +0.9% | 1.1% / 0.1% | +1.2% | 1.3% / 0.2% |
| 3,000 samples | +0.2% | 0.7% / 0.1% | +0.2% | 0.0% / 0.1% |
| 2 alleles | +0.9% | 1.1% / 0.2% | +1.7% | 0.4% / 0.0% |
| 4 alleles | +5.9% | 1.8% / 0.3% | +8.6% | 0.1% / 0.7% |
| 6 alleles | +25.5% | 0.9% / 0.8% | +22.6% | 0.0% / 0.0% |

What to read from these:

- **In the joint fit, 32 samples is the cheapest case of each half on both platforms**, +7 to +18%.
  The dearest is the SNP/indel fit of one sample on Linux, +46%. The SNP/indel half costs more on
  Linux than on macOS in every case, and the repeat-tract half less.
- **Site quality's cost grows with allele count and shrinks with sample count**: +23 to +26% at six
  alleles, under 2% from 1,000 samples up.
- **Three unchanged passes on macOS were noisy, and there the table understates the cost** — site
  quality at 63 samples (passes 10.4% apart), and the SNP/indel fit at 8 and 32 samples (7.5% and
  6.2%). In each, one unchanged pass sits near A3's figure for the same case and the other above it;
  the slow pass raises the unchanged mean. Measured against the pass near A3's figure, the cost is
  +16.9% (reported +10.8%), +18.1% (+13.7%) and +11.8% (+8.3%). The repeat-tract fit at 32 samples
  (passes 2.5% apart) reads +10.1% the same way. So the low end of macOS's joint-fit range is nearer
  +10% than +8%. The two change passes agree to within 1.7% in every case on both platforms.
- **Two Linux figures A3 could not trust are replaced.** A3's libm build measured site quality at
  four and six alleles as +13% and +56%, from a pass whose confidence intervals were wide. The change,
  with both trees in one session, measures +8.6% and +22.6%.
- **Against A3's libm build, the change costs no more on macOS in any of the 18 cases** — for example
  the SNP/indel fit at 8 samples, +13.7% (or +18.1% by the correction above) against +23%, and site
  quality at 63 samples, +10.8% (+16.9%) against +21%. On Linux the two agree within about 3
  percentage points, except at six alleles. A3's figures came from another session, so this hints at
  a gain from inlining without showing one.

**Do the two trees time the same work?** The B5 review expected not: the joint-fit bench draws its
simulated evidence with `ln`, `sin` and `cos`, which the change converted, and the fit's iteration
count depends on the draw. A throwaway program built in both trees, on macOS, drew every bench case's
evidence and hashed it (`tmp/ab2/b7_draw_probe.rs`, output `tmp/ab2/probe_{std,float}_macos.txt`).
**All 11 draws hash identically in the two trees.** So on macOS the conversion's rounding differences
did not change the drawn evidence; the draws on Linux were not checked. Every SNP/indel fit runs its
full eight passes in both trees — eight is the bench's cap, and none converged earlier — so the
number of passes is the same, but that says nothing about the work inside a pass. The fitted values
differ in their last bits, as expected; the repeat-tract fit's number of objective evaluations is
not exposed and was not counted.

**The other three benches** — the repeat-tract delimiter, the SNP/indel pileup and the psp reader and
writer — ran one pass of the change, compared with A3's two passes of the unchanged tree from
2026-09-14.

- **macOS: all 27 cases are within 2.6% of A3's mean.** Seven moved by more than A3's own two passes
  were apart, none of them by more than 1.4%; the three largest differences (−2.6%, −1.9%, −1.8%) are
  psp cases where A3's two passes were themselves 3.8% to 5.5% apart. Table in
  `tmp/ab2/benches_macos2/table.md`.
- **Linux: no measurable change either, but that took a second run.** Against A3's passes, the
  delimiter's five depth cases came out 1.5% to 3.7% slower where A3's two passes had been 0.3% to
  3.2% apart; the pileup and psp cases were within 2.0%. The delimiter runs once for every read that
  crosses a repeat tract, so it was run again in the same session on both trees, in the order
  unchanged, change, change, unchanged (`tmp/ab2/delim_linux/`). There its 13 cases moved between
  −2.3% and +3.5%, and every case beyond ±1.1% is one where a tree's own two passes were 2.1% to 3.4%
  apart; the depth cases that had looked slower now read −2.3%, −1.9%, +0.4%, −0.5% and +0.5%. The
  first reading was the difference between two sessions, not the change.

**Two review notes these benches answer.**

- The hidden-duplication filter adds two log-probabilities with a cubic series below a threshold and
  `ln_1p` above it (`src/paralog/locus_score.rs`, `log_add_exp`); the threshold was chosen from
  glibc's `log1p` cost, and libm's `ln_1p` is a different price (B4 review). **No bench exercises that
  filter**, so the threshold was not re-tuned. It runs inside `call-from-psps` and
  `call-from-alignments`, whose whole cost is at most 2.0% (§2.2).
- `crate::float`'s functions are marked `#[inline]`, but the libm crate's own function bodies are not,
  and they can be inlined into the caller only where link-time optimisation runs. The release profile
  sets `lto = "fat"` and `codegen-units = 1`; the `profiling` and `soak` profiles set `lto = false`,
  so a profile taken with either should be expected to show more time in libm than a release build
  spends there (B3 review); that was not measured. Every figure in this report is from the release or
  bench profile.

### 2.4 Peak memory

Not re-measured in the interleaved session. The change's release binaries were run earlier in this
step with `scripts/portable_float_baseline.sh` (`tmp/measure_B7/{linux,macos}_{runs,runs160,fit20}/`),
and their peak resident memory shows **no consistent direction** against A3's figures for the
unchanged build. In 13 of the 14 command-and-platform cells the change's range is not inside A3's,
but it is above in some cells and below in others. The largest gap is Linux `generate-psps` over 160
regions, whose highest run, 491.4 MiB, is 12% above A3's highest, 438.8 MiB; yet that command's own
three runs spread from 409.7 to 491.4 MiB, 20%. Two more examples: `call-from-alignments` over 160
regions on macOS, 1,378–1,451 MiB against A3's 1,339–1,447, and the fit on macOS, 240–243 MiB against
248–249. The two builds were never measured side by side for memory, so a difference of a few per cent
either way is not excluded. Those runs' **times** were taken under the load that §2.1 discards and are
not quoted.

## 3. Output

### 3.1 Calling with default parameters does not move

- **The identity oracle's five checksums are unchanged on both platforms** (`b74f3edf…` twice,
  `5e6d874d…`, `25d388c8…`, `259baf04…`; `tmp/measure_B7/oracle_{linux,macos}.log`).
- **Over all 160 regions, the two builds' VCFs are the same record for record** — 120,538 records,
  from `call-from-alignments` and from `call-from-psps`, on macOS and on Linux. Each build's records
  are also the same on the two platforms. (These are the outputs of the timing runs in §2.2,
  `tmp/ab2/{linux,macos}/{cfa,cfp}_{std,float}.vcf`, compared field by field with `##` headers left
  out.) The B2 review asked for a check of this kind. The repeat-tract aligner has two algorithms —
  a flat-gap one and the unit-slip one (`SsrUnitRobustAligner`) the caller ships — and its
  recorded-answer test, which compares against answers the deleted older caller wrote down, covers
  only the flat-gap one. The 160-region VCFs exercise the shipped one.
- **`estimate-contamination` could not be compared.** The two builds wrote identical JSON over the 20
  regions (`tmp/measure_B7/outputs/contamination_{std,float}.json`), but neither file holds an
  estimate: with four samples, all four were refused ("the panel is too small or too uniform"), over
  0 varying positions, and every contamination value is null. The files could only have differed in
  their position counts. The B5 review's question — does the conversion move this command's
  estimates — is still open (§5).

### 3.2 The fit moves, once, and then agrees across platforms

`estimate-parameters` on the four 20-region psps written on macOS:

- **The change writes the same parameters file on macOS and on Linux** — MD5 prefix `0fb3e50d`, on
  all six runs earlier in this step, and again in the last fit round of each platform in §2.2. The unchanged build's two
  platforms still differ from each other (`dfdefb5d` on Linux, `9d8c5539` on macOS; A3 found 7 of
  574 lines apart).
- **Against the unchanged build, 56 of 574 lines differ on Linux and 57 on macOS**
  (`tmp/measure_B7/outputs/{std,float}.parameters.toml`; `tmp/ab2/macos/fit_{std,float}.toml`). A3
  counted 57 for the libm build on Linux and
  described them: rows of the repeat-tract slippage model, per-repeat-length concentrations, one read
  group's error-probability multiplier, the four inbreeding coefficients and the two genome-wide
  concentrations, the largest relative change among values above 1e-4 being 0.29%.
- **Against A3's libm build, 8 lines differ**, and one constant explains all of them. The joint fit's
  `ln(3)` used to be computed while building, with the build machine's library, which gives the
  correctly rounded `0x3ff193ea7aad030b`; through libm it is one step lower, `…030a`. (A *step* here
  is the gap between one `f64` and the next representable one.) With the old bits
  put back for one run, the change's macOS fit reproduced the libm build's file (`3f5e5e7f`) exactly
  (`tmp/parity/B5_probe_ln3_macos/`). libm's value is kept: the plan writes a constant out as its old
  bits only where converting it moves a recorded answer, and none moved.

### 3.3 What the moved fit does to calls

Linux, the macOS psps, each build calling `call-from-psps` with its own fit
(`tmp/measure_B7/outputs/calls_{std,float}.vcf`):

| records | differ | fields that differ | genotypes, QUAL, FILTER or sites that differ |
|---:|---:|---|---:|
| 6,735 | 30 | AF 9, PARALOG_POST 19, PARALOG_LR 2 | 0 |

AF is the allele frequency written to the VCF; PARALOG_LR and PARALOG_POST are the hidden-duplication
filter's likelihood ratio and posterior. Each differs in its last printed digit, for example
`AF=0.099730` against `AF=0.099731`. **No genotype moved, so the plan's condition for running the
GIAB benchmark was not met.** For scale, A3 measured the unchanged build calling with the macOS fit
against the Linux fit: 46 records differ, one of them in QUAL.

And the other direction: called with the libm build's fitted parameters, the change writes that
build's VCF records byte for byte on both platforms (`1e40bd1e`; checked at B4 and B5).

### 3.4 Recorded answers that moved in B, and how

Two kinds, both in steps B2 and B4, both recorded in those commits:

- **One tie flipped, and it was undone.** The flat-gap aligner's gap-open cost for a repeat tract is
  `ln(0.01)`, which sits almost exactly halfway between two `f64`s. glibc and Apple give
  `0xc0126bb1bbb55515`, libm `…516`. That single unit changed where one generated read's tract began
  in `alignment::delimit_parity` (seed `0x5eed0001`, case 28). The cost is now written out as its old
  bits, and a test ties it to `float::ln` within one step.
- **The hidden-duplication filter's recorded-answer tests now compare probabilities to within 64
  representable numbers**, instead of bit for bit. Over those tests 17 of 465 probabilities moved,
  identically on macOS and Linux — 10 by one step, 5 by two, 1 by eight, and the fitted duplication
  rate by 22 steps, a relative 3.6e-15. Every decision the tests check — which record is dropped, the
  EM's convergence, the counts, the cut — is still compared exactly.

## 4. Tests

`cargo test --lib --tests --all-features`:

| platform | passed | failed | ignored |
|---|---:|---:|---:|
| Linux (container) | 4,796 | 0 | 4 |
| macOS | 4,795 | 0 | 4 |

Run on `cd056a3f` with this report uncommitted, 2026-09-15. Before the work began, at `75722336`,
the counts were 4,787 and 4,786 passed, so B1–B6 added nine tests on each platform. macOS runs one
test fewer because `psp::writer`'s failed-flush check is Linux-only. The container's other gates
passed on the same tree: `cargo fmt --all -- --check`, `cargo clippy --all-targets --all-features --
-D warnings` and `cargo doc --no-deps --all-features` with warnings as errors.

**Three example tests fail, on both trees**, under `cargo test --examples` on macOS: two in
`ng_generic_loci_dump` (`only_the_rows_that_departed_from_the_reference_carry_a_chain_id`,
`a_deletion_across_a_region_boundary_keeps_the_support_a_single_region_walk_gives`) and one in
`ng_ssr_loci_dump` (`a_cap_above_the_depth_is_invisible_and_a_cap_below_it_bites`). The unchanged tree
fails the same three (`tmp/ab2/example_tests_{std,float}_macos.log`), so the conversion did not cause
them. The unchanged tree's examples are those of `main`. CI runs `cargo test --lib --tests --examples`
on x86_64 Linux, where this was not run; if the three fail there as they do here, **that CI step fails
on `main` today**.

## 5. What was not covered

- **Linux on x86_64**, which CI and the Linux dev box run. glibc there picks one of two builds of
  `exp` and `log` at start-up, with and without fused multiply-add, so the unchanged build's output
  can differ between two x86_64 machines. The change's output should not, since libm is the same Rust
  on every target; nothing here has shown it. Speed on x86_64 is unmeasured.
- **GIAB.** Not run, because no genotype moved (§3.3). Four samples at 3× is one corner of the range;
  a decimal that moves here could move a genotype near a threshold in a deeper or larger cohort.
- **`estimate-contamination`, speed and output.** It was not timed in the quiet session. Its output
  could not be compared, because on four samples it refuses every sample (§3.1); whether the
  conversion moves its estimates needs a panel large enough to produce one.
- **Cohorts of more than four samples, deeper than about 3×**, except as the benches sweep them in
  miniature.

## 6. Deviations from the plan

- **The plan said to repeat A3's runs with `scripts/portable_float_baseline.sh`.** That script times
  one build at a time, and on this host that design did not survive the noise (§2.1). The timings
  here come from a script that interleaves the two builds (`tmp/ab2/ab2.sh`, `gen.sh`); the baseline
  script's runs supply peak memory only (§2.4).
- **A3's libm build was left out of the timing.** Its question — does routing through libm cost
  about 30% on the fit — is answered, and one Linux fit takes 7.4 to 9.9 minutes per build.
- **The 20-region calling commands were not re-timed** in the quiet session; only the 160 regions,
  where a small per-region cost is less hidden by start-up time.
- **The parity oracle was not re-run at `cd056a3f`.** §3.3 cites its result from B4 and B5. B6, the
  only commit since, changed nothing in `src/` that computes: it rewrapped a comment and added a
  lint `allow` to a test.
- **The joint-fit and site-quality benches were compared in the same session against the unchanged
  tree**, rather than against A3's numbers; the other three once against A3's, as planned. On Linux the delimiter bench was also
  run in the same session against the unchanged tree, after its comparison with A3 looked slower
  (§2.3).
- **`generate-psps` was timed in its own interleaved runs**, after the benches, because the script
  that timed the other commands did not include it.
- **The bench script failed once** on a variable-name clash after its first build, and was rerun into
  a fresh directory (`benches_macos2`); no timing was lost.

## Appendix A. Binaries and commands

| build | platform | SHA-256 prefix | source |
|---|---|---|---|
| unchanged | macOS | `5252cff30207d56c` | `../pop_var_caller-portable-float-std`, `c6a4394b`; rebuilt during this step to the same checksum |
| unchanged | Linux | `5fd016bfbb484c6d` | the same worktree, container build |
| change | macOS | `68b3954875b1517a` | this worktree, `cd056a3f` |
| change | Linux | `a197b7f30030ab5b` | this worktree, container build |

Copies in `tmp/ab/bin/{std,float}_{macos,linux}`, which is what every timing ran.

```sh
# commands and fit, interleaved: Linux then macOS
DEV_EXTRA_MOUNT=/Users/jose/devel/pop_var_caller/benchmarks ./scripts/dev.sh sh tmp/ab2/ab2.sh linux 5 3 tmp/ab2/linux
sh tmp/ab2/ab2.sh macos 5 4 tmp/ab2/macos
DEV_EXTRA_MOUNT=/Users/jose/devel/pop_var_caller/benchmarks ./scripts/dev.sh sh tmp/ab2/gen.sh linux 3 tmp/ab2/gen_linux
sh tmp/ab2/gen.sh macos 3 tmp/ab2/gen_macos

# host load, alongside
sh tmp/ab2/loadlog.sh tmp/ab2/hostload.tsv

# benches: both trees, into target-bench / target-container-bench
sh tmp/ab2/benches.sh macos tmp/ab2/benches_macos2
sh tmp/ab2/benches.sh linux tmp/ab2/benches_linux
uv run --no-project python tmp/ab2/benchtable.py <platform> <dir> tmp/baseline_A3/<platform>_benches

# VCF comparison, field by field
uv run --no-project python tmp/ab2/vcfdiff.py <a.vcf> <b.vcf>
```

Every command's inputs are in `ab2.sh`: the reference `S_lycopersicum_chromosomes.4.00.fa`, the
catalog `S_lycopersicum_chromosomes.4.00.repeats.parquet`, the four CRAMs under
`benchmarks/tomato2/bams/`, the 160-region psps in `tmp/baseline_A3/redo/macos_runs160/work/psps/`
and the 20-region psps in `tmp/baseline_A3/macos_runs/work/psps/`.
