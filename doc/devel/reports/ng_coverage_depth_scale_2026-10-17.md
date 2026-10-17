# The coverage model refused at high depth: cause, fix and measurements

2026-10-17, branch `coverage-depth-scale` (from main 354754b8).

## 1. Summary

A colleague ran pvc 0e905bbb on GIAB HG002, HG003 and HG004, chromosome 22 only, one sample a
run. At 201–223 reads a position every sample's coverage model was refused, with 0.346 of its
windows past the top of the depth axis. Without a coverage model the hidden-duplication filter
scores nothing. At 5, 10 and 30 reads a position the same samples were accepted.

- **The defect is real, and it is in how the depth axis is scaled.** The axis was fitted from the
  sample's first 10,000 windows. A synthetic input whose first 19 kb are twenty times shallower
  than the rest reproduces the refusal exactly: 0.949 of the windows overflow, and the fit refuses
  the sample. The mirror case, a first 19 kb ten times deeper than the rest, is worse because it is
  silent: the fit accepts the sample, but squeezes its whole one-copy peak into the axis's first five
  bins and reports a scatter that is only the axis's own resolution.
- **The chromosome-22 refusal itself was not reproduced here**, because its input files are not on
  this machine. Every local GIAB input is a panel of 100 regions a sample, and on those the fit is
  accepted at every depth. §3.3 says what the chromosome-22 runs would show if the same cause is
  behind them, and which files would settle it.
- **The fix fits the axis from every window, not from the first ones**, in one pass and with no
  extra memory (§4). On the panel, at 5 to 300 reads a position, the fitted one-copy depth moves by
  at most 2.2% and its scatter by at most 0.007. Every written record is unchanged apart from its
  two filter scores, except one more record dropped at HG004 300x, which the truth set confirms is
  a false positive. Both synthetic inputs are fixed.

## 2. What the depth axis is, and why it is scaled per sample

The hidden-duplication filter compares each locus's read depth against what one copy of the genome
gives in that sample. It learns "one copy" from a histogram of every 500-base window's mean depth,
split by GC content, and the coverage fit takes the histogram's main peak as the one-copy depth.
The fit's spread around that peak is the scatter, σ₀.

The histogram has 400 depth bins plus an overflow bin. A fixed bin width cannot serve both three
reads a position and three hundred, so the width is scaled to the sample: the 400 bins span ten
times the sample's median window depth. **If more than a fifth of the windows land in the overflow
bin, the fit refuses the sample**, because the one-copy peak itself has then left the axis
(`src/paralog/coverage_model.rs`, `DEFAULT_MAX_OVERFLOW_FRACTION`).

Until this change, the median the width was scaled to was taken over **the sample's first 10,000
windows** (`DEPTH_SCALE_WINDOWS`), the windows were held back until then, and the width was never
refitted. There is one window per covered position, so those are the first 10 kb or so of covered
ground. The code's own comment recorded that this prefix median came out low in 8 of 10 stores, by
up to 34%, and that ten medians of headroom absorbed it.

## 3. Is it real?

### 3.1 The GIAB per-sample panel: accepted at every depth

`benchmarks/giab/per_sample/`: each sample has its own 100 random regions, 572 kb in all for HG002.
The 300x inputs are HG002's "reads selected" CRAM and HG003's and HG004's BAMs; the lower depths are
the `seed42` downsamples. Each sample was run alone through `generate-psps` and then
`call-from-psps --defaults`, as the colleague ran it.

Every one of the 12 runs fitted a coverage model. The prefix median, read with
`ng_window_coverage_probe --bin-scheme` on the unchanged code, sat between 20% below and 6% above
the median over every window:

| depth | sample | median of the first 10,000 windows | median of every window |
|---|---|---|---|
| 5x | HG002 / HG003 / HG004 | 4.48 / 4.70 / 4.48 | 4.95 / 4.88 / 4.96 |
| 30x | HG002 / HG003 / HG004 | 24.09 / 29.15 / 30.28 | 30.14 / 29.86 / 30.08 |
| 300x | HG002 / HG003 / HG004 | 245.7 / 315.7 / 352.4 | 304.6 / 308.7 / 333.8 |

An axis 20% short still reaches eight medians, so no window overflowed in any of these runs. The
panel does not reproduce the refusal. It shows that **the prefix error is in both directions on
real data** (HG004 300x is 6% high).

### 3.2 Synthetic inputs: the prefix alone causes both failures

Two inputs built from HG003's 300x BAM with `samtools` (`tmp/depth_scale/make_synthetic.sh`), run
over HG003's own 100 regions. The first three regions, 19.3 kb, hold the first 10,000 windows.

| input | first 19 kb | the rest | prefix median | median of every window | before the fix |
|---|---|---|---|---|---|
| shallow prefix | 5% of reads | all reads | 14.6 | 305.2 | **refused**: 0.949 of windows over the top |
| deep prefix | all reads | 10% of reads | 315.7 | 30.9 | accepted, one copy 30.01, scatter 0.132 |

The deep-prefix case is accepted, and its one-copy depth is right. Its bins are 7.9 reads wide,
though, so the whole one-copy peak sits in the first five bins. The fit cannot report a scatter
finer than half a bin, which here is 0.131 of the one-copy depth, and the scatter it reported is
0.132. The same sample's real 30x run reports 0.140 to 0.143. The filter dropped 25 records on this input;
with the fix it drops 19, and the 6 it now keeps are all variants in HG003's truth set (five
homozygous ones inside the deep first 19 kb, and chr17:35,690,325).

### 3.3 The chromosome-22 runs: what would confirm the cause

A prefix that is a fixed fraction of the sample's depth would refuse the sample at every depth, and
the colleague's runs were refused only at 300x. So if the prefix is the cause there, **which ground
makes up the prefix must change with depth**. Two mechanisms in the code can do that. Neither has
been tested.

- **The window floor.** A window is used only if 50 of its 500 positions have reads. On sparse,
  poorly mapped ground near the start of chromosome 22, windows that fall under the floor at 30x
  can pass it at 300x, so the prefix at 300x can be made of thin ground the 30x prefix skips.
- **The depth ceiling.** A record over 1,000 reads in one read group is not measured at all
  (`run::depth_ceiling`). At 300x a collapsed satellite at four or more copies passes 1,000 reads and
  drops out, while at 30x it stays in the prefix and holds its median up.

Either way the fix below removes the dependence, because no part of the run is privileged. To
confirm the cause rather than just the cure, one of these:

- the colleague's 300x and 30x psps for any one sample, run through the **unchanged** probe:
  `ng_window_coverage_probe --reference <GRCh38> --bin-scheme <sample>.psp`. Its
  `as-the-run-fits-it` row prints the prefix median, and its `scale-sample=every-window` row the
  whole-run median. The prefix median should be about a tenth of the whole-run median at 300x, and
  not at 30x. That probe binary is kept at `tmp/depth_scale/probe_before`;
- or the 300x chromosome-22 CRAM or BAM for any one sample, from which the psp can be made here.

### 3.4 Can the prefix come out too high on real data?

Yes, by 6% on HG004 300x and 2% on HG003 300x in §3.1. No real input here is ten times too high,
but nothing prevents one: a run whose first 10 kb is a collapsed repeat, such as a panel that
starts in a segmental duplication, is the §3.2 deep-prefix case. That case is never refused; it
only loses resolution, so it would never have been reported.

## 4. The fix

**Every window is counted, as it is finalised, on a fixed grid whose bins are equal steps of
relative depth. When the pass ends, the median is read off that grid, and the grid is re-cut into
the linear axis the fit has always read.** The coverage fit is not changed.

- **The grid.** 50 GC bins and 400 depth bins, with edges 2.33% apart from 1 read a position to
  10,000, and an overflow bin above. One read is the floor because a covered position has at least
  one read. 10,000 is ten times the default depth ceiling, so a default run never reaches it. The
  grid has as many cells as the histogram, 80.2 kB a sample.
- **The median** is the upper middle window, located inside its grid bin by interpolation, so it is
  good to within one grid step.
- **The re-cut.** Each grid bin's count is shared out among the linear bins it overlaps, in
  proportion to the overlap. Each GC row is then rounded back to whole counts so that its total is
  exact (the largest remainders get the missing counts, ties to the lower bin). This blurs every
  window by at most one grid step, 2.33% of its depth. That widens a one-copy scatter of 8% (the
  narrowest measured, at 300x) by 0.03%.
- **Memory falls.** The 10,000 held-back windows cost 262 kB a sample until the axis was fitted. They
  are gone, so a sample costs 106.4 kB for the whole pass instead of 368 kB at its peak: 319 MB
  instead of 1.1 GB at three thousand samples.
- **Every log and exponential goes through `crate::float`**, so the grid bins and the median are the
  same bits on macOS and Linux.

The code is `src/window_coverage/depth_grid.rs` (new) and `accumulator.rs`. `mod.rs` replaces
`DEPTH_SCALE_WINDOWS` with the three grid constants, each documented with the reason for its value.
`window_coverage.md` §3.4, §4, §5, §6 and §9 and `hidden_paralog_filter.md` §4 are updated.

### 4.1 The alternatives weighed

| option | refused sample (prefix too low) | squashed sample (prefix too high) | cost |
|---|---|---|---|
| **fixed grid, re-cut at the end (chosen)** | fixed | fixed | each window blurred by ≤ 2.33%; memory falls by 262 kB a sample |
| widen the axis in place: when overflow passes a threshold, merge bins in pairs and double the width | fixed | **not fixed**: bins cannot be split again | windows already in the overflow bin before each doubling stay there; keeps the 262 kB |
| scale from windows spread over the run (a reservoir sample) | fixed | fixed | needs every window kept until the end, or a second pass |
| hand the fit a log-scale axis | fixed | fixed | rewrites the fit's mode, median and scatter estimates, which assume equal bins, and changes the fitted one-copy depth at low depth (a peak on a log axis sits higher by about σ₀², 9% at σ₀ = 0.3) |

## 5. Before and after

Fitted one-copy depth (reads a window) and scatter, and records the filter dropped. Per-sample panel
of §3.1 and the synthetic inputs of §3.2. Raw logs: `tmp/depth_scale/{before,after2}/`.

| input | before: one copy, scatter, dropped | after: one copy, scatter, dropped |
|---|---|---|
| 5x HG002 | 4.77, 0.256, 0 | 4.81, 0.261, 0 |
| 5x HG003 | 4.39, 0.237, 0 | 4.43, 0.244, 0 |
| 5x HG004 | 4.61, 0.263, 0 | 4.59, 0.256, 0 |
| 10x HG002 | 9.59, 0.191, 0 | 9.77, 0.188, 0 |
| 10x HG003 | 8.27, 0.188, 0 | 8.15, 0.193, 0 |
| 10x HG004 | 10.45, 0.214, 0 | 10.46, 0.215, 0 |
| 30x HG002 | 29.64, 0.136, 0 | 29.57, 0.137, 0 |
| 30x HG003 | 29.55, 0.143, 0 | 30.19, 0.140, 0 |
| 30x HG004 | 28.38, 0.136, 0 | 28.34, 0.136, 0 |
| 300x HG002 | 299.55, 0.092, 7 | 302.07, 0.092, 7 |
| 300x HG003 | 305.36, 0.084, 0 | 306.62, 0.084, 0 |
| 300x HG004 | 310.71, 0.089, 1 | 310.33, 0.091, 2 |
| synthetic, shallow prefix | **refused** | 307.87, 0.086, 0 |
| synthetic, deep prefix | 30.01, 0.132, 25 | 29.38, 0.147, 19 |

- **Every sample has a coverage model after the fix.** The one refusal, the shallow-prefix input,
  now fits at 307.9 reads with scatter 0.086, beside 306.6 and 0.084 from the same sample's own 300x
  run.
- **At 5x to 30x the fitted numbers stay where they were**: the one-copy depth moves by at most 2.2%
  (HG003 30x, 29.55 to 30.19) and the scatter by at most 0.007 (HG003 5x). These moves come mostly
  from the bins being cut at the whole-run median rather than the prefix median, a different but
  equally valid set of bin edges.
- **The written calls do not change.** In all 12 panel runs, every record written before is written
  after with the same position, alleles, quality, genotype and other fields; only `PARALOG_LR` and
  `PARALOG_POST` differ. One record is dropped that was written before: HG004 300x, chr13:45,646,237
  A>G at QUAL 17.6, which is not in HG004's truth set. It is under the benchmark's QUAL-30 gate, so
  it moves no headline count.
- **What a working filter removes at 300x** was measured only where it had been refused, the
  shallow-prefix input: it scored 724 records and dropped none, as the same sample's real 300x run
  does. The colleague's segmental-duplication false positives cannot be measured until the
  chromosome-22 runs are repeated with the fix.

## 6. How a run reports a sample with no coverage model

Before this change the reason went to the log at debug level ("sample HG002 has no coverage model:
…"), and the VCF header said `samples_with_coverage_model=0/1`. The terminal showed only a
consequence, and only when no sample was scored at all: "the hidden-duplication rate could not be
fitted … the filter used the fallback rate". That warning does not say the filter scored nothing,
and a run where 3 of 1,000 samples were refused printed nothing on the terminal.

**Now one terminal warning is printed whenever any sample lacks a model** (decided by the owner,
2026-10-17). It gives the count, names up to three samples with their reasons, and says how many
more are named in the log:

```
warning: 5 of 6 sample(s) have no coverage model, so the hidden-duplication filter used no
coverage evidence from them: sample A (reason); sample B (reason); sample C (reason); and 2 more,
each named with its reason in the log
```

The per-sample lines stay in the log as before. The line is built in
`run::paralog_filter::finish::what_to_tell_the_operator`; the report already puts every line that
starts `warning:` on the terminal. Three tests cover it: two refused samples named with their
reasons, five refused of which three are named, and no warning when every sample has a model.

## 7. The oracle, the checksums and the tests

**The identity oracle** (`scripts/promote_ng_oracle.sh`: four tomato2 accessions over 20 regions,
about 5 reads a window; called with default parameters, fitted, and called with the fit) was run
on the macOS host with the change, and again with the pre-change binary, which reproduced every
recorded checksum. Its two calling routes agree, and the windows are identical bit for bit.

- **Unchanged:** the windows, and both parameters files.
- **Moved, as expected:** the histograms; with them the four samples' one-copy depths, by at most
  2.8% (4.66 to 4.53 reads), and scatters, by at most 0.006. The fitted duplication rates barely
  move: 0.031618 to 0.031693 with default parameters, 0.050664 to 0.051914 with the fit.
- **The calls:** apart from `PARALOG_LR` and `PARALOG_POST`, every record is unchanged except a
  few the filter now decides the other way. With default parameters (7,508 records to 7,507), 2
  more are dropped and 1 is kept. With the fit (6,680 to 6,682), 1 more is dropped and 3 are kept.
  The baseline is re-recorded, with these positions listed in its comment.

**The three cross-platform checksum tests** (`cli::cross_platform_digests`) are re-recorded. Their
parameters files did not move. In each VCF only `PARALOG_LR` and `PARALOG_POST` moved; with those
taken out, the old and new files are identical line for line, and the old code's files reproduce
the old checksums. Those fixtures hold fewer than 10,000 windows, so the old axis was already
fitted from every window, and the move is the grid's rounding alone.

**Tests.**

- `window_coverage`: 70 pass, seven of them new in `depth_grid` (the grid's step, the median at 3,
  30 and 300 reads, the re-cut keeping every window, overflow, rounding). Seven accumulator tests
  that pinned the prefix are removed: the size of the scale sample, held-back windows folding the
  same as live ones, and the prefix median. Two replace them. One shows that a sample's axis is
  fitted over all its contigs. The other is the regression test: deep windows first, then the real
  depth, and the axis comes from the real depth. Four tests that pinned an exact width or cell now
  allow one grid step. The test of a depth on the top edge moved to `depth_grid`, where the re-cut
  decides it.
- The probe's tests: 21 pass. Its existing test of this very failure, 10,300 windows at 1 read then
  19,700 at 100, now asserts that the run's own setting accepts the sample. One test fixture
  predated two fields `SequenceObservation` gained on main; it is fixed.
- The whole library passes, except three should-panic tests of `debug_assert!`s, which a release
  build turns off; they test alignment and `lgamma` code this change does not touch.
- Clippy is clean on every changed file.

**Still open:** the colleague's chromosome-22 runs, to confirm the cause there (§3.3) and to
measure what a working filter removes at 300x.
