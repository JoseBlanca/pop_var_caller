# Parameter estimation that fits in memory

*Status: plan, 2026-09-26. Nothing here is built yet.*

This plan is the build order for three changes to `estimate-parameters`, and a final measurement
that checks them, decided by the owner on 2026-09-26 after the command ran out of memory on a real
cohort. It is for the people implementing them. Each issue says what to change, where, how to know it
is right, and what it must not change.

**None of the three changes a fitted number.** Every one either frees memory earlier, changes how
the work is scheduled, or rearranges the same integer sums. The cross-platform checksum test
(`src/cli/cross_platform_digests.rs`) pins the parameters file and the VCF called with it to exact
bytes; **it must pass unchanged after each issue.** If it moves, the change did more than it was
meant to, and the implementer stops and finds out why.

## 1. What happened

On 2026-09-25 the reads pipeline ran `estimate-parameters` over 2,169 tomato psps (whole genome,
SL4.00) on kimura, a machine with 96 cores and 125 GiB of memory. The kernel killed it 9 h 22 min in,
holding 119 GiB.

A run of the same command over a random 100 of those psps, with every progress line printing the
process's memory, showed where the memory goes:

| step | resident at its end | most so far |
|---|---|---|
| reading the censuses and the reference | 0.3 GiB | 0.4 GiB |
| SNP/indel fit, throughout | 0.7–0.8 GiB | 0.8 GiB |
| contamination, done | 2.0 GiB | 2.5 GiB |
| repeat-tract evidence gathered (86,688 tracts in 141 strata) | 2.7 GiB | 2.7 GiB |
| repeat-tract fit, one minute in | 11.8 GiB | 13.2 GiB |
| repeat-tract fit, later | 12–13 GiB | 14.5 GiB |

Every figure except the first grows with the number of samples, so at 2,169 samples the repeat-tract
fit alone would need several hundred GiB.

## 2. What the command does, and what each step holds

`estimate-parameters` runs five steps, strictly one after another. The first two form one group and
the last three another; **the only thing that crosses from the first group to the second is each
sample's inbreeding coefficient, one number per sample.**

1. **The SNP/indel fit** (`fit_jointly`, `src/parameter_estimation/joint/fit.rs`). One statistical
   model fitted by repeated passes over the data: sequencing error rates, each sample's inbreeding
   coefficient, the population's allele-frequency distribution, and two classes of positions
   (mismapped and duplicated) that keep those numbers honest. It reads every sample's census
   evidence at the 2 million ordinary positions the census selected.
2. **Contamination** (`fit_contamination_over`, `contamination.rs`). Runs on the same evidence, using
   the error rates, inbreeding coefficients and per-position mismapping probabilities from step 1.
3. **Gathering the repeat-tract evidence** (`gather_strata`, `ssr_fit.rs`). Reads every sample's
   census evidence at the selected repeat tracts and arranges it by *stratum* — all tracts of one
   motif length whose reference carries the same number of copies.
4. **The per-stratum fit** (`fit_strata`, `ssr_fit.rs`). Each stratum's slippage and allele-length
   spectrum, from three starting points. Uses the inbreeding coefficients from step 1.
5. **Curves across strata**, and values for strata too thin to fit on their own. A quick closed-form
   step over step 4's results.

What each holds, with S the number of samples:

| step | held in memory | size |
|---|---|---|
| 1 | every sample's ordinary-position evidence, decoded together | about 2.5 MB per sample and read group; 7.5 GB at 2,169 samples (measured) |
| 2 | step 1's evidence, plus **two arrays of S × 2 million counts** | the arrays are 8 bytes per sample per position: 1.6 GB at 100 samples, **35 GB at 2,169** |
| 3–5 | the repeat-tract evidence (`Vec<StratumEvidence>`) | 0.7 GiB at 100 samples; about 15 GiB at 2,169 |
| 4 | a likelihood table per stratum being fitted | tracts × samples with reads × 91 genotypes × 8 bytes: up to 0.36 GB at 100 samples, **7.9 GB at 2,169**, for a 5,000-tract stratum |

Step 4 runs one stratum per thread at once, so its tables add up across threads. A fix on
2026-09-25 (commit `b435adfe`) made the fit switch to one stratum at a time when that would not fit
in half the available memory. Issue 2 replaces that guess with the user's explicit choice.

## 3. The order

1. ✅ **Issue 1: release the repeat-tract evidence as soon as the fit no longer needs it.** (§4)
2. ☐ **Issue 2: the user says how many strata are fitted at once.** (§5)
3. ☐ **Issue 3: build contamination's markers in two passes, without the S × 2 million arrays.** (§6)
4. ☐ **Final step: simulate cohorts of growing size and check the memory against what §7 expects.** (§7)

The owner's order. Each issue is one branch, merged before the next starts, so a checksum that moves
is traced to one change. The final step runs after all three are merged.

**Checkpoints — stop and hand back to the owner:**

- **after issue 3**, before the simulation starts, with the three merged changes summarised;
- **after the simulator is checked against kimura** (§7, "Check the simulator first"), before any
  cohort larger than 100 samples is run, with the 100-sample comparison;
- **at the end**, with the report.

## 4. Issue 1: release the repeat-tract evidence after the fit

**What it does.** After the per-stratum fit, the command keeps the whole repeat-tract evidence
alive until the parameters file is written, for one purpose: to read each stratum's substitution
rate, which is two counts per stratum. This issue keeps those counts and drops the evidence.

**What it does not do.** It does not lower the command's peak. The per-stratum fit reads the
evidence throughout, so the evidence and the fit's tables are resident together, and that moment is
the peak whether or not the evidence outlives the fit. What this issue frees is about 15 GiB at
2,169 samples during the assembly of the parameters file, and it removes a large value from
`CohortFit`, which tests and examples hold. **Making the evidence itself smaller is what would
lower the peak**; that is listed in §7 as the next step, not part of this plan.

**Where.**

- `CohortFit::tract_evidence: Vec<StratumEvidence>` (`src/run/census_fit.rs`, the struct near the top
  of the file) is filled in `fit_a_cohort` and read only in `parameters_from_the_fit`, by the loop
  that builds `ssr_substitution_rate` from `stratum.substitution_rate()` and
  `stratum.bases_compared`.
- Replace the field with a per-stratum list of exactly what that loop reads: the stratum, its
  `mismatching_bases` and its `bases_compared` (or the rate `substitution_rate()` returns, plus
  `bases_compared`). Name the new type for what it holds, e.g. `StratumSubstitutionCounts`.
- In `fit_a_cohort`, build that list from `evidence` right after `ssr_fit::fit_strata` returns, then
  let `evidence` go out of scope.

**Checks.**

- The loop in `parameters_from_the_fit` produces the same `ssr_substitution_rate` map. The existing
  tests of `parameters_from_the_fit` and `parameters_file_of`, and the checksum test, cover this.
- Search for other readers of `tract_evidence` (tests, `examples/`) and move them to the new field.

## 5. Issue 2: the user chooses how many strata are fitted at once

**Why explicit.** The owner's rule: when a setting decides whether a run fits in memory, the user
states it; the program does not guess it from what the machine reports. Commit `b435adfe` guesses —
it reads `MemAvailable` and switches arm when 96 tables would not fit in half of it. That guess is
removed.

**The two schedules** the fit already has (`WhereTheThreadsGo`, `ssr_fit.rs`), which return the
same bits by contract and by test:

- **One stratum at a time** (`AcrossTheTractsOfOneStratum`): every thread works on one stratum's
  tracts. Memory: one table.
- **Several strata at once** (`AcrossStrata`): each stratum's fit runs on one thread. Memory: one
  table per stratum running.

**The option.** `estimate-parameters --str-param-estimates-at-once <N>`:

- **N = 1, the default**: one stratum at a time, with every thread on its tracts. Always the smallest
  memory, whatever the cohort.
- **N ≥ 2**: N stratum fits at once, each on one thread, on a pool of N threads built for this step
  (`rayon::ThreadPoolBuilder::new().num_threads(N).build()` and `install`). The other cores are idle
  during this step when N is below the core count; that is the cost of the simple scheme, stated in
  the option's help.

The default is 1 rather than "every core" because the safe value is the one to fall back on: a run
that has not been told anything should not be able to run out of memory in this step.

**Tell the user what N costs.** The stage's first progress line already prints which schedule runs.
Make it also print the largest stratum's table size and what N of them would take, so the next run
can choose N from a number instead of a guess:

```text
estimating: repeat-tract fit: fitting 141 strata from 3 starting point(s) each, one stratum at a time; the largest stratum's table is 7.9 GiB, so --str-param-estimates-at-once N needs about N × 7.9 GiB
```

**Make the table size true.** When a walk moves to new slippage values, `Scorer::refresh`
(`ssr_fit.rs`) builds the new tables while the old ones are still held in `self.prepared`, so a walk
briefly holds two tables. Set `self.prepared = None` before building the new ones. This halves the
walk's peak and makes the printed size the real one. It changes no number: the old tables are not
read while the new ones are built.

**Where.**

- `EstimateParametersArgs` (`src/cli/estimate_parameters.rs`): the new option, a `NonZeroUsize`.
  Its help says what one stratum's table costs and that 1 is the smallest memory.
- `SsrFitConfig::threads` (`ssr_fit.rs`) becomes the count of strata at once, or carries it; the CLI
  passes it where it builds `SsrFitConfig::default()` in `fit_and_assemble`.
- `fit_strata` (`ssr_fit.rs`): delete `arm_that_fits` and its call, and the `the_arm_that_fits` test
  module. Keep `largest_table_bytes` for the progress line. Delete `memory_available`
  (`src/parameter_estimation/progress.rs`) if nothing else uses it.
- Run the several-at-once arm inside the N-thread pool.

**Checks.**

- The checksum test passes unchanged at the default.
- A test fits one small cohort's strata at N = 1, 2 and 4 and requires identical outcomes. The
  existing parity tests between the two arms stay.
- A test of `Scorer::refresh` that the tables it leaves are the ones a fresh `Scorer` builds.

## 6. Issue 3: contamination's markers in two passes

**What the arrays are for.** `markers` (`src/parameter_estimation/joint/contamination.rs`) chooses
the positions contamination is measured at — *markers*, positions where the cohort's common
alternative allele is at a usable frequency — and estimates each sample's genotype at each one. To
do that it fills two arrays of S × 2 million: each sample's reads on the common alternative allele,
and each sample's depth, pooled over its read groups. At 2,169 samples they take 35 GB.

**Neither question needs them.**

- **Choosing the markers** reads, at each position, only sums across samples: how many samples have
  any depth there (`covered`), the total alternative reads and the total depth. Those can be added
  up one sample at a time into three arrays the size of the position list.
- **Each sample's genotype** is needed only at the markers. On the 63-sample tomato cohort that was
  52,525 markers out of 2 million positions, one in 38.

**The new shape** — the same arithmetic in a different order:

1. Compute `major` (the common alternative allele at each position) exactly as now.
2. **First pass, one sample at a time.** Fill two scratch arrays of 2 million with that sample's
   pooled alternative reads and depth, exactly as the current code fills `alternative[s]` and
   `depth[s]`, including the `saturating_add` over read groups. Add them into `covered` (count of
   samples with depth above zero), `total_alternative` and `total_depth`, as `u64` sums in sample
   order. Clear the scratch arrays and move to the next sample — one pair of buffers for the whole
   pass, not one a sample.
3. **Choose the markers** from those totals and `noisy_posterior`, with the same tests in the same
   order as the current loop: mismapping posterior, `MIN_SAMPLES_WITH_DATA`, zero depth, then
   `MIN_FREQUENCY`. This gives `pooled` per marker and `marker_at`.
4. **Second pass, one sample at a time.** Refill the scratch arrays for that sample and, for each
   marker, compute that sample's dosage with the same formula as now, from the same two counts and
   the marker's `pooled`. Write it into `marker.dosage[s]`, allocated with S entries when the marker
   was created.
5. The existing per-read-group pass that fills each marker's `alternative`, `depth`, `depth_low`
   and `depth_high` is unchanged.

**Why the answer cannot move.** The counts are integers, so a sum is the same whatever order it is
taken in, and each sample's pooled counts are computed by the same code as now. The dosage for
(sample, marker) is computed by the same formula from the same numbers; only the loop order across
samples and markers changes, and no floating-point value is summed across that loop.

**Memory.** The two S × 2 million arrays (35 GB at 2,169 samples) become two scratch arrays of 2
million plus three totals, about 40 MB whatever S is. What remains is what the fit itself needs: the
dosages, S × markers × 8 bytes, about 0.9 GB at 2,169 samples and 52,525 markers.

**Cost in time.** The per-sample pooling runs twice instead of once. It is a walk over each
sample's depth codes, far cheaper than one pass of the SNP/indel fit, which runs hundreds of times.

**Checks.**

- **A differential test**: keep the current dense `markers` as a test-only function and require the
  new one to return identical markers — every field, compared exactly — on drawn cohorts. Include one
  with more than one read group in a sample, one where a position is covered by exactly
  `MIN_SAMPLES_WITH_DATA` samples, and one with a mismapped position.
- The existing contamination tests and the checksum test pass unchanged.

## 7. Final step: simulated cohorts of growing size

**What it answers.** Whether the three changes give the memory this plan predicts, and **how many
samples fit in 40 GB** on the development Mac, which has 64 GB. 40 GB leaves the machine room to
work while a run is going.

**What is expected.** After issues 1–3, with one stratum fitted at a time, the peak is in the
repeat-tract fit: the repeat-tract evidence plus one stratum's table. Per sample, from the kimura
measurement at 100 samples and three reads a position:

| held at the peak | per sample | where the figure comes from |
|---|---|---|
| repeat-tract evidence | about 7 MB | 0.7 GiB added by gathering, at 100 samples |
| the largest stratum's table (5,000 tracts) | about 3.6 MB | 5,000 × 91 × 8 bytes, one row per sample with reads |
| fixed cost (reference, selection, census directories) | — | about 0.3 GiB, before any fit |

So the peak should grow by about **10.6 MB a sample**, and 40 GB should hold **about 3,700 samples**.
The SNP/indel fit and contamination should stay below it, at about 5 MB a sample (0.8 GiB at 100
samples, less the fixed cost) plus contamination's dosages, 0.4 MB a sample at 50,000 markers. A
simulation that disagrees with these figures by more than about a fifth has found something this
plan did not account for, and that is the result to report, not a number to tune away.

**The cohorts.** Build them on disk as the real command would read them, because a file-backed
census and an in-memory one hold different things: a real run reads each sample's census from its
psp file on demand (`Sections::Backed`), while a census built in memory keeps every section decoded
(`Sections::Resident`). So each simulated sample is drawn in memory, written with
`census_file::write_census`, dropped, and reopened with `census_file::open_census`.

The drawing code already exists in `examples/dhat_ng_joint_fit.rs` (`drawn_generic_cohort`,
`drawn_ssr_cohort`, `draw_stratum`), at sizes too small for this. The new example draws at the
kimura cohort's shape, so its measurement can be checked against the real one:

- 2 million ordinary positions, at three reads a position;
- 86,688 repeat tracts in 141 strata, with the largest stratum at the census cap of 5,000 tracts;
- one read group per sample, and a second in one sample in ten.

**Check the simulator first.** At 100 samples the simulation must reproduce the kimura figures in
§1 — 0.8 GiB during the SNP/indel fit, 0.7 GiB added by the gathering — within about a fifth. If it
does not, the simulated samples are not the size of real ones, and the simulator is fixed before any
larger cohort is run.

**The runs.** Call the same functions `estimate-parameters` calls, in the same order: `fit_a_cohort`
(`src/run/census_fit.rs`) where a matching `CensusLoci` can be built for the drawn censuses, or else
its three parts, `fit_jointly`, `gather_strata` and `fit_strata`, which is all it adds to the digest
check. Print the memory at each step with the progress lines' own figures (`VmRSS` and `VmHWM`), so
the numbers are read the same way as kimura's.

**Short runs, because the peaks do not depend on run length.** Each sample's evidence is decoded
once, and each stratum's table is built in its first round, so one cycle of the SNP/indel fit and one
round of one starting point in the repeat-tract fit reach every peak. Set `max_passes` to 3 (one
accelerated cycle) and the repeat-tract fit to one starting point and one round. At full length a
cohort of thousands would take days on eight cores.

**Where to run.** Memory is read from `/proc`, so the runs go in the Linux dev container, started
with `DEV_MEM=48g DEV_CPUS=8 ./scripts/dev.sh …`. The container's 48 GB is also the backstop: a run
that overshoots is stopped by the container, not by the Mac.

**The series.**

1. Double the cohort from 100 samples — 100, 200, 400, 800, 1,600, 3,200 — until the peak passes
   40 GB or the next doubling would exceed the container.
2. Between the last size under 40 GB and the first over it, halve the interval until the largest
   cohort under 40 GB is known to within 100 samples.
3. **One depth point.** Repeat the 400-sample cohort at 30 reads a position. The caller has to work
   from a few reads a position to several hundred, and the evidence per sample grows with depth
   (more non-reference observations, more reads at each tract), so the per-sample cost at three
   reads a position is not the whole answer.

**What to report**, in `doc/devel/reports/implementations/`:

- a table of cohort size against peak memory at each step (SNP/indel fit, contamination, gathering,
  repeat-tract fit), with the expected figure beside each measured one;
- the per-sample cost of each step, from a straight-line fit across the sizes, and the fixed cost;
- the largest cohort under 40 GB, at three reads a position, and what the 30-reads point says about
  deeper data;
- anything that grew faster than its expected figure, with what grew.

## 8. After this plan

**Measure on kimura.** Rerun the 100-sample and 400-sample subsets with the merged build and compare
the memory on each progress line with §1's table:

- contamination should add almost nothing where it added 1.7 GiB at 100 samples;
- the repeat-tract fit at the default of one stratum at a time should add one table, not tens of GiB;
- the memory after the repeat-tract fit should drop by the size of the evidence.

Then the full 2,169-sample run. Expected peak, from the arithmetic above: the repeat-tract evidence
(about 15 GiB) plus one 7.9 GB table, with the SNP/indel fit and contamination below it at about
8–9 GB.

**Next, not in this plan: a smaller repeat-tract evidence.** The evidence is the largest thing left
at the peak. `StratumEvidence` holds, for every tract and every sample with reads there, a separately
allocated vector per read group of 17 counters, about 150 bytes per tract per sample at three reads a
position where most counters are zero. A flat layout, like the one `TractLikelihoods` already uses
for its rows, would cut that several-fold. Gathering one band of strata at a time, as the census
specification describes, would also stop the decoded census sections and the arranged evidence from
peaking together. Both are separate designs.
