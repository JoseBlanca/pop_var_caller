# The repeat-tract fit by penalised EM — implementation plan

**Draft, 2026-10-14; amended the same day after an adversarial review. Amended 2026-10-15: step
C8, the stutter prior, and step C9, the stated F's median over classes with 8 or more loci, both
approved by the owner 2026-10-15. Amended again 2026-10-15 with the owner's answers to C8's review:
the slip rate's error rule (C8), the one-sample cost recorded and measured at D1, the weight kept at
one read, and step C10, the cohort in one fixed order.** Turns the settled design —
specs [`../spec/str_penalised_fit.md`](../spec/str_penalised_fit.md) and
[`../spec/str_fit_to_calling.md`](../spec/str_fit_to_calling.md), architecture
[`../arch/str_penalised_fit.md`](../arch/str_penalised_fit.md) — into build order. **It invents no
design**: a question that comes up here goes back to the specs or the architecture.

## Scope

**In:** the penalised fit as `estimate-parameters`' repeat-tract fit; parameters-file format 3 with
reading of older files; calling's new fallback step; the checks of spec §4; deleting the climb and
what goes with it (arch §4.1).

**Out, each with its home:**
- stutter per slippage group and library multipliers → `spec/str_library_stutter.md`, rebuilt on this
  fit;
- reading the census one motif length at a time → the census-read step (spec §6);
- the full-scale run on kimura's real data → when kimura is back (§D3 measures it on a simulated
  cohort here meanwhile);
- re-recording main's own stale oracle baseline → done on this branch by the merge (`e0b41ad4`).

## Principles (how the order was chosen)

- **Types first, then behaviour.** Milestone A changes shapes only; every number stays the same except
  the one step that says otherwise (A5).
- **The hand-over before the fit.** Calling must hold an absent F (milestone B) before any fit can
  produce one (milestone C); in the other order calling panics (spec §3).
- **The prototype is the oracle, recorded before it moves.** The research bench's `D_tract_prior`
  gives the numbers the production fit must reproduce; C0 writes them down before C1 renames
  anything and C2 removes the bench methods that measured alternatives.
- **Isolate the silent steps.** Eight steps can change numbers without failing anything: the stated
  F's median (A5), F's newtype and conversions (B1), the keys holding F (B2), making the fit
  deterministic (C3), switching `estimate-parameters` to it (C7), the stutter prior (C8), the
  classes the stated F's median counts (C9) and the cohort's fixed order (C10). Each is its own
  commit, with its oracle run before and after.
- **Delete last.** The climb goes only after the new fit has passed on real data (D2) and in calling
  (D5).

## Preconditions (already in place)

- Branch `str-penalised-fit` with these amended documents: main merged in (`e0b41ad4`), specs
  approved.
- The research code on the branch: `ssr_fit/penalised_fit.rs`, `composite_fit.rs`,
  `counting_fit_errors.rs`, `bench_fixtures/{fit_bench,cohort_draw}.rs`,
  `examples/ng_str_fit_bench.rs`.
- **tomato2**: the psps in format 2.0, `/Users/jose/devel/pop_var_caller-str-guard/tmp/tomato2/psps`
  (59 accessions, written 2026-10-14 by main's `pvc`), and the reference and catalog they were
  written against, which `estimate-parameters` checks (`refuse_another_selection`,
  `run/census_fit.rs:233`): `$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa`
  and the main checkout's `benchmarks/tomato1/crams/S_lycopersicum_chromosomes.4.00.repeats.parquet`
  (from that directory's `generate.sh`). **The spec's 5.8 s was measured on a different tomato2 set**
  (`../pop_var_caller-stutter-multiplier/tmp/route2/psps_tomato2_32mb`, `cost_2026-10-14.md`), so a
  time on this set is not directly comparable; say which set a timing used.
- **HG002 and the GIAB trio**: the data sit in the main checkout's `benchmarks/` (`ssr_hg002/`,
  `giab/`); this worktree holds only their scripts. Mount them with `DEV_EXTRA_MOUNT`.
- **For D3, kimura's loci per class and samples with reads per tract.** The repository holds only
  kimura's totals (86,688 tracts in 141 classes, 107 with fewer than 8 tracts with reads, 66.2 million
  pairs with reads; `kimura_shapes_span9_2026-10-14.md`). Copy the per-class counts off kimura before
  D3, or D3's report states the per-class spread it assumed.
- Known test failures before step 1: three `cli::cross_platform_digests` tests (stale, spec §3) and
  three debug-only panic tests in release.

## The steps

### Milestone A — shapes only

- ✅ **A1** `StratumFit::concentration` becomes `Option<f64>`; the climb always fills `Some`. Every
  reader updated. *Depends:* —. *Source:* arch §1.2.
- ✅ **A2** The climb's own fields move into `climb: Option<ClimbRecord>`; its log summaries read
  them there. Updates the four examples that read those fields
  (`examples/ng_joint_records_walk.rs:951-957, 1843-1847`, `examples/dhat_ng_joint_fit.rs:382`,
  `examples/ng_str_neighbour_start.rs:133-136, 185`, `examples/ng_ssr_thin_stratum_gate.rs:146`), the
  16 `StratumFit { … }` literals in 9 files (`ssr_fit.rs`, `stratum_fits.rs`,
  `bench_fixtures/fit_bench.rs`, `calling/mod.rs`, `calling/inference/summarise_condition.rs`,
  `calling/parameters_file/{from,to}_run_parameters.rs`, `tests/ng_calling_loop_calls_genotypes.rs`,
  `examples/shared/stutter_table.rs`) and `benches/ng_joint_fit_perf.rs`; `cargo clippy
  --all-targets` passes. *Depends:* A1. *Source:* arch §1.2. *Done:* 10 literals in 8 files —
  `examples/shared/stutter_table.rs` defines its own `StratumFit`, so it was not touched.
- ✅ **A3** `LevelSource::PenalisedCurve`, `ShareSource::PenalisedCurve`, `PenalisedPoint`,
  `StratumError::TooFewLoci` — added, unused. `curve_weight()` returns `Option<f64>`; `warrant_of`
  matches the source exhaustively. *Depends:* —. *Source:* arch §1.3–1.4.
- ✅ **A4** `FittedAlleleLengthVariation::concentration` becomes `Option<f64>`, with
  `concentration_standard_error`; `AlleleLengthVariationRung::StratumsOwnSharesStatedSpread` and
  `FittedFrom::ThisStratumsSharesStatedSpread` added; `allele_length_variation_at`,
  `checked_allele_length_variation` and `fitted_allele_length_variation_of_each_stratum` handle
  `None`; the rung match at `to_run_parameters.rs:2016-2020` updated. Tests for the new rung with a
  hand-built table, and for its rank under `BundlePrior::rung`. *Depends:* A1. *Source:* arch §1.5;
  companion §3.2.
- ✅ **A5** **Own commit, do not bundle.** The stated concentration is the median taken in F; its
  warrant is `Defaulted` when no class has an F. *Oracle:* the oracle script before and after; with
  the climb every fitted class has a concentration, so a move can come only from an even count of
  fitted classes, and is explained. A test at one sample pins the warrant. *Depends:* A4. *Source:*
  arch §1.5; companion §3.5.

> **Checkpoint A:** the oracle (`scripts/promote_ng_oracle.sh`) gives the same checksums as at the
> precondition commit, or A5's move, explained; the touched tests pass. Pause for review.

### Milestone B — the parameters file, format 3

- ✅ **B1** **Own commit, do not bundle.** `BetweenLocusShare` and `between_locus_share_of`,
  `concentration_of`, `between_locus_share_error_of`, with tests: the round trip, F of 0 and 1 refused,
  an F below about 5.6×10⁻³⁰⁹ refused, a concentration below about 1.1×10⁻¹⁶ refused. No key changes.
  *Oracle:* the oracle's checksums unchanged. *Depends:* A4. *Source:* arch §1.6, §2.3.
  *Also (owner, 2026-10-14, from milestone A's review):* the stated value's local F conversion
  (`median_concentration`, A5) moves onto these helpers, keeping the even-count closed form
  (c1+c2+2c1c2)/(2+c1+c2), which cannot cancel to 0.
  *Done:* the median keeps its own closed form — `parameter_estimation` imports nothing from
  `calling` (`parameters_file/mod.rs:36-37`) — and a test beside the helpers pins the round trip
  through them equal to it.
- ✅ **B2** **Own commit, do not bundle.** The rows' and fallback's keys renamed and holding F;
  `between_locus_share_standard_error` written from the concentration's error;
  `FORMAT_VERSION = 3`; validation on `BetweenLocusShare`. Golden files regenerated by both ignored
  tests (`regenerate_the_golden_file`, `regenerate_the_written_golden_file`), diffs checked;
  `validate.rs:1737` flipped (version 3 read, 4 refused);
  `benchmarks/ssr_hg002/src/sweep_tract_parameters.sh:136` writes `fallback_between_locus_share`.
  *Oracle:* as checkpoint B. *Depends:* B1. *Source:* arch §1.6; companion §2.1, §4.
  *Also (owner, 2026-10-14, from milestone A's review):* validation refuses a period row
  (`allele_length_variation_by_period`) with no F — nothing produces one and calling has no rung for
  it (companion §3.4); a present F error is range-checked; `fallback_between_locus_share` marked
  `fitted_here` is refused when no stratum row carries an F (today's rule, "no stratum rows", at
  `validate.rs:752`, misses a one-sample file); the file's reader-facing note "this run fitted no
  stratum" (`to_toml.rs:347`, `:927-930`) reworded to "no stratum has an F".
  *Done:* `between_locus_share_of` and the error conversion pick, among F's nearest `f64`s, one
  that converts back exactly, so a file written, read and written again is byte-identical (every one
  of 100,000 draws; 7 million in the review); added helper `concentration_error_of`; a row's F is
  checked as the file is parsed; the writer stops, naming the class, on a concentration or error with
  no valid F. Oracle: only the two parameters files moved, re-recorded in the baseline.
- ✅ **B3** `migrate_to_format_3` on the parsed TOML; a version-3 file still deserialised from its
  text; `version_1_as_written.toml` restored to its old keys (from `aa42d629~1`); `from_toml.rs:367`
  expects version 3. Tests: a format-1 file and both format-2 spellings read back with F = 1/(1+c); an
  old file's error names its key path. *Depends:* B2. *Source:* arch §2.3; companion §2.2.
- ✅ **B4** `LevelSmoothing::PenalisedCurve` and `ShareSmoothing::PenalisedCurve`, written, read and
  validated (`strength` finite and positive); `of_the_numbers_emitted` (`stratum_fits.rs:107-126`)
  and `an_own_fit_error`'s callers (`validate.rs:810-817`, `:851-865`) made exhaustive matches that
  keep and accept an error beside `PenalisedCurve`. Test: an error from a `PenalisedCurve` outcome
  reaches the file and reads back. *Depends:* A3, B2. *Source:* arch §1.3, §1.6; companion §2.3.

> **Checkpoint B:** a format-3 file round-trips, F absent included. **What the oracle should show:**
> the VCFs called with default parameters unchanged; `from_alignments.parameters.toml`'s checksum
> moved (format version, renamed keys, the fallback as F); `fitted.parameters.toml`'s moved (the same,
> plus every row's F); the VCF called with the fitted file unchanged, or its move traced to the
> concentration read back as 1/F − 1 of a rounded F. Pause for review.

### Milestone C — the fit

- ✅ **C0** The prototype's reference numbers recorded as a fixture beside the bench: `D_tract_prior`
  at the cells of `compare_D_prior_2026-10-14.md` with the production settings (spec §2.7), every
  number and every class's rounds. *Depends:* —. *Source:* spec §4.2; arch §6.
  *Done:* `bench_fixtures/d_tract_prior_reference.jsonl.gz`, 153 runs, compared by
  `scripts/research/str_fit_fixture.py` (exact = bits; tolerance = counts equal, 1e-6 relative + 1e-12).
  The prototype fits every class of a motif length in one loop, so the rounds recorded are the loop's,
  not a class's. Found: its numbers move by up to 1.3e-8 relative with the thread count — C3's target.
- ✅ **C1** `composite_fit.rs` renamed `expected_counts.rs`; no code change. *Depends:* C0.
  *Source:* arch, Module home.
- ✅ **C2** `PenalisedFitConfig::production()`, `allele_span` and the named constants; the functions on
  the fit's path take the window from `PenalisedFitConfig`; the research-only fields removed, and with
  them the bench's `D_prior_*` methods and the clustering argument of `str_fit_penalised_em`; the
  bench's `D_tract_prior` built from `production()`. *Oracle:* C0's fixture reproduced. *Depends:* C1.
  *Source:* arch §1.1; spec §2.7.
- ✅ **C3** **Own commit, do not bundle.** Every float sum on the fit's path taken in a fixed order,
  in every module the fit calls, `SamplePairs::log_likelihood` (`composite_fit.rs:1583-1590`)
  included; a test pins two runs and 1 against 18 threads. *Oracle:* the bench's numbers equal to C0's
  within 1e-6 relative, and every class's rounds equal, before and after; a discrete decision that
  flips (spec §3) is investigated, not tolerated. *Depends:* C2. *Source:* spec §3; arch §2.1.
  *Done:* two parallel sums in `expected_counts.rs` (`SamplePairs::log_likelihood` and the counting
  fit's spread objective) now add per-locus values in locus order; the full grid is bit-identical at 1
  and 18 threads; against C0 every count equal, largest move 1.14e-7 relative (a standard error),
  estimates ≤ 2.75e-13. The fixture is re-recorded at this step, so C4 compares bit for bit.
- ✅ **C4** Memory bounded: classes one at a time within a period, periods one at a time, the
  parallel loops over classes (`penalised_fit.rs:588`, `:852`) made sequential. *Oracle:* bench
  numbers bit-identical before and after; peak memory at 2,000 samples × 30 reads below 2 × the
  counting fit's 1.1 GB. *Depends:* C3. *Source:* spec §2.2, §3; arch §2.1.
  *Done:* the loops over classes are sequential; peak at 2,000 × 30 fell from 3,263 to 1,085 MiB
  (one-base) and 2,201 to 634 MiB (two-base); grid bit-identical; 1.3–1.45× the wall time at 2,000
  samples. Within a class, `DistinctRows::of` blocks by ~8,192 rows and the errors build a class's
  table once. **"Periods one at a time" lands with C5**, which owns the loop over motif lengths.
- ✅ **C5** `fit_repeat_tracts`: per-period grouping, one slippage group asserted, the
  `StratumOutcome` conversion, F absent at an end of its search, the error mapping with
  `LOCI_FOR_AN_ERROR` and `LOCI_FOR_A_SHARE_ERROR`. Tested on drawn cohorts: outcomes in input order;
  `Fitted` for a class with reads, `Refused { NoSpanningReads }` for one without, never `Derived`;
  provenance `PenalisedCurve` with `reach` inside; `samples_fitted_on` `None`; F `None` at one sample;
  no `Estimated(0)`. *Depends:* C4, A2, A3. *Source:* arch §1.4, §2.1.
- ✅ **C6** The run's log of spec §2.9: per motif length, rounds, the five strengths, eased bends,
  time in each step and peak memory (the progress lines' `VmHWM`, `parameter_estimation/progress.rs:79-110`,
  Linux only; nothing is printed on macOS); per class, rows with reads, distinct rows, F absent; the
  count of withheld errors. *Depends:* C5. *Source:* spec §2.9.
- ✅ **C7** **Own commit, do not bundle.** `fit_a_cohort` calls `fit_repeat_tracts` and takes
  `&PenalisedFitConfig`; `repeat_tract_config` returns `production()`;
  `--str-param-estimates-at-once` removed with an error, and the `EstimateParametersArgs` literals
  (`cli/cross_platform_digests.rs:371`, `cli/estimate_parameters/tests.rs:98`) updated. *Oracle:* the
  oracle script, run before and after; the move is explained (expected: repeat-tract rows of the
  fitted file, and the calls made with it). *Depends:* C6, B4. *Source:* arch §2.2, §2.4.

> **Checkpoint C:** `estimate-parameters` runs the new fit end to end on tomato2; the oracle's move
> is written down. Pause for review.

**Before C8: the milestone C review's two defects in the errors are fixed in a follow-up commit**
(approved by the owner 2026-10-15): standard errors far too small at motif lengths with stiff
curves, yet written as estimated (finding B1, the review's blocker), and a failed factorisation in
`newton` that returns the starting values silently, which the choice of strength then scores
(finding M2). They land before C8 so that C8's oracle move is the prior's alone. The review is
`doc/devel/reports/reviews/str_penalised_fit_milestone_C_2026-10-15.md` on this branch, not committed.

- ✅ **C8** *(approved by the owner 2026-10-15, spec §2.10 and its answers in §2.10.9.)* **Own commit,
  do not bundle.** The stutter prior: `stutter_prior_reads` and `stutter_prior_centre` in
  `PenalisedFitConfig`, `STUTTER_PRIOR_READS` and `STUTTER_PRIOR_CENTRE`; the prior added in the
  M-step's three stutter binomials at points with a class; its curvature in the errors' penalty; the
  two direction-trial fields on `PenalisedClass` and the new `NotPlaced` row; the run log's lines;
  the sentence in the file's note.
  *Depends:* C7, and the follow-up commit fixing the milestone C review's findings B1 and M2 (above). *Source:* spec §2.10; arch §1.1a, §1.4.
  *Tests (spec §4.1):*
  - a one-class motif length with no slipped read returns 0.1/(n + 1), 0.5 and 0.05, its shorter
    share and fall-off errors `NotPlaced` (today the ridge gives about 10⁻⁶);
  - one slipped read, one repeat shorter, in a one-class motif length: shorter share 0.75, fall-off
    0.025;
  - *(owner, 2026-10-15, after C8's review; spec §2.10.4, arch §1.4)* **the slip rate of a class with
    fewer expected slipped reads than `stutter_prior_reads` is `NotPlaced`**, by the same check as
    the direction numbers, on `expected_slipped_reads`; the lone-class test with no slipped read
    (10 loci × 3 reads, where the sandwich gave a logit error of 0.65) asserts it. The doc of
    `errors_of_a_class` drops "the slip rate needs no such rule", and the file note's sentence says
    the slip rate of such a row carries no error either (spec §2.10.5);
  - the calling-side test that `STUTTER_PRIOR_CENTRE` converts to `hipstr_shipped`'s whole-repeat
    shares;
  - weight 0 reproduces the fit before this step bit for bit.
  *Simulator checks* (the bench's drawn cohorts, `bench_fixtures/cohort_draw.rs`):
  - **genuinely low stutter on many reads keeps its rate**: the `pcr_free` truth (every slip rate a
    tenth of tomato's, 5×10⁻⁴ at one-base 8 repeats) at 2,000 samples × 30 reads; every class with at
    least 8 loci holds the truth within two of its errors, and its rate moves from the weight-0 fit by
    less than one part in a thousand;
  - **a thin class no longer drifts to zero**: a drawn motif length with a tail of classes of 1 to 3
    loci, at 1 sample × 3 reads and 4 samples × 5 reads (the oracle's shape); no stutter logit within
    one unit of ±40, and, since every true rate there is above 0.005, no class's slip rate below a tenth of `0.1/(n + 1)` for its n reads.
  - **the weight (spec §2.10.9, question 1)**: held-out prediction with `fit_bench.rs`'s scorer at
    weights 0, 1 and 5 on the cells 1 × 3, 10 × 3, 60 × 3 and 500 × 30 reads, one seed each. The
    owner chose one read; this sweep checks it. If 0 or 5 predicts better, the result goes back to
    the owner and the constant is not changed in this step. *Answered 2026-10-15:* the review re-ran
    the six thin cells at seeds 2 and 3; weight 0 beats 1 on two-base cohorts of 10 and 60 samples by
    0.2–0.5 log units at every seed, weight 1 beats 0 at one sample by 0.3–1.4, weight 5 loses
    3.2–6.1 on two-base runs at 1 and 10 samples. **The owner kept one read** (spec §2.10.2).
  *Oracle:*
  - **C0's fixture is re-recorded at this step**, because every stutter number moves: alone, a class
    moves by `a·(p₀ − p)/(n + a)` (spec §2.10.2), and the allele lengths and F follow through step 1.
    Before re-recording, compare the old fixture to the new one cell by cell and write the largest
    relative move of each number. Expected: below 10⁻³ at cells of 60 samples or more; largest at
    1 sample × 3 reads, and largest in the thin classes there. A move against that pattern is
    investigated, not re-recorded.
  - **The oracle script (`scripts/promote_ng_oracle.sh`), run before and after.** Expected: the
    fitted file's stutter rows move; no slip rate is below 10⁻⁴ unless its class has thousands of
    reads; the three-base slip rate no longer runs to logit −40; the shorter share and fall-off at
    three bases × 6 repeats and four bases × 7 come off 0.99998 / 3.5×10⁻⁵ to about 0.75 / 0.025; the
    one-base classes with 8 or more loci move by under 1%, while the thin tails of the one- and two-base curves rise (the review found them resting on classes with 0.002 to 1.5 expected slipped reads); the VCF called with the fitted file changes QUAL in
    three- to six-base records, and a genotype only where a homozygote has 4 to 6 reads with one read
    a repeat off (spec §2.10.7). Write the before-and-after table of the near-zero rows that the
    milestone C review listed.
  *Two criteria above missed, and are restated in units of each number's own error* (owner,
  2026-10-15, approved with the step; C8's review, `tmp/C8_review.md` on this branch):
  - **moves at 60 samples or more**: the relative moves reached 2.8×10⁻³ on the slip rate and
    5.8×10⁻³ on the direction numbers, against the 10⁻³ written above, which is the wrong scale for a
    rate of 5×10⁻⁴. Measured against each number's error, over every cell of 60 samples or more and
    every seed, the largest move is 0.023 of an error for the slip rate, 0.23 for the shorter share
    and 0.086 for the fall-off. The criterion becomes: under a tenth of the number's own error for the
    slip rate, a quarter for the direction numbers;
  - **low stutter, "every class within two errors"**: 3 of 35 classes fell outside, at weight 0 and
    weight 1 alike, and the prior moved those three by 1.2–2.7 parts in a million. Honest errors
    put about 1.6 of 35 outside two errors, and 3 or more happens about 22 times in 100. The criterion
    becomes: about 1 class in 20 outside two errors, and the prior moves no class by more than 10⁻³.

  **What C8 changes in milestone D.** D1 compares with `compare_D_prior_2026-10-14.md`, measured
  without the prior: the comparison is against C8's re-recorded fixture, and a cell that differs from
  the report beyond run-to-run noise is explained by the prior or investigated. D2, D4, D5 and D6 run
  after C8, so they measure the fit with the prior; D4 re-records the oracle from C10's state, which includes C8's prior. D6's
  answer is not expected to move: at 500 and 2,000 samples × 3 reads the prior shifts a rich class's
  slip rate by about one part in 100,000. D3 (cost) is unaffected.

- ✅ **C9** *(owner, 2026-10-15.)* **Own commit, do not bundle; oracle before and after.** The stated
  F's median counts only classes with an F and at least `LOCI_FOR_AN_ERROR` (8) loci with reads; if
  no class has 8, every class with an F, as before. `StratumFits::over` collects each fitted class's
  `StratumFit::tracts_of_its_own` beside its concentration and hands both to `median_concentration`
  (`stratum_fits.rs:1250`); the warrant stays `FittedHere` whenever any class has an F and
  `Defaulted` otherwise, its code unchanged. *Depends:* C7. *Source:* companion §3.5; arch §1.5.
  *Tests* (in `stratum_fits.rs`, beside `the_stated_concentration_is_the_median_of_the_runs_own_fitted_ones`):
  - a cohort where thin classes would move the median: three classes of 8 or more loci and four
    classes of fewer than 8 with higher F; the median is the three thick classes' alone, and the
    warrant `FittedHere`;
  - no class with 8 loci: the median is taken over every class with an F, `FittedHere`;
  - one sample (no class has an F): 0.5, `Defaulted`, as before.
  The existing median tests build their classes with the `fitted` helper, 40 loci each
  (`stratum_fits.rs:1833`, `:1853`), and should hold unchanged.
  *Oracle:* the oracle script, run before and after. Expected: only the parameters file's
  `fallback_between_locus_share` moves (on tomato2's oracle the milestone C review computed 0.615
  over the 15 classes of 8 or more loci, against 0.723 over all 46; C8 may move both, so the exact
  figure is the one after C8), and with it calls at classes without an F and at tracts with no row,
  which take the fallback. Any other move is investigated before re-recording.

- ✅ **C10** *(owner, 2026-10-15.)* **Own commit, do not bundle; oracle before and after.** The cohort
  in one fixed order: in `open_the_censuses` (`src/cli/estimate_parameters.rs`), before anything is
  summed, samples sorted by the sample name in each psp's header, and read groups within a sample by
  their declared ID. The fitted file then lists samples in name order, which `parameters_file.md` §6
  allows. *Depends:* C9. *Source:* spec §3 (the order trap); arch §2.4a;
  `parameter_prepass_cohort.md` test 7 ("The gather is order-independent"); `parameters_file.md` §6.
  *Why:* the same four tomato2 psps fitted in the oracle's (VCF) order and in alphabetical order give
  different files, md5 `4c4efa9…` against `5c53b9dc…`; the first difference is in the SNP/indel fit,
  which stops at its pass limit; inbreeding 0.96645 against 0.96617; a two-base fall-off strength of
  5.7×10⁹ against 1.2. Thread count is not the cause (`tmp/thread_determinism_report.md`).
  **The defect is on main too**, so it is fixed on main when this branch merges.
  *Test:* one cohort fitted from two orderings of its psps gives identical bytes.
  *Oracle:* the oracle script before and after. Expected: every fitted number of the fitted file
  moves, since the oracle names its psps in the VCF's order and the fit now sums in name order, and
  with them the calls made with that file; the VCFs called with default parameters do not move. The
  move is re-recorded at D4.

### Milestone D — checks (spec §4)

Every D step runs after C10 *(owner, 2026-10-15)*, so each measures the fit as it will merge.

- ✅ **D1** The bench reproduces `compare_D_prior_2026-10-14.md`'s `D_tract_prior` rows with the
  production fit, within run-to-run noise, and adds the cells spec §1.2 names as not run: 500 × 30,
  2,000 × 10, and three seeds at every 2,000-sample cell. It measures error coverage on classes with 1
  to 7 loci, for the owner's decision on `LOCI_FOR_AN_ERROR`. **It also counts, at one sample with the
  stutter prior on, the two-base classes whose direction numbers lie more than two errors from the
  truth** (spec §2.6; C8's review found 31/63, 28/126 and 21/126 against 18/62, 18/126 and 12/126
  without the prior, where honest errors miss about 1 in 22), for the owner. *Depends:* C7.
  *Source:* spec §1.2, §2.6, §4.2.
  *Done (2026-10-16):* `reports/implementations/str_penalised_fit_D1_2026-10-15.md`. The production
  fit predicts held-out loci better than the counting fit in all 48 cells (+0.009 to +0.162 log units
  a locus; 167 of 168 seed-runs). Recommendation for checkpoint D: stutter-number errors from 3 loci
  with reads, F and shares keep 8, `LOCI_FOR_A_SHARE_ERROR` from 3 to 1. The runs that lost every
  stutter error are diagnosed separately (`str_penalised_fit_thin_errors_diagnosis_2026-10-15.md`).
- ✅ **D2a** The bench gains a real-psp mode: it reads a psp set's censuses, folds loci by whole locus
  over the real evidence, fits on the kept folds and scores the held-out ones with `fit_bench.rs`'s
  scorer. *Depends:* C7. *Source:* spec §4.4.
  *Done (2026-10-16):* `--psp`, `--reference`, `--catalog`, `[--parameters]`; the evidence comes from
  `estimate-parameters`' own path (`open_the_censuses` … `the_tract_strata_of_a_cohort`), checked
  bit-identical to the oracle's fitted slip numbers on its 4 psps. The held-out loop moved into
  `fit_bench.rs` with a test that every locus is scored by the fit made without it. No full-cohort
  fit in this mode (no truth). Review: `reports/reviews/str_penalised_fit_D2a_2026-10-16.md`.
- ☐ **D2** tomato2, held-out by whole locus: the production fit against the counting fit, with D2a's
  mode. *Depends:* D2a. *Source:* spec §4.4.
- ✅ **D3** Full scale, simulated: 2,169 samples with kimura's loci per class (about 66 million
  pairs with reads), the production fit's wall time, thread scaling and peak memory on this Mac.
  Within one hour, with a memory watchdog. *Depends:* C7. *Source:* spec §2.2, §4.5 (kimura's real
  run deferred).
  *Done (2026-10-16):* `reports/implementations/str_penalised_fit_D3_2026-10-16.md`. 83 to 173 s at
  18 threads, 151 to 252 s at 6, peak 10,801 MiB, at 66.2 million pairs; kimura's per-class spread
  assumed from the tomato catalog and kimura's totals (two variants).
- ☐ **D4** The oracle and `cli::cross_platform_digests` re-recorded, each move explained, from
  C10's state (the cohort in name order). *Depends:* C7, C10. *Source:* spec §4.6.
- ✅ **D5** Calling on HG002 and the GIAB trio, **run once with a parameters file fitted by the climb
  and once with one fitted by the new fit**, on the same psps; nothing else differs.
  - HG002: `scripts/ng_str_accuracy/hg002_two_binaries.sh`, scored by
    `benchmarks/lib/tract_qual_experiment.py` against GIAB's tandem-repeat truth; the number is right
    genotypes by period class and depth (`genotype.tsv`), 30× and 50×, as in
    `doc/devel/reports/implementations/str_partial_read_phase_2026-10-12.md` §4. The script calls with
    `--defaults`; each arm passes `--parameters` instead.
  - The trio: `benchmarks/giab/src/run_ng_per_sample.sh` and `score_ng_recall.sh`, as in
    `doc/devel/reports/ng_str_routing_recovery_2026-09-02.md`.
  - Data: the main checkout's `benchmarks/` (preconditions).
  *Depends:* C7. *Source:* spec §4.7.
  *Done (2026-10-16), for the owner to judge:* `reports/implementations/str_penalised_fit_D5_2026-10-16.md`.
  HG002: within 10 tracts of the climb's right genotypes in every class and depth; trio: SNPs identical,
  one indel fewer at each depth. Six one-sample fits per arm (the trio's HG002 reused the tier's fit).
- ✅ **D6** Does a slip rate 0.5% of its value low change calls? On the bench's read counts, at 500 and
  2,000 samples × 3 reads: each sample's genotype call from its genotype weights under the fitted
  stutter against under the true stutter. **A proxy for calling**: nothing in the project simulates
  reads for that many samples. *Depends:* C7. *Source:* spec §5.1.
  *Done (2026-10-16):* `reports/implementations/str_penalised_fit_D6_2026-10-16.md`. A slip rate 0.5%
  low changes about 4 calls in 100,000 (500 and 2,000 × 3), both directions alike, against 5–6 in 100
  already wrong; the fitted stutter changes 2 to 5 in 10,000.

> **Checkpoint D:** a short report under `doc/devel/reports/implementations/` with every D result.
> Pause for review: the owner decides whether D2 and D5 pass before the climb is deleted, and whether
> D1's thin-class coverage lifts `LOCI_FOR_AN_ERROR`.

### Milestone E — delete the climb

- ☐ **E1** Everything arch §4.1 marks deleted, in one commit: the climb's code, tests and examples;
  `fit_period_allele_length_variations` (tests that used its output build their tables by hand);
  `ssr_fit/library_stutter.rs`, `genotype_rows.rs`, `neighbour_start.rs` and their examples; the
  bench's methods `A` and `B_curves`; `SsrFitConfig` (examples that name it only for the window
  ported); `StratumRefusal::BelowTheFloor`. Kept: what older files need to read. *Depends:* D2, D5 and
  the owner's go-ahead. *Source:* arch §4.1; spec §7.
- ☐ **E2** The specs the climb built on get a dated note that the fit they describe was replaced
  (`parameter_prepass_ssr.md`, `str_slippage_level_curve.md`), pointing at the new spec.
  *Depends:* E1. *Source:* spec-authoring's rule for decisions that were acted on.

> **Checkpoint E:** full suite and code review. Pause for review before any merge into main.

## Verification summary

| milestone | proven by |
|---|---|
| A | oracle checksums unchanged (A5's move explained); new-rung and warrant tests |
| B | conversion, round-trip and migration tests; golden-file diffs; oracle as checkpoint B states |
| C | C0's fixture reproduced (C2, C3, C4); determinism tests, threads and psp order (C3, C10); peak memory; oracle move explained (C7, C8, C9, C10) |
| D | the bench, tomato2 held-out, simulated full scale, HG002 and GIAB calling, the D6 proxy |
| E | full suite and code review |

## Out of scope (next plans)

Slippage groups and library multipliers (`str_library_stutter.md`); per-period census reading;
kimura's full-scale run; the read-crossing bias for long alleles (`str_allele_window.md` §5.3).
