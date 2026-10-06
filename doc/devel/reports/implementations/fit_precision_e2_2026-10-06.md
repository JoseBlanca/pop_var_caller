# Fit precision, step E2 — the parameters file, version 2: each number's standard error in the file

**Date:** 2026-10-06. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step E2.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §5.2–5.3, §8 question 1. **Branch:** `fit-precision`.
Review: [fit_precision_e2_2026-10-06.md](../reviews/fit_precision_e2_2026-10-06.md), fixes:
[fixes_applied_fit_precision_e2_2026-10-06.md](../reviews/fixes_applied_fit_precision_e2_2026-10-06.md).

## 1. What was built

The parameters file a run writes goes from format version 1 to version 2. Every key version 2 adds is optional, so
this build reads a version-1 file as the same numbers without errors; a version-1 build refuses a version-2 file with
its existing message.

| what the file now says | where | from |
|---|---|---|
| `standard_error` beside `value` | every warranted number: a library's base-quality multiplier, a sample's inbreeding coefficient, a repeat-tract substitution rate | the multiplier: the library's rate error over its mean claimed error; the coefficient: the SNP/indel fit's error on the sample's homozygote excess; the substitution rate: the binomial error of the count (step E1) |
| `own_fit_standard_error` | a slippage row's level origin, **only where the level is the stratum's own fit** | the stratum's own fit (plan step C1) |
| `shorter_share_own_fit_standard_error`, `fall_off_own_fit_standard_error` | a slippage row's shares origin, each only where that share is the stratum's own fit | the same |
| `samples_fitted_on` | a slippage row, where a large cohort's stratum was fitted on a subset of samples | `StratumFit::samples_fitted_on` (plan step D2) |
| `snp_indel_fit_starts` | `[fitted_from]`: each start's number, how it ended (`converged`, `at_the_pass_limit`, `agreed_with_an_earlier_start` with which), its passes | `JointFit::starts` (plan step B3) |

- **Absent is not zero**: a number with no error writes no key. A `defaulted` number never carries one, by the rule
  it carries no count; a `supplied` number keeps the error it came with, as it keeps its count.
- **The file's notes say what each key is**, and, as the owner asked before this milestone, that where a row
  carries `samples_fitted_on` its errors describe the subset's own fit, not its distance from a fit on every sample;
  and that its `expected_slipped_reads` are among the samples the stratum was fitted on.
- **`validate` checks the new keys**: an error is finite and above zero; none sits beside a `defaulted` number; an
  own-fit error sits only beside a number that is the stratum's own fit; `samples_fitted_on` is at least one; the
  starts are numbered from one in order, and only a start that agreed names an earlier one.
- **The errors survive a round trip.** `ReadGroupCalibration` carries the multiplier's error
  (`scale_standard_error`); `FittedSlippage` carries each cell's own-fit errors and samples, built in
  `StratumFits::over` from the own fit's errors, kept only for numbers whose provenance says the stratum's own; a
  calling run carries the read file's starts back into the file it writes (`TheRunsNumbers::snp_indel_fit_starts`,
  `ParametersFile::with_snp_indel_fit_starts`).

### A number with no information is `defaulted` (spec question 1; the owner at checkpoint A′)

- **A read group whose ordinary-position rate has no information** — a library with no reads, or one whose sample's
  census holds repeat tracts only, which the fit lists with nothing scored — comes out `Defaulted` at the stated rate
  (one error in a thousand), from no observation, with no error. Before, it was `fitted_here` at the start the fit
  never moved, counted over the run's every position.
- **A sample whose homozygote excess has no information** comes out `Defaulted` at no excess, from no observation; its
  inbreeding coefficient is then written `defaulted` at 0.0. The default value is written rather than the fit's: no
  read informs the excess, yet the fit's accelerated steps can move it off its start.
- **The excess is fitted only across samples that inform it** (found by the review). The rule "two samples or more"
  counted every sample, so one sample with reads beside a tracts-only sample was fitted as a two-sample cohort and its
  excess written `fitted_here` with an error, where alone it is held. It now counts samples with a read at an ordinary
  position, at the fit, and samples whose own excess carries information, in the error code — the six places that ask.

### A defect this exposed, fixed here

**A cohort holding a sample whose census has repeat tracts only stopped `estimate-parameters` with a panic.** The
census sums what a library's qualities claimed at ordinary positions only
([calibration.rs:287](../../../../src/parameter_estimation/calibration.rs#L287)), so such a read group has no totals;
the fit still lists it with a rate, and assembly refuses a rate without totals
([run_parameters.rs:697](../../../../src/calling/run_parameters.rs#L697)). `parameters_from_the_fit` now gives it
empty totals, which have no mean, so it takes the defaulted multiplier of one — the documented "accumulator saw no
read" case. **Only a read group whose rate is defaulted from no observation gets them**: a fitted rate without totals
is a lost read-group axis, and assembly still refuses it. A read group with totals and a defaulted rate is unchanged: by the owner's ruling of 2026-08-31 it is
charged one error in a thousand over its own claimed mean, marked `defaulted`.

### Departures from the spec's wording, recorded

- **The shares origin has two own-fit keys**, not one `own_fit_standard_error`: it covers two numbers, smoothed
  separately, and either can be the stratum's own while the other comes from a curve.
- **`samples_fitted_on` is on the slippage row**, beside its two origin blocks, whose counts it qualifies.
- **The starts are written by a builder beside `ParametersFile::of_run`**, not by an eighth argument to it: only the
  fit that made them and a calling run writing back a file it read have any to say. Both calling commands write their
  file through one door, `TheRunsNumbers::parameters_file`, which adds them, so neither can forget.
- **The file's type for how a start ended is `StartOutcome`**, a noun distinct from the fit's `StartEnding`; the file
  spells its three values `converged`, `at_the_pass_limit`, `agreed_with_an_earlier_start`.

## 2. What it does to the files, measured

**The four-accession oracle cohort** (`scripts/promote_ng_oracle.sh`, 20 regions; `tmp/fit_precision/oracle_e2b/`, on
the code as committed, against `oracle_d3/`, the run before this milestone, which step E1 did not change; the run on
the reviewed code, `oracle_e2/`, differs from it only in note lines):

- every call checksum is unchanged, the calls made with the fit included;
- the two parameters files moved, and **with every key and note version 2 added taken out each is the old file line
  for line** (`tmp/fit_precision/e2/strip_v2.py`: 340 of 340 lines in the fitted file, 116 of 116 in the defaults
  run's);
- the fitted file gains 144 standard errors — 4 multipliers, 4 inbreeding coefficients, 136 substitution rates; the
  other 52 substitution rates are exactly zero and carry none — one own-fit error, on the one stratum of 37 whose level
  is its own fit (the other 36 levels and all 74 shares come from a curve or a blend), and three starts, all
  `at_the_pass_limit` (198, 198 and 200 passes). Four samples draw no subset, so no row carries `samples_fitted_on`.

**The cross-platform fixture** ([cross_platform_digests.rs](../../../../src/cli/cross_platform_digests.rs)): the same
taking-out turns its new fitted file back into the bytes of the recorded `22dda760d260572b786958e69b15b2b2`; the calls
checksum does not move. The file as committed differs from that one in 19 note lines and nothing else. Re-recorded
`a5e87936e968f9b855f91f88db828083`, and the oracle baseline's two parameters-file lines (`70a98909…`, `351f244e…`).
**These pin the format and the notes, not the new errors**: the fixture's rates sit on their bounds and its tracts
never mismatch, so its file carries no `standard_error` (the review's M3, for the owner).

## 3. Tests

| test | what it shows |
|---|---|
| `a_sample_with_no_ordinary_positions_and_its_read_group_come_out_defaulted` (fit) | four drawn samples and a fifth holding one repeat-tract section: its read group `Defaulted` at 0.001 from 0 observations with no error, its excess `Defaulted` at 0 from 0; the four keep fitted numbers with errors |
| `its_read_group_and_its_coefficient_are_written_defaulted` (run) | the same cohort, the four drawn samples given claimed-error totals, through assembly and the file: no panic; the four multipliers fitted with errors, the fifth `defaulted` 1.0 with no count or error; its coefficient `defaulted` 0.0; the others `fitted_here` with errors; the file reads back and validates |
| `one_sample_with_reads_beside_a_tracts_only_sample_holds_its_excess` | one drawn sample alone and beside a tracts-only sample: its excess held, `defaulted` with no error, both times |
| `one_sample_informing_its_excess_holds_it` | in the error code, two samples of which one informs its excess: that excess is held |
| `a_library_with_no_reads_beside_one_with_reads_comes_out_defaulted` | three samples of two libraries, one library drawn at no depth: it alone is `defaulted`; its sibling and its sample keep fitted numbers with errors |
| `each_ending_is_written_with_its_number_and_passes` | the fit's three endings, `Agreed` with its start, as the file writes them |
| `a_run_writes_back_the_starts_its_file_carried`, `a_calling_run_writes_back_the_fit_starts_it_read` | a calling run handed a file with starts writes them back, through the whole psp command and through the one door |
| `own_fit_errors_are_kept_only_for_the_numbers_the_stratum_emits_as_its_own` | of an own fit's three errors per group, only the cell's level and the stratum's own fall-off are kept; a blended level and curve shares keep none; a derived stratum has none and no sample count; the subset's 300 samples travel with each group |
| `the_multipliers_error_is_the_rates_over_the_claimed_mean` | the multiplier's error is the rate's over the mean, the same share of the multiplier as of the rate; none from a rate with none |
| `a_defaulted_number_is_written_with_no_error_and_a_supplied_one_keeps_its_own` | the writer's rule |
| `every_key_version_two_added_is_refused_where_it_means_nothing` | each check in `validate`, by key — and no refusal message carries a run of spaces |
| `a_file_a_version_one_build_wrote_reads_back_without_the_new_keys` | the last version-1 golden file, kept as `testdata/version_1_as_written.toml`, reads as today's fixture without the new keys, and validates |
| `a_version_one_file_is_accepted`, `a_version_this_build_does_not_read_is_refused` | 1 and 2 read, 0 and 3 refused |
| `the_fits_defaulted_excess_is_the_runs_default_coefficient` | the fit's 0.0 and calling's 0.0 are one number |
| the golden files and every round trip (extended) | each new key written, read and written back; the fixture carries one of each, and absent ones where a number has none |

Five mutations, each run against its test and restored (`tmp/fit_precision/e2/mut/`), all killed: the empty totals
reverted (the run panics in assembly), the no-information test forced true, the level's own-fit filter removed,
`validate`'s own-fit check disabled, the writer's defaulted-error rule removed.

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- Full suite (`cargo test --all-targets --all-features --no-fail-fast`): 5,040 passed, 3 failed (the pre-existing
  `examples/ng_generic_loci_dump.rs` × 2 and `examples/ng_ssr_loci_dump.rs` × 1), 9 ignored.
