# Fit precision, step A6 (second part) — each library's rates get their own standard errors

**Date:** 2026-09-28. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step A6,
second part; the first part is [fit_precision_a6a](fit_precision_a6a_2026-09-28.md). **Source:**
checkpoint A's decision 1; spec [§3.2](../../ng/spec/fit_precision.md)'s per-library block.
**Branch:** `fit-precision`.

## 1. What was built

The first part made the likelihood read each library under its own two error rates. This part gives
those rates their slopes, their share of the information and their standard errors, so a sample read
from `k` libraries now carries `1 + 2k` parameters of its own in the errors, as spec §3.2 sizes it.
**No fitted number moves**: the scores, the information and the errors are computed from the fit and
never fed back into it.

- **The slope of a library's rate at a position** is its own reads' slope in that rate, weighted by
  the posterior of the genotype its sample's libraries share (`fill_read_slopes` now fills one entry
  a library, and the three branches' scorers credit each library's slot).
- **A sample's row of its own parameters** is laid out as its first library's two rates, its
  homozygote excess, then each further library's two rates (`information::sample`,
  `sample::rate(library, class)`). A sample of one library keeps exactly the three slots it had, so
  its scores, blocks and errors are computed with the same operations as before.
- **The information blocks and the errors hold a row of each sample's own length**
  (`PositionScores::samples`, `InformationSums::sample_blocks` and `sample_cohort_blocks`,
  `StandardErrors::samples` are vectors of rows). The arrow inversion was already written over the
  parameters a block keeps, so it needed only the block's size read from the sums.
- **The interim of the first part is gone**: the reason `AwaitingEachLibrarysScores`, the call that
  marked every error of a cohort holding a sample of several libraries, and its log clause.
- **The log line counts the rates a read group each**: "error rates at ordinary positions (N read
  groups)", where it said "(N samples)".
- **The trace gives each read group its own two errors**, from its own slots.

Memory: a sample's blocks are `8(n² + 8n)` bytes for `n` own parameters — 264 bytes for one library,
as before, 520 for two.

## 2. What moves

Nothing fitted. For a cohort of one library a sample, the standard errors are the same bits as the
first part's: see §6, the oracle. For a cohort holding a sample of several libraries, errors now
exist where the first part reported none.

## 3. Deviations and assumptions

1. **The row layout puts the first library before the excess** and further libraries after it,
   rather than all rates before the excess. It keeps a one-library sample's row and every existing
   test's slot constants unchanged; `sample::rate` hides the order from callers.
2. **A sample holding no library keeps a row of three**, its two rate slots at zero, so the excess
   stays in its slot; its rates say they have no information.
3. **One test's parameters were moved off a flat spot.** In
   `a_samples_libraries_are_scored_each_under_its_own_rates`, sample 2's mismapped rate sat where
   its slope is −0.038, and a central difference at a step of 10⁻⁵ of the rate cannot resolve that
   from the likelihood's rounding (disagreement 5.7 × 10⁻⁵; 5.1 × 10⁻⁶ at a step of 10⁻⁴, as rounding
   predicts; 9.2 × 10⁻⁴ at 10⁻³, where the step's own curvature takes over). At 0.03 the slope is
   112 and agrees to 1.6 × 10⁻⁸ (`tmp/fit_precision/a6b_tests.log` for the first run,
   `a6b_fd_step.log`, `a6b_fd_step3.log`, `a6b_fd_final.log`).
4. **The cost, measured in review**: the final pass, which collects the information, went from
   1.53–1.56 to 1.69–1.71 times a plain pass for one library a sample (about 10% longer), and from
   1.37–1.45 to 1.71–1.79 times for two (the first part computed no slopes for them). The final pass
   runs three times a fit, against about 600 plain passes on the oracle.
5. **The first review commit of this step held `information.rs` double-encoded** (every —, §, ×
   stored garbled, 117 lines) by a perl edit that wrote a wide character; all three reviewers found
   it. The file was decoded once, keeping the one character written correctly after the damage.

## 4. Changes

- [fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs): the layout
  (`ONE_LIBRARY_SAMPLE_PARAMETERS`, `own_parameters`, `sample::rate`, `sample::class_of`);
  `PositionScores` over rows and a read-slope table over libraries; `fill_read_slopes` and the
  three branch scorers per library; `InformationSums` over rows; the module doc; tests.
- [fit/standard_errors.rs](../../../../src/parameter_estimation/joint/fit/standard_errors.rs): rows
  of each sample's own length; bounds by slot (`bounds_of_own`); `named` and `described` per read
  group; the interim reason and marking removed; the module doc; tests.
- [fit.rs](../../../../src/parameter_estimation/joint/fit.rs): the pass builds the rows from
  `group_index`; `fit_jointly` no longer marks the errors; `JointFit::standard_errors`' doc; the
  two-library whole-fit test.

## 5. Tests

| test | what it shows |
|---|---|
| `a_samples_libraries_are_scored_each_under_its_own_rates` (information.rs), the finite-difference test per library | Sample 0's reads split over two libraries at rates unlike each other's: every library's two rate slopes and every sample's excess slope against a central difference of the log-likelihood — 14 slopes, the largest disagreement 6.8 × 10⁻⁷ relative (a homozygote excess), every rate within 2.4 × 10⁻⁸. |
| `each_librarys_slopes_hold_where_its_depths_are_ranges` (information.rs, new) | Two samples, each of a 40-read and a 132-read library under a cap of 140: every library's two slopes and each excess against a central difference, the largest disagreement 1.6 × 10⁻⁷. |
| `the_blocks_give_the_errors_of_an_arrow_of_mixed_libraries` (standard_errors.rs, new) | Samples of one, two, three and one library: every error equals the dense inverse's to 10⁻¹² relative, and each sample's row has its own length (3, 5, 7, 3). |
| `a_samples_two_libraries_come_back_at_their_own_rates` (fit.rs) | Now: each of the twelve libraries' clean rates lies within three of its own standard errors of the rate it was drawn at (from −1.91 to +1.27 errors); and, from review, within each sample the library drawn at the higher rate has the larger error, for both rates — which an error read from a sibling's slot cannot satisfy. |
| `a_sample_of_three_libraries_has_every_slope_right` (information.rs, from review) | One sample of one library beside one of three, the duplicated class fitted: every library's two slopes and each excess against a central difference, the largest disagreement 4.4 × 10⁻⁸ (3.1 × 10⁻⁸ on the three-library sample). Dropping the third library from the segregating branch fails it. |
| `a_sample_holding_no_library_carries_no_information` (information.rs, from review) | A sample with no library: the likelihood bit-for-bit the one where its library holds no read, its row zeros, its three errors "no information", and the log counting three read groups. |
| `each_slot_of_a_samples_row_is_its_own_parameter_with_its_own_bounds` (standard_errors.rs, from review), and the width check added to the mixed-arrow test | Samples of one to four libraries: every rate's slot is distinct, names its own class and carries its class's bounds; through the inversion, a second library's rate whose error is 0.3 is an error as a mismapped rate and wider than its range as a clean one. |
| `a_pass_sums_the_scores_multiplied_pairwise` (information.rs), extended in review | Also on samples of one, two and three libraries (blocks 3, 5, 7): a pass's sums against a by-hand account, the worst entry within 1.2 × 10⁻¹⁵ of its bound. |
| `the_logged_line_summarises_each_kind`, `the_trace_names_each_error_after_its_value` (standard_errors.rs) | Rewritten: a second library's two errors come from its own slots, into the log's counts a read group each and into the trace under its read group. |

`cargo test --release --lib parameter_estimation::joint::fit`: 68 passed, 1 ignored (63 at the first
part; five added). The coverage of each library's errors was measured in review over 240 drawn
cohorts: 65 to 69 in 100 within one error and 93 to 96 within two (review report §2).

## 6. Validation

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- `cargo test --all-targets --all-features --no-fail-fast` (`tmp/fit_precision/suite_a6b.log`):
  4,923 passed, 3 failed, 5 ignored; the three failures are the pre-existing ones
  (`examples/ng_generic_loci_dump.rs` × 2, `examples/ng_ssr_loci_dump.rs` × 1). Both
  `cli::cross_platform_digests` tests pass unchanged.
- **The oracle** (four tomato accessions, one library each; `tmp/fit_precision/oracle_a6b/`): all
  seven checksums equal the first part's, and the logged line of standard errors is the first part's
  to every printed digit — the only difference is the rates' count reading "(4 read groups)" where
  it read "(4 samples)".
