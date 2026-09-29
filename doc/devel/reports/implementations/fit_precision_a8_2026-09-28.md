# Fit precision, step A8 — the whole information matrix at small cohorts

**Date:** 2026-09-28. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step A8.
**Source:** checkpoint A's decision 3. **Branch:** `fit-precision`. **Review:**
[fit_precision_a8_2026-09-28.md](../reviews/fit_precision_a8_2026-09-28.md); **fixes:**
[fixes_applied_fit_precision_a8_2026-09-28.md](../reviews/fixes_applied_fit_precision_a8_2026-09-28.md).

## 1. What was built

The SNP/indel fit's standard errors came from the information summed in blocks: the cohort's eight
parameters with each other, each sample's own parameters with each other, and each sample's with the
cohort's. **What the blocks leave out is the pairing of two samples' parameters**, and at small cohorts
that matters: at 4 samples and 3 reads a position the estimates scattered 1.23 to 1.67 times the
blocks' errors, and 0.89 to 1.02 times the errors of the whole matrix
([A4 report](fit_precision_a4_2026-09-27.md)).

**For a cohort of at most 20 samples (`FULL_MATRIX_SAMPLES`) the fit's final pass now also sums the
whole matrix** — every parameter's score times every other's, over the positions — and the errors come
from that matrix inverted at once. Above 20 samples, nothing changes: the blocks alone, as before.

- **`FullInformation`** (`fit/information.rs`) holds the matrix's upper triangle, row after row, in the
  layout the errors are reported in: the cohort's eight, then each sample's own parameters in turn — a
  sample of k libraries has 1 + 2k (step A6). A position adds its scores' products to it, skipping a
  row whose score is zero; chunks' sums are added in the pass's fixed order, like every other count.
  `InformationSums` carries it as `full`, set by `InformationSums::for_a_cohort_of`, which is what the
  pass asks for; `InformationSums::new` stays blocks-only.
- **`StandardErrors::of`** reads the whole matrix when the sums hold it (`of_the_whole_matrix`) and
  the blocks otherwise (`of_the_blocks`, the former body). The whole matrix is inverted over every
  parameter with information that the fit moves, **taken in the order the blocks take them** — each
  sample's own, the samples in turn, then the cohort's — so a parameter the others mimic is dropped by
  the same rule, and on an arrow (no two samples paired) the two ways give the same errors. The
  reasons an error is absent are the blocks' four: no information, held fixed, not identified, wider
  than the parameter's range.

**No fitted number moves.** Only the final pass, which already summed the blocks, sums more; the
parameters it reports are unchanged. The cross-platform checksums pass unchanged.

## 2. What was measured

### 2.1 The errors mean what they say, now at 2 and 4 samples too

A4's coverage test, 200 drawn cohorts a regime, now reading the errors the fit itself reports — the
whole matrix's in every regime, since none has more than 20 samples — beside the blocks' errors for the
same cohorts (`tmp/fit_precision/a8_coverage_full.log`, 1,843 s). The test also sums the whole matrix
itself and requires every reported error to agree with it to 10⁻⁶ relative; it did, for every estimate.
A share within one error of about 0.68, within two of about 0.95, and a spread of about 1 mean the
error is right:

| regime | kind | within one error: fit / blocks | within two | spread in errors: fit / blocks |
|---|---|---|---|---|
| 2 samples, 3 reads | clean error rates | 0.823 / 0.530 | 0.940 / 0.838 | 0.84 / 1.39 |
| | mismapped error rates | 0.892 / 0.641 | 0.993 / 0.862 | 0.63 / 1.44 |
| | mismapped share | 0.804 / 0.529 | 0.920 / 0.834 | 1.44 / 1.77 |
| 4 samples, 3 reads | clean error rates | 0.741 / 0.591 | 0.978 / 0.899 | 0.89 / 1.23 |
| | mismapped error rates | 0.738 / 0.590 | 0.960 / 0.877 | 0.94 / 1.45 |
| | mismapped share | 0.775 / 0.455 | 0.950 / 0.790 | 1.04 / 1.67 |
| 20 samples, 3 reads | mismapped error rates | 0.729 / 0.678 | 0.955 / 0.939 | 0.96 / 1.06 |
| | mismapped share | 0.650 / 0.620 | 0.970 / 0.940 | 1.00 / 1.10 |

**At 4 samples and 3 reads the errors the fit reports now mean what they say** for every kind: 0.655
to 0.810 within one error and 0.950 to 0.985 within two, where the blocks' put the mismapped share
within one only 455 times in 1,000. **At 2 samples they err the other way**: the error rates' errors
are now too wide (spreads of 0.63 and 0.84), where the blocks' were too narrow (1.39 and 1.44) — an
over-cautious error in place of a confident wrong one. At 4 samples and 30 reads and at 20 samples
the two differ by at most 0.051 in the share within one error (the mismapped rates at 20 samples and
3 reads). The homozygote excess, the density's shapes, and the invariant and fixed shares barely
move: at 2 samples 0.618 against 0.623 for the excess, and at 4 and 20 samples the shares and shapes
by at most 0.005.

Every regime's other figures — fits converged, the likelihood's flat point, the shapes' distance
from the truth — are the same as step A7's, as they must be: no fitted number moved.

### 2.2 What the whole matrix costs

Single thread, median of 9 interleaved runs of each of three passes over the same positions: one that
sums nothing for the errors, one that also scores each position and sums the blocks (the final pass
before this step), and one that sums the blocks and the whole matrix (the final pass now, up to 20
samples); `tmp/fit_precision/a8_timing.log`:

| cohort | plain pass | with the blocks | with the blocks and the whole matrix |
|---|---|---|---|
| 20 samples × 20,000 positions, 3 reads (68 parameters) | 98.1 ms | 154.3 ms (1.57 ×) | 160.8 ms (1.64 ×; 1.04 × the blocks') |
| 20 × 20,000, 30 reads | 104.5 ms | 206.9 ms (1.98 ×) | 209.2 ms (2.00 ×; 1.01 ×) |
| 4 × 60,000, 3 reads (20 parameters) | 87.1 ms | 128.7 ms (1.48 ×) | 130.8 ms (1.50 ×; 1.02 ×) |
| 21 × 20,000, 3 reads (blocks only) | 95.4 ms | 160.3 ms (1.68 ×) | 160.2 ms (1.68 ×; 1.00 ×) |

**The whole matrix adds 1 to 4% to the final pass at 20 samples of one library.** The pass is
dominated by computing each position's scores, not by summing their products. The pass holds the
matrix once for each of its chunks and once for their total: 129 × 2,346 numbers, 2.4 MB.

**With several libraries a sample it grows as the square of the parameters.** The correctness
reviewer's measurement, same method, on a loaded host where two identical passes differed by 1.04 to
1.05 times, and the final pass's peak memory at 8 threads:

| 20 samples of … libraries each | parameters | whole matrices, 129 copies | time × the blocks' pass |
|---|---|---|---|
| 1 | 68 | 2.4 MB | 1.04 to 1.23 at 3 reads, 1.06 at 30 |
| 4 | 188 | 18.3 MB | 1.10 to 1.12 |
| 16 | 668 | 231 MB (the pass's peak rose by 261 MB, the blocks' 28 MB included) | 1.33 to 1.64 |
| 50 | 2,028 | 2.1 GB, not run | not run |

Only the one final pass a start sums the information, so its time matters little; its memory has no
bound as libraries a sample grow (§3 item 1).

## 3. Deviations and assumptions

1. **The threshold is on samples**, as the plan states, not on parameters — so samples of many
   libraries make the matrix, and the final pass's memory, much larger (design review, Major; the
   pass holds 129 copies):

   | 20 samples of … libraries each | parameters | products a position | the pass's whole matrices |
   |---|---|---|---|
   | 1 | 68 | 2,346 | 2.4 MB |
   | 5 | 228 | 26,106 | 27 MB |
   | 16 (lane-level read groups) | 668 | 223,446 | 231 MB |
   | 50 | 2,028 | 2,057,406 | 2.1 GB |

   **Decided at checkpoint A′ (owner, 2026-09-29): the parameters are capped too**, at
   `FULL_MATRIX_PARAMETERS` = 188 — 20 samples of four libraries, 18.3 MB and 1.1 times the blocks'
   pass, measured — so a cohort of many-library samples takes the blocks rather than gigabytes.
   Built in its own commit after this step's.
2. **The blocks are still summed below the threshold**, beside the whole matrix. They cost about a
   third of the whole matrix's products at 20 one-library samples (724 against 2,346), keep the two
   paths separate, and let the tests compare the blocks' errors with the whole matrix's on the same
   pass.
3. **The coverage test now reads the fit's own errors.** Every regime of it has at most 20 samples, so
   what it tallies as the fit's coverage is the whole matrix's; the blocks' errors, which a larger
   cohort is given, are tallied beside them. The test also sums the whole matrix itself and requires
   every error the fit reports to agree with it to 10⁻⁶ relative.

## 4. Changes

- [fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs):
  `FULL_MATRIX_SAMPLES`, `FullInformation`, `InformationSums::full` and `for_a_cohort_of`; the coverage
  test reads the fit's errors (§3 item 3).
- [fit/standard_errors.rs](../../../../src/parameter_estimation/joint/fit/standard_errors.rs):
  `StandardErrors::of` chooses; `of_the_whole_matrix`, over `InvertedParameter`s; `of_the_blocks`; the
  module doc describes both ways and the order parameters are judged in.
- After review: two chunks that disagree on keeping the whole matrix are refused in every build, not
  only in debug builds; three tests from the reviewers (§5).
- [fit.rs](../../../../src/parameter_estimation/joint/fit.rs): the pass asks for
  `InformationSums::for_a_cohort_of`.

## 5. Tests

| test | what it shows |
|---|---|
| `the_whole_matrix_gives_the_errors_of_its_dense_inverse` (standard_errors.rs) | A random dense matrix, every pair informed, over samples of 1, 2, 3 and 1 libraries (26 parameters): every error equals the dense inverse's to 10⁻¹². |
| `on_an_arrow_the_whole_matrix_and_the_blocks_agree` (standard_errors.rs) | On an arrow the whole matrix and the blocks give the same errors to 10⁻¹². |
| `the_whole_matrix_says_why_a_parameter_has_no_error` (standard_errors.rs) | A sample with empty rows has no information; a share made twice another's row is dropped alone as not identified; at one sample the excess is held fixed; every other error equals the dense inverse's without those. |
| `the_whole_matrix_is_kept_up_to_twenty_samples` (information.rs) | Sums for 20 samples hold it, sized for their own parameters (a two-library sample counts five); for 21 they do not. |
| `a_pass_sums_the_scores_multiplied_pairwise` (information.rs), extended | A pass's whole matrix equals every position's scores multiplied pairwise and summed by hand, to 10⁻¹² of each entry's bound, on four one-library samples, one sample, and samples of one, two and three libraries, over three chunks. |
| `a_fit_carries_the_errors_at_the_parameters_it_returns` (information.rs), extended | A four-sample fit's errors are the whole matrix's, and differ from the blocks'. |
| `the_information_is_the_same_bits_at_any_pool_width` (information.rs), unchanged | Its three-sample cohort now keeps the whole matrix, and it prints alike at 1, 4 and 8 threads. |
| `on_an_arrow_the_whole_matrix_and_the_blocks_drop_the_same_parameters` (standard_errors.rs, from the review) | A cohort parameter made to copy a sample's, and a sample's second library's mismapped rate made to copy its clean rate: both ways drop the same two, keep the rest, and agree to 10⁻⁹. Taking the cohort's parameters first fails it. |
| `the_whole_matrix_judges_each_slot_against_its_own_interval` (standard_errors.rs, from the review) | A second library's mismapped rate whose error is 0.3 keeps it — inside its own interval, outside the clean rate's. |
| `a_sample_without_reads_has_no_row_in_the_whole_matrix` (information.rs, from the review) | Through the pass's own sums, a sample with no reads has zero rows and columns, every other entry equals the products summed by hand, and its errors say no information. |
| `the_errors_mean_what_they_say` (information.rs, ignored), rewritten | §2.1; every reported error equals the test's own whole-matrix inverse to 10⁻⁶. |

`cargo test --release --lib parameter_estimation::joint::fit`: 79 passed, 1 ignored.
`cli::cross_platform_digests`: both pass with the checksums unchanged.

## 6. Validation

- **No fitted number moves**: both `cli::cross_platform_digests` tests pass unchanged; the correctness
  reviewer's three drawn fits that take the whole-matrix path give identical digests of everything but
  the errors on this step and its parent.
- **Same bits at 1, 4 and 8 threads**, the whole matrix included: `the_information_is_the_same_bits_at_any_pool_width`,
  and the correctness reviewer's probe fitting 20 samples of one to four libraries over 128 chunks,
  with and without the duplicated class.
- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- `cargo test --all-targets --all-features --no-fail-fast` (`tmp/fit_precision/suite_a8.log`): 4,937
  passed, 3 failed, 5 ignored; the three are the pre-existing failures
  (`examples/ng_generic_loci_dump.rs` × 2, `examples/ng_ssr_loci_dump.rs` × 1).
- The coverage run at its full count passed every recorded check (§2.1).
