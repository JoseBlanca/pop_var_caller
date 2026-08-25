# Fixes applied — ng calling loop, C3a

**Date:** 2026-08-25
**Review:** [ng_calling_loop_c3a_2026-08-25.md](ng_calling_loop_c3a_2026-08-25.md)
**Branch:** `ng-calling-loop`

## The Blocker

**`next.fill(0.0)` beside `current.fill(0.0)`** in the count-axis fold — one line. The comment
above it now records what the alternative costs, because it is a silent wrong answer rather than
a crash: 46.3 Phred against production's 733.7 at 63 diploid samples.

Verified singly against the fixed tree: deleting it again fails three tests
(`the_exact_zero_term_is_what_keeps_a_fifty_sample_cohort_off_the_ceiling`,
`the_collapse_strides_by_the_allele_count_and_not_by_the_ploidy`,
`a_finite_quality_above_the_ceiling_is_capped`).

## Code fixes

| change | why |
|---|---|
| `genotype_quality`'s check moved from the winner to the **row's total**, computed in the same walk | a single `NaN` beside real values never wins the fold — `[0.7, NaN, 0.1]` returned an ordinary 5.2288 Phred. `NaN` survives addition, so the total sees it; and the same check catches a row that does not sum to one |
| an assertion in the collapse: **every sample must have one possible genotype** | an all-`−∞` row drives `−∞ − −∞` to a `NaN` that surfaced at the end as a complaint about the normalisation, two steps from the cause |
| the final `panic!`'s message names both causes | it named only the above-zero one |
| `cohort_summing_buffers_mut` given its doc comment back | the new accessor's doc had been inserted inside it |

## Tests added — six, each against a mutation that had survived or a check nothing reached

| test | the mutation or check it covers |
|---|---|
| `the_collapse_strides_by_the_allele_count_and_not_by_the_ploidy` | the stride, invisible at diploid biallelic loci where both numbers are 2 |
| `shifting_every_likelihood_by_a_constant_leaves_the_quality_alone` | `log_scale += largest`, inert on every zero-peaked fixture |
| `a_finite_quality_above_the_ceiling_is_capped` | the finite-above-ceiling arm |
| `a_single_nan_beside_real_probabilities_is_refused` | the `NaN` hole above |
| `a_sample_with_no_possible_genotype_is_refused` | the new collapse assertion |
| `a_winning_probability_above_one_is_refused_even_where_the_row_totals_one`, `a_copy_count_table_of_the_wrong_shape_is_refused`, `a_log_axis_of_the_wrong_length_is_refused` | three release-held checks the first battery showed were unreached |

Two tests were **rewritten rather than kept**, because they were passing on the defect:
`a_confident_cohort_survives_the_fold_rather_than_hitting_the_ceiling` became
`the_exact_zero_term_is_what_keeps_a_fifty_sample_cohort_off_the_ceiling`, and
`reference_looking_samples_do_not_inflate_the_site_quality` became
`a_locus_nobody_carries_stays_at_nothing_however_many_samples_look_at_it` — see §8 of the review
for why the property the first one claimed is not one this module can demonstrate.

## Documentation corrections

- The paragraph saying the exact zero-term override never matters, and that *"a later step
  trimming it should know that the tests will not object"*, is replaced by the measurement on the
  fixed fold: nothing at 20 samples, 4295.97 against the ceiling at 50, 6899.70 against the
  ceiling at 80.
- The clamp test's claim that removing the clamp causes a panic — it does not; `f32::clamp` maps
  `+∞` to the ceiling. The comment now says the clamp is currently untested.
- The `NaN` assertion's message, which claimed to catch a `NaN` anywhere in the row.
- The bystander fixture's two contradictory descriptions ("firmly say reference" and "thin on
  purpose").

## Renames

`kernel` is gone — `copy_count_log_likelihoods`, `copy_counts_per_sample`, `copy_counts`;
`count_axis_current/next/log` → `allele_count_distribution`, `allele_count_distribution_next`,
`log_allele_count_distribution`; `fold_count_axis` → `fold_samples_into_allele_counts`;
`genotype_quality` → `score_best_genotype`; `site_quality_baseline` →
`score_uncorrected_site_quality`.

## Not applied

**`ArtifactTestCounts` keeps its name.** The reviewer is right that `primary_alternative` is an
`AlleleId` and not a count, but the name is `spec/calling_quality.md` §10's. Renaming a type the
design document names is a design change, not a naming fix.

## Validation after the fixes

- `cargo fmt --all -- --check` — exit 0.
- `cargo clippy --all-targets --all-features -- -D warnings` — exit 0.
- `cargo test --lib` — `4643 passed; 0 failed; 14 ignored` (4,635 before the fixes; 4,620 before
  C3a).
- `cargo test --release --lib ng::calling --all-features` — `597 passed; 0 failed; 3 ignored`.
- The release-held-assertion battery, all nine downgraded together under `--release`:
  `586 passed; 11 failed`, every check reached.
