# Fix Application Report: fit_precision_d2_2026-10-02.md

**Date:** 2026-10-02
**Source review:** `doc/devel/reports/reviews/fit_precision_d2_2026-10-02.md`
**Source state reviewed against:** `3e5b9dff` (a review object of the working tree, parent `cbb27a9c`)
**Execution mode:** interactive — six design questions put to the owner, all accepted (2026-10-02)
**Overall status:** Completed

---

## 1. Executive summary

- Review: 1 Blocker, 6 Majors, 11 Minors, Nits grouped.
- Applied: B1, M1–M4, M6, Mi1–Mi3, Mi5–Mi9, Mi11. Applied in part: M5, Mi4, the Nits. Deferred: Mi10 (accepted
  as recorded).
- Validation: fmt and clippy clean; module 71 passed; full suite 5,025 passed and the 3 pre-existing failures; the oracle's and the cross-platform checksums unchanged.
- Performance check: not applicable — no file under `benches/` reaches the subset path; its cost is measured in
  the implementation report §3.

## 2. Findings table

| ID | Severity | Title | Final status | Files changed |
|---|---|---|---|---|
| B1 | Blocker | the floor on the subset path untested | Applied | `ssr_fit.rs` test |
| M1 | Major | one walk from the last answer settles low | Applied (owner decision 1) | `ssr_fit.rs`, spec §4.4 |
| M2 | Major | subsets do not nest | Applied | `ssr_fit.rs` |
| M3 | Major | answered from fewer tracts than the floor | Applied (owner decision 5) | `ssr_fit.rs`, spec §4.4 |
| M4 | Major | the ≤ 256 test cannot fail | Applied | `ssr_fit.rs` test |
| M5 | Major | warm start, doubling, stop untested | Applied in part | `ssr_fit.rs` tests |
| M6 | Major | evidence counts and the file | Applied (owner decision 6) | `ssr_fit.rs` docs, plan E2 |
| Mi1 | Minor | detached doc comment, negated comparison | Applied | `ssr_fit.rs` |
| Mi2 | Minor | an unreachable path described and "tested" | Applied | `ssr_fit.rs` |
| Mi3 | Minor | a full copy at every sample | Applied | `ssr_fit.rs` |
| Mi4 | Minor | order unchecked; `subsets` ignored by `fit_strata` | Applied in part | `ssr_fit.rs` |
| Mi5 | Minor | a reader at low depth untested | Applied | `ssr_fit.rs` test |
| Mi6 | Minor | the summary untested | Applied | `ssr_fit.rs` test |
| Mi7 | Minor | top-up only when none | Applied (owner decision 2) | `ssr_fit.rs`, spec §4.4 |
| Mi8 | Minor | target on groups that cannot grow | Applied (owner decision 3) | `ssr_fit.rs`, spec §4.4 |
| Mi9 | Minor | the 2,048 step | Applied (owner decision 4) | `ssr_fit.rs`, spec §4.4 |
| Mi10 | Minor | one thread a stratum under several at once | Deferred — recorded in spec §4.4 | — |
| Mi11 | Minor | the curves' weight change unrecorded | Applied | spec §4.4 |
| Nits | Nit | | Applied in part | `ssr_fit.rs` |

## 3. Questions asked and answers

1. The six design choices of the review (the answer's subset from every starting point and the last answer;
   top-up below 8 readers; the target only on groups still growing; every sample past three quarters; a thin
   subset grows; the parameters file to record the subset size). — **Answer (owner, 2026-10-02):** yes to all.

## 4. Per-finding log

### B1, M3 — the floor
- `the_floor_is_judged_on_the_whole_stratum_and_a_thin_subset_grows`: forty tracts each read by one sample and a
  first subset of 4 are fitted on 8 samples (8 tracts), not refused; six such tracts are refused below the floor.
  `fit_on_growing_subsets` fits no subset with fewer tracts with reads than the floor, unless it holds every
  sample.

### M1 — the answer's subset
- When a subset other than the first stops the growth, it is fitted from every starting point and the walk from
  the last answer is added, the best winning (`the_best_of_the_walks`, the rule `the_better_walk` keeps).
  `a_subset_grows_…` asserts four walks and a score at least the three-start fit's on the same samples.

### M2 — nesting
- `grow_subset` adds to the previous subset's flags instead of drawing them afresh;
  `a_subset_tops_up_a_thin_group_and_keeps_its_samples_as_it_grows` checks every earlier sample stays in.

### M4 — the boundary
- `a_cohort_no_larger_than_the_first_subset_is_fitted_whole` runs a first subset of 20 at 20 samples (no subset)
  and at 21 (a subset).

### M5 — the schedule
- **Applied:** `next_subset_size` is its own function, tested at 2,169 and 64 samples; the final refit is asserted
  by its walk count. **Not applied:** a test that the decision walks start from the last answer rather than a
  starting point — they now only decide whether to grow, and the answer's subset is fitted from every start
  either way.

### M6 — what the counts mean
- `tracts_of_its_own`, `reads_crossing` and `samples_fitted_on` say the counts are among the samples the stratum
  was fitted on, and when `samples_fitted_on` is `None`. Plan step E2 now writes it in the parameters file and
  says so in the file's text for `expected_slipped_reads` (owner).

### Mi1, Mi2, Mi3
- `refused_before_any_walk`'s doc is back on it and names the subset path; the comparison uses `partial_cmp`. The
  "passed over" sentence is gone: a subset is now skipped only when thinner than the floor, which is tested. A
  stratum that grows to every sample is fitted on its own evidence, not a copy.

### Mi4
- **Applied:** an assertion that the order ranks as many samples as the homozygote excess has, and `fit_strata`'s
  doc says it ignores `config.subsets`. **Deferred:** a `SampleOrder` type guaranteeing a permutation — the one
  caller builds it with `sample_order`.

### Mi5, Mi6
- `one_read_makes_a_sample_a_reader_of_its_group`; `the_subsets_summary_counts_the_sizes`.

### Mi7, Mi8, Mi9
- `grow_subset` tops up every group below 8 readers; `level_is_measured_to` asks only groups with readers outside
  the subset (`the_level_is_judged_only_in_groups_that_can_still_grow`); `next_subset_size` takes every sample
  past three quarters.

### Mi10
- **Deferred**, accepted and recorded in the spec amendment: the default schedule is one stratum at a time.

### Mi11
- The spec amendment records the curve weight a subset-fitted stratum carries.

### Nits
- **Applied:** `level_is_measured_to`; `fit_from` and `parameters_of_fit` removed; plan labels out of the doc
  comments. **Deferred:** `SampleSubsets::first`'s name, `first.max(1)`, logging the subset rule, the median,
  the parallel arm's progress line, the unread guard counts, `walks` of the answer's subset only.

## 5. Deferred findings to carry forward
- Mi10 — the parallel arm runs a stratum's subsets on one thread (recorded).
- Mi4 in part — a `SampleOrder` type.
- The Nits listed above.

## 6. Disputed findings
None.

## 7. Mutations of the final code

Seven mutations of the final code, each run against the module's tests and restored (`tmp/fit_precision/d2_mut/`):

| mutation | result |
|---|---|
| the floor not judged on the whole stratum | killed (`the_floor_is_judged_…`) |
| a subset thinner than the floor fitted anyway | survived at first — the test's top-up already reached 8 tracts; with the test topping up to one reader, killed |
| `<` for `<=` at the first subset's size | killed (`a_cohort_no_larger_than_the_first_subset_…`) |
| a group topped up only when it holds none | killed (two tests) |
| each subset drawn afresh, not added to | survived — and changes nothing found: with every group topped up to 8, a group's kept readers are always the earliest of its readers in the order, so the subsets nest either way; the cumulative flags are a guard |
| the answer's subset not refitted from every start | killed (`a_subset_grows_…`) |
| the target asked of every live group | killed (`the_level_is_judged_only_in_groups_that_can_still_grow`) |

## 8. Commands run
- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::ssr_fit`
- `tmp/fit_precision/d2_mut/run.sh`
- the oracle (`scripts/promote_ng_oracle.sh`), clippy and the full suite (§1)
