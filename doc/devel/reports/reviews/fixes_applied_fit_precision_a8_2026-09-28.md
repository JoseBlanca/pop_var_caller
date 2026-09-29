# Fix Application Report: fit_precision_a8_2026-09-28.md

**Date:** 2026-09-29
**Source review:** `doc/devel/reports/reviews/fit_precision_a8_2026-09-28.md`
**Source state reviewed against:** `c6ec817d` and `e3641146` (review-only commits; parent `11da7a62`)
**Execution mode:** non-interactive
**Overall status:** Completed (one design limit carried to the owner at checkpoint A′)

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 3
- Minors: 10
- Nits: 1 (five parts)

### Outcome totals
- Applied: 10
- Applied with adaptation: 2
- Deferred: 2
- Awaiting user answer: 0 (M1 goes to the owner at the checkpoint with a recommendation)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --all-targets --all-features --no-fail-fast` → 101: 4,937 passed, 3 failed, 5 ignored; the three are the pre-existing failures (`examples/ng_generic_loci_dump.rs` × 2, `examples/ng_ssr_loci_dump.rs` × 1). `tmp/fit_precision/suite_a8.log`.
- `cargo test --release --lib parameter_estimation::joint::fit` → 0, 79 passed, 1 ignored (last run before commit)
- `cargo test --release --lib cli::cross_platform_digests` → 0, both pass, checksums unchanged
- The coverage test at 200 cohorts a regime → 0, every recorded check holds (`tmp/fit_precision/a8_coverage_full.log`)
- `cargo doc`, `cargo audit` → not run

### Unresolved high-priority findings
- M1 — the limit on samples; the owner decides whether to cap parameters too.

## 2. Findings table

| ID | Severity | Title | Final status | Files changed |
|---|---|---|---|---|
| M1 | Major | the limit is on samples, not size | Deferred (owner) | report §3 item 1, the constant's doc |
| M2 | Major | the inversion order is untested | Applied | `standard_errors.rs` |
| M3 | Major | a sample slot's own interval is untested | Applied | `standard_errors.rs` |
| m1 | Minor | module doc contradicts the code | Applied | `standard_errors.rs` |
| m2 | Minor | spec §3.2 | Deferred (owner) | None |
| m3 | Minor | `absorb` checks only in debug builds | Applied | `information.rs` |
| m4 | Minor | bare tuple | Applied | `standard_errors.rs` |
| m5 | Minor | the 0.08 comparison over different sets; bounds widened blind | Applied | `information.rs` |
| m6 | Minor | misplaced doc comment | Applied | `information.rs` |
| m7 | Minor | the coverage assertion cannot see a wrong reason | Applied with adaptation | `standard_errors.rs` |
| m8 | Minor | no no-reads test through the whole matrix | Applied | `information.rs` |
| m9 | Minor | two docs call the blocks "what the fit reports" | Applied | `information.rs` |
| m10 | Minor | the report's §2, §5 stale | Applied | report |
| n1 | Nit | names, scratch, a tally's doc | Applied with adaptation | `information.rs` |

## 3. Per-finding log

- **M1 — Deferred to the owner; decided at checkpoint A′ (2026-09-29): cap at 188 parameters, built in its own commit.** The plan fixes the limit on samples. The report (§2.2, §3 item 1)
  and the constant's doc give the sizes: 20 samples of 16 libraries have 668 parameters and 223,446
  products a position; the correctness reviewer measured 231 MB of whole matrices in the final pass.
  **Recommendation for checkpoint A′:** cap the parameters at 188 as well (20 samples of four
  libraries, 18.3 MB, 1.1 times the blocks' pass), falling back to the blocks above it.
- **M2 — Applied.** `on_an_arrow_the_whole_matrix_and_the_blocks_drop_the_same_parameters`: an arrow
  where the density's first shape copies sample 0's clean rate and sample 1's second library's
  mismapped rate copies its clean rate; both ways drop the shape and the mismapped rate, keep the
  sample's clean rate, and agree on every other error to 10⁻⁹. Adapted from the design reviewer's
  probe, with `blocks_of_dense` factored out of `an_arrow_of` so the matrix can be edited before the
  blocks are cut from it. The design reviewer showed the cohort-first order survives the committed
  tests at `c6ec817d`; this test's fixture is the case that order decides.
- **M3 — Applied.** `the_whole_matrix_judges_each_slot_against_its_own_interval`, from the
  reliability reviewer's probe, which failed under the mutation that judges every slot by the clean
  rate's interval.
- **m1 — Applied.** The module doc describes both ways, the order ("each sample's own first, the
  samples in turn, then the cohort's"), that the drops match only on an arrow, and what the whole
  matrix does at two samples.
- **m2 — Deferred:** the spec is the owner's; offered at checkpoint A′ with the other §3.2 amendments.
- **m3 — Applied.** `InformationSums::absorb` refuses a mismatch with `unreachable!`, as
  `Statistics::absorb` does; `FullInformation::absorb` asserts the size in every build.
- **m4 — Applied.** `InvertedParameter { owner: Owner, slot, row }`, with `Owner::Cohort` or
  `Owner::Sample(s)`.
- **m5 — Applied.** The blocks are tallied only for estimates that also have a reported error, so the
  0.08 comparison is over the same estimates (at 4 and 20 samples every estimate has both, so no
  figure moved). The within-one band is back to 0.58–0.78, which the 200-cohort run's 0.650 to 0.741
  fits.
- **m6, m9, m10 — Applied.**
- **m7 — Applied with adaptation:** not in the coverage test, whose cohorts drop nothing at 4 and 20
  samples; the reasons are held by the order test (M2) and `the_whole_matrix_says_why_a_parameter_has_no_error`.
- **m8 — Applied.** `a_sample_without_reads_has_no_row_in_the_whole_matrix`, from the reliability
  reviewer's probe (it failed under 10 of its mutations).
- **n1 — Applied with adaptation:** `sample_parameters` renamed `rows_of_sample`; the `blocks` tally's
  doc rewritten. Kept: `FULL_MATRIX_SAMPLES` (the plan's name), "full" in type names beside "whole
  matrix" in prose, and the scratch row with its hand-written `Debug` (the reliability reviewer found
  the determinism test still sees every summed number).

## 4. Deferred findings to carry forward
- M1 — the cap on parameters, the owner's.
- m2 — spec §3.2, the owner's.

## 5. Performance check
No criterion comparison: the change touches only the final pass. Measured directly (implementation
report §2.2): the whole matrix adds 1 to 4% to that pass at 20 one-library samples; 1.10 to 1.12 times
at four libraries a sample; 1.33 to 1.64 at sixteen.

## 6. Notes
- All three reviewers' worktrees were removed after their findings and probes were copied to
  `tmp/review_2026-09-28_fit_precision_a8/evidence_*`.
