# Fixes applied — ng calling loop, C2

**Date:** 2026-08-25
**Review:** [ng_calling_loop_c2_2026-08-25.md](ng_calling_loop_c2_2026-08-25.md)
**Branch:** `ng-calling-loop`

**All applied but one, and that one is recorded as a decision rather than a skip.**

## Signature changes

| change | why |
|---|---|
| `cohort_copies_have_settled` → `cohort_expected_copies_have_settled`, `previous_copies`/`current_copies` → `previous_/current_expected_copies` | the crate's word for the quantity is *expected copies*, and the call site already spells `previous_cohort_expected_copies()` |
| the stopping rule takes `ploidy: Ploidy, sample_count: usize` in place of `cohort_chromosomes: f64` | two adjacent `f64` parameters were transposable and every assertion admitted the transposition; the new pair cannot be swapped, and an infinite or `NaN` chromosome count is no longer expressible — **one release-held check deleted rather than tested** |
| `run_frequency_loop` no longer takes a `Ploidy` | it already takes the `GenotypeTableView`, which was *built* for a `(ploidy, allele count)` shape. Two sources, nothing comparing them: measured, a `ploidy` of 64 against a diploid table returns `passes: 2, converged: true` where the truth is `passes: 4`, identically in debug and release |
| `#[must_use]` on `FrequencyLoopOutcome` and on the stopping rule | with no `Result` in this module the outcome is the only carrier of *did not settle*, and §6 requires the flag to reach the output; a bare-statement call compiled under `-D warnings` |
| the row-shape assertion split into two, each with its own message | the two guard tests were matching two substrings of one message, so neither could say which condition fired |

## Tests added — seven, each against a mutation that had survived

| test | the mutation it kills |
|---|---|
| `each_sample_is_scored_against_its_own_inbreeding_coefficient` | every sample scored with `inbreeding_by_sample[0]` |
| `a_fall_larger_than_every_rise_has_not_settled` | `.abs()` deleted |
| `run_frequency_loop_reports_converged_when_the_last_allowed_pass_settles` | the cap test moved ahead of the settled test |
| `more_inbreeding_coefficients_than_samples_are_refused` | `assert_eq!` weakened to `assert!(len >= sample_count)` |
| `an_infinite_threshold_is_refused` | `threshold.is_finite() &&` deleted |
| `a_movement_exactly_at_the_threshold_has_not_settled` | `<` relaxed to `<=` |
| `the_loop_settles_at_three_alleles_and_the_copies_still_sum_to_the_chromosomes` | no mutation — it closes the gap that **no C2 test ran the loop past two alleles at all** |

**Each was re-run singly against the fixed tree**: one failing test per mutation, the named one,
with no other test disturbed.

## Renames and doc corrections

- `never_settles_before` → `settles_only_at_a_bitwise_fixed_point`. The old name was false of
  every fixture: a threshold of `1e-300` is met the moment two passes agree bitwise, and
  measured, `three_samples_pulling_apart` reaches that at **pass 29** and reports
  `converged: true`.
- `loop_over`'s `seed` → `seed_concentration` — `seed` reads as an RNG seed where the crate's
  term is the prior's concentration.
- The stopping rule's doc claimed the `NaN` sentinel row is *"the row this function is first
  handed"*. **It is not** — the prior-free initialisation's M-step writes finite copies before
  the first swap. Reworded: the guarantee is for C3's final pass and D1's outer rounds.
- *"`u8` and `usize` both widen to `f64` exactly at every value either can hold"* — `usize` does
  not: `9007199254740993usize as f64 == 9007199254740992.0`.
- *"Both sides are given a threshold of `1e-300`, so neither can stop early"* — the hand-driven
  side has no stopping rule at all, which is a better reason for the test than the one written.
- The guard test's comment said a short coefficient slice *"would panic on an index"*; the pass
  walks the coefficients, so it silently scores fewer samples. Corrected, and the fixture now
  carries real likelihoods so the only reachable panic is the one it names.
- `is step D1's` → `step D1 of the calling-loop plan`, matching the `expect(dead_code)` strings.

## Not applied

**`FrequencyLoopOutcome::passes` stays a `u32`.** The review is right that the field's doc calls
zero impossible while the type admits it. But its consumer is `LocusInference::passes`, itself a
`u32` whose constructor asserts `passes > 0` for callers other than this loop. Retyping only
this end leaves two adjacent types in the same plan disagreeing and removes no check. **C3 is
where the two meet**, and the question is recorded there.

## Validation after the fixes

- `cargo fmt --all -- --check` — exit 0.
- `cargo clippy --all-targets --all-features -- -D warnings` — exit 0.
- `cargo test --lib` — `4620 passed; 0 failed; 14 ignored` (4,613 before the fixes; 4,603
  before C2).
- `cargo test --release --lib ng::calling --all-features` — `574 passed; 0 failed; 3 ignored`.
- The release-held-assertion battery, all five downgraded together under `--release`:
  `567 passed; 7 failed`, every check reached.
