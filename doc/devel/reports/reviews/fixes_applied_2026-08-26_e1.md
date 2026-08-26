# Fixes applied — ng calling loop E1 review

**Review:** [`ng_calling_loop_e1_2026-08-26.md`](ng_calling_loop_e1_2026-08-26.md).
**Date:** 2026-08-26. **Everything raised was applied.**

| finding | fix |
|---|---|
| Blocker — the per-sample leftover's reset untested | `a_sample_ruled_uncallable_at_one_locus_is_callable_at_the_next`, plus a destructured reset a sixth field cannot slip past |
| M1 — the sort's justification backwards | the module doc, the inline comment and the test's own doc now say the sort is insurance against a contract change, not a description of today's selection |
| M2 — the remapping's reset untested | `the_allele_mapping_is_this_locus_own_and_not_the_last_ones` |
| M3 — the cross-locus guard compared a count | the locus's region, and then the join itself at every sample; two tests |
| M4 — the sort unpinned past the first sample | `a_shifted_covering_sample_with_two_alternatives_is_sorted_like_the_first` |
| M5 — a comment naming an unreachable hazard | corrected, and `rows_left_for` added so the test can see what it claims |
| Minors | the module and scratch renamed; three parallel `Vec`s folded into one; the tolerance made relative and named; the §7 citation corrected; a dead `truncate` deleted; `#[must_use]`; a zero-sample refusal; the loop-shape test |

**Six tests added**, 4,710 → 4,716, each closing something measured to survive the first suite.

## Not changed

`GenericObservation::fill_from_supported_alleles`' doc tells a caller to add its return value to
selection's pool, which doubles the leftover. In `likelihood/`, whose plan has merged; raised for
its owner in the review.

## Validation

- `cargo fmt --all -- --check` — exit 0; `cargo clippy --all-targets --all-features -- -D warnings`
  — exit 0.
- `cargo test --lib` — `4716 passed; 0 failed; 14 ignored`.
- `cargo test --release --lib ng::calling --all-features` — `670 passed; 0 failed; 3 ignored`.
- The eight release-held checks downgraded to `debug_assert` in one run: **9 failures**, every
  check reached; the file restored from a byte-identical backup and the suites re-run green.
