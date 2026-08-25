# Fixes applied — ng calling loop D1 review

**Review:** [`ng_calling_loop_d1_2026-08-26.md`](ng_calling_loop_d1_2026-08-26.md).
**Date:** 2026-08-26. **Everything raised was applied.** Nothing was deferred; one finding
(the module split for a second seam arm) was accepted as a note for the plan that adds that arm.

## Code

| finding | fix |
|---|---|
| Blocker — the row-to-sample join untested, and the mutation silent | a fixture with the uncallable sample **first**, where a row's index and its sample's differ |
| M1 — the warrant stamped `FittedHere` over a defaulted calibration | derived: `weakest_warrant_at_the_locus` folds the calibrations of the read groups whose reads reached the locus with `Provenance::weaker_of`. `LocusInference::weakest_provenance`'s stale "no ordering" paragraph corrected |
| M2 — the map derived twice, compared by count | the final pass reads the stored map and asserts agreement at every sample; the count check kept for rows claimed past the end of the run |
| M3 — three `sample_count()`s, one changed meaning | `CallingScratch::sample_count` → `row_count`, field, parameter and messages with it |
| M4 — two doc comments stating the order the code refuses | corrected, with the measured consequence ("0 of them were claimed") |
| M5–M7 — three tests that could not fail | the pass cap, the per-sample inbreeding coefficient, and three emission-count fixtures |
| M8 — the refusals fired after the scratch was prepared | both moved to `call_locus`'s front door; the inner restatements are now `debug_assert` and `unreachable!` |
| Minors | `// PANIC-FREE:` markers, `checked_mul` for the error-spread table, `#[expect]` over `#[allow]`, `EmissionCost` to `pub(crate)`, four renames, `GenericLocusSample::is_callable`, a mangled panic message |

## Prose

**Ten claims corrected**, all mechanisms — the four inverted row-map stories, the outer loops'
"cannot be reached", the `never_loop` reason's mutation claim, the claim/prepare consequence, where
the length check fires, and two counts. The implementation report carried three of them and is
corrected too.

## Tests

**Eleven added**, 4,680 → 4,691. Each closes something the review measured to survive: the
gap-first join, the derived warrant and its identity case, the pass cap, the per-sample
coefficient, three emission counts, a four-pass build, a cohort of one, the repeat-tract callable
ruling as a unit, and rows claimed past the end of the run.

## Validation

- `cargo fmt --all -- --check` — exit 0; `cargo clippy --all-targets --all-features -- -D warnings`
  — exit 0.
- `cargo test --lib` — `4691 passed; 0 failed; 14 ignored`.
- `cargo test --release --lib ng::calling --all-features` — `645 passed; 0 failed; 3 ignored`.
- The six release-held checks downgraded to `debug_assert` in one run: `638 passed; 7 failed`,
  every check reached. Files restored from byte-identical backups and the suites re-run green
  before the commit.
