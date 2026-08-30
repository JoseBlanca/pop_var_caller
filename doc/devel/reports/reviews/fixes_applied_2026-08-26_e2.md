# Fixes applied — ng calling loop E2 review

**Review:** [`ng_calling_loop_e2_2026-08-26.md`](ng_calling_loop_e2_2026-08-26.md).
**Date:** 2026-08-26. **Everything raised was applied**, except one item deferred with a
recommendation.

| finding | fix |
|---|---|
| Blocker — the gap's mechanism wrong in six places | all six corrected: nothing slides, the highest read group is dropped, and the failure lands at a locus rather than at assembly — which is why the check stays |
| M1 — the contamination map's keys unchecked, and its test unable to see it | an axis check, a fourth read group with no estimate in the fixture, and `a_contamination_estimate_off_the_read_group_axis_is_refused` |
| M2 — an unmeasured view's fabricated `source` | documented at the constant, pinned by a test; making it unrepresentable is `ContaminationView`'s owner's |
| M3 — one direction of the pairing rule untested | `a_minted_total_without_its_fitted_rate_is_refused` |
| M4 — a broken intra-doc link the crate denies | the full path, which needs no import |
| M5 — the ploidy (and the read group) in the lookup key unpinned | a haploid fixture and a two-library fixture |
| Minors | the unreachable `.zip()` arm's comment; the palindrome fixture; `checked_read_group_count_of`; two stale counts; a stray import |

**Seven tests added**, 4,726 → 4,733, each closing something measured to survive the first suite.

## Deferred, with a recommendation

The substitution lookup returns `Option` where `StratumFits::at` returns a typed error that can say
*this read group is not in the fit*. Both fill the same scoring context. **Mirror the error type
when the tract's context assembly is built** — the step that first holds both — rather than
inventing the shape here against no caller.

## Validation

- `cargo fmt --all -- --check` — exit 0; `cargo clippy --all-targets --all-features -- -D warnings`
  — exit 0; `cargo doc --no-deps --lib` — nothing from this module.
- `cargo test --lib` — `4733 passed; 0 failed; 14 ignored`.
- `cargo test --release --lib ng::calling --all-features` — `687 passed; 0 failed; 3 ignored`.
- The five release-held checks downgraded to `debug_assert` in one run: **6 failures**, every check
  reached; the file restored from a byte-identical backup and the suites re-run green.
