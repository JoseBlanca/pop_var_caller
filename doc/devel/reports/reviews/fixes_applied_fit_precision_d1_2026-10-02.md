# Fix Application Report: fit_precision_d1_2026-10-02.md

**Date:** 2026-10-02
**Source review:** `doc/devel/reports/reviews/fit_precision_d1_2026-10-02.md`
**Source state reviewed against:** `c15c624f` (a review object of the working tree, parent `9854283e`)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

- Review: 0 Blockers, 0 Majors, 3 Minors, Nits grouped.
- Applied: Mi1, Mi2, the four Nits. Deferred: Mi3, to step D2, which it is about.
- Validation: `cargo fmt --check` → 0; `cargo clippy --all-targets --all-features -- -D warnings` → 0; module
  (`cargo test --release --lib parameter_estimation::joint::sample_order`) → 4 passed; full suite → 5,015 passed, 3 failed (the pre-existing `examples/ng_generic_loci_dump.rs` × 2 and `examples/ng_ssr_loci_dump.rs` × 1), 8 ignored.
- Performance check: not applicable (no caller yet).

## 2. Findings table

| ID | Severity | Title | Final status | Files changed |
|---|---|---|---|---|
| Mi1 | Minor | the name reaches xxh3 through `str`'s unstable `Hash` layout | Applied | `sample_order.rs` |
| Mi2 | Minor | the doc's "over the name's bytes" | Applied | `sample_order.rs` |
| Mi3 | Minor | what D2 must pass in | Deferred to D2 | — |
| Nits | Nit | | Applied | `sample_order.rs` |

## 3. Questions asked and answers

None.

## 4. Per-finding log

### Mi1, Mi2
- `hash_sample` is `xxh3_64_with_seed(name.as_bytes(), seed)`; the module doc says the hash is over the name's
  UTF-8 bytes and nothing else, and why it does not go through `Hash` as the census does. The pinned order changed,
  as the review predicted, and is re-recorded from the dev container.

### Mi3
- **Deferred to D2**, which builds the caller: the order computed once in `census_fit.rs` from
  `cohort.sample_names()` and handed to the fit with each sample's rank.

### Nits
- The pinned test's doc says the order must be the same on every platform; "four tomato accessions" replaces the
  internal label; the function's doc says a run never holds two samples of one name; the duplicate-name test is
  removed with that promise.

## 5. Deferred findings to carry forward
- Mi3 — to D2.

## 6. Disputed findings to return to reviewer
None.

## 7–8. Failed validation; blocked by context mismatch
None.

## 9. Performance check
Skipped — nothing calls the function yet.
