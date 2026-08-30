# Fixes applied — ng calling loop D2 review

**Review:** [`ng_calling_loop_d2_2026-08-26.md`](ng_calling_loop_d2_2026-08-26.md).
**Date:** 2026-08-26. **Everything raised was applied.**

| finding | fix |
|---|---|
| M1 — the allocation count was declared out of reach on wrong grounds | `tests/ng_calling_loop_allocation.rs`: dhat's counting allocator in a test binary of its own, `total_blocks` around two runs at two and four passes — 8 blocks each, and 8 against 10 with one `Vec::with_capacity` added to the seeded pass. The fingerprint tests stay as the cheap guard that runs without the feature |
| M2 — the emission counter's own reset untested | one scratch across both loci, as a worker has: with `prepare_for_locus`'s reset deleted the second locus reports `2 / 6 / 36` |
| `3 × 3 × 3 = 27` belonged to the other fixture | replaced by the two shapes this fixture separates, measured: 9 and 54 |
| "no sample is certain" — inverted | the most evenly shared sample is called at GQ 54.7 and the single-row one at 12.3 |
| "every buffer's bytes" — seventeen of twenty | the doc now names what it covers, and why the two row scratches are excluded: they size themselves per sample inside the table build, which runs once, outside the loop |

**One thing the new test found about itself:** its first fixture cloned one call's allele table
inside the measured region and moved the other's in, reporting 10 blocks against 6 and blaming the
loop. Both are cloned before either reading now.

## Not changed

`benches/psp_writer_perf.rs` panics under `cargo test --all-targets`, indexing one past the end of
its own record fixture. Pre-existing, untouched by this branch, and raised for its owner in the
review.

## Validation

- `cargo fmt --all -- --check` — exit 0; `cargo clippy --all-targets --all-features -- -D warnings`
  — exit 0.
- `cargo test --lib` — `4694 passed; 0 failed; 14 ignored`.
- `cargo test --release --lib ng::calling --all-features` — `648 passed; 0 failed; 3 ignored`.
- `cargo test --test ng_calling_loop_allocation --features dhat-heap` — `1 passed; 0 failed`.
