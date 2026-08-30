# Code review — ng calling loop D2: the two cost invariants

**Scope:** the working-tree diff of step D2 of
[`calling_loop.md`](../../ng/impl_plan/calling_loop.md), on top of `e26c62f2` — three tests and
one test-only accessor.
**Date:** 2026-08-26. **Verdict: request changes** — **2 Majors, and 4 of 20 claims wrong**. All
applied; see [the fix report](fixes_applied_2026-08-26_d2.md).

**One agent, in its own worktree**, carrying the reliability checklist and step 8a together,
because the diff is three tests. **8 mutations run, 2 survived, 0 changed no behaviour.**

---

## M1 — the surrogate was honestly bounded and wrongly justified, and the real measurement was
twenty lines away

The diff measured the loop's zero-allocation invariant with a **surrogate**: the data pointer and
length of the scratch's per-locus buffers, compared across two runs at different pass counts. Its
doc said an allocation *count* was unavailable because a counting `#[global_allocator]` needs
`unsafe impl GlobalAlloc`, which `src/lib.rs`'s `#![forbid(unsafe_code)]` refuses.

**That reasoning is wrong, and the reviewer demonstrated it rather than argued it.**
`#[global_allocator]` is a safe *attribute*; the `unsafe impl` behind `dhat::Alloc` belongs to
dhat, which is already this repository's heap-profiling dependency. Declared inside the forbid
scope, it compiles and runs with no `unsafe` written anywhere — and against a per-pass temporary
`Vec` in the frequency loop it fails where the fingerprint test does not.

**The real obstacle is narrower than the one claimed**: a global allocator counts the whole
process, and the library suite runs its tests in parallel. That argues for a test binary of its
own, not for abandoning the measurement.

**Fixed.** `tests/ng_calling_loop_allocation.rs` counts `total_blocks` around two runs of one
locus at two passes and at four: **8 blocks each**, and one `Vec::with_capacity` in the seeded pass
takes them apart at 8 against 10 — one block per extra pass. The fingerprint tests stay as the
cheap guard that runs without the feature.

## M2 — the counter's own reset was untested

`the_emission_count_…` took a **fresh** scratch inside its `for config` loop. Deleting
`self.emission_cost = EmissionCost::default();` from `prepare_for_locus` left the whole suite
green — and it is not a no-op: two loci on one scratch then report `table_builds: 2, row_builds: 6,
emission_evaluations: 36` where the clean tree reports `1 / 3 / 18`. **One scratch across many loci
is the only shape a real worker has.** Fixed by hoisting the scratch out of the loop, which
strengthens the test for free.

## Four wrong claims of twenty

- **"A three-way product would give `3 × 3 × 3 = 27`" — wrong on this fixture.** Measured, the
  first-row-count bug gives **9** and a pooled-total bug gives **54**; 27 is the *correct* `Σ_s` of
  the *equal-count* fixture beside it. The idea was right and the number belonged to the other
  fixture. Corrected to name both wrong shapes with their measured values.
- **"the reads shared out so no sample is certain" — wrong, and inverted.** The most evenly shared
  sample is called `0/1` at **GQ 54.7**; the single-row sample is the least certain at **12.3**.
  Spreading the reads is what makes a heterozygote's two alleles both visible. Corrected.
- **"Where every buffer's bytes are" — wrong.** Seventeen of twenty: the two row scratches'
  buffers are not fingerprinted. **And they should not be** — `GenericRowScratch` sizes itself per
  *sample*, inside the table build, so it legitimately grows within a locus; the table build
  happens once, outside the frequency loop, so a pass still allocates nothing. The doc now says
  which buffers it covers and why the others are excluded.
- The fingerprint doc's `unsafe` justification, above.

**Everything else checked out**: `3 × (1 + 2 + 3) = 18`, `builds = 1`, `row_builds = 3`, the pass
counts of 2 and 4, the unequal `[1, 2, 3]` counts, the samples not being identical, and a wider
locus moving eleven of the seventeen fingerprint slots (`genotype_likelihoods` 3 → 10,
`error_spreads` 6 → 40).

## Mutations

**Caught:** a per-pass table rebuild (`1/3/18` → `2/6/36`); the first row's count charged for every
row (18 → 9); a constant fingerprint; per-pass buffer growth; an ignored pass cap. **Survived:**
the per-pass temporary — the disclosed blind spot, now caught by the counted test — and the deleted
counter reset, filed as M2. Dropping partials from the charge survives D2's three tests, whose
fixture has no partials, but D1's own test catches it; that is a note on D2's reach rather than a
gap in the suite.

## One thing found while fixing, worth recording

**The counted test caught a flaw in its own first fixture.** It cloned one call's allele table
*inside* the measured region and moved the other's in — 10 blocks against 6, with the failure
message blaming the loop. Cloning both before either reading is what makes the two runs comparable.

## Out of scope, and pre-existing

`cargo test --all-targets --all-features` fails in `benches/psp_writer_perf.rs`, which runs its
criterion benches once in test mode and indexes one past the end of its own record fixture
(`index out of bounds: the len is 3300000 but the index is 3300000`). Nothing on this branch
touches the PSP writer or its benches. **Raised for its owner**; this step validated with the
branch's established gate.

## Verification

- `cargo fmt --all -- --check` — exit 0; `cargo clippy --all-targets --all-features -- -D warnings`
  — exit 0.
- `cargo test --lib` — `4694 passed; 0 failed; 14 ignored`.
- `cargo test --release --lib ng::calling --all-features` — `648 passed; 0 failed; 3 ignored`.
- `cargo test --test ng_calling_loop_allocation --features dhat-heap` — `1 passed; 0 failed`.
