# ng calling loop — D2: what a locus pays for, and what a pass does not

**Step:** D2 of [`calling_loop.md`](../../ng/impl_plan/calling_loop.md) — the two cost invariants.
**Design authority:** [`spec/calling_em_loop.md`](../../ng/spec/calling_em_loop.md) §13 tests 5
and 7; [`arch/calling_em_loop.md`](../../ng/arch/calling_em_loop.md), *Test & bench shape*.
**Date:** 2026-08-26. **Branch:** `ng-calling-loop`.

---

## 1. Why this is its own commit

**Both failures it guards against are invisible in the output.** A genotype-likelihood table
rebuilt on every pass of the frequency loop returns exactly the same genotypes as one built once;
a `Vec` allocated inside a pass changes no number at all. Neither shows up as a wrong call, a
panic, or a failing assertion — only as a run that costs more than the design says it does. So the
instrument *is* the test, and the plan asks for it in a commit of its own so that the counter and
the code it counts cannot be adjusted together.

## 2. The first invariant: `candidates × Σ_s (observations in sample s) × builds`

`the_emission_count_is_candidates_times_the_sum_over_samples_at_every_pass_count` asserts the
formula rather than a literal, over a fixture whose per-sample observation counts are
**deliberately unequal**: three samples showing one, two and three observations over three
candidate alleles.

- The right answer is `3 × (1 + 2 + 3)` = **18**.
- **The two wrong shapes it separates, both measured on this fixture:** charging the first row's
  count for every row gives `3 × 1 × 3` = **9**, and charging the locus's pooled total for every
  row gives `3 × 6 × 3` = **54**. A fixture whose three samples showed three observations each
  would report **27** under all three — which is why the spec says `Σ_s`, and why a fixture built
  so the samples match is the one shape that hides the bug.
- **The same 18 at two passes and at four.** The test runs the same locus at a cap of two and at
  the shipped default, asserts the pass counts really are 2 and 4, and asserts the identical
  `EmissionCost` from both. Without that second run the invariant is a fact about one pass count.
- **One scratch across both, which is the only shape a real worker has.** The review measured what
  a fresh scratch per locus hides: with `prepare_for_locus`'s counter reset deleted, the second
  locus reports `table_builds: 2, row_builds: 6, emission_evaluations: 36` and every test still
  passes.

`EmissionCost` itself landed with D1; what D2 adds is the fixture shape that can fail.

## 3. The second invariant: nothing allocates inside a pass

**It is measured twice, and the first draft of this report was wrong about why.** That draft said
the allocation count was out of reach because a counting `#[global_allocator]` needs
`unsafe impl GlobalAlloc` and `src/lib.rs` forbids `unsafe`. **The review showed otherwise in
twenty lines**: `#[global_allocator]` is a safe *attribute*, the `unsafe impl` behind `dhat::Alloc`
is dhat's own, and dhat is already this repository's heap-profiling dependency. The forbid stands
untouched. The real obstacle is narrower — a global allocator counts the whole process and the
library suite runs its tests in parallel — and its answer is a test binary of its own.

**So `tests/ng_calling_loop_allocation.rs` counts it.** It installs dhat's allocator, warms the
scratch and the genotype-table cache on the locus's shape, then reads `total_blocks` around two
runs of the same locus at two passes and at four. **Measured: 8 blocks each**, and one
`Vec::with_capacity` added to the seeded pass takes them apart — 8 against 10, one block per extra
pass. It runs under `--features dhat-heap`; without the feature the file compiles to nothing.

**The test caught a flaw in its own first fixture, which is worth recording.** Its first draft
cloned one call's allele table *inside* the measured region and moved the other's in — 10 blocks
against 6, with the message blaming the loop. Cloning both before either reading is what makes the
two runs comparable.

**And the cheap half still runs on every build.**
`CallingScratch::buffer_fingerprints()` (test-only) returns the data pointer and length of the
seventeen per-locus buffers. A `Vec` that **grew** during the loop moves its bytes, so its pointer changes;
one refilled in place does not. Two tests use it:

- `no_buffer_of_the_scratch_moves_or_grows_however_many_passes_the_loop_takes` calls the same
  locus twice on **one** scratch, at two passes and at four, and asserts every fingerprint is
  identical.
- `a_wider_locus_than_the_worker_has_met_does_grow_the_scratch` asserts the opposite at a locus of
  four alleles after one of two — **without it the first test passes against an implementation
  whose buffers never change because nothing ever prepares them.**

**What the surrogate cannot see is exactly what the counted test does:** a temporary allocated and
dropped inside a pass leaves no trace in any buffer. The two together are the invariant.

**The two row scratches are deliberately outside the fingerprint**, and their absence is the
invariant rather than a gap: `GenericRowScratch` sizes itself per *sample*, inside the table build,
so it legitimately grows within a locus when a wider sample arrives. The table build happens once,
outside the frequency loop, so a pass still allocates nothing — which is what test 7 claims.
Fingerprinting them would fail on correct code.

## 4. Tests

**Three in the library suite**, 4,691 → 4,694, plus **one in a test binary of its own** under
`--features dhat-heap`. **No new release-held assertions**: D2 adds observables and tests and
asserts nothing in shipped code, so the downgrade battery has nothing new to reach.

## 5. Validation

- `cargo fmt --all -- --check` — exit 0; `cargo clippy --all-targets --all-features -- -D warnings`
  — exit 0.
- `cargo test --lib` — `4694 passed; 0 failed; 14 ignored`. Before D2: **4,691**.
- `cargo test --release --lib ng::calling --all-features` — `648 passed; 0 failed; 3 ignored`.
  Before D2: **645**.
- `cargo test --test ng_calling_loop_allocation --features dhat-heap` — `1 passed; 0 failed`, and
  it fails as `8 against 10` with one `Vec::with_capacity` added to the seeded pass.

## 6. Follow-up

- **The counted test runs only under a feature, so CI's `--lib` steps do not reach it.** The commit
  gate (`cargo test --all-targets --all-features`) does. Whoever tunes the CI matrix should know
  that `--features dhat-heap` is what makes this file exist.
