# Code Review: estimation_memory_issue3
**Date:** 2026-09-26
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** issue 3 of `doc/devel/implementation_plans/estimation_memory.md` — contamination's markers chosen in two passes, without the samples × positions arrays
**Status:** Approve-with-changes

---

## 1. Scope

- **What was reviewed:** commit `dbdba774` against `f009bebe` (main), branch
  `contamination-two-pass`.
- **In-scope files:** [contamination.rs](../../../../src/parameter_estimation/joint/contamination.rs)
  (the new `markers` and its helpers, the dense oracle and the differential tests), the
  implementation report, the status block.
- **Categories dispatched:** two agents in their own worktrees — reliability + errors +
  refactor_safety; naming + idiomatic + smells + float_portability + module_structure + extras.
  `unsafe_concurrency`, `defaults` and `tooling` were not dispatched: the change adds no
  concurrency, no option and no dependency, in one file.

## 2. Verdict

Approve-with-changes. **Both agents compared the new `markers` with the dense one line by line
and found no semantic difference**: the filters run in the same order, each sample's counts are
pooled with the same saturating `u32` additions and widened to `u64` in sample order, and an empty
cohort, a cohort under the coverage floor and a mismapping list shorter than the positions all give
what they gave. The findings are about what the tests reach and what the report says.

## 3. Execution status

- clippy `-D warnings` clean; `cargo test --lib -- contamination cross_platform` 37 passed;
  full suite 4,873 passed, 3 failed (the pre-existing probe examples).
- Reliability agent: 9 mutations, 5 caught, 4 survived, and none of the 4 changed a result on the
  test panels — two because they are equivalent in production (swapping the floor and zero-total
  checks, whose second can never fire; `+` for `saturating_add` at depths far below overflow), two
  because the fixtures never reach the branch (M1).
- The naming agent confirmed the report's fixture claims by printing them: 2,816 and 2,855 markers
  from 3,000 positions; 631 positions covered by exactly 8 samples and 568 by 7; markers 3,767
  without the mismapping posteriors and 3,390 with.

## 4. Open questions and assumptions

None.

## 5. Top 3 priorities

1. **M1** — the oracle's panels never reach the common-allele filter, the depth cap or the dosage
   fallback.
2. **Mi1** — the report undercounts what still grows with the cohort.
3. **Mi3** — the comparison skips any field later added to `Marker`.

## 6. Findings

### Major

**M1: contamination.rs, the differential tests — three branches of the pooling are unreachable
from their fixtures.** **Categories:** reliability. Confidence: High. Every panel draws its
alternative reads on one allele (`C`) at 0.33 to 6 reads a position. So making the
common-allele filter always true, and dropping the depth cap's clamp, both left every test green,
and the `2.0 * pooled` fallback in `dosage_of` (every genotype's weight underflowing, which needs
about 1,074 reads) is never taken. The code is correct on all three — it is the old text — but the
test that guards bit-identity cannot fail there. Fix: a panel at 150 reads a position with some
alternative reads on a second allele, asserting it reaches both branches; a unit test of the
fallback.

### Minor

- **Mi1** the report counts only the dosages (0.9 GB) as what still grows; each marker also holds
  four per-library `u32` lists, 16 bytes × units × markers — 1.8 GB at 2,169 single-library samples
  and 52,525 markers. And 52,525 is the tomato count from before mismapped positions were refused;
  with the current default the count is 31,758. Plan §6 has the same gap.
- **Mi2** the report points at the commit message for the full-suite figures.
- **Mi3** `the_same_markers` reads fields one by one; a field added to `Marker` would be skipped
  silently. Destructure it.
- **Mi4** no differential test has fewer than 30 samples; one sample and five are not pinned.
- **Mi5** `pool_one_sample` takes two `&mut [u32]` in a row and `dosage_of` two `u32` then two
  `f64`; swapping either pair compiles.

### Nits

- `fill_each_librarys_reads` (possessive) and its parameter `out`; `positions` holds a count.
- The `#[allow(clippy::too_many_arguments)]` on `markers_dense` does nothing.
- The report's 64 MB differs from the plan's 40 MB (64 MB is right) without saying so.
- A plain `+=` in `major_alleles` (would need about 16.8 million observations at one position), an
  `.expect` without a comment, a `..Default::default()` in a test helper — all carried over.
- The status block's "Open:" line is stale.

## 7. Out of scope observations

- **Performance** (not measured): the first pass moves about 130 MB a sample, roughly 15–30 s
  single-threaded at 2,169 samples — no slower than the dense loop over 35 GB. The totals cannot be
  folded into `pool_one_sample` bit-identically, since `covered` needs a sample's depth after all
  its read groups are pooled. Two savings that keep the output: `mem::take` instead of `fill(0)`,
  and computing each sample's counts at the markers inside the per-library pass instead of pooling
  every position a second time.

## 8. Missing tests to add now

- A deep, two-allele differential test (M1); a unit test of the dosage fallback (M1); a panel under
  the coverage floor (Mi4).

## 9. What's good

- The pooling is one function both passes call, so the two cannot disagree about a sample's counts.
- The dense version is kept verbatim as the oracle, inside the test module where the library cannot
  reach it.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo test --lib -- contamination cross_platform`
