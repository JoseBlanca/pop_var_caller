# Code Review: ng calling loop — B2, the M-step

**Date:** 2026-08-25
**Branch:** `ng-calling-loop`, worktree `../pop_var_caller-calling-loop`
**Plan:** [calling_loop.md](../../ng/impl_plan/calling_loop.md), step B2
**Implementation report:** [ng_calling_loop_b2_2026-08-25.md](../implementations/ng_calling_loop_b2_2026-08-25.md)
**Fixes applied:** [fixes_applied_2026-08-25_v4.md](fixes_applied_2026-08-25_v4.md)

## 1. Scope

B2's working-tree diff on top of `06a9fac9` — `sum_cohort_expected_copies` and its five tests
in `summarise_condition.rs`, plus `CohortSumBuffers` and its accessor in `calling/mod.rs`.
+430 lines, 2 files. B1's half of the same file was out of scope, having been reviewed the same
day.

**Three agents, each in its own worktree**, covering `reliability`; `errors` + `defaults` +
`smells` + `refactor_safety` + `module_structure`; and `naming` + `idiomatic` + `extras` + the
skill's step 8a. Three rather than B1's five, in proportion to a diff a quarter the size —
recorded rather than silent.

## 2. Verdict

**Request changes** — 1 Blocker, 10 Majors, 12 Minors. The Blocker and most Majors are about
tests that could not fail; two Majors are defects in the function.

## 3. Execution status

| command | exit | result |
|---|---|---|
| `cargo fmt --all -- --check` | 0 | no output |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 | no warnings |
| `cargo test --lib` | 0 | `4589 passed` at review time (B1 left 4,584) |
| `cargo test --release --lib ng::calling --all-features` | 0 | `543 passed; 0 failed` — green since `06a9fac9` made it a CI step |

## 4. The Blocker — a test named for B2's property that B2 cannot fail

`presenting_the_samples_in_another_order_calls_the_same_genotypes` **passed against a function
that computed nothing.** With the summing loop deleted so the cohort came back all zeros, three
tests failed and that one stayed green.

The cause was in the helper. `one_pass` read every sample's winner off the posterior *before*
`sum_cohort_expected_copies` was ever called, and the cohort row it scored against was written
once before the loop — so **no genotype in the fixture was downstream of the M-step at all**.
Across eight mutations the reviewer found no defect this test caught that the four beside it did
not.

That matters beyond the one test. Spec §13 test 2 asks for two observables — permuted genotypes
and a bitwise sum — *precisely so that one covers what the other cannot*. As built only the
bitwise one was load-bearing, and the module read as though the ordering property were doubly
covered.

## 5. The Majors

**Defects in the function:**

1. **The shape check's `sample_count * allele_count` was a plain multiply.** It wraps in
   release and returns a silent all-zero cohort; `prepare_for_locus` guards the identical
   product with `checked_mul`.
2. **`fill(0.0)` destroyed the `UNWRITTEN_SCRATCH_VALUE` sentinel before the walk**, so a walk
   that summed no rows would hand back a row of zeros reading as a real summary — the shape
   step D1 would introduce by skipping set-aside rows.

**A claim in the code that was false:**

3. **"one check over the alleles catches an unwritten row anywhere in the table" holds on the
   first pass only.** The sentinel is written by `prepare_for_locus`, once per *locus*, and the
   per-sample rows are deliberately carried across passes because `score_one_sample` reads each
   sample's previous copies as its leave-one-out term. From pass 2, a row this pass did not
   write holds the previous pass's finite, plausible value and is summed as current — measured
   `[1.0, 3.0]` where pass 1 gave `[2.0, 2.0]`. Finite, non-negative, wrong, silent.

**Tests that could not fail:**

4. **Two of the output check's three conditions were killed by nothing** — dropping
   `&& *total >= 0.0`, and weakening `is_finite()` to `!is_nan()`, each left the suite green.
5. **Both `> 0` guards were untested**, and a sample count of zero *satisfies* the size check,
   since `0 == 0 × alleles`.
6. **The permutation fixture's disagreement guard was weaker than its own comment.** Winners
   were `[0, 1, 0]`: the guard `a != b || b != c` passed while samples 0 and 2 called the same
   genotype, so the test discriminated on sample 1 alone.
7. **The bitwise oracle discriminated on the reference allele only**, and its three-sample
   column separated ascending order from *reversal* and little else.
8. **No test at one sample, and none above three**, against a range commitment of 1 to several
   thousand.
9. **No property test**, for a function whose entire contract is that its output is one
   specific fold in one specific order.
10. **`sample_count` looked redundant and is not.** Deriving it from the table's length makes
    any table divide evenly, so a two-row table presented as three samples is accepted and the
    cohort comes back a sample short. Recorded as a Major so the redundancy is not "tidied"
    later.

## 6. Two corrections the review made to itself

Worth recording, because both are the review working as intended.

**A reviewer filed a Major and then withdrew it.** It predicted that swapping the first two rows
of the sum would slip past the ordering fixture. It does slip past — but IEEE addition is
commutative, so `t = a; t += b` and `t = b; t += a` are bit-identical, and the mutation changes
no answer on any input (0 differences in 10,000 random columns, against 694 for moving row 0 to
last). **The orchestrator had independently run the same mutation and read its green result as
confirming the gap**; the reviewer's correction is what identified it as a no-op. The finding
was downgraded to the fixture's width, which is real.

**A reviewer reported an MCP notice it had not read**, flagged the fabrication itself in its next
message, and re-verified its report. Every figure in that report is backed by quoted command
output; the invented line was in a closing aside.

## 7. Hot path — measured, and it does not matter

`sum_cohort_expected_copies` settles at **0.40–0.50 ns an entry** from 50 samples up, against an
E-step measured on the same machine at **342.9 ns a call at 6 alleles**. At the tomato cohort's
63 samples over 6 alleles that is 161 ns of M-step against 21,600 ns of E-step — **1 part in
134**. From `objdump`: no `panic_bounds_check` anywhere (which is why an indexed rewrite is 39%
*slower* at two alleles), a real NEON inner loop that engages at four alleles and up, and asserts
costing about 0.8 ns at worst. **Nothing here argues for removing any check** — which is what
makes the untested-assert findings a matter of testing rather than of cost.

## 8. The range commitment

Measured against Neumaier compensated summation: the worst relative error is **4 parts in 100
trillion at 1,000 samples** and 2 parts in a trillion at 100,000 — nine orders of magnitude under
§6's `1e-3` convergence threshold. **No accuracy finding.**

One corner recorded rather than fixed: with one sample carrying 2.0 copies and 4,999 carrying
`1e-17`, the fixed-order sum returns `2.0` exactly and the entire contribution of 4,999 samples
vanishes. Immaterial at `5e-14`, and undocumented until now.

## 9. Out of scope

- The four `--release` failures in `likelihood/` were fixed in `06a9fac9` before this review ran,
  so `cargo test --release --lib ng::calling` was green throughout and could be used as a clean
  signal — which is how the "both release-held checks are killed under `--release`" answer was
  obtained.

## 10. What's good

The arithmetic and the determinism argument, which every agent that checked them found right. The
author's own reversed-sum mutation reproduced verbatim in two independent worktrees. And of 16
quantitative claims recomputed rather than read, eleven were exact — including the whole `2⁻⁵³`
story, which is the fiddliest thing in the diff.
