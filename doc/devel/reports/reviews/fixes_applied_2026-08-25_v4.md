# Fix Application Report: ng_calling_loop_b2_2026-08-25.md

**Date:** 2026-08-25
**Branch:** `ng-calling-loop`, worktree `../pop_var_caller-calling-loop`
**Review:** [ng_calling_loop_b2_2026-08-25.md](ng_calling_loop_b2_2026-08-25.md)

## 1. Executive summary

**1 Blocker, 10 Majors, 12 Minors. The Blocker and all 10 Majors are fixed.** Tests went from
five to **thirteen**; the library target went 4,589 → **4,596**.

### What changed in the code, as against in the tests

Three behaviour changes:

1. **The shape product is `checked_mul`**, matching `prepare_for_locus`'s guard on the identical
   product. A plain multiply wraps in release, and a wrapped product that happened to equal the
   table's length would let the sum run over a shape nobody asked for.
2. **The cohort row is seeded from the first sample's row, not from `fill(0.0)`.** Bit-identical
   — `0.0 + x` is exactly `x` for every value here — but it makes an empty walk
   unrepresentable, where the zeroed version would hand back a plausible row of zeros. That is
   the shape D1 would introduce by skipping spec §5.0's set-aside rows.
3. **The `# Panics` block documents all four release-held checks**, where it listed two. The
   two it omitted are the `> 0` guards, and one of them — `sample_count > 0` — is the only
   thing between a zero count and an all-zero cohort row, because `0 == 0 × alleles` satisfies
   the size check.

And one correction to a claim that was false:

4. **The `NaN` sentinel's guarantee holds on the first pass at a locus only.** The comment said
   one check over the alleles catches an unwritten row "anywhere in the table". It does — until
   pass 2, after which a row this pass did not write holds the previous pass's finite value and
   is summed as current. The comment now says so, names the measurement, and names D1 as owing
   the per-pass written mask a skipping loop would need.

### The Blocker's fix, and how it was verified

`one_pass` became `run_passes(n, …)`, and the permutation test takes its winners after **two**
passes, so the genotypes it asserts on are scored against the row the first pass's M-step
produced.

**Verified by the mutation that exposed the defect.** With the summing loop deleted:

- before the fix, `presenting_the_samples_in_another_order_calls_the_same_genotypes` **passed**
  while three other tests failed;
- after it, that test **fails**, and nine tests fail in total.

### Validation

| command | result |
|---|---|
| `cargo fmt --all -- --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo test --lib` | `4596 passed; 0 failed; 14 ignored` |
| `cargo test --release --lib ng::calling --all-features` | `550 passed; 0 failed` |

**All four release-held checks are reached by a test that fails under `--release` without
them** — measured with all four downgraded together: `14 passed; 6 failed`, the six being the
`should_panic` tests across the four checks.

## 2. Per-finding log

### BL1 — the permutation test could not fail — **Fixed**

Above. The helper's doc now carries the measurement, so a later reader cannot quietly turn it
back into a one-pass helper.

### M1 — unchecked multiply — **Fixed**

`checked_mul` with a message naming the shape, and a comment saying `prepare_for_locus` guards
the same product the same way.

### M2 — `fill(0.0)` before the walk — **Fixed**

Seeded from the first row via an iterator whose `next()` the shape check guarantees. The comment
records that the two are bit-identical and that the change is about representability, not
arithmetic.

### M3 — the sentinel's first-pass-only guarantee — **Fixed (comment + test)**

`from_the_second_pass_an_unwritten_row_is_summed_as_current` pins the behaviour, asserting the
measured `[1.0, 3.0]` against pass 1's `[2.0, 2.0]`. Its doc says plainly that this is not a bug
in the M-step, and that if a later step makes mid-run skipping impossible the test should start
failing and be deleted with a note.

### M4 — two output-check conditions untested — **Fixed**

`an_infinite_total_is_refused` reaches the `is_finite` arm (via `f64::MAX + f64::MAX`), and
`a_negative_total_is_refused` reaches the `>= 0.0` arm. The latter's doc also carries the
check's **limit**: it catches a negative *total*, not negative inputs that cancel — two rows of
`-1.0` and `3.0` sum to an acceptable `2.0`. Checking every input would cost
`samples × alleles` on the hot path to catch a state no producer can reach.

### M5 — both `> 0` guards untested — **Fixed**

`a_cohort_of_no_samples_is_refused` and `a_locus_of_no_alleles_is_refused`. Independently
confirmed before fixing: with both guards deleted, `ng::calling` was `556 passed; 0 failed`.

### M6 — the disagreement guard was weaker than its comment — **Fixed**

The guard now requires **all three** samples to call different genotypes, and the likelihoods
were retuned until they do. The old fixture called `[0, 1, 0]` and passed a guard allowing it,
so the test discriminated on one sample.

### M7 — the bitwise oracle was narrow — **Fixed**

Four samples, reference column `[2⁻⁵³, 1.0, 2⁻⁵³, 2⁻⁵²]`, chosen by searching the space of
power-of-two columns. **The fixture now certifies itself**: the test builds reversal and every
adjacent transposition and asserts they differ from ascending.

**With one exception that is a fact rather than a gap.** Swapping the *first two* samples is
bit-identical on every input, IEEE addition being commutative, so no fixture of any shape can
separate it and no implementation can be wrong by making it. The test asserts that
equality explicitly, so the exception is documented in executable form rather than in prose.

**And the residual limit is stated:** 3 of the 23 non-identity permutations of four samples sum
bit-identically to ascending. No column separates every permutation.

### M8 — no test at one sample — **Fixed**

`the_cohort_of_one_sample_is_that_samples_own_row_bit_for_bit`. At one sample the M-step is not
an approximation of a sum, it *is* the row, so the assertion is on the bits.

### M9 — no property test — **Fixed**

A `proptest` over `1..24` samples × `1..8` alleles with copies in `0.0..=2.0`, comparing on the
bits against an ascending fold written independently in the test. `proptest` was already a
dev-dependency.

### M10 — `sample_count` looks redundant — **Fixed (doc)**

Kept, with the field's doc now saying why: derive it and any table divides evenly, so a two-row
table presented as three samples is accepted and the cohort comes back a sample short. The
reviewer measured that — with the count derived, the wrong-size test starts accepting the table
and returning `[2.0, 2.0]`.

### Minors applied

- **The "21 genotypes" attribution was wrong.** The paragraph said "measured on this module's
  own fixture" and then quoted 21; no fixture here has 21 genotypes — `generic_locus` tops out
  at 4 alleles → 10, and the B2 tests use 3 or 6. 21 is the spec's six-allele example. The
  sentence now names both.
- **The ploidy test asserted only the grand total**, which any allele-permuting M-step passes.
  It now asserts per-allele properties too.
- **`CohortSumBuffers` → `CohortSummingBuffers`** and `cohort_sum_buffers_mut` →
  `cohort_summing_buffers_mut`, so the pair scans beside B1's `SampleScoringBuffers` /
  `sample_scoring_buffers_mut`.
- **The M-step is now defined in `calling/mod.rs`**, where the type's doc uses the term to
  carry an argument; it was defined only in the other file.
- **`one_pass` → `run_passes`** — a verb, like the module's other functions.
- **`forward`/`backward` → `ascending` and an explicit fold**, because they held the *reference
  allele's* total and the doc called it "the total".
- **`straight`/`rotated` → `called_in_order`/`called_rotated`** and the cohort rows likewise.

### Not applied

- **A negative-input check on every entry.** It would cost `samples × alleles` on the loop's hot
  path to catch a state no producer can reach. The output check's limit is documented instead,
  on the test that reaches it.

## 3. Carried forward

1. **Spec §5.0's set-aside samples** — D1's, and now with a second reason attached: a loop that
   skips a sample mid-run needs a **per-pass written mask**, because the scratch's sentinel is
   armed once per locus and not once per pass.
2. **C2 needs no third bundle** — both cohort accessors are `&self`, so the convergence delta
   compiles. But **three of four natural spellings of the max-delta falsely converge on pass 1**
   against the `NaN` sentinel, because `f64::max` discards `NaN`; only `.all(|d| d < tol)`
   behaves as `advance_cohort_expected_copies`'s doc promises. C2 should start from that.
3. **The vanishing-tail corner**: one sample at 2.0 copies and 4,999 at `1e-17` returns `2.0`
   exactly, the 4,999 contributing nothing. Immaterial at `5e-14` against a `1e-3` threshold,
   recorded so it is not rediscovered as a bug.

## 4. One artifact, and where it came from

The new property test wrote
`proptest-regressions/ng/calling/inference/summarise_condition.txt`, holding a seed that shrank
to `sample_count = 1, allele_count = 1`. **It records a failure under the deliberate
summing-loop mutation, not a defect** — the test passes on every tree that was ever committed,
and the seed was saved while that mutation was in the working tree.

It is committed, because the repository already tracks seven such files and proptest's own
header asks for it, and because the case it names is a legitimate one to re-run first. **This
paragraph exists so the next reader does not take it for evidence of a past real failure.** The
same case has an explicit test of its own in
`the_cohort_of_one_sample_is_that_samples_own_row_bit_for_bit`.
