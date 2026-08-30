# ng calling loop — B2: the M-step, in the run's fixed sample order

**Date:** 2026-08-25
**Plan:** [calling_loop.md](../../ng/impl_plan/calling_loop.md), step B2
**Design authority:** [spec/calling_em_loop.md](../../ng/spec/calling_em_loop.md) §2, §5, §5.0,
§8, §13 test 2; [arch/calling_em_loop.md](../../ng/arch/calling_em_loop.md) §1, §4
**Branch:** `ng-calling-loop`, worktree `../pop_var_caller-calling-loop`

> **Read this against the review that followed it.** Three category reviews raised one Blocker
> and ten Majors, and the fixes changed what §3 and §4 describe: the shape product became
> `checked_mul`, the cohort row is seeded from the first sample rather than zeroed, the tests
> went from five to thirteen, and the bitwise oracle was rebuilt on four samples because the
> three-sample one separated ascending order from reversal and little else. **Two claims in
> this report were wrong**: the `NaN` sentinel's guarantee holds on the first pass at a locus
> only, and the permutation test could not fail on any defect in the function.
> [The review](../reviews/ng_calling_loop_b2_2026-08-25.md), and
> [what was done about it](../reviews/fixes_applied_2026-08-25_v4.md).

---

## 1. Plan

The other half of one pass. B1 turned the cohort's summary into each sample's genotype
probabilities; `sum_cohort_expected_copies` turns those back into the summary the next pass
conditions on. It is the whole of the M-step — no second quantity moves while the frequency
loop runs (spec §5's table).

**Its own commit, because the defect it can carry is silent.** A sum in the wrong order gives
a different answer in the last bits at a different worker count, and never crashes.

## 2. Assumptions and deviations

### 2.1 The function takes the whole table, not a row at a time

The plan says *"a sum over samples in the run's fixed sample order"*. A signature taking one
row per call would put the order in the caller's loop, where a later step that parallelised
over samples would change it without touching this function. So `CohortSumBuffers` carries the
whole `samples × alleles` table and the walk is here.

### 2.2 A second borrow bundle, and the review predicted it

`CohortSumBuffers` is `SampleScoringBuffers`' counterpart for the same reason: the per-sample
copies and the cohort row are two fields of one `CallingScratch`, and each per-buffer accessor
borrows the whole of it. B1's review measured this in advance — an M-step written against the
per-buffer accessors gives two `error[E0502]`. **That there are now two bundles is a fact about
A1's private-fields decision, and it is recorded in `PROJECT_STATUS.md` as something
`arch/calling_em_loop.md` §2 should absorb.**

### 2.3 Spec §5.0's set-aside sample is *not* implemented here, and fails loudly

A sample whose own reads earned an allele the candidate cap cut must contribute nothing to this
sum: its posterior would sit over the wrong allele set and would pull the locus's frequencies
toward the reference by exactly the samples carrying the rarest alleles. `LocusEvidence::Generic`
carries the flag and keeps such a sample's **index**, so the table has a row for it.

What happens today is that the row is never written, so it still holds the scratch's `NaN`
sentinel, and the release-held finiteness check refuses the locus. **That is the right failure
while the exclusion is unbuilt** — the alternative is a quietly wrong cohort frequency — but it
is not the answer. Step D1 assembles the loop and owns the choice between skipping those rows
here and never giving them one. The function's doc says so, and
`the_cohort_carries_the_ploidy_for_every_sample_it_summed` is written so that whoever builds it
has to edit a test that states the old identity.

## 3. Changes made

- **[summarise_condition.rs](../../../../src/ng/calling/inference/summarise_condition.rs)** —
  `sum_cohort_expected_copies`, and the `one_pass` test helper that runs a whole pass.
- **[calling/mod.rs](../../../../src/ng/calling/mod.rs)** — `CohortSumBuffers` and
  `CallingScratch::cohort_sum_buffers_mut`.

**Four checks held in release** *(amended after the review, which found the report and the
`# Panics` block both listing two)*: the cohort row names at least one allele; the sample count
is at least one; the table is `sample_count × alleles`, with the product taken by `checked_mul`;
and every cohort entry comes out finite and non-negative.

The last is the load-bearing one and it is cheap for a reason worth stating: `NaN` survives every
addition, so one pass over the alleles catches an unwritten sample row anywhere in a table of any
size — `alleles` work, not `samples × alleles`. **On the first pass at a locus only**, which this
report originally missed: the sentinel is armed by `prepare_for_locus`, once per *locus*, and the
per-sample rows are deliberately carried across passes because the E-step reads each sample's
previous copies as its leave-one-out term. From pass 2 an unwritten row holds the previous pass's
finite value and is summed as current — measured `[1.0, 3.0]` where pass 1 gave `[2.0, 2.0]`.

The two `> 0` guards are not belt-and-braces: a sample count of zero **satisfies** the size
check, since `0 == 0 × alleles`.

## 4. Tests added

Five, plus the `one_pass` helper.

| test | what it pins |
|---|---|
| `presenting_the_samples_in_another_order_calls_the_same_genotypes` | the semantic half of spec §13 test 2 — which sample sits at which index is an accident of assembly and no call may depend on it |
| `the_sum_runs_in_ascending_sample_order` | the **bitwise** half, and the mutation oracle |
| `the_cohort_carries_the_ploidy_for_every_sample_it_summed` | a dropped or double-counted sample |
| `a_sample_row_that_was_never_written_is_refused` | the sentinel reaching the cohort summary |
| `a_per_sample_table_of_the_wrong_size_is_refused` | a table one row short summing `n − 1` and reporting it as `n` |

**The mutation, and it is the plan's own point demonstrated.** With the sum reversed
(`.chunks_exact(allele_count).rev()`), **only** `the_sum_runs_in_ascending_sample_order` fails —
`left: 4607182418800017409  right: 4607182418800017408`, one bit — while both
genotype-observable tests pass. That is why the plan insists the mutation check is on the
summed copies compared bitwise and not on the argmax.

**The bitwise fixture computes both orders rather than quoting them**, so the claim that they
differ cannot go stale: three samples carrying `1.0`, `2⁻⁵³` and `2⁻⁵³` copies of the reference
allele sum forward to `1.0`, because each tiny addend rounds away against the one, and backward
to `1.0000000000000002`, because the two tiny ones meet each other first.

## 5. Validation

- `cargo fmt --all -- --check` — exit 0.
- `cargo clippy --all-targets --all-features -- -D warnings` — exit 0.
- `cargo test --lib` — `4596 passed; 0 failed; 14 ignored`. Before B2: **4,584**, so B2 adds
  twelve: five as written and seven the review added.
- `cargo test --release --lib ng::calling --all-features` — `550 passed; 0 failed`, the run that
  can tell `assert!` from `debug_assert!` and a CI step since `06a9fac9`.

The two aggregate gates that are red are the standing pre-existing ones
(`benches/psp_writer_perf.rs:386` and `--example ng_generic_loci_dump`), both verified against
an unpatched tree during B1.

## 6. Trade-offs and follow-ups

- **§5.0's exclusion, above** — D1's.
- **The convergence test (C2) needs a third view of the scratch**: the current and previous
  cohort rows together. Whether that is a third bundle or a method on the scratch is C2's to
  settle; two bundles is not yet a pattern.
- **Nothing here is exercised above 3 samples.** The sum is linear and the fixed order is what
  it is at any size, but the accuracy of a 1,000-addend sum is not measured.
